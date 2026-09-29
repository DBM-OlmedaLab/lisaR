# Concurrency, process-lifetime and export regression tests.
# These use synthetic job records without rendering figures.

explore_repair_ws <- function(name) {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), paste0("lisa-explore-repair-", name))
  unlink(root, recursive = TRUE, force = TRUE)
  lisa_explore_open(run, root)
}

# --- finding 4: process identity, lock ownership ----------------------------

test_that("process identity carries a start time and detects PID reuse", {
  self <- lisaR:::lisa_explore_process_identity()
  expect_identical(self$pid, as.integer(Sys.getpid()))
  skip_if(is.na(self$start_time), "no process start time available on this host")

  expect_true(lisaR:::lisa_explore_identity_alive(self))

  # The same live PID with a different start time is a REUSED number: the
  # process that recorded it is gone, so the identity must read as dead.
  expect_false(lisaR:::lisa_explore_identity_alive(
    list(host = self$host, pid = self$pid, start_time = "1")))

  # An unknown start time degrades to the old PID-only answer rather than
  # guessing.
  expect_true(lisaR:::lisa_explore_identity_alive(
    list(host = self$host, pid = self$pid, start_time = NA_character_)))

  dead <- 4000000L
  expect_false(lisaR:::lisa_explore_identity_alive(
    list(host = self$host, pid = dead, start_time = "1")))

  # A foreign host cannot be judged, so it is reported alive (fail-closed): we
  # never steal a lock or declare a job dead on evidence we do not have.
  expect_true(lisaR:::lisa_explore_identity_alive(
    list(host = "some-other-host", pid = dead, start_time = "1")))
})

test_that("a lock is released only by its owner", {
  ws <- explore_repair_ws("lock")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  lock <- lisaR:::lisa_explore_lock_path(ws)

  token <- lisaR:::lisa_explore_acquire_lock(ws)
  owner <- lisaR:::lisa_explore_read_json(file.path(lock, "owner.json"))
  expect_true(all(c("host", "pid", "start_time", "token") %in% names(owner)))
  expect_identical(as.character(owner$token), as.character(token))

  # Releasing with someone else's token must not hand the index to two writers.
  lisaR:::lisa_explore_release_lock(ws, token = "a-different-token")
  expect_true(dir.exists(lock))

  lisaR:::lisa_explore_release_lock(ws, token = token)
  expect_false(dir.exists(lock))
})

test_that("a stale lock is reclaimed and leaves nothing behind", {
  ws <- explore_repair_ws("reclaim")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  lock <- lisaR:::lisa_explore_lock_path(ws)

  dir.create(lock, showWarnings = FALSE)
  lisaR:::lisa_explore_write_json(
    list(host = lisaR:::lisa_explore_host(), pid = 4000000L, start_time = "1",
         token = "stale", acquired_at = lisaR:::lisa_explore_now()),
    file.path(lock, "owner.json"), run_root = ws$root)

  token <- lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 5)
  expect_true(nzchar(token))
  lisaR:::lisa_explore_release_lock(ws, token)
  expect_identical(
    list.files(ws$root, pattern = "^[.]index-lock[.]stale-", all.files = TRUE),
    character(0))
})

test_that("a live lock is respected", {
  ws <- explore_repair_ws("livelock")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  token <- lisaR:::lisa_explore_acquire_lock(ws)
  on.exit(lisaR:::lisa_explore_release_lock(ws, token), add = TRUE)
  expect_error(lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0.2),
               "LISA-EXPLORE-015")
})

explore_write_partial_owner <- function(ws, owner, suffix = "fixture") {
  lock <- lisaR:::lisa_explore_lock_path(ws)
  dir.create(lock, showWarnings = FALSE, mode = "0700")
  path <- file.path(lock, paste0(".lisa-write-", suffix, ".json"))
  lisaR:::lisa_explore_write_json(owner, path, run_root = ws$root)
  path
}

test_that("an interrupted complete owner promotion is recovered only when dead", {
  ws <- explore_repair_ws("partial-dead")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  partial <- explore_write_partial_owner(ws, list(
    host = lisaR:::lisa_explore_host(), pid = 4000000L, start_time = "1",
    token = "partial-dead", acquired_at = lisaR:::lisa_explore_now()))
  expect_true(file.exists(partial))
  expect_false(file.exists(file.path(dirname(partial), "owner.json")))

  token <- lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0.2)
  expect_true(nzchar(token))
  expect_false(identical(token, "partial-dead"))
  expect_true(lisaR:::lisa_explore_release_lock(ws, token))
  expect_false(dir.exists(lisaR:::lisa_explore_lock_path(ws)))
})

test_that("partial live, foreign and unknown ownership is never stolen", {
  identities <- list(
    live = lisaR:::lisa_explore_process_identity(),
    foreign = list(host = "foreign-host", pid = 4000000L, start_time = "1"))
  for (kind in names(identities)) {
    ws <- explore_repair_ws(paste0("partial-", kind))
    lock <- lisaR:::lisa_explore_lock_path(ws)
    owner <- c(identities[[kind]],
               list(token = paste0("partial-", kind),
                    acquired_at = lisaR:::lisa_explore_now()))
    partial <- explore_write_partial_owner(ws, owner, suffix = kind)
    expect_error(lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0),
                 "LISA-EXPLORE-015.*recovery was refused")
    expect_true(dir.exists(lock))
    expect_true(file.exists(partial))
    expect_false(lisaR:::lisa_explore_release_lock(ws, "not-the-owner"))
    expect_true(dir.exists(lock))
    unlink(ws$root, recursive = TRUE, force = TRUE)
  }

  ws <- explore_repair_ws("partial-malformed")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  lock <- lisaR:::lisa_explore_lock_path(ws)
  dir.create(lock, mode = "0700")
  malformed <- file.path(lock, ".lisa-write-malformed.json")
  writeLines("{", malformed, useBytes = TRUE)
  expect_error(lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0),
               "LISA-EXPLORE-015.*malformed JSON")
  expect_true(file.exists(malformed))
})

test_that("a reused PID in a complete partial owner record is reclaimed", {
  self <- lisaR:::lisa_explore_process_identity()
  skip_if(is.na(self$start_time), "no process start time available on this host")
  ws <- explore_repair_ws("partial-reused-pid")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  explore_write_partial_owner(ws, c(self[c("host", "pid")],
    list(start_time = "1", token = "reused-pid",
         acquired_at = lisaR:::lisa_explore_now())))
  token <- lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0.2)
  expect_true(lisaR:::lisa_explore_release_lock(ws, token))
})

test_that("token races restore changed ownership instead of deleting it", {
  ws <- explore_repair_ws("lock-token-race")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  lock <- lisaR:::lisa_explore_lock_path(ws)
  explore_write_partial_owner(ws, list(
    host = lisaR:::lisa_explore_host(), pid = 4000000L, start_time = "1",
    token = "replacement", acquired_at = lisaR:::lisa_explore_now()))
  expect_false(lisaR:::lisa_explore_reclaim_stale_lock(lock, "inspected-earlier"))
  evidence <- lisaR:::lisa_explore_lock_evidence(lock)
  expect_true(evidence$verified)
  expect_identical(as.character(evidence$owner$token), "replacement")
})

test_that("failed owner initialization removes only its own claimed lock", {
  ws <- explore_repair_ws("lock-init-failure")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  original_write <- lisaR:::lisa_explore_write_json
  testthat::local_mocked_bindings(
    lisa_explore_write_json = function(object, path, run_root) {
      if (identical(basename(path), "owner.json")) {
        original_write(object,
          file.path(dirname(path), ".lisa-write-simulated.json"), run_root)
        stop("simulated owner promotion failure", call. = FALSE)
      }
      original_write(object, path, run_root)
    },
    .package = "lisaR")
  expect_error(lisaR:::lisa_explore_acquire_lock(ws),
               "simulated owner promotion failure")
  expect_false(dir.exists(lisaR:::lisa_explore_lock_path(ws)))
  expect_identical(
    list.files(ws$root, pattern = "^[.]index-lock[.](released|stale)-",
               all.files = TRUE),
    character(0))
})

test_that("an acquired owner with unavailable start time can still release itself", {
  ws <- explore_repair_ws("lock-unknown-start")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  testthat::local_mocked_bindings(
    lisa_explore_process_identity = function(pid = Sys.getpid())
      list(host = lisaR:::lisa_explore_host(), pid = as.integer(pid),
           start_time = NA_character_),
    .package = "lisaR")
  token <- lisaR:::lisa_explore_acquire_lock(ws)
  expect_true(lisaR:::lisa_explore_release_lock(ws, token))
  expect_false(dir.exists(lisaR:::lisa_explore_lock_path(ws)))
})

# --- finding 2: one costly worker across distinct requests ------------------

`%||%` <- function(x, y) if (is.null(x)) y else x

explore_fake_job <- function(ws, job_id, state, launched, background = TRUE,
                             pid = Sys.getpid(), created_at = NULL, ...) {
  # modifyList, not c(): concatenating would leave two entries for any field the
  # caller overrides, and `job$field` would silently return the default one.
  job <- list(
    job_id = job_id, key = paste0("key-", job_id), request_id = job_id,
    request = list(), source_run = ws$source_run, root = ws$root,
    state = state, stage = state, cancel_requested = FALSE,
    background = background, launched = launched,
    host = lisaR:::lisa_explore_host(), pid = pid,
    start_time = lisaR:::lisa_explore_proc_start_time(pid),
    created_at = created_at %||% lisaR:::lisa_explore_now(), message = "")
  overrides <- list(...)
  if (length(overrides)) job <- utils::modifyList(job, overrides)
  lisaR:::lisa_explore_write_job(ws, job)
}

test_that("the worker slot is held by a launched, live job only", {
  ws <- explore_repair_ws("slot")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)

  expect_null(lisaR:::lisa_explore_slot_holder(ws))

  # Queued but never dispatched: it is waiting, not occupying.
  explore_fake_job(ws, "job-waiting", "queued", launched = FALSE)
  expect_null(lisaR:::lisa_explore_slot_holder(ws))

  # Running under a live process: it holds the slot.
  explore_fake_job(ws, "job-live", "running", launched = TRUE)
  holder <- lisaR:::lisa_explore_slot_holder(ws)
  expect_identical(as.character(holder$job_id), "job-live")

  # A dead owner does not keep the slot hostage.
  explore_fake_job(ws, "job-live", "running", launched = TRUE, pid = 4000000L)
  expect_null(lisaR:::lisa_explore_slot_holder(ws))
})

test_that("the slot is claimed exclusively and the queue is fair", {
  ws <- explore_repair_ws("claim")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)

  explore_fake_job(ws, "job-b", "queued", launched = FALSE,
                   created_at = "2026-01-01T00:00:02+0000")
  explore_fake_job(ws, "job-a", "queued", launched = FALSE,
                   created_at = "2026-01-01T00:00:01+0000")

  # Oldest first, regardless of file order.
  expect_identical(
    as.character(lisaR:::lisa_explore_next_dispatchable(ws)$job_id), "job-a")

  expect_true(lisaR:::lisa_explore_claim_slot(ws, "job-a")$claimed)
  expect_true(isTRUE(lisaR:::lisa_explore_read_job(ws, "job-a")$launched))

  # With the slot taken, a second claim is refused and says by whom.
  second <- lisaR:::lisa_explore_claim_slot(ws, "job-b")
  expect_false(second$claimed)
  expect_identical(as.character(second$blocked_by), "job-a")
  expect_false(isTRUE(lisaR:::lisa_explore_read_job(ws, "job-b")$launched))
})

test_that("a cancelled job is never dispatched", {
  ws <- explore_repair_ws("cancelled")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  explore_fake_job(ws, "job-x", "queued", launched = FALSE,
                   cancel_requested = TRUE)
  expect_null(lisaR:::lisa_explore_next_dispatchable(ws))
  expect_false(lisaR:::lisa_explore_claim_slot(ws, "job-x")$claimed)
  expect_null(lisa_explore_pump(ws))
})

test_that("a foreground submit is refused while another worker holds the slot", {
  ws <- explore_repair_ws("foreground")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  explore_fake_job(ws, "job-busy", "running", launched = TRUE)

  request <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOMF",
    category_id = "SYN_SIGNAL", product = "volcano")
  traced <- explore_with_generator_trace(
    expect_error(lisa_explore_submit(ws, request, background = FALSE),
                 "LISA-EXPLORE-036"))
  # Refused, not run: no generator was reached.
  expect_identical(length(traced$post), 0L)
})

test_that("an undispatched job has no owning process to orphan", {
  ws <- explore_repair_ws("undispatched")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  explore_fake_job(ws, "job-waiting", "queued", launched = FALSE,
                   pid = NA_integer_)
  job <- lisaR:::lisa_explore_read_job(ws, "job-waiting")
  reconciled <- lisaR:::lisa_explore_reconcile_job(ws, job)
  expect_identical(reconciled$state, "queued")
  expect_match(reconciled$message, "waiting for the single worker slot")
})

test_that("a job whose dispatching session died is failed, freeing the slot", {
  ws <- explore_repair_ws("dispatcher")
  on.exit(unlink(ws$root, recursive = TRUE, force = TRUE), add = TRUE)
  explore_fake_job(
    ws, "job-lost", "queued", launched = TRUE, pid = NA_integer_,
    dispatcher = list(host = lisaR:::lisa_explore_host(), pid = 4000000L,
                      start_time = "1"))
  job <- lisaR:::lisa_explore_read_job(ws, "job-lost")
  expect_identical(lisaR:::lisa_explore_reconcile_job(ws, job)$state, "failed")
  expect_null(lisaR:::lisa_explore_slot_holder(ws))
})

# --- finding 3: the background handle must survive garbage collection -------

test_that("the process registry keeps a strong reference and reaps finished workers", {
  skip_if_not_installed("callr")
  marker <- tempfile(fileext = ".done")
  process <- callr::r_bg(
    func = function(marker) { Sys.sleep(3); file.create(marker) },
    args = list(marker = marker), supervise = FALSE)
  lisaR:::lisa_explore_register_process("job-gc", process)
  rm(process)

  # Dropping every local reference and forcing a collection must NOT kill the
  # worker. Without the registry, processx's cleanup finaliser does exactly that.
  gc(full = TRUE); gc(full = TRUE)
  expect_false(is.null(lisaR:::lisa_explore_process_handle("job-gc")))
  lisaR:::lisa_explore_process_handle("job-gc")$wait(timeout = 30000)
  expect_true(file.exists(marker))

  # Finished handles are released, so a long session does not accumulate them.
  lisaR:::lisa_explore_reap_processes()
  expect_null(lisaR:::lisa_explore_process_handle("job-gc"))
  unlink(marker)
})

test_that("callr children are registered for cleanup, which is why the registry exists", {
  skip_if_not_installed("callr")
  # Documents the measured reason for the repair: `cleanup` cannot be reached
  # through r_bg, so every child is GC-killable and the handle must be retained.
  expect_false("cleanup" %in% names(formals(callr::r_bg)))
  expect_false("cleanup" %in% names(callr::r_process_options()))
  expect_true(isTRUE(formals(
    processx::process$public_methods$initialize)$cleanup))
})

# --- finding 5: stale entries are not current -------------------------------

test_that("an entry from a changed source run is not current", {
  shared <- explore_shared()
  entries <- lisaR:::lisa_explore_read_index(shared$ws)
  expect_gte(length(entries), 1L)
  entry <- entries[[1L]]
  expect_true(lisaR:::lisa_explore_entry_is_current(shared$ws, entry))

  stale <- entry
  stale$source_manifest_hash <- strrep("0", 64)
  currency <- lisaR:::lisa_explore_entry_currency(shared$ws, stale)
  expect_false(currency$current)
  expect_match(currency$reason, "source run changed")
})

test_that("an entry whose recorded key no longer matches today's key is not current", {
  shared <- explore_shared()
  entry <- lisaR:::lisa_explore_read_index(shared$ws)[[1L]]
  stale <- entry
  stale$key <- paste0("fig-", strrep("a", 24))
  currency <- lisaR:::lisa_explore_entry_currency(shared$ws, stale)
  expect_false(currency$current)
  expect_match(currency$reason, "engine or resources changed")
})

test_that("an entry for a figure the run no longer offers is not current", {
  shared <- explore_shared()
  entry <- lisaR:::lisa_explore_read_index(shared$ws)[[1L]]
  stale <- entry
  stale$request$category_id <- "NO_SUCH_CATEGORY"
  currency <- lisaR:::lisa_explore_entry_currency(shared$ws, stale)
  expect_false(currency$current)
  expect_match(currency$reason, "no longer offers")
})

test_that("file hashes alone do not make an entry current", {
  # The point of finding 5: lisa_explore_validate_artifact() can pass while the
  # entry is stale, because it only proves the files are unchanged since they
  # were recorded -- not that the engine that made them is still the current one.
  shared <- explore_shared()
  entry <- lisaR:::lisa_explore_read_index(shared$ws)[[1L]]
  stale <- entry
  stale$source_manifest_hash <- strrep("0", 64)
  expect_true(lisaR:::lisa_explore_validate_artifact(shared$ws, stale)$ok)
  expect_false(lisaR:::lisa_explore_entry_is_current(shared$ws, stale))
})

# --- finding 6: the reuse key is cheaper but unchanged in value -------------

test_that("the memoised key context produces identical keys", {
  shared <- explore_shared()
  ws <- shared$ws
  catalog <- lisa_explore_catalog(ws)
  rows <- seq_len(min(12L, nrow(catalog)))
  context <- lisaR:::lisa_explore_key_context(ws)
  for (index in rows) {
    row <- catalog[index, , drop = FALSE]
    request <- lisaR:::lisa_explore_request_from_row(row)
    expect_identical(lisaR:::lisa_explore_key(ws, request, row, context = context),
                     lisaR:::lisa_explore_key(ws, request, row),
                     info = as.character(row$category_id))
  }
})

test_that("the engine digest is computed once per distinct planned-script set", {
  shared <- explore_shared()
  ws <- shared$ws
  catalog <- lisa_explore_catalog(ws)
  volcano <- catalog[catalog$product == "volcano", , drop = FALSE]
  skip_if(nrow(volcano) < 3L, "not enough volcano rows in the fixture")
  context <- lisaR:::lisa_explore_key_context(ws)
  for (index in seq_len(3L)) {
    lisaR:::lisa_explore_context_engine_digest(context,
                                               volcano[index, , drop = FALSE])
  }
  # Three rows of the same product share one script set, so one cached digest.
  expect_identical(length(ls(envir = context$engine_digests)), 1L)
})


test_that("the portability report finds absolute references and symlinks", {
  root <- file.path(tempdir(), "lisa-explore-portability")
  unlink(root, recursive = TRUE, force = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(root, recursive = TRUE)

  writeLines('<html><body><a href="page2.html">ok</a></body></html>',
             file.path(root, "clean.html"))
  clean <- lisaR:::lisa_explore_portability_report(root)
  expect_identical(clean$absolute_reference_files, character(0))
  expect_identical(clean$symlinks, character(0))

  writeLines('<html><body><img src="/external/figure.png"></body></html>',
             file.path(root, "bad.html"))
  bad <- lisaR:::lisa_explore_portability_report(root)
  expect_identical(bad$absolute_reference_files, "bad.html")

  skip_on_os("windows")
  link_ok <- suppressWarnings(file.symlink(file.path(root, "clean.html"),
                                           file.path(root, "link.html")))
  skip_if_not(link_ok, "symlinks are not supported here")
  linked <- lisaR:::lisa_explore_portability_report(root)
  expect_true("link.html" %in% linked$symlinks)
})

test_that("a non-relocatable bundle is refused rather than promoted", {
  built <- explore_shared_bundle()
  # Reuse the already exported bundle as a source tree that DOES contain an
  # absolute reference, so the refusal is exercised without a further export of
  # the full fixture.
  root <- file.path(tempdir(), "lisa-explore-refusal")
  unlink(root, recursive = TRUE, force = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(root, recursive = TRUE)
  writeLines('<html><body><a href="/external/elsewhere.html">x</a></body></html>',
             file.path(root, "page.html"))

  report <- lisaR:::lisa_explore_portability_report(root)
  expect_identical(report$absolute_reference_files, "page.html")

  # And the exported shared bundle, which is clean, reports clean.
  expect_identical(built$bundle$html_files_with_absolute_references, character(0))
  expect_identical(built$bundle$symlinks, character(0))
})
