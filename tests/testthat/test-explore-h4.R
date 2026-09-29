
explore_h4_workspace <- function(name) {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), paste0("lisa-explore-h4-", name))
  unlink(root, recursive = TRUE, force = TRUE)
  withr::defer(unlink(root, recursive = TRUE, force = TRUE),
               envir = parent.frame())
  lisa_explore_open(run, root)
}

explore_h4_request <- function(category = "SYN_SIGNAL") {
  lisa_figure_request(unit_type = "single_de", analysis_id = "response_a",
    collection = "GOBP-C2", category_id = category, product = "volcano")
}

explore_h4_add_artifact <- function(ws, request) {
  plan <- lisa_explore_plan(ws, request)
  stopifnot(isTRUE(plan$applicable))
  artifact <- file.path(ws$root, "extensions", plan$key)
  product <- file.path(artifact, "artifacts", request$analysis_id,
                       request$collection, request$product)
  dir.create(product, recursive = TRUE, showWarnings = FALSE)
  writeLines("png", file.path(product, "figure.png"), useBytes = TRUE)
  writeLines("pdf", file.path(product, "figure.pdf"), useBytes = TRUE)
  writeLines("symbol\tlog2FC\nGENE\t1", file.path(product, "figure_source.tsv"),
             useBytes = TRUE)
  writeLines("# reproducible fixture recipe", file.path(product, "figure_recipe.R"),
             useBytes = TRUE)
  files <- lisaR:::lisa_explore_collect_artifact_files(artifact)
  entries <- lisaR:::lisa_explore_index_upsert(
    lisaR:::lisa_explore_read_index(ws), list(
      key = plan$key, request_id = plan$request_id,
      request = as.list(lisaR:::lisa_explore_request_vector(request)),
      associations = as.list(request$category_id), state = "available",
      job_id = paste0("fixture-", plan$request_id),
      source_manifest_hash = ws$source_manifest_hash,
      updated_at = lisaR:::lisa_explore_now(), message = "",
      artifact_dir = plan$key, files = files))
  lisaR:::lisa_explore_write_index(ws, entries)
  plan$key
}

explore_h4_hold_slot <- function(ws) {
  job <- list(kind = "figure", job_id = "job-h4-holder",
    key = "fig-h4-holder", request_id = "req-h4-holder", request = list(),
    source_run = ws$source_run, root = ws$root, state = "queued",
    stage = "queued", cancel_requested = FALSE, background = TRUE,
    launched = TRUE, host = lisaR:::lisa_explore_host(), pid = NA_integer_,
    start_time = NA_character_, dispatcher = lisaR:::lisa_explore_process_identity(),
    created_at = lisaR:::lisa_explore_now(), message = "")
  lisaR:::lisa_explore_write_job(ws, job)
  job
}

test_that("wrong-token release never moves or exposes a live owner's lock", {
  ws <- explore_h4_workspace("wrong-release")
  token <- lisaR:::lisa_explore_acquire_lock(ws)
  lock <- lisaR:::lisa_explore_lock_path(ws)
  before <- lisaR:::lisa_explore_lock_evidence(lock)

  expect_false(lisaR:::lisa_explore_release_lock(ws, "not-the-owner"))
  expect_true(dir.exists(lock))
  after <- lisaR:::lisa_explore_lock_evidence(lock)
  expect_true(before$verified)
  expect_true(after$verified)
  expect_identical(as.character(after$owner$token), token)
  expect_error(lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0),
               "LISA-EXPLORE-015")
  expect_true(lisaR:::lisa_explore_release_lock(ws, token))
})

test_that("an overlapping acquirer waits for a slow normal release", {
  skip_if_not_installed("callr")
  ws <- explore_h4_workspace("slow-release")
  token <- lisaR:::lisa_explore_acquire_lock(ws)
  release_started <- tempfile("lisa-h4-release-started-")
  acquire_started <- tempfile("lisa-h4-acquire-started-")
  acquired <- tempfile("lisa-h4-acquired-")
  child <- callr::r_bg(function(workspace, release_started, acquire_started,
                               acquired_marker) {
    suppressPackageStartupMessages(library(lisaR))
    deadline <- Sys.time() + 5
    while (!file.exists(release_started) && Sys.time() < deadline) Sys.sleep(0.01)
    if (!file.exists(release_started)) stop("release did not reach its transition")
    writeLines("attempting", acquire_started, useBytes = TRUE)
    acquired <- lisaR:::lisa_explore_acquire_lock(workspace, timeout_seconds = 3)
    writeLines("acquired", acquired_marker, useBytes = TRUE)
    released <- lisaR:::lisa_explore_release_lock(workspace, acquired)
    list(token = acquired, released = released)
  }, args = list(workspace = ws, release_started = release_started,
                 acquire_started = acquire_started,
                 acquired_marker = acquired),
     libpath = .libPaths())

  original_evidence <- lisaR:::lisa_explore_lock_evidence
  delayed <- FALSE
  acquired_during_transition <- NA
  testthat::local_mocked_bindings(
    lisa_explore_lock_evidence = function(lock, allow_marker_only = FALSE) {
      if (!delayed && grepl(".index-lock.releasing-", basename(lock),
                            fixed = TRUE)) {
        delayed <<- TRUE
        writeLines("releasing", release_started, useBytes = TRUE)
        deadline <- Sys.time() + 2
        while (!file.exists(acquire_started) && Sys.time() < deadline) {
          Sys.sleep(0.01)
        }
        if (!file.exists(acquire_started)) {
          stop("overlapping acquirer did not start deterministically")
        }
        Sys.sleep(0.15)
        acquired_during_transition <<- file.exists(acquired)
        Sys.sleep(0.05)
      }
      original_evidence(lock, allow_marker_only = allow_marker_only)
    }, .package = "lisaR")

  expect_true(lisaR:::lisa_explore_release_lock(ws, token))
  child$wait(timeout = 10000)
  expect_false(child$is_alive())
  result <- child$get_result()
  expect_true(nzchar(result$token))
  expect_true(result$released)
  expect_true(delayed)
  expect_false(acquired_during_transition)
  expect_true(file.exists(acquired))
  expect_identical(lisaR:::lisa_explore_lock_quarantines(ws), character(0))
  expect_identical(lisaR:::lisa_explore_lock_transitions(ws), character(0))
})

test_that("reclaim prechecks replacement tokens and quarantine fails closed", {
  ws <- explore_h4_workspace("reclaim-precheck")
  token <- lisaR:::lisa_explore_acquire_lock(ws)
  lock <- lisaR:::lisa_explore_lock_path(ws)
  expect_false(lisaR:::lisa_explore_reclaim_stale_lock(lock, "stale-inspection"))
  expect_identical(
    as.character(lisaR:::lisa_explore_lock_evidence(lock)$owner$token), token)
  expect_true(lisaR:::lisa_explore_release_lock(ws, token))

  quarantine <- paste0(lock, ".stale-", strrep("a", 24L))
  dir.create(quarantine, mode = "0700")
  expect_error(lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0),
               "LISA-EXPLORE-075.*No new writer was admitted")
  expect_false(dir.exists(lock))
  expect_true(dir.exists(quarantine))
})

test_that("a post-create quarantine check relinquishes the new token", {
  ws <- explore_h4_workspace("quarantine-after-create")
  calls <- 0L
  testthat::local_mocked_bindings(
    lisa_explore_assert_no_lock_quarantine = function(ws) {
      calls <<- calls + 1L
      if (calls >= 2L) stop("LISA-EXPLORE-075 simulated retained evidence",
                            call. = FALSE)
      invisible(TRUE)
    }, .package = "lisaR")
  expect_error(lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0),
               "LISA-EXPLORE-075")
  expect_false(dir.exists(lisaR:::lisa_explore_lock_path(ws)))
})

test_that("expected owner marker residue stays verifiable and recoverable", {
  ws <- explore_h4_workspace("marker-residue")
  lock <- lisaR:::lisa_explore_lock_path(ws)
  dir.create(lock, mode = "0700")
  old_token <- strrep("b", 24L)
  owner <- list(host = lisaR:::lisa_explore_host(), pid = 4000000L,
    start_time = "1", token = old_token,
    acquired_at = lisaR:::lisa_explore_now())
  lisaR:::lisa_explore_write_json(owner, file.path(lock, "owner.json"),
                                  run_root = ws$root)
  dir.create(lisaR:::lisa_explore_lock_marker(lock, old_token), mode = "0700")
  expect_true(lisaR:::lisa_explore_lock_evidence(lock)$verified)
  token <- lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0.2)
  expect_false(identical(token, old_token))
  expect_true(lisaR:::lisa_explore_release_lock(ws, token))
})

test_that("numeric-looking catalog identities remain exact text", {
  path <- tempfile(fileext = ".tsv")
  writeLines("analysis_id\tcollection\tcategory_id\trank\n001\t0002\t03\t7", path)
  catalog <- lisaR:::lisa_explore_read_catalog(path)
  expect_identical(catalog$analysis_id, "001")
  expect_identical(catalog$collection, "0002")
  expect_identical(catalog$category_id, "03")
  expect_identical(catalog$rank, "7")
})

test_that("export snapshots queue behind renders, deduplicate and cancel cleanly", {
  skip_if_not_installed("callr")
  ws <- explore_h4_workspace("export-queue")
  first_key <- explore_h4_add_artifact(ws, explore_h4_request("SYN_SIGNAL"))
  holder <- explore_h4_hold_slot(ws)
  destination <- file.path(tempdir(), "lisa-h4-queued-export")
  unlink(destination, recursive = TRUE, force = TRUE)
  withr::defer(unlink(destination, recursive = TRUE, force = TRUE),
               envir = parent.frame())

  first <- lisaR:::lisa_explore_submit_export(ws, destination, background = TRUE)
  expect_identical(first$action, "wait")
  job <- lisaR:::lisa_explore_read_job(ws, first$job_id)
  expect_identical(lisaR:::lisa_explore_job_kind(job), "export")
  expect_false(isTRUE(job$launched))
  expect_identical(as.character(unlist(job$snapshot_keys)), first_key)
  expect_false(dir.exists(destination))

  duplicate <- lisaR:::lisa_explore_submit_export(ws, destination,
                                                   background = TRUE)
  expect_identical(duplicate$action, "join")
  expect_identical(duplicate$job_id, first$job_id)
  expect_identical(sum(vapply(lisaR:::lisa_explore_all_jobs(ws), function(item)
    identical(lisaR:::lisa_explore_job_kind(item), "export"), logical(1L))), 1L)

  # A later artifact cannot enter the already-persisted snapshot.
  second_key <- explore_h4_add_artifact(ws, explore_h4_request("SYN_MATRIX"))
  expect_false(identical(first_key, second_key))
  unchanged <- lisaR:::lisa_explore_read_job(ws, first$job_id)
  expect_identical(as.character(unlist(unchanged$snapshot_keys)), first_key)

  other <- file.path(tempdir(), "lisa-h4-second-export")
  unlink(other, recursive = TRUE, force = TRUE)
  expect_error(lisaR:::lisa_explore_submit_export(ws, other, background = TRUE),
               "LISA-EXPLORE-083")

  before_index <- lisaR:::lisa_explore_read_index(ws)
  cancelled <- lisa_explore_cancel(ws, first$job_id)
  expect_true(cancelled$cancelled)
  expect_identical(cancelled$state, "failed")
  expect_identical(lisaR:::lisa_explore_read_index(ws), before_index)
  expect_false(dir.exists(destination))
  holder$state <- "failed"
  holder$launched <- FALSE
  lisaR:::lisa_explore_write_job(ws, holder)
})

test_that("the Shiny export action queues promptly and remains interactive", {
  skip_if_not_installed("callr")
  skip_if_not_installed("shiny")
  ws <- explore_h4_workspace("shiny-export")
  destination <- file.path(tempdir(), "lisa-h4-shiny-export")
  unlink(destination, recursive = TRUE, force = TRUE)
  testthat::local_mocked_bindings(
    lisa_explore_start_process = function(...)
      list(get_pid = function() Sys.getpid(), is_alive = function() TRUE),
    .package = "lisaR")

  shiny::testServer(lisaR:::lisa_explore_shiny_server(ws, "wsprefix"), {
    session$setInputs(explore_export = 1, explore_export_path = destination)
    expect_match(export_status(), "Export (queued|running)")
    job_id <- export_job_id()
    expect_true(nzchar(job_id))

    # A second click joins the same persisted export and an unrelated input can
    # still be processed while its fake background owner holds the slot.
    session$setInputs(explore_export = 2, explore_export_path = destination)
    expect_identical(export_job_id(), job_id)
    session$setInputs(explore_export_path = paste0(destination, " later"))
    session$elapse(1000)
    expect_match(export_status(), "Export (queued|running)")
    expect_identical(sum(vapply(lisaR:::lisa_explore_all_jobs(ws), function(job)
      identical(lisaR:::lisa_explore_job_kind(job), "export"), logical(1L))), 1L)
    session$setInputs(explore_export_cancel = 1)
    expect_match(export_status(), "cancelled before starting")
    lisaR:::lisa_explore_forget_process(job_id)
  })
  # The behavioral proof above covers prompt queueing, deduplication, polling,
  # unrelated input processing, and cancellation. Cold testServer startup is
  # intentionally not treated as a submit-latency measurement.
  expect_false(dir.exists(destination))
})

test_that("a refused export remains visible over an older successful job", {
  skip_if_not_installed("callr")
  skip_if_not_installed("shiny")
  ws <- explore_h4_workspace("export-feedback")
  old_destination <- file.path(tempdir(), "lisa-h4-old-success")
  completed <- list(kind = "export", job_id = "export-old-success",
    key = "export-old-success", request_id = "export-old-success",
    request = list(), source_run = ws$source_run, root = ws$root,
    destination = old_destination, state = "completed", stage = "completed",
    message = "", cancel_requested = FALSE, background = TRUE,
    launched = TRUE, host = lisaR:::lisa_explore_host(), pid = Sys.getpid(),
    start_time = lisaR:::lisa_explore_process_identity()$start_time,
    created_at = "2026-09-13T00:00:00Z")
  lisaR:::lisa_explore_write_job(ws, completed)
  refused_destination <- file.path(tempdir(), "lisa-h4-refused-existing")
  unlink(refused_destination, recursive = TRUE, force = TRUE)
  dir.create(refused_destination, mode = "0700")
  withr::defer(unlink(refused_destination, recursive = TRUE, force = TRUE),
               envir = parent.frame())
  accepted_destination <- file.path(tempdir(), "lisa-h4-feedback-next")
  unlink(accepted_destination, recursive = TRUE, force = TRUE)
  testthat::local_mocked_bindings(
    lisa_explore_start_process = function(...)
      list(get_pid = function() Sys.getpid(), is_alive = function() TRUE),
    .package = "lisaR")

  shiny::testServer(lisaR:::lisa_explore_shiny_server(ws, "wsprefix"), {
    expect_match(export_status(), "Saved expanded report to")
    session$setInputs(explore_export = 1,
                      explore_export_path = refused_destination)
    refused <- export_status()
    expect_match(refused, "LISA-EXPLORE-031.*already exists")
    session$elapse(1000)
    expect_identical(export_status(), refused)
    expect_identical(export_feedback()$kind, "error")

    # A later successful request explicitly replaces the refused-action
    # feedback and becomes the job that polling follows.
    session$setInputs(explore_export = 2,
                      explore_export_path = accepted_destination)
    expect_null(export_feedback())
    expect_match(export_status(), "Export (queued|running)")
    lisaR:::lisa_explore_forget_process(export_job_id())
  })
})

test_that("running export cancellation stays pending until terminal", {
  skip_if_not_installed("callr")
  skip_if_not_installed("shiny")
  ws <- explore_h4_workspace("cancel-feedback")
  running <- list(kind = "export", job_id = "export-cancel-running",
    key = "export-cancel-running", request_id = "export-cancel-running",
    request = list(), source_run = ws$source_run, root = ws$root,
    destination = file.path(tempdir(), "lisa-h4-cancel-running"),
    state = "running", stage = "copying STANDARD", message = "copied 128 files",
    cancel_requested = FALSE, background = TRUE, launched = TRUE,
    host = lisaR:::lisa_explore_host(), pid = Sys.getpid(),
    start_time = lisaR:::lisa_explore_process_identity()$start_time,
    created_at = "2026-09-13T00:00:01Z")
  lisaR:::lisa_explore_write_job(ws, running)

  shiny::testServer(lisaR:::lisa_explore_shiny_server(ws, "wsprefix"), {
    session$setInputs(explore_export_cancel = 1)
    expect_match(export_status(), "cancellation requested")
    expect_identical(export_feedback()$kind, "cancel_pending")
    session$elapse(1000)
    expect_match(export_status(), "cancellation requested.*next safe point")
    expect_identical(export_feedback()$kind, "cancel_pending")

    terminal <- lisaR:::lisa_explore_read_job(ws, running$job_id)
    terminal$state <- "failed"
    terminal$stage <- "failed"
    terminal$message <- "cancelled at a safe point"
    lisaR:::lisa_explore_write_job(ws, terminal)
    session$elapse(1000)
    expect_null(export_feedback())
    expect_match(export_status(), "Export failed: cancelled at a safe point")
  })
})

test_that("export progress does not invalidate the figure-control payload", {
  ws <- explore_h4_workspace("export-fingerprint")
  baseline <- lisaR:::lisa_explore_control_fingerprint(ws)
  export_job <- list(kind = "export", job_id = "export-progress", key = "export-x",
    request_id = "export-x", request = list(), source_run = ws$source_run,
    root = ws$root, state = "running", stage = "copying STANDARD",
    message = "copied 128 files", launched = TRUE, background = TRUE,
    host = lisaR:::lisa_explore_host(), pid = Sys.getpid(),
    start_time = lisaR:::lisa_explore_process_identity()$start_time,
    created_at = lisaR:::lisa_explore_now())
  lisaR:::lisa_explore_write_job(ws, export_job)
  after_export <- lisaR:::lisa_explore_control_fingerprint(ws)
  expect_identical(after_export, baseline)
  export_job$stage <- "building manifest"
  export_job$message <- "hashing staged bundle"
  lisaR:::lisa_explore_write_job(ws, export_job)
  expect_identical(lisaR:::lisa_explore_control_fingerprint(ws), baseline)

  figure_job <- export_job
  figure_job$kind <- "figure"
  figure_job$job_id <- "figure-progress"
  figure_job$key <- "fig-progress"
  figure_job$request_id <- "req-progress"
  lisaR:::lisa_explore_write_job(ws, figure_job)
  expect_false(identical(lisaR:::lisa_explore_control_fingerprint(ws), baseline))
})

test_that("orphaned export staging is reclaimed only from a proven dead owner", {
  ws <- explore_h4_workspace("export-recovery")
  destination <- file.path(tempdir(), "lisa-h4-recovered-export")
  staging <- lisaR:::lisa_explore_export_staging_path(destination)
  owner_path <- lisaR:::lisa_explore_export_owner_path(destination)
  unlink(staging, recursive = TRUE, force = TRUE)
  unlink(owner_path, force = TRUE)
  withr::defer({
    unlink(staging, recursive = TRUE, force = TRUE)
    unlink(owner_path, force = TRUE)
  }, envir = parent.frame())
  dir.create(staging, recursive = TRUE, mode = "0700")
  writeLines("partial", file.path(staging, "partial.txt"))
  lisaR:::lisa_explore_write_json(list(
    host = lisaR:::lisa_explore_host(), pid = 4000000L, start_time = "1",
    token = "old", job_id = "old-job", destination = destination,
    claimed_at = lisaR:::lisa_explore_now()), owner_path, run_root = NULL)
  token <- lisaR:::lisa_explore_claim_export_staging(
    ws, list(job_id = "new-job"), destination)
  expect_false(dir.exists(staging))
  claimed_staging <- lisaR:::lisa_explore_export_staging_path(destination, token)
  claimed_owner <- lisaR:::lisa_explore_export_owner_path(destination, token)
  expect_true(dir.exists(claimed_staging))
  expect_identical(as.character(
    lisaR:::lisa_explore_read_json(claimed_owner)$job_id), "new-job")
  expect_true(lisaR:::lisa_explore_release_export_staging(destination, token))
  unlink(claimed_staging, recursive = TRUE, force = TRUE)

  dir.create(staging, recursive = TRUE, mode = "0700")
  expect_error(lisaR:::lisa_explore_claim_export_staging(
    ws, list(job_id = "ambiguous"), destination), "LISA-EXPLORE-080")
  expect_true(dir.exists(staging))
})

test_that("KEGG action labels expose readable title and exact ID", {
  expect_identical(lisaR:::lisa_explore_action_label(
    "kegg_pathway_map", "hsa04110", "", "Cell cycle"),
    "Generate KEGG pathway map (Cell cycle (hsa04110))")
  expect_identical(lisaR:::lisa_explore_action_label(
    "contrast_kegg_map", "hsa05200", "", "Pathways in cancer"),
    "Generate contrast KEGG map (Pathways in cancer (hsa05200))")
})
