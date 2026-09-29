
explore_h5_workspace <- function(name) {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), paste0("lisa-explore-h5-", name))
  unlink(root, recursive = TRUE, force = TRUE)
  withr::defer(unlink(root, recursive = TRUE, force = TRUE),
               envir = parent.frame())
  lisa_explore_open(run, root)
}

explore_h5_job <- function(ws, job_id, kind = "export", state = "running") {
  identity <- lisaR:::lisa_explore_process_identity()
  list(kind = kind, job_id = job_id, key = paste0("key-", job_id),
       request_id = paste0("request-", job_id), request = list(),
       source_run = ws$source_run, root = ws$root,
       destination = file.path(tempdir(), paste0("destination-", job_id)),
       state = state, stage = state, message = "",
       cancel_requested = FALSE, background = TRUE, launched = TRUE,
       host = identity$host, pid = identity$pid,
       start_time = identity$start_time,
       created_at = lisaR:::lisa_explore_now())
}

test_that("locked progress waits and preserves a concurrent cancellation", {
  skip_if_not_installed("callr")
  ws <- explore_h5_workspace("progress-cancel")
  job <- explore_h5_job(ws, "h5-progress-cancel")
  job$stage <- "copying STANDARD"
  lisaR:::lisa_explore_write_job(ws, job)
  started_marker <- tempfile("lisa-h5-progress-started-")

  token <- lisaR:::lisa_explore_acquire_lock(ws)
  child <- callr::r_bg(function(workspace, job_id, marker) {
    suppressPackageStartupMessages(library(lisaR))
    writeLines("started", marker, useBytes = TRUE)
    started <- Sys.time()
    outcome <- tryCatch({
      lisaR:::lisa_explore_job_progress(
        workspace, job_id, "building manifest", "hashing staged bundle")
      "updated"
    }, error = conditionMessage)
    list(outcome = outcome,
         elapsed = as.numeric(difftime(Sys.time(), started, units = "secs")))
  }, args = list(workspace = ws, job_id = job$job_id,
                 marker = started_marker), libpath = .libPaths())
  released <- FALSE
  withr::defer({
    if (!released) lisaR:::lisa_explore_release_lock(ws, token)
    if (child$is_alive()) child$wait(timeout = 35000)
  }, envir = parent.frame())

  deadline <- Sys.time() + 5
  while (!file.exists(started_marker) && Sys.time() < deadline) Sys.sleep(0.01)
  expect_true(file.exists(started_marker))
  Sys.sleep(0.2)
  cancelled <- lisaR:::lisa_explore_read_job(ws, job$job_id)
  cancelled$cancel_requested <- TRUE
  lisaR:::lisa_explore_write_job(ws, cancelled)
  released <- lisaR:::lisa_explore_release_lock(ws, token)
  expect_true(released)

  child$wait(timeout = 10000)
  result <- child$get_result()
  expect_match(result$outcome, "cancelled at the safe point")
  expect_gte(result$elapsed, 0.15)
  persisted <- lisaR:::lisa_explore_read_job(ws, job$job_id)
  expect_true(isTRUE(persisted$cancel_requested))
  expect_identical(as.character(persisted$stage), "copying STANDARD")
})

test_that("foreground export refusal leaves no undispatchable job", {
  ws <- explore_h5_workspace("foreground-refusal")
  holder <- explore_h5_job(ws, "h5-slot-holder", kind = "figure",
                           state = "queued")
  holder$dispatcher <- lisaR:::lisa_explore_process_identity()
  lisaR:::lisa_explore_write_job(ws, holder)
  destination <- file.path(tempdir(), "lisa-h5-foreground-export")
  unlink(destination, recursive = TRUE, force = TRUE)

  expect_error(
    lisaR:::lisa_explore_submit_export(ws, destination, background = FALSE),
    "LISA-EXPLORE-084.*No export was queued")
  jobs <- lisaR:::lisa_explore_all_jobs(ws)
  expect_identical(length(jobs), 1L)
  expect_identical(as.character(jobs[[1L]]$job_id), holder$job_id)
  expect_false(any(vapply(jobs, function(item)
    identical(lisaR:::lisa_explore_job_kind(item), "export"), logical(1L))))
  expect_null(lisaR:::lisa_explore_next_dispatchable(ws))
  expect_false(dir.exists(destination))
})

test_that("export reports final verification before atomic publication", {
  source <- file.path(tempdir(), "lisa-h5-minimal-standard")
  workspace <- file.path(tempdir(), "lisa-h5-publication-workspace")
  destination <- file.path(tempdir(), "lisa-h5-publication-stage")
  unlink(c(source, workspace, destination), recursive = TRUE, force = TRUE)
  dir.create(source, mode = "0700")
  writeLines("<html><body>minimal STANDARD</body></html>",
             file.path(source, "report_index.html"), useBytes = TRUE)
  writeLines(paste("path", "sha256", sep = "\t"),
             file.path(source, "run_manifest.tsv"), useBytes = TRUE)
  withr::defer(unlink(c(source, workspace, destination), recursive = TRUE,
                      force = TRUE), envir = parent.frame())
  ws <- lisa_explore_open(source, workspace)
  testthat::local_mocked_bindings(
    lisa_explore_artifacts = function(ws) data.frame(), .package = "lisaR")
  events <- list()
  lisa_explore_export(ws, destination, snapshot_artifacts = data.frame(),
    progress = function(stage, message) {
      events[[length(events) + 1L]] <<- c(stage = stage, message = message)
    })
  stages <- vapply(events, `[[`, character(1L), "stage")
  messages <- vapply(events, `[[`, character(1L), "message")
  expect_true("final verification" %in% stages)
  expect_match(messages[match("final verification", stages)],
               "before atomic publication")
  expect_true(file.exists(file.path(destination,
                                    "explore_export_manifest.tsv")))
})

test_that("two workspaces isolate staging and the winning destination", {
  first_ws <- explore_h5_workspace("destination-first")
  second_ws <- explore_h5_workspace("destination-second")
  destination <- file.path(tempdir(), "lisa-h5-shared-destination")
  unlink(destination, recursive = TRUE, force = TRUE)
  withr::defer({
    unlink(destination, recursive = TRUE, force = TRUE)
    candidates <- list.files(dirname(destination), all.files = TRUE,
      full.names = TRUE, pattern = paste0("^[.]", basename(destination),
                                         "[.]export-(staging|owner)"))
    unlink(candidates, recursive = TRUE, force = TRUE)
  }, envir = parent.frame())

  first_token <- lisaR:::lisa_explore_claim_export_staging(
    first_ws, list(job_id = "first-export"), destination)
  second_token <- lisaR:::lisa_explore_claim_export_staging(
    second_ws, list(job_id = "second-export"), destination)
  expect_false(identical(first_token, second_token))
  first_staging <- lisaR:::lisa_explore_export_staging_path(
    destination, first_token)
  second_staging <- lisaR:::lisa_explore_export_staging_path(
    destination, second_token)
  expect_true(dir.exists(first_staging))
  expect_true(dir.exists(second_staging))
  writeLines("winner", file.path(first_staging, "identity.txt"), useBytes = TRUE)
  writeLines("loser", file.path(second_staging, "identity.txt"), useBytes = TRUE)

  expect_true(file.rename(first_staging, destination))
  expect_true(lisaR:::lisa_explore_release_export_staging(
    destination, first_token))
  expect_error(
    lisaR:::lisa_promote_managed_directory(second_staging, destination),
    "Refusing to overwrite an existing managed directory")
  expect_identical(readLines(file.path(destination, "identity.txt")), "winner")
  expect_true(file.exists(file.path(second_staging, "identity.txt")))
  expect_true(lisaR:::lisa_explore_release_export_staging(
    destination, second_token))
})

test_that("refusal clears when the followed export changes to terminal failure", {
  skip_if_not_installed("shiny")
  ws <- explore_h5_workspace("feedback-transition")
  running <- explore_h5_job(ws, "h5-followed-export")
  running$stage <- "copying STANDARD"
  lisaR:::lisa_explore_write_job(ws, running)
  other <- file.path(tempdir(), "lisa-h5-refused-while-running")
  unlink(other, recursive = TRUE, force = TRUE)

  shiny::testServer(lisaR:::lisa_explore_shiny_server(ws, "wsprefix"), {
    session$setInputs(explore_export = 1, explore_export_path = other)
    expect_match(export_status(), "LISA-EXPLORE-083")
    expect_identical(export_feedback()$followed_state, "running")

    failed <- lisaR:::lisa_explore_read_job(ws, running$job_id)
    failed$state <- "failed"
    failed$stage <- "failed"
    failed$message <- "source changed during export"
    failed$launched <- FALSE
    lisaR:::lisa_explore_write_job(ws, failed)
    session$elapse(1000)
    expect_null(export_feedback())
    expect_match(export_status(),
                 "Export failed: source changed during export")
  })
})

test_that("reopen derives persisted cancellation status without pinned feedback", {
  skip_if_not_installed("shiny")
  ws <- explore_h5_workspace("reopen-cancellation")
  running <- explore_h5_job(ws, "h5-reopen-cancel")
  running$stage <- "copying STANDARD"
  running$cancel_requested <- TRUE
  lisaR:::lisa_explore_write_job(ws, running)

  shiny::testServer(lisaR:::lisa_explore_shiny_server(ws, "wsprefix"), {
    expect_null(export_feedback())
    expect_match(export_status(),
                 "cancellation requested; stopping at the next safe point")
  })
})

test_that("post-create transition makes an acquirer relinquish and retry", {
  ws <- explore_h5_workspace("post-create-yield")
  calls <- 0L
  testthat::local_mocked_bindings(
    lisa_explore_recover_orphaned_lock_transitions = function(...) invisible(0L),
    lisa_explore_lock_transitions = function(ws) {
      calls <<- calls + 1L
      if (calls == 2L) "synthetic-transition" else character(0)
    }, .package = "lisaR")
  token <- lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 3)
  expect_gte(calls, 4L)
  expect_true(dir.exists(lisaR:::lisa_explore_lock_path(ws)))
  expect_true(lisaR:::lisa_explore_release_lock(ws, token))
  expect_identical(lisaR:::lisa_explore_lock_transitions(ws), character(0))
})

test_that("orphaned transitions recover only from verified dead ownership", {
  ws <- explore_h5_workspace("orphan-transitions")
  lock <- lisaR:::lisa_explore_lock_path(ws)
  make_transition <- function(kind, suffix, owner = NULL) {
    path <- paste0(lock, ".", kind, "-", suffix)
    dir.create(path, mode = "0700")
    if (!is.null(owner)) {
      lisaR:::lisa_explore_write_json(owner, file.path(path, "owner.json"),
                                      run_root = ws$root)
    }
    path
  }
  owner <- function(token, host = lisaR:::lisa_explore_host(),
                    pid = 4000000L, start_time = "1") {
    list(host = host, pid = pid, start_time = start_time, token = token,
         acquired_at = lisaR:::lisa_explore_now())
  }

  dead <- make_transition("releasing", strrep("a", 24L),
                          owner(strrep("b", 24L)))
  token <- lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0.5)
  expect_false(dir.exists(dead))
  expect_true(lisaR:::lisa_explore_release_lock(ws, token))

  live_identity <- lisaR:::lisa_explore_process_identity()
  live <- make_transition("reclaiming", strrep("c", 24L),
    c(live_identity, list(token = strrep("d", 24L),
                          acquired_at = lisaR:::lisa_explore_now())))
  expect_error(lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0),
               "LISA-EXPLORE-015.*verified owner is live")
  expect_true(dir.exists(live))
  unlink(live, recursive = TRUE, force = TRUE)

  foreign <- make_transition("releasing", strrep("e", 24L),
    owner(strrep("f", 24L), host = "foreign.example", pid = 1L))
  expect_error(lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0),
               "LISA-EXPLORE-015.*another host")
  expect_true(dir.exists(foreign))
  unlink(foreign, recursive = TRUE, force = TRUE)

  ambiguous <- make_transition("reclaiming", strrep("0", 24L))
  expect_error(lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0),
               "LISA-EXPLORE-075.*could not be verified.*retained")
  expect_true(dir.exists(ambiguous))
  unlink(ambiguous, recursive = TRUE, force = TRUE)

  marker_only <- make_transition("releasing", strrep("1", 24L))
  marker_token <- strrep("2", 24L)
  dir.create(lisaR:::lisa_explore_lock_marker(marker_only, marker_token),
             mode = "0700")
  expect_error(lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0),
               "LISA-EXPLORE-075.*complete process identity is absent")
  expect_true(dir.exists(marker_only))
})
