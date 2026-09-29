
explore_g1_request <- function(category_id, collection = "GOBP-C2",
                               analysis_id = "response_a") {
  lisa_figure_request(unit_type = "single_de", analysis_id = analysis_id,
                      collection = collection, category_id = category_id,
                      product = "volcano")
}

# A fresh workspace over the saved fixture run. Opening is a read; nothing is
# rendered here.
explore_g1_workspace <- function(name) {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), paste0("lisa-explore-g1-", name))
  unlink(root, recursive = TRUE, force = TRUE)
  withr::defer(unlink(root, recursive = TRUE, force = TRUE),
               envir = parent.frame())
  lisa_explore_open(run, root)
}

# A job record that occupies the single worker slot without any worker existing:
# queued, dispatched, and owned by a process that really is alive (this one).
explore_g1_hold_slot <- function(ws) {
  job <- list(job_id = "job-g1-holder", key = "fig-g1-holder",
              request_id = "req-g1-holder", request = list(),
              source_run = ws$source_run, root = ws$root,
              state = "queued", stage = "queued", cancel_requested = FALSE,
              background = TRUE, launched = TRUE,
              host = lisaR:::lisa_explore_host(), pid = NA_integer_,
              start_time = NA_character_,
              dispatcher = lisaR:::lisa_explore_process_identity(),
              created_at = lisaR:::lisa_explore_now(), message = "")
  lisaR:::lisa_explore_write_job(ws, job)
  job
}

explore_g1_release_slot <- function(ws, job) {
  job$state <- "failed"
  job$launched <- FALSE
  lisaR:::lisa_explore_write_job(ws, job)
  invisible(NULL)
}

# A process handle that behaves like a live callr child without being one.
explore_g1_fake_process <- function() {
  list(get_pid = function() Sys.getpid(), is_alive = function() TRUE)
}


test_that("a refused foreground submit leaves no record and does not poison the queue", {
  ws <- explore_g1_workspace("refuse")
  request <- explore_g1_request("SYN_MATRIX")
  plan <- lisa_explore_plan(ws, request)
  expect_true(plan$applicable)

  holder <- explore_g1_hold_slot(ws)
  expect_identical(as.character(lisaR:::lisa_explore_slot_holder(ws)$job_id),
                   "job-g1-holder")

  traced <- explore_with_generator_trace(
    expect_error(lisa_explore_submit(ws, request, background = FALSE),
                 "LISA-EXPLORE-036"))
  # A refusal is a refusal: nothing was rendered to produce it.
  expect_identical(length(traced$post), 0L)

  # The defect was that the refusal left a `queued` job and index entry which the
  # pump would never dispatch (it only dispatches background jobs) and which
  # every later submit -- background or not -- would join instead of starting.
  # Nothing at all may survive the refusal.
  expect_null(lisaR:::lisa_explore_index_find(lisaR:::lisa_explore_read_index(ws),
                                              plan$key))
  jobs <- lisaR:::lisa_explore_all_jobs(ws)
  expect_identical(length(jobs), 1L)
  expect_identical(as.character(jobs[[1L]]$job_id), "job-g1-holder")
  expect_identical(lisa_explore_plan(ws, request)$state, "ungenerated")

  # There is no undispatchable active job left behind, by the engine's own test.
  expect_null(lisaR:::lisa_explore_next_dispatchable(ws))

  # Once the slot is free the very same figure submits successfully and is
  # actually dispatched -- it is not silently joined to a ghost.
  explore_g1_release_slot(ws, holder)
  expect_null(lisaR:::lisa_explore_slot_holder(ws))
  testthat::local_mocked_bindings(
    lisa_explore_start_process = function(...) explore_g1_fake_process(),
    .package = "lisaR")
  withr::defer(lisaR:::lisa_explore_forget_process(result$job_id))
  result <- lisa_explore_submit(ws, request, background = TRUE)
  expect_false(result$reused)
  expect_identical(result$key, plan$key)
  job <- lisaR:::lisa_explore_read_job(ws, result$job_id)
  expect_true(isTRUE(job$launched))
  expect_identical(as.character(job$state), "queued")
})


test_that("a launch failure fails the job, frees the slot, and lets the next job run", {
  ws <- explore_g1_workspace("launch")
  first <- explore_g1_request("SYN_MATRIX")
  second <- explore_g1_request("SYN_STRESS")

  testthat::local_mocked_bindings(
    lisa_explore_start_process = function(...)
      stop("callr refused to start a child", call. = FALSE),
    .package = "lisaR")
  traced <- explore_with_generator_trace(
    lisa_explore_submit(ws, first, background = TRUE))
  # The failure is injected at launch, so nothing was ever rendered.
  expect_identical(length(traced$post), 0L)
  failed <- traced$value

  expect_identical(failed$state, "failed")
  job <- lisaR:::lisa_explore_read_job(ws, failed$job_id)
  expect_identical(as.character(job$state), "failed")
  expect_false(isTRUE(job$launched))
  expect_match(as.character(job$message), "could not be started")

  # The index agrees, so the UI cannot keep showing the figure as queued.
  entry <- lisaR:::lisa_explore_index_find(lisaR:::lisa_explore_read_index(ws),
                                           failed$key)
  expect_identical(as.character(entry$state), "failed")

  # The slot is free again without restarting the session.
  expect_null(lisaR:::lisa_explore_slot_holder(ws))

  testthat::local_mocked_bindings(
    lisa_explore_start_process = function(...) explore_g1_fake_process(),
    .package = "lisaR")
  next_result <- lisa_explore_submit(ws, second, background = TRUE)
  withr::defer(lisaR:::lisa_explore_forget_process(next_result$job_id))
  expect_true(isTRUE(lisaR:::lisa_explore_read_job(ws, next_result$job_id)$launched))
  expect_identical(as.character(lisaR:::lisa_explore_slot_holder(ws)$job_id),
                   next_result$job_id)
})


explore_g1_make_staging <- function(ws, key, suffix = "deadbeef") {
  path <- file.path(ws$root, "extensions",
                    paste0(lisaR:::lisa_explore_staging_prefix(key), suffix))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  writeLines("half-written", file.path(path, "leftover.txt"))
  path
}

test_that("orphaned deterministic staging is recovered so a retry can run", {
  ws <- explore_g1_workspace("staging")
  key <- "fig-g1staging0000000000"
  other_key <- "fig-g1other00000000000"
  job <- list(job_id = "job-g1-retry")

  stale <- explore_g1_make_staging(ws, key)
  untouched <- explore_g1_make_staging(ws, other_key)
  artifact <- file.path(ws$root, "extensions", key)
  dir.create(artifact, recursive = TRUE, showWarnings = FALSE)
  writeLines("previous", file.path(artifact, "kept.txt"))

  expect_identical(lisaR:::lisa_explore_recover_staging(ws, job, key), 1L)

  # The tombstone that made LISA-EXTENSION-022 permanent is gone...
  expect_false(dir.exists(stale))
  # ...another request's staging is not, and neither is a valid previous artifact.
  expect_true(dir.exists(untouched))
  expect_true(file.exists(file.path(artifact, "kept.txt")))

  # Ownership is now recorded, so a second worker can tell whether it may clean.
  owner <- lisaR:::lisa_explore_read_json(
    lisaR:::lisa_explore_staging_owner_path(ws, key))
  expect_identical(as.character(owner$job_id), "job-g1-retry")
  lisaR:::lisa_explore_clear_staging_owner(ws, key)
  expect_false(file.exists(lisaR:::lisa_explore_staging_owner_path(ws, key)))
})

test_that("a live owner's staging is refused, a dead owner's is reclaimed", {
  ws <- explore_g1_workspace("staging-owner")
  key <- "fig-g1owner00000000000"
  live <- explore_g1_make_staging(ws, key)

  # A marker naming a different, genuinely live process.
  lisaR:::lisa_explore_write_json(
    c(lisaR:::lisa_explore_process_identity(),
      list(job_id = "job-other", key = key)),
    lisaR:::lisa_explore_staging_owner_path(ws, key), run_root = ws$root)
  expect_error(
    lisaR:::lisa_explore_recover_staging(ws, list(job_id = "job-mine"), key),
    "LISA-EXPLORE-039")
  # Refusing means refusing: the other worker's staging is still there.
  expect_true(dir.exists(live))

  # Same PID, different process start time: the recorded owner is gone and its
  # number has been reused. That identity is dead and its debris is reclaimable.
  identity <- lisaR:::lisa_explore_process_identity()
  identity$start_time <- "0"
  lisaR:::lisa_explore_write_json(
    c(identity, list(job_id = "job-other", key = key)),
    lisaR:::lisa_explore_staging_owner_path(ws, key), run_root = ws$root)
  expect_identical(
    lisaR:::lisa_explore_recover_staging(ws, list(job_id = "job-mine"), key), 1L)
  expect_false(dir.exists(live))
})


explore_g1_portability_tree <- function(name, files) {
  root <- file.path(tempdir(), paste0("lisa-explore-g1-port-", name))
  unlink(root, recursive = TRUE, force = TRUE)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  for (relative in names(files)) {
    target <- file.path(root, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    writeLines(files[[relative]], target)
  }
  root
}

test_that("portability rejects every absolute local reference, not only /home", {
  offenders <- list(
    unix_home   = '<a href="/home/fixture-user/a.png">x</a>',
    unix_tmp    = '<img src="/tmp/a.png">',
    unix_srv    = '<a href="/srv/reports/a.html">x</a>',
    unix_users  = '<img src="/Users/fixture-user/a.png">',
    unix_mnt    = '<a href="/mnt/share/a.html">x</a>',
    root_rel    = '<link rel="stylesheet" href="/assets/report.css">',
    file_url    = '<a href="file:///data/a.html">x</a>',
    windows_drv = '<img src="C:/reports/a.png">',
    windows_bs  = '<img src="C:\\reports\\a.png">',
    unc         = '<a href="\\\\server\\share\\a.html">x</a>',
    protocol_rel = '<script src="//cdn.example.org/a.js"></script>')
  for (name in names(offenders)) {
    root <- explore_g1_portability_tree(name, stats::setNames(
      list(paste0("<html><body>", offenders[[name]], "</body></html>")),
      "page.html"))
    report <- lisaR:::lisa_explore_portability_report(root)
    expect_identical(report$absolute_reference_files, "page.html", info = name)
    unlink(root, recursive = TRUE, force = TRUE)
  }
})

test_that("portability accepts the relative links, query routes and anchors the report uses", {
  root <- explore_g1_portability_tree("clean", list(
    "page.html" = paste0(
      '<html><head><link rel="stylesheet" href="../../report_assets/lisa-explore.css">',
      '<script defer src="./assets/category-evidence.js"></script></head><body>',
      # The report's real category routing: a relative path plus a query route
      # plus a fragment. None of this breaks relocation and none may be flagged.
      '<a href="../../../evidence/response_a/GOBP-C2/index.html?category=SYN_SIGNAL',
      '&amp;analysis_id=response_a">open</a>',
      '<a href="#category-SYN_SIGNAL">jump</a>',
      '<a href="mailto:someone@example.org">mail</a>',
      '<img src="figures/volcano.png"></body></html>'),
    "report_assets/style.css" = paste0(
      "@import url('theme.css');\n",
      ".x{background:url(\"../figures/a.png\")}\n",
      ".y{background:url(data:image/png;base64,AAAA)}")))
  report <- lisaR:::lisa_explore_portability_report(root)
  expect_identical(report$absolute_reference_files, character(0))
  expect_identical(report$symlinks, character(0))
  unlink(root, recursive = TRUE, force = TRUE)
})

test_that("an absolute reference inside CSS is caught, and external URLs are recorded not refused", {
  root <- explore_g1_portability_tree("css", list(
    "a.css" = ".x{background:url(/assets/a.png)}",
    "b.html" = '<html><body><a href="https://example.org/x">out</a></body></html>'))
  report <- lisaR:::lisa_explore_portability_report(root)
  expect_identical(report$absolute_reference_files, "a.css")
  # An external URL is a network dependency, not a relocation failure: it is
  # reported so the receipt is honest, and it does not make the bundle fail.
  expect_identical(report$external_reference_files, "b.html")
  unlink(root, recursive = TRUE, force = TRUE)
})

test_that("the symlink check is real and runs before the guarded scan would refuse", {
  skip_on_os("windows")
  root <- explore_g1_portability_tree("symlink", list(
    "page.html" = "<html><body><a href=\"a.png\">x</a></body></html>"))
  expect_true(file.symlink(file.path(root, "page.html"),
                           file.path(root, "link.html")))
  report <- lisaR:::lisa_explore_portability_report(root)
  expect_identical(report$symlinks, "link.html")
  # The previous ordering called the guarded scan first, which *throws* on the
  # first symlink, so this branch could never be reached and the claim that
  # symlinks were enforced held only because the guarded copy cannot make one.
  expect_error(lisaR:::lisa_explore_relative_paths(root))
  unlink(root, recursive = TRUE, force = TRUE)
})

test_that("the ZIP bundle is produced when a zip program exists on PATH", {
  skip_if(!nzchar(Sys.which("zip")), "no zip program on PATH")
  withr::local_envvar(R_ZIPCMD = "")

  root <- file.path(tempdir(), "lisa-explore-g1-zip")
  unlink(root, recursive = TRUE, force = TRUE)
  unlink(paste0(root, ".zip"), force = TRUE)
  withr::defer(unlink(c(root, paste0(root, ".zip")), recursive = TRUE, force = TRUE))
  dir.create(file.path(root, "report_pages"), recursive = TRUE, showWarnings = FALSE)
  writeLines("<html><body>x</body></html>", file.path(root, "report_index.html"))

  archive <- lisaR:::lisa_explore_zip_bundle(root)
  expect_true(file.exists(archive))
  expect_true(file.info(archive)$size > 0)
  expect_true("report_index.html" %in%
                basename(utils::unzip(archive, list = TRUE)$Name))
  # An existing archive is never silently overwritten.
  expect_error(lisaR:::lisa_explore_zip_bundle(root), "LISA-EXPLORE-034")
})


test_that("the default workspace is user-owned and outside the run's results tree", {
  run <- explore_fixture_run()
  home <- file.path(tempdir(), "lisa-explore-g1-userdata")
  unlink(home, recursive = TRUE, force = TRUE)
  withr::local_envvar(R_USER_DATA_DIR = home)
  withr::defer(unlink(home, recursive = TRUE, force = TRUE))

  root <- lisaR:::lisa_explore_default_root(run)
  # Not beside the run, and not anywhere inside the tree that holds it.
  expect_false(lisaR:::lisa_path_within(root, dirname(run)))
  expect_false(lisaR:::lisa_path_within(root, run))
  expect_false(is.null(lisaR:::lisa_path_relative_within(root, home)))
  # Public package code must not carry a project-private path. Checked against
  # the deparsed function bodies rather than the file text, because the comments
  # legitimately *name* the accepted run whose tree the old default intruded on.
  code <- paste(c(deparse(lisaR:::lisa_explore_default_root),
                  deparse(lisaR:::lisa_explore_user_data_dir),
                  deparse(lisaR:::lisa_explore_assert_default_root_usable),
                  deparse(lisaR:::lisa_explore_open)), collapse = "\n")
  expect_false(grepl("/home/|/Users/|[A-Za-z]:/", code))

  # Opening with the default writes nothing into the run or its parent tree.
  before_run <- explore_tree_digest(run)
  before_parent <- sort(list.files(dirname(run), all.files = TRUE, no.. = TRUE))
  ws <- lisa_explore_open(run)
  expect_false(is.null(lisaR:::lisa_path_relative_within(ws$root, home)))
  expect_identical(explore_tree_digest(run), before_run)
  expect_identical(sort(list.files(dirname(run), all.files = TRUE, no.. = TRUE)),
                   before_parent)
})

test_that("an unusable default asks for an explicit extension_root instead of guessing", {
  run <- explore_fixture_run()
  blocker <- tempfile("lisa-explore-g1-blocked")
  writeLines("not a directory", blocker)
  withr::defer(unlink(blocker, force = TRUE))
  # R_user_dir() resolves under this path, whose parent is a regular file, so no
  # safe default can be created.
  withr::local_envvar(R_USER_DATA_DIR = file.path(blocker, "data"))
  expect_error(lisa_explore_open(run), "LISA-EXPLORE-016")
  expect_error(lisa_explore_open(run), "extension_root")
})


explore_g1_presentation <- function(categories, route) {
  do.call(rbind, lapply(categories, function(category) data.frame(
    key = paste0("fig-", tolower(category)), category_id = category,
    product = "volcano", attach_route = route,
    attach_anchor = paste0("category-", category),
    attach_query = paste0("category=", category),
    png = paste0("explore_extensions/fig-", tolower(category), "/a.png"),
    pdf = "", source_data = "", recipe = "",
    stringsAsFactors = FALSE)))
}

test_that("a real category sheet gets the block inside the open category's section", {
  run <- explore_fixture_run()
  route <- explore_attach_route()
  page <- file.path(run, route)
  expect_true(file.exists(page))
  html <- paste(readLines(page, warn = FALSE), collapse = "\n")

  presentation <- explore_g1_presentation(c("SYN_SIGNAL", "SYN_MATRIX"), route)
  attached <- lisaR:::lisa_explore_presentation_attach_html(html, presentation, route)

  # Placement, on the page the reader actually opens: inside `main`, after the
  # category's own heading block, and before the analytical panels -- not at the
  # foot of the document, which is where the audited build put it.
  at <- regexpr('data-lisa-explore-attachment="true"', attached, fixed = TRUE)[[1L]]
  expect_gt(at, regexpr('id="category-info"', attached, fixed = TRUE)[[1L]])
  first_panel <- regexpr(lisaR:::lisa_explore_presentation_panel_pattern(),
                         attached, perl = TRUE)[[1L]]
  expect_lt(at, first_panel)
  expect_lt(at, regexpr("(?i)</main\\s*>", attached, perl = TRUE)[[1L]])

  # Category scoping: each block declares its own category and starts hidden, so
  # the sheet for SYN_SIGNAL cannot present SYN_MATRIX's figure as its own.
  expect_match(attached, 'data-lisa-explore-category="SYN_SIGNAL"', fixed = TRUE)
  expect_match(attached, 'data-lisa-explore-category="SYN_MATRIX"', fixed = TRUE)
  blocks <- regmatches(attached, gregexpr(
    '<section class="lisa-explore-attachment"[^>]*>', attached, perl = TRUE))[[1L]]
  expect_identical(length(blocks), 2L)
  expect_true(all(grepl(" hidden>", blocks, fixed = TRUE)))

  # The sheet is left otherwise intact and still has exactly one main element.
  expect_identical(length(gregexpr("(?i)</main\\s*>", attached,
                                   perl = TRUE)[[1L]]), 1L)
  expect_match(attached, "leading-edge gene connections", fixed = TRUE)
  # And re-attaching is still a no-op, so repeated exports cannot stack blocks.
  expect_identical(
    lisaR:::lisa_explore_presentation_attach_html(attached, presentation, route),
    attached)
})

test_that("the collection navigator fallback shows every category and hides none", {
  run <- explore_fixture_run()
  route <- "report_pages/category_navigation/response_a/GOBP-C2/index.html"
  html <- paste(readLines(file.path(run, route), warn = FALSE), collapse = "\n")
  expect_false(lisaR:::lisa_explore_presentation_scoped(html))

  attached <- lisaR:::lisa_explore_presentation_attach_html(
    html, explore_g1_presentation(c("SYN_SIGNAL", "SYN_MATRIX"), route), route)
  blocks <- regmatches(attached, gregexpr(
    '<section class="lisa-explore-attachment"[^>]*>', attached, perl = TRUE))[[1L]]
  expect_identical(length(blocks), 2L)
  # This page legitimately shows the whole collection at once, so nothing that
  # belongs on it may be hidden by the category-scoping rule.
  expect_false(any(grepl(" hidden>", blocks, fixed = TRUE)))
})

test_that("the attachment route is the category sheet, with the navigator as fallback", {
  run <- explore_fixture_run()
  request <- explore_g1_request("SYN_SIGNAL")
  resolved <- lisaR:::lisa_explore_attachment(request, run)
  expect_identical(resolved$route, explore_attach_route())
  expect_identical(resolved$surface, "evidence")
  expect_true(file.exists(file.path(run, resolved$route)))

  # A run without that sheet is not given a route that is not in its bundle.
  empty <- file.path(tempdir(), "lisa-explore-g1-norun")
  unlink(empty, recursive = TRUE, force = TRUE)
  dir.create(empty, recursive = TRUE, showWarnings = FALSE)
  withr::defer(unlink(empty, recursive = TRUE, force = TRUE))
  fallback <- lisaR:::lisa_explore_attachment(request, empty)
  expect_identical(fallback$surface, "category_navigation")
  expect_match(fallback$route, "category_navigation", fixed = TRUE)

  # A contrast product resolves to the contrast sheet, not the single-DE one.
  contrast <- lisaR:::lisa_explore_attachment(lisa_figure_request(
    unit_type = "contrast", contrast_id = "c1", collection = "H",
    category_id = "C", product = "contrast_profile"))
  expect_identical(contrast$route, "report_pages/contrast_evidence/c1/H/index.html")
})


test_that("a refused submission is reported in the category context and can be retried", {
  skip_if_not_installed("shiny")
  ws <- explore_g1_workspace("shiny")

  # A legacy catalogue product the exact engine does not generate
  # (LISA-EXPLORE-021). It is a real refusal
  # reachable from a real control payload, and before the repair it propagated
  # out of the observeEvent and closed the Shiny session outright.
  refused <- list(unit_type = "single_de", analysis_id = "response_a",
                  contrast_id = "", collection = "GOBP-C2",
                  category_id = "SYN_MATRIX", product = "kegg")
  refused_id <- lisa_request_id(lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a",
    collection = "GOBP-C2", category_id = "SYN_MATRIX", product = "kegg"))
  accepted <- list(unit_type = "single_de", analysis_id = "response_a",
                   contrast_id = "", collection = "GOBP-C2",
                   category_id = "SYN_MATRIX", product = "volcano")
  accepted_id <- lisa_request_id(explore_g1_request("SYN_MATRIX"))

  testthat::local_mocked_bindings(
    lisa_explore_start_process = function(...) explore_g1_fake_process(),
    .package = "lisaR")

  shiny::testServer(lisaR:::lisa_explore_shiny_server(ws, "wsprefix"), {
    session$setInputs(explore_generate = c(refused, list(nonce = 1)))
    notes <- feedback()
    expect_true(refused_id %in% names(notes))
    expect_false(isTRUE(notes[[refused_id]]$ok))
    expect_match(notes[[refused_id]]$message, "LISA-EXPLORE-021")

    # A malformed payload is refused the same way rather than ending the session.
    session$setInputs(explore_generate = list(nonce = 2))
    expect_true(length(feedback()) >= 1L)

    # The session is still alive and still usable: a different action is
    # accepted and answered.
    session$setInputs(explore_export = 1, explore_export_path = "")
    expect_match(export_status(), "Choose a new destination")

    # And a valid figure can still be submitted afterwards: valid work is
    # preserved and the refusal is not sticky.
    session$setInputs(explore_generate = c(accepted, list(nonce = 3)))
    expect_false(accepted_id %in% names(feedback()))
    expect_true(refused_id %in% names(feedback()))
  })

  jobs <- lisaR:::lisa_explore_all_jobs(ws)
  expect_identical(length(jobs), 1L)
  expect_true(isTRUE(jobs[[1L]]$launched))
  expect_identical(as.character(jobs[[1L]]$state), "queued")
  lisaR:::lisa_explore_forget_process(as.character(jobs[[1L]]$job_id))
})

test_that("a failing engine call does not take the polling observer down with it", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("callr")
  ws <- explore_g1_workspace("shiny-pump")
  testthat::local_mocked_bindings(
    lisa_explore_pump = function(...) stop("LISA-EXPLORE-015 lock held", call. = FALSE),
    .package = "lisaR")
  shiny::testServer(lisaR:::lisa_explore_shiny_server(ws, "wsprefix"), {
    session$elapse(1000)
    # One normal short hold is suppressed; persistence makes it actionable.
    expect_identical(engine_status(), "")
    session$elapse(1000)
    expect_match(engine_status(), "LISA-EXPLORE-015")
  })
})

test_that("one healthy transient lock hold never flashes as a toolbar error", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("callr")
  ws <- explore_g1_workspace("shiny-transient-pump")
  attempts <- 0L
  testthat::local_mocked_bindings(
    lisa_explore_pump = function(...) {
      attempts <<- attempts + 1L
      if (attempts == 1L) stop("LISA-EXPLORE-015 short healthy hold",
                               call. = FALSE)
      invisible(NULL)
    }, .package = "lisaR")
  shiny::testServer(lisaR:::lisa_explore_shiny_server(ws, "wsprefix"), {
    session$elapse(1000)
    expect_identical(engine_status(), "")
    session$elapse(1000)
    expect_identical(engine_status(), "")
  })
})

test_that("an unverifiable lock does not block the Shiny polling event loop", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("callr")
  ws <- explore_g1_workspace("shiny-incomplete-lock")
  lock <- lisaR:::lisa_explore_lock_path(ws)
  dir.create(lock, mode = "0700")
  writeLines("{", file.path(lock, ".lisa-write-malformed.json"), useBytes = TRUE)

  started <- proc.time()[["elapsed"]]
  shiny::testServer(lisaR:::lisa_explore_shiny_server(ws, "wsprefix"), {
    session$elapse(1000)
    expect_identical(engine_status(), "")
    # A normal independent interaction is answered while polling keeps the
    # last good payload and retries on the next tick.
    session$setInputs(explore_export = 1, explore_export_path = "")
    expect_match(export_status(), "Choose a new destination")
    session$elapse(1000)
    expect_match(engine_status(), "LISA-EXPLORE-015.*malformed JSON")
  })
  expect_lt(proc.time()[["elapsed"]] - started, 2)
})
