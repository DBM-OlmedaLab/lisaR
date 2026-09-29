test_that("opening a workspace never writes inside the source run", {
  run <- explore_fixture_run()
  before <- explore_tree_digest(run)
  root <- file.path(tempdir(), "lisa-explore-open-ws")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)

  ws <- lisa_explore_open(run, root)
  expect_s3_class(ws, "lisa_explore_workspace")
  expect_true(dir.exists(file.path(root, "extensions")))
  expect_true(file.exists(file.path(root, "session.json")))
  expect_identical(explore_tree_digest(run), before)

  # A workspace inside the immutable run is refused before anything is created.
  expect_error(lisa_explore_open(run, file.path(run, "explore")),
               "LISA-EXPLORE-011")
})

test_that("navigating and planning produce zero generator calls", {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), "lisa-explore-plan-ws")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)
  ws <- lisa_explore_open(run, root)
  before <- explore_tree_digest(run)

  traced <- explore_with_generator_trace({
    catalog <- lisa_explore_catalog(ws)
    status <- lisa_explore_status(ws)
    plans <- lapply(seq_len(min(8L, nrow(catalog))), function(index) {
      lisa_explore_plan(ws, lisaR:::lisa_explore_request_from_row(
        catalog[index, , drop = FALSE]))
    })
    list(catalog = catalog, status = status, plans = plans)
  })

  # Browsing the whole catalog, computing every product state and planning
  # several figures must not reach a generator even once.
  expect_identical(length(traced$post), 0L)
  expect_identical(traced$pathway, 0L)
  expect_gt(nrow(traced$value$catalog), 0L)
  expect_true(all(traced$value$status$state %in% lisaR:::lisa_explore_states()))
  expect_true(all(vapply(traced$value$plans,
                         function(plan) plan$state, character(1)) ==
                    "ungenerated"))
  expect_identical(explore_tree_digest(run), before)
})

test_that("an unavailable product is reported not applicable, not rendered", {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), "lisa-explore-na-ws")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)
  ws <- lisa_explore_open(run, root)

  request <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOBP-C2",
    category_id = "NO_SUCH_CATEGORY", product = "volcano")
  plan <- lisa_explore_plan(ws, request)
  expect_false(plan$applicable)
  expect_identical(plan$state, "not_applicable")
  expect_match(plan$reason, "no volcano")
  expect_error(lisa_explore_submit(ws, request, background = FALSE),
               "LISA-EXPLORE-022")
})

test_that("non-generatable products are refused rather than approximated", {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), "lisa-explore-h1-ws")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)
  ws <- lisa_explore_open(run, root)

  request <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOBP-C2",
    category_id = "SYN_SIGNAL", product = "kegg")
  traced <- explore_with_generator_trace(
    expect_error(lisa_explore_submit(ws, request, background = FALSE),
                 "LISA-EXPLORE-021"))
  expect_identical(length(traced$post), 0L)
})

test_that("one volcano request renders exactly one volcano and no sibling", {
  shared <- explore_shared()
  trace <- shared$trace

  expect_identical(trace$value$state, "completed")
  expect_false(trace$value$reused)

  # Exactly one native generator call, for the volcano builder, carrying one
  # category. No gene-set or KEGG gene-set generator ran at all.
  scripts <- explore_generator_scripts(trace)
  expect_identical(length(scripts), 1L)
  expect_identical(scripts, "build_single_de_category_volcano_overlays.R")
  expect_identical(trace$pathway, 0L)

  args <- explore_generator_args(trace)
  expect_match(args, "--categories SYN_SIGNAL(\\s|$)")
  expect_false(grepl("SYN_STRESS|SYN_MATRIX|SYN_PROLIF", args))

  # The figure is delivered as PNG and PDF with its source table and recipe,
  # and nothing else was drawn.
  artifacts <- lisa_explore_artifacts(shared$ws)
  expect_identical(nrow(artifacts), 1L)
  expect_identical(artifacts$product, "volcano")
  expect_identical(artifacts$category_id, "SYN_SIGNAL")

  artifact_dir <- file.path(shared$ws$root, "extensions", artifacts$artifact_dir)
  files <- list.files(file.path(artifact_dir, "artifacts"), recursive = TRUE)
  expect_identical(sum(grepl("[.]png$", files)), 1L)
  expect_identical(sum(grepl("[.]pdf$", files)), 1L)
  expect_identical(sum(grepl("_source[.]tsv$", files)), 1L)
  expect_identical(sum(grepl("_recipe[.]R$", files)), 1L)
  expect_true(all(grepl("^response_a/GOBP-C2/volcano/", files)))
})

test_that("rendering leaves the immutable source run unchanged", {
  shared <- explore_shared()
  expect_identical(explore_tree_digest(shared$run), shared$before)
})

test_that("two requests in different contexts create no cross-product", {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), "lisa-explore-cross-ws")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)
  ws <- lisa_explore_open(run, root)
  catalog <- lisa_explore_catalog(ws)

  # The existing selection API applies categories and products as independent
  # vectors, so asking for two categories and two products expands to four
  # units. That behaviour is unchanged and is exactly why the engine never
  # passes a vector selection.
  legacy <- lisaR:::lisa_extension_filter(
    catalog,
    lisaR:::lisa_extension_read_selection(list(selection = list(
      mode = "selected", analyses = "response_a", collections = "GOBP-C2",
      categories = c("SYN_SIGNAL", "SYN_STRESS"),
      products = c("volcano", "heatmap")))))
  expect_identical(nrow(legacy), 6L)

  # The engine asks for "volcano of SYN_SIGNAL" and "heatmap of SYN_STRESS".
  # Each resolves to exactly one unit, and the union is two units, not four.
  volcano_a <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOBP-C2",
    category_id = "SYN_SIGNAL", product = "volcano")
  heatmap_b <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOBP-C2",
    category_id = "SYN_STRESS", product = "heatmap", variant = "zscore")

  resolved <- lapply(list(volcano_a, heatmap_b), function(request) {
    lisaR:::lisa_extension_filter(
      catalog,
      lisaR:::lisa_extension_read_selection(
        lisaR:::lisa_explore_selection(request)))
  })
  expect_identical(vapply(resolved, nrow, integer(1)), c(1L, 1L))

  union <- unique(do.call(rbind, resolved))
  expect_identical(nrow(union), 2L)
  expect_setequal(paste(union$category_id, union$product),
                  c("SYN_SIGNAL volcano", "SYN_STRESS heatmap"))
  # Neither of the two cross terms exists.
  expect_false(any(paste(union$category_id, union$product) %in%
                     c("SYN_STRESS volcano", "SYN_SIGNAL heatmap")))
})

test_that("submitting the same figure twice reuses it instead of rendering again", {
  shared <- explore_shared()

  traced <- explore_with_generator_trace(
    lisa_explore_submit(shared$ws, shared$request, background = FALSE))

  # The second submission must reach no generator at all.
  expect_identical(length(traced$post), 0L)
  expect_identical(traced$pathway, 0L)
  expect_true(traced$value$reused)
  expect_identical(traced$value$state, "available")

  # And it must not have created a second job.
  jobs <- lisa_explore_poll(shared$ws)
  expect_identical(nrow(jobs), 1L)
})

test_that("an in-flight job is joined rather than duplicated", {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), "lisa-explore-dedupe-ws")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)
  ws <- lisa_explore_open(run, root)
  request <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOCC",
    category_id = "SYN_SIGNAL", product = "volcano")
  plan <- lisa_explore_plan(ws, request)

  # Simulate the first submission having reached `running` in its own process:
  # a live PID (this one) owns it. A second submission must join that job.
  job_id <- "job-inflight-0001"
  lisaR:::lisa_explore_write_job(ws, list(
    job_id = job_id, key = plan$key, request_id = plan$request_id,
    request = as.list(lisaR:::lisa_explore_request_vector(request)),
    source_run = ws$source_run, root = ws$root, state = "running",
    stage = "rendering", cancel_requested = FALSE,
    host = lisaR:::lisa_explore_host(), pid = Sys.getpid(),
    created_at = lisaR:::lisa_explore_now(), message = ""))
  lisaR:::lisa_explore_write_index(ws, list(list(
    key = plan$key, request_id = plan$request_id,
    request = as.list(lisaR:::lisa_explore_request_vector(request)),
    state = "running", job_id = job_id,
    source_manifest_hash = ws$source_manifest_hash,
    updated_at = lisaR:::lisa_explore_now(), message = "",
    artifact_dir = NULL, files = list())))

  traced <- explore_with_generator_trace(
    lisa_explore_submit(ws, request, background = FALSE))
  expect_identical(length(traced$post), 0L)
  expect_true(traced$value$reused)
  expect_identical(traced$value$job_id, job_id)
  expect_identical(traced$value$state, "running")
  expect_identical(nrow(lisa_explore_poll(ws)), 1L)
})

test_that("a job whose worker died is reported failed, never available", {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), "lisa-explore-orphan-ws")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)
  ws <- lisa_explore_open(run, root)
  request <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOMF",
    category_id = "SYN_SIGNAL", product = "volcano")
  plan <- lisa_explore_plan(ws, request)

  dead_pid <- 4000000L  # above any plausible live PID on this host
  job_id <- "job-orphan-0001"
  lisaR:::lisa_explore_write_job(ws, list(
    job_id = job_id, key = plan$key, request_id = plan$request_id,
    request = as.list(lisaR:::lisa_explore_request_vector(request)),
    source_run = ws$source_run, root = ws$root, state = "running",
    stage = "rendering", cancel_requested = FALSE,
    host = lisaR:::lisa_explore_host(), pid = dead_pid,
    created_at = lisaR:::lisa_explore_now(), message = ""))
  lisaR:::lisa_explore_write_index(ws, list(list(
    key = plan$key, request_id = plan$request_id,
    request = as.list(lisaR:::lisa_explore_request_vector(request)),
    state = "running", job_id = job_id,
    source_manifest_hash = ws$source_manifest_hash,
    updated_at = lisaR:::lisa_explore_now(), message = "",
    artifact_dir = NULL, files = list())))

  status <- lisa_explore_poll(ws, job_id)
  expect_identical(status$state, "failed")
  expect_match(status$message, "ended without completing")
  expect_identical(lisa_explore_plan(ws, request)$state, "failed")
  expect_identical(nrow(lisa_explore_artifacts(ws)), 0L)
})

test_that("a corrupted artifact is not reused", {
  shared <- explore_shared()
  artifacts <- lisa_explore_artifacts(shared$ws)
  artifact_dir <- file.path(shared$ws$root, "extensions", artifacts$artifact_dir)
  png <- list.files(file.path(artifact_dir, "artifacts"), pattern = "[.]png$",
                    recursive = TRUE, full.names = TRUE)[[1]]
  original <- readBin(png, "raw", file.info(png)$size)
  on.exit(writeBin(original, png), add = TRUE)

  expect_identical(lisa_explore_plan(shared$ws, shared$request)$state, "available")

  # Content identity fails closed: the recorded hash no longer matches.
  writeBin(c(original, as.raw(0L)), png)
  plan <- lisa_explore_plan(shared$ws, shared$request)
  expect_identical(plan$state, "ungenerated")
  expect_match(plan$reason, "changed on disk")
  expect_identical(nrow(lisa_explore_artifacts(shared$ws)), 0L)

  # A missing file is equally disqualifying.
  file.remove(png)
  plan <- lisa_explore_plan(shared$ws, shared$request)
  expect_identical(plan$state, "ungenerated")
  expect_match(plan$reason, "missing")

  writeBin(original, png)
  expect_identical(lisa_explore_plan(shared$ws, shared$request)$state, "available")
})

test_that("an incompatible identity does not reuse the old artifact", {
  shared <- explore_shared()
  ws <- shared$ws
  catalog <- lisa_explore_catalog(ws)
  row <- lisaR:::lisa_explore_catalog_row(ws, shared$request, catalog)
  key <- lisaR:::lisa_explore_key(ws, shared$request, row)

  # A different presentation policy is a different product identity.
  other_policy <- lisaR:::lisa_explore_report_policy()
  other_policy$formats$pdf <- FALSE
  expect_false(identical(key, lisaR:::lisa_explore_key(ws, shared$request, row,
                                                       other_policy)))

  # A changed source run is a different identity too, so the stored artifact
  # can never be served for it: the key simply does not match any index entry.
  moved <- ws
  moved$source_manifest_hash <- paste0(substr(ws$source_manifest_hash, 1L, 60L),
                                       "dead")
  changed_key <- lisaR:::lisa_explore_key(moved, shared$request, row)
  expect_false(identical(key, changed_key))
  expect_null(lisaR:::lisa_explore_index_find(
    lisaR:::lisa_explore_read_index(ws), changed_key))
})

test_that("closing and re-opening the workspace preserves validated work", {
  shared <- explore_shared()
  root <- shared$ws$root

  # A brand new workspace handle over the same root: nothing is carried over in
  # memory, everything must come from the persisted index.
  reopened <- lisa_explore_open(shared$run, root)
  traced <- explore_with_generator_trace({
    artifacts <- lisa_explore_artifacts(reopened)
    plan <- lisa_explore_plan(reopened, shared$request)
    submit <- lisa_explore_submit(reopened, shared$request, background = FALSE)
    list(artifacts = artifacts, plan = plan, submit = submit)
  })

  expect_identical(length(traced$post), 0L)
  expect_identical(traced$pathway, 0L)
  expect_identical(nrow(traced$value$artifacts), 1L)
  expect_identical(traced$value$plan$state, "available")
  expect_true(traced$value$submit$reused)
})

test_that("a second session does not become a second index owner", {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), "lisa-explore-lock-ws")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)
  ws <- lisa_explore_open(run, root)

  token <- lisaR:::lisa_explore_acquire_lock(ws)
  on.exit(lisaR:::lisa_explore_release_lock(ws, token), add = TRUE)
  # A live owner is respected rather than overridden.
  expect_error(lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 0.2),
               "LISA-EXPLORE-015")

  # A lock left behind by a process that no longer exists is reclaimed.
  lisaR:::lisa_explore_write_json(
    list(host = lisaR:::lisa_explore_host(), pid = 4000000L,
         start_time = "1", token = "stale-engine-test",
         acquired_at = lisaR:::lisa_explore_now()),
    file.path(lisaR:::lisa_explore_lock_path(ws), "owner.json"), run_root = root)
  replacement <- lisaR:::lisa_explore_acquire_lock(ws, timeout_seconds = 2)
  expect_true(lisaR:::lisa_explore_release_lock(ws, replacement))
})

test_that("cancelling a queued job is immediate and stops no foreign process", {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), "lisa-explore-cancel-ws")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  unlink(root, recursive = TRUE, force = TRUE)
  ws <- lisa_explore_open(run, root)
  request <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOCC",
    category_id = "SYN_STRESS", product = "volcano")
  plan <- lisa_explore_plan(ws, request)

  job_id <- "job-queued-0001"
  lisaR:::lisa_explore_write_job(ws, list(
    job_id = job_id, key = plan$key, request_id = plan$request_id,
    request = as.list(lisaR:::lisa_explore_request_vector(request)),
    source_run = ws$source_run, root = ws$root, state = "queued",
    stage = "queued", cancel_requested = FALSE,
    host = lisaR:::lisa_explore_host(), pid = NA_integer_,
    created_at = lisaR:::lisa_explore_now(), message = ""))

  outcome <- lisa_explore_cancel(ws, job_id)
  expect_true(outcome$cancelled)
  expect_identical(outcome$state, "failed")
  expect_identical(lisa_explore_poll(ws, job_id)$state, "failed")

  # A running job is asked to stop; it is never killed.
  lisaR:::lisa_explore_write_job(ws, list(
    job_id = "job-running-0001", key = plan$key, request_id = plan$request_id,
    request = as.list(lisaR:::lisa_explore_request_vector(request)),
    source_run = ws$source_run, root = ws$root, state = "running",
    stage = "rendering", cancel_requested = FALSE,
    host = lisaR:::lisa_explore_host(), pid = Sys.getpid(),
    created_at = lisaR:::lisa_explore_now(), message = ""))
  running <- lisa_explore_cancel(ws, "job-running-0001")
  expect_false(running$cancelled)
  expect_match(running$message, "safe point")
  expect_true(isTRUE(lisaR:::lisa_explore_read_job(ws, "job-running-0001")$cancel_requested))
})
