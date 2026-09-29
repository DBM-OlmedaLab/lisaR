test_that("pipeline.workers has a strict Linux-multicore contract", {
  policies <- list(
    de_table_duplicate_policy = "error",
    matrix_duplicate_policy = "error",
    mapped_id_collision_policy = "error"
  )
  config <- function(workers = NULL) {
    pipeline <- list(
      schema_version = "1.0.0", profile = "targeted",
      evidence_mode = "full_de", duplicate_policies = policies
    )
    if (!is.null(workers)) pipeline$workers <- workers
    list(pipeline = pipeline)
  }
  expect_identical(lisaR:::lisa_validate_pipeline_config(config())$workers, 4L)
  expect_identical(lisaR:::lisa_validate_pipeline_config(config(10))$workers, 10L)
  invalid <- list(0, -1, 1.5, TRUE, "10", NA_real_, Inf, .Machine$integer.max + 1)
  for (value in invalid) {
    expect_error(lisaR:::lisa_validate_pipeline_config(config(value)), "pipeline.workers")
  }
})

test_that("parallel policy is Linux multicore and serial on Windows and macOS", {
  linux <- lisaR:::lisa_parallel_plan(
    4L, 10L, platform = "linux", check_limit = "false",
    allocation = list(workers = 8L, reason = "visible_cores")
  )
  windows <- lisaR:::lisa_parallel_plan(
    4L, 10L, platform = "windows", check_limit = "false",
    allocation = list(workers = 8L, reason = "visible_cores")
  )
  macos <- lisaR:::lisa_parallel_plan(
    4L, 10L, platform = "macos", check_limit = "false",
    allocation = list(workers = 8L, reason = "visible_cores")
  )

  expect_identical(linux$backend_requested, "auto")
  expect_identical(linux$backend_effective, "multicore")
  expect_identical(linux$workers_effective, 4L)
  for (plan in list(windows, macos)) {
    expect_identical(plan$backend_effective, "sequential")
    expect_identical(plan$workers_effective, 1L)
    expect_match(plan$cap_reason, "unsupported_platform_serial", fixed = TRUE)
  }
  one_task_windows <- lisaR:::lisa_parallel_plan(
    4L, 1L, platform = "windows", check_limit = "false",
    allocation = list(workers = 8L, reason = "visible_cores")
  )
  expect_match(one_task_windows$cap_reason, "task_count", fixed = TRUE)
  expect_match(one_task_windows$cap_reason, "unsupported_platform_serial", fixed = TRUE)

  task_ids <- c("platform-a", "platform-b")
  serial_reference <- lisaR:::lisa_pipeline_map(task_ids, identity, workers = 1L, platform = "linux")
  expect_identical(
    lisaR:::lisa_pipeline_map(task_ids, identity, workers = 4L, platform = "windows"),
    serial_reference
  )
  expect_identical(
    lisaR:::lisa_pipeline_map(task_ids, identity, workers = 4L, platform = "macos"),
    serial_reference
  )
  expect_error(
    lisaR:::lisa_map_with_task_seeds(task_ids, identity, backend = "parallel", workers = 2L),
    "should be one of"
  )
})

test_that("legacy backend options cannot select PSOCK", {
  skip_if_not(identical(lisaR:::lisa_runtime_platform(), "linux"), "Linux multicore contract only.")
  skip_if_not(parallel::detectCores(logical = TRUE) > 1L, "Multicore verification requires at least two cores.")
  old_backend <- getOption("lisaR.parallel_backend", NULL)
  on.exit(options(lisaR.parallel_backend = old_backend), add = TRUE)
  options(lisaR.parallel_backend = "parallel")
  observed <- lisaR:::lisa_pipeline_map(
    c("backend-a", "backend-b"),
    function(task_id) getOption("lisaR.outer_parallel_context")$backend_effective,
    workers = 2L
  )
  expect_identical(unname(unlist(observed)), c("multicore", "multicore"))
})

test_that("task count, R check, and visible allocation cap effective workers", {
  task_cap <- lisaR:::lisa_parallel_plan(
    10L, 3L, platform = "linux", check_limit = "false",
    allocation = list(workers = 24L, reason = "visible_cores")
  )
  check_cap <- lisaR:::lisa_parallel_plan(
    10L, 20L, platform = "linux", check_limit = "TRUE",
    allocation = list(workers = 24L, reason = "visible_cores")
  )
  allocation_cap <- lisaR:::lisa_parallel_plan(
    10L, 20L, platform = "linux", check_limit = "false",
    allocation = list(workers = 3L, reason = "scheduler_allocation")
  )

  expect_identical(task_cap$workers_effective, 3L)
  expect_match(task_cap$cap_reason, "task_count", fixed = TRUE)
  expect_identical(check_cap$workers_effective, 2L)
  expect_match(check_cap$cap_reason, "r_check_limit_cores", fixed = TRUE)
  expect_identical(allocation_cap$workers_effective, 3L)
  expect_match(allocation_cap$cap_reason, "scheduler_allocation", fixed = TRUE)
  expect_identical(
    lisaR:::lisa_visible_worker_limit(
      detected_cores = 16L,
      scheduler = c(SLURM_CPUS_PER_TASK = "6", PBS_NP = "")
    ),
    list(workers = 6L, reason = "scheduler_allocation")
  )
  expect_identical(lisaR:::lisa_positive_worker_limit("1.5"), NA_integer_)
})

test_that("workers 1, 2, and 4 preserve task RNG, order, and caller RNG on Linux", {
  skip_if_not(identical(lisaR:::lisa_runtime_platform(), "linux"), "Linux multicore contract only.")
  old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit(if (is.na(old_limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  task_ids <- sprintf("task-%02d", seq_len(10L))
  fun <- function(task_id) c(task_id = task_id, value = sprintf("%.17g", stats::runif(1)))
  set.seed(4401)
  before <- .Random.seed
  serial <- lisaR:::lisa_pipeline_map(task_ids, fun, workers = 1L)
  expect_identical(.Random.seed, before)
  for (workers in c(2L, 4L)) {
    multicore <- lisaR:::lisa_pipeline_map(task_ids, fun, workers = workers)
    expect_identical(.Random.seed, before)
    expect_identical(multicore, serial)
    expect_identical(names(multicore), task_ids)
  }
  expect_named(lisaR:::lisa_pipeline_map(character(), identity, workers = 4L), character())
  expect_error(lisaR:::lisa_pipeline_map(c("same", "same"), identity, workers = 4L), "unique")
})

test_that("typed task envelopes fail closed on worker errors and NULL values", {
  skip_if_not(identical(lisaR:::lisa_runtime_platform(), "linux"), "Linux multicore contract only.")
  old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit(if (is.na(old_limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  task_ids <- c("task-ok", "task-bad")

  expect_error(
    lisaR:::lisa_pipeline_map(
      task_ids,
      function(task_id) if (identical(task_id, "task-bad")) stop("intentional worker error") else task_id,
      workers = 2L
    ),
    "LISA-PARALLEL-004.*task-bad.*intentional worker error"
  )
  expect_error(
    lisaR:::lisa_pipeline_map(
      task_ids,
      function(task_id) if (identical(task_id, "task-bad")) NULL else task_id,
      workers = 2L
    ),
    "LISA-PARALLEL-004.*task-bad.*returned NULL"
  )
})

test_that("parent reconciles envelope cardinality, names, and schema", {
  skip_if_not(identical(lisaR:::lisa_runtime_platform(), "linux"), "Linux multicore contract only.")
  old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit(if (is.na(old_limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  task_ids <- c("task-a", "task-b")
  valid <- function(task_id) list(task_id = task_id, ok = TRUE, value = task_id, error = NULL)

  expect_error(
    testthat::with_mocked_bindings(
      lisaR:::lisa_pipeline_map(task_ids, identity, workers = 2L),
      lisa_multicore_apply = function(task_ids, worker, workers, preschedule = FALSE) {
        stats::setNames(list(valid("task-a")), "task-a")
      },
      .package = "lisaR"
    ),
    "LISA-PARALLEL-001"
  )
  expect_error(
    testthat::with_mocked_bindings(
      lisaR:::lisa_pipeline_map(task_ids, identity, workers = 2L),
      lisa_multicore_apply = function(task_ids, worker, workers, preschedule = FALSE) {
        stats::setNames(list(valid("task-a"), valid("task-b")), c("wrong-a", "wrong-b"))
      },
      .package = "lisaR"
    ),
    "LISA-PARALLEL-002"
  )
  expect_error(
    testthat::with_mocked_bindings(
      lisaR:::lisa_pipeline_map(task_ids, identity, workers = 2L),
      lisa_multicore_apply = function(task_ids, worker, workers, preschedule = FALSE) {
        stats::setNames(
          list(valid("task-a"), list(task_id = "task-b", ok = TRUE, value = "task-b")),
          task_ids
        )
      },
      .package = "lisaR"
    ),
    "LISA-PARALLEL-003.*task-b"
  )
})

test_that("outer worker context and thread controls are observable and recorded", {
  root <- tempfile("lisaR-parallel-receipt-")
  dir.create(root)
  old_root <- getOption("lisaR.run_root", NULL)
  on.exit(options(lisaR.run_root = old_root), add = TRUE)
  options(lisaR.run_root = root)
  observed <- lisaR:::lisa_pipeline_map(
    c("single_de:A:PATHWAYS", "single_de:A:HALLMARKS"),
    function(task_id) {
      list(
        context = getOption("lisaR.outer_parallel_context"),
        threads = lisaR:::lisa_thread_environment()
      )
    },
    workers = 4L,
    platform = "windows"
  )

  expect_true(all(vapply(observed, function(value) {
    identical(value$context$backend_effective, "sequential") &&
      identical(value$context$workers_effective, 1L) &&
      identical(unname(value$threads), rep("1", length(value$threads)))
  }, logical(1))))
  execution <- lisaR:::read_lisa_tsv(file.path(root, "parallel_execution.tsv"))
  seeds <- lisaR:::read_lisa_tsv(file.path(root, "parallel_task_seeds.tsv"))
  expect_identical(execution$backend_requested, "auto")
  expect_identical(execution$backend_effective, "sequential")
  expect_identical(as.integer(execution$workers_requested), 4L)
  expect_identical(as.integer(execution$workers_effective), 1L)
  expect_match(execution$cap_reason, "unsupported_platform_serial", fixed = TRUE)
  expect_identical(seeds$task_id, c("single_de:A:PATHWAYS", "single_de:A:HALLMARKS"))
  expect_equal(length(unique(seeds$effective_seed)), 2L)
})

test_that("parallel SHA-256 uses coarse chunks without changing results", {
  skip_if_not(identical(lisaR:::lisa_runtime_platform(), "linux"), "Linux multicore contract only.")
  old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit(if (is.na(old_limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  root <- tempfile("lisaR-worker-hash-")
  dir.create(root)
  files <- file.path(root, sprintf("file-%02d.txt", seq_len(31L)))
  for (i in seq_along(files)) writeLines(rep(sprintf("row-%02d", i), i), files[[i]])
  serial <- lisaR:::lisa_hash_files(files, workers = 1L)
  for (workers in c(2L, 4L)) {
    multicore <- lisaR:::lisa_hash_files(files, workers = workers)
    expect_identical(multicore, serial)
    expect_length(multicore, length(files))
  }
})

test_that("parallel SHA-256 rejects short and misnamed worker results", {
  skip_if_not(identical(lisaR:::lisa_runtime_platform(), "linux"), "Linux multicore contract only.")
  skip_if(lisaR:::lisa_parallel_plan(2L, 2L)$workers_effective < 2L, "Two effective workers required.")
  old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit(if (is.na(old_limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  root <- tempfile("lisaR-worker-invalid-hash-")
  dir.create(root)
  files <- file.path(root, sprintf("file-%02d.txt", seq_len(4L)))
  for (path in files) writeLines("artifact", path)
  valid <- paste(rep("0", 64L), collapse = "")

  testthat::with_mocked_bindings(
    expect_error(lisaR:::lisa_hash_files(files, workers = 2L), "LISA-SHA256-004"),
    lisa_sha256_files_batch = function(files, ...) {
      files <- normalizePath(files, winslash = "/", mustWork = TRUE)
      stats::setNames(rep(valid, length(files) - 1L), files[-length(files)])
    },
    .package = "lisaR"
  )
  testthat::with_mocked_bindings(
    expect_error(lisaR:::lisa_hash_files(files, workers = 2L), "LISA-SHA256-004"),
    lisa_sha256_files_batch = function(files, ...) {
      files <- normalizePath(files, winslash = "/", mustWork = TRUE)
      stats::setNames(rep(valid, length(files)), rev(files))
    },
    .package = "lisaR"
  )
})

test_that("parallel SHA-256 worker failures retain a typed integrity condition", {
  skip_if_not(identical(lisaR:::lisa_runtime_platform(), "linux"), "Linux multicore contract only.")
  skip_if(lisaR:::lisa_parallel_plan(2L, 2L)$workers_effective < 2L, "Two effective workers required.")
  old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit(if (is.na(old_limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  root <- tempfile("lisaR-worker-hash-error-")
  dir.create(root)
  files <- file.path(root, c("a.txt", "b.txt"))
  writeLines("a", files[[1L]])
  writeLines("b", files[[2L]])

  condition <- tryCatch(
    testthat::with_mocked_bindings(
      lisaR:::lisa_hash_files(files, workers = 2L),
      lisa_sha256_files_batch = function(...) stop("intentional hash worker failure"),
      .package = "lisaR"
    ),
    error = identity
  )
  expect_s3_class(condition, "lisa_sha256_error")
  expect_identical(condition$code, "LISA-SHA256-006")
  expect_match(conditionMessage(condition), "LISA-PARALLEL-004", fixed = TRUE)
  expect_match(conditionMessage(condition), "hash:0001", fixed = TRUE)
  expect_match(conditionMessage(condition), "intentional hash worker failure", fixed = TRUE)
})

test_that("parallel SHA-256 batches avoid shared subprocess diagnostics under stress", {
  skip_on_os("windows")
  old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit(if (is.na(old_limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  root <- tempfile("lisaR-worker-batch-")
  dir.create(root)
  files <- file.path(root, sprintf("file-%03d.txt", seq_len(100L)))
  for (i in seq_along(files)) writeLines(rep(sprintf("row-%03d", i), i %% 7L + 1L), files[[i]])
  old_root <- getOption("lisaR.run_root", NULL)
  old_batch <- getOption("lisaR.sha256_batch_size", NULL)
  options(lisaR.run_root = root)
  options(lisaR.sha256_batch_size = 7L)
  on.exit({
    options(lisaR.run_root = old_root)
    options(lisaR.sha256_batch_size = old_batch)
  }, add = TRUE)
  testthat::local_mocked_bindings(
    lisa_run_subprocess = function(...) stop("batch hashing must not create subprocess diagnostics"),
    .package = "lisaR"
  )
  serial <- lisaR:::lisa_hash_files(files, workers = 1L)
  for (i in seq_len(8L)) {
    expect_identical(lisaR:::lisa_hash_files(files, workers = 10L), serial)
  }
  expect_false(dir.exists(file.path(root, ".lisa_subprocess")))
})

test_that("parallel SHA-256 rejects a malformed chunk result", {
  skip_if_not(identical(lisaR:::lisa_runtime_platform(), "linux"), "Linux multicore contract only.")
  skip_if(lisaR:::lisa_parallel_plan(2L, 2L)$workers_effective < 2L, "Two effective workers required.")
  old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit(if (is.na(old_limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  root <- tempfile("lisaR-worker-hash-malformed-")
  dir.create(root)
  files <- file.path(root, c("a.txt", "b.txt"))
  writeLines("a", files[[1L]])
  writeLines("b", files[[2L]])

  expect_error(
    testthat::with_mocked_bindings(
      lisaR:::lisa_hash_files(files, workers = 2L),
      lisa_sha256_files_batch = function(files, run_root = NULL) character(),
      .package = "lisaR"
    ),
    "LISA-SHA256-004"
  )
})

test_that("SHA-256 digest batches reject invalid size options", {
  root <- tempfile("lisaR-hash-batch-option-")
  dir.create(root)
  path <- file.path(root, "artifact.txt")
  writeLines("artifact", path)
  old_batch <- getOption("lisaR.sha256_batch_size", NULL)
  on.exit(options(lisaR.sha256_batch_size = old_batch), add = TRUE)
  for (value in list(0L, -1L, 1.5, NA_integer_, "128")) {
    options(lisaR.sha256_batch_size = value)
    expect_error(lisaR:::lisa_sha256_files_batch(path), "one positive integer")
  }
})

test_that("run manifests exclude shared subprocess diagnostics", {
  root <- tempfile("lisaR-manifest-subprocess-")
  dir.create(root)
  writeLines("artifact", file.path(root, "artifact.txt"))
  dir.create(file.path(root, ".lisa_subprocess"))
  writeLines("diagnostic", file.path(root, ".lisa_subprocess", "stdout.txt"))
  manifest <- lisaR:::lisa_run_manifest(root)
  expect_identical(manifest$path, "artifact.txt")
})

test_that("single-DE orchestration is identical with workers 1, 2, and 4 on Linux", {
  skip_if_not(identical(lisaR:::lisa_runtime_platform(), "linux"), "Linux multicore contract only.")
  old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit(if (is.na(old_limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  testthat::local_mocked_bindings(
    run_LISA_DE = function(...) {
      args <- list(...)
      if (!dir.exists(dirname(args$output_dir))) {
        stop("single-DE collection parent was not created before dispatch")
      }
      dir.create(args$output_dir, recursive = FALSE, showWarnings = FALSE)
      write.table(
        data.frame(analysis_id = args$comparison_name, collection = args$universes),
        file.path(args$output_dir, "scientific.tsv"), sep = "\t",
        row.names = FALSE, quote = FALSE
      )
      list(status = "completed")
    },
    .package = "lisaR"
  )
  de <- data.frame(
    analysis_id = c("A", "B"), de_path = c("A.tsv", "B.tsv"),
    species = "Homo sapiens", stringsAsFactors = FALSE
  )
  collections <- c("PATHWAYS", "HALLMARKS", "GOBP-C2", "GOMF", "GOCC")
  registry <- data.frame(
    analysis_collection = collections,
    dictionary_id = paste0("DICT_", collections),
    run_ora = FALSE, run_category_pathway_plots = FALSE,
    post_lisa_profile = "full", allow_missing_dictionary = FALSE,
    stringsAsFactors = FALSE
  )
  run_one <- function(root, workers) lisaR:::run_lisa_pipeline_single_de(
    de, registry, getwd(), root, "term2gene.tsv", "dictionary.tsv",
    "category-map.tsv", NULL, "core", "semantic", "png", "tsv",
    workers = workers
  )
  serial_root <- tempfile("lisaR-serial-")
  serial <- run_one(serial_root, 1L)
  stable_columns <- c("analysis_id", "collection", "dictionary_id", "status", "message")
  inventory <- function(root) {
    paths <- list.files(root, recursive = TRUE, full.names = TRUE)
    relative <- substring(paths, nchar(root) + 2L)
    stats::setNames(vapply(paths, lisaR:::lisa_sha256_file, character(1)), relative)
  }
  for (workers in c(2L, 4L)) {
    multicore_root <- tempfile(sprintf("lisaR-multicore-%d-", workers))
    multicore <- run_one(multicore_root, workers)
    expect_identical(multicore[stable_columns], serial[stable_columns])
    expect_identical(inventory(multicore_root), inventory(serial_root))
  }
  expect_identical(serial$analysis_id, rep(c("A", "B"), each = length(collections)))
  expect_identical(serial$collection, rep(collections, times = 2L))
})

test_that("post-lisa units finish before the serial root report", {
  skip_if_not(identical(lisaR:::lisa_runtime_platform(), "linux"), "Linux multicore contract only.")
  old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit(if (is.na(old_limit)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit), add = TRUE)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "false")
  arg_value <- function(args, key) {
    index <- match(key, args)
    if (is.na(index) || index == length(args)) "" else args[[index + 1L]]
  }
  testthat::local_mocked_bindings(
    # This fixture tests task joining, not saved scientific inputs. The real
    # serial inference phase has its own numerical and integration tests.
    lisa_run_category_inference = function(...) data.frame(),
    lisa_resolve_package_dir = function(...) normalizePath(getwd()),
    lisa_write_code_identity_ledger = function(package_dir, script_names, path,
                                               executable_recipes = character()) {
      structure(
        list(path = path, sha256 = "test", entries = data.frame()),
        class = "lisa_code_identity_ledger"
      )
    },
    lisa_build_single_gene_level_tables = function(analysis_id, collection, outputs_root, file_label_prefix) {
      marker <- file.path(outputs_root, "markers", paste(analysis_id, collection, "tables", sep = "-"))
      dir.create(dirname(marker), recursive = TRUE, showWarnings = FALSE)
      writeLines("complete", marker)
      list(status = "completed", output_path = marker, message = "")
    },
    lisa_run_post_script = function(package_dir, script_name, args,
                                    trusted_run_root = NULL, code_ledger = NULL) {
      project <- arg_value(args, "--project-dir")
      marker_dir <- file.path(project, "outputs", "markers")
      dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)
      if (identical(script_name, "build_LISA_report.R")) {
        ready <- list.files(marker_dir, recursive = TRUE)
        status <- if (length(ready) >= 40L) "completed" else "failed"
        return(list(status = status, output_path = "", message = paste("markers", length(ready))))
      }
      identifier <- arg_value(args, "--analysis-id")
      if (!nzchar(identifier)) identifier <- arg_value(args, "--contrast-id")
      collection <- arg_value(args, "--universe")
      marker <- file.path(marker_dir, paste(script_name, identifier, collection, sep = "-"))
      writeLines("complete", marker)
      list(status = "completed", output_path = marker, message = "")
    },
    .package = "lisaR"
  )
  collections <- c("PATHWAYS", "HALLMARKS", "GOBP-C2", "GOMF", "GOCC")
  registry <- data.frame(
    analysis_collection = collections, run_enrichmentmap = TRUE,
    run_gene_cards = TRUE, run_volcano_overlays = TRUE,
    run_recurrent_gene_screen = TRUE, run_kegg_layers = FALSE,
    stringsAsFactors = FALSE
  )
  single <- expand.grid(collection = collections, analysis_id = c("A", "B"),
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  single$status <- "completed"
  contrasts <- data.frame(
    contrast_id = "A_vs_B", collection = collections, status = "completed",
    stringsAsFactors = FALSE
  )
  contrast_index <- data.frame(
    contrast_id = "A_vs_B", analysis_a = "A", analysis_b = "B",
    stringsAsFactors = FALSE
  )
  run_one <- function(root, workers) {
    dir.create(file.path(root, "outputs"), recursive = TRUE, showWarnings = FALSE)
    lisaR:::run_lisa_pipeline_post_lisa(
      de_index = data.frame(), contrast_index = contrast_index, registry = registry,
      output_dir = root, outputs_root = file.path(root, "outputs"),
      single_status = single, contrast_status = contrasts,
      file_label_prefix = "semantic",
      run_gene_level = TRUE, run_reports = TRUE, run_kegg = FALSE,
      workers = workers
    )
  }
  serial <- run_one(tempfile("lisaR-post-serial-"), 1L)
  parallel <- run_one(tempfile("lisaR-post-parallel-"), 4L)
  stable <- function(value) {
    value[c("timestamp", "output_path")] <- NULL
    value
  }
  expect_identical(stable(parallel), stable(serial))
  expect_identical(tail(parallel$stage, 1L), "root_html_report")
  expect_identical(tail(parallel$status, 1L), "completed")
})
