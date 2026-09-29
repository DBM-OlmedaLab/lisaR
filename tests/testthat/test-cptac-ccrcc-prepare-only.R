# Focused, installed tests for the CPTAC ccRCC `LISAR_CPTAC_CCRCC_PREPARE_ONLY`
# mode. These never run DE/limma, GSEA, reports or KEGG. The first test runs
# the real installed example script (validate_lisa_config()/plan_lisa_outputs()
# for real) against an explicitly supplied bundle and resource cache,
# and only structurally proves run_lisa() was skipped (no runs/ directory, no
# run receipt). The second test additionally proves this at the function-call
# level with a mocked lisaR in an isolated subprocess, and confirms the
# unset-flag default still dispatches the full run.

cptac_prepare_only_paths <- function() {
  list(
    bundle = Sys.getenv("LISAR_TEST_CPTAC_BUNDLE", ""),
    cache = Sys.getenv("LISAR_TEST_RESOURCE_CACHE", "")
  )
}

test_that("installed CPTAC prepare-only mode really validates, copies and promotes, then stops before run_lisa", {
  paths <- cptac_prepare_only_paths()
  skip_if_not(
    dir.exists(paths$bundle) && dir.exists(paths$cache),
    "set LISAR_TEST_CPTAC_BUNDLE and LISAR_TEST_RESOURCE_CACHE for this integration test"
  )
  script <- system.file(
    "examples/cptac-ccrcc/prepare_and_run_cptac_ccrcc.R", package = "lisaR", mustWork = TRUE
  )
  root <- withr::local_tempdir("cptac-prepare-only-")
  project <- file.path(root, "prepared-project")
  rscript <- file.path(R.home("bin"), "Rscript")
  run_env <- c(
    paste0("LISAR_CPTAC_CCRCC_BUNDLE=", paths$bundle),
    paste0("LISAR_CPTAC_CCRCC_PROJECT_DIR=", project),
    paste0("LISAR_CPTAC_CCRCC_RESOURCE_CACHE=", paths$cache),
    "LISAR_CPTAC_CCRCC_PREPARE_ONLY=1",
    paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep))
  )

  out <- suppressWarnings(system2(
    rscript, c("--vanilla", script), env = run_env, stdout = TRUE, stderr = TRUE
  ))
  status <- attr(out, "status")
  expect_true(is.null(status) || identical(status, 0L))
  expect_true(any(grepl("^CPTAC_CCRCC_LOCAL_ROUTE=PREPARE_ONLY_PASS$", out)))
  expect_true(any(out == paste0("PROJECT_DIR=", project)))
  next_line <- grep("^NEXT_RUN_COMMAND=", out, value = TRUE)
  expect_length(next_line, 1L)
  expect_match(next_line, "lisaR::run_lisa", fixed = TRUE)
  expect_match(next_line, "study.yml", fixed = TRUE)

  # Real preparation artefacts exist: the bundle was validated and copied,
  # normal resources resolved, and the project promoted.
  expect_true(dir.exists(project))
  expect_true(file.exists(file.path(project, "study.yml")))
  expect_true(file.exists(file.path(project, "data", "tumor_nat.tsv")))
  expect_true(file.exists(file.path(project, "logs", "validation_readiness.tsv")))
  expect_true(file.exists(file.path(project, "logs", "output_plan_summary.tsv")))
  expect_identical(
    unname(tools::md5sum(file.path(project, "data", "tumor_nat.tsv"))),
    unname(tools::md5sum(file.path(paths$bundle, "de", "de_primary_conservative_80pairs.tsv")))
  )

  # run_lisa()/verify_lisa_run() were never reached: no runs/ directory and no
  # run receipt, even though the promoted project and config exist.
  expect_false(dir.exists(file.path(project, "out")))
  expect_false(file.exists(file.path(project, "logs", "run_verification.tsv")))
  expect_false(file.exists(file.path(project, "logs", "cptac_ccrcc_prepared_receipt.tsv")))

  # Preservation check: the existing-project guard is untouched under
  # prepare-only. Repeating the same call must fail closed, not silently
  # overwrite or re-promote the already-prepared project.
  out2 <- suppressWarnings(system2(
    rscript, c("--vanilla", script), env = run_env, stdout = TRUE, stderr = TRUE
  ))
  status2 <- attr(out2, "status")
  expect_true(!is.null(status2) && !identical(status2, 0L))
  expect_true(any(grepl("already exists; preserve", out2, fixed = TRUE)))
  expect_false(dir.exists(file.path(project, "out")))
  expect_true(file.exists(file.path(project, "study.yml")))
})

test_that("prepare-only never calls run_lisa()/verify_lisa_run(), and the unset default still does", {
  paths <- cptac_prepare_only_paths()
  skip_if_not(dir.exists(paths$bundle), "set LISAR_TEST_CPTAC_BUNDLE for this integration test")
  script <- system.file(
    "examples/cptac-ccrcc/prepare_and_run_cptac_ccrcc.R", package = "lisaR", mustWork = TRUE
  )

  driver_root <- withr::local_tempdir("cptac-prepare-only-mock-")
  mock_src <- file.path(driver_root, "mock-lisaR")
  dir.create(file.path(mock_src, "R"), recursive = TRUE)
  writeLines(c(
    "Package: lisaR", "Type: Package", "Title: prepare-only mock", "Version: 0.0.0",
    "Authors@R: person('T', 'A', email = 'test@example.invalid', role = c('aut', 'cre'))",
    "Description: Test-only mock for the CPTAC prepare-only dispatch boundary.",
    "License: GPL-3", "Encoding: UTF-8"
  ), file.path(mock_src, "DESCRIPTION"))
  writeLines(
    "export(lisa_sha256_file, lisa_dictionary_cache_root, validate_lisa_config, plan_lisa_outputs, run_lisa, verify_lisa_run)",
    file.path(mock_src, "NAMESPACE")
  )
  writeLines(c(
    ".mock_state <- new.env(parent = emptyenv())", ".mock_state$events <- character()",
    "lisa_sha256_file <- function(path) sub(' .*$', '', system2('sha256sum', shQuote(path), stdout = TRUE))",
    "lisa_dictionary_cache_root <- function() tempdir()",
    "lisa_compact_example_project <- function(project, example, config_path) {",
    "  cfg <- yaml::read_yaml(config_path); cfg$pipeline$output_dir <- 'out'",
    "  yaml::write_yaml(cfg, file.path(project, 'study.yml')); invisible(NULL)",
    "}",
    "validate_lisa_config <- function(path, check_files = TRUE, strict = TRUE) {",
    "  .mock_state$events <- c(.mock_state$events, 'validate'); cfg <- yaml::read_yaml(path)",
    "  list(readiness = data.frame(gate = 'PASS'), resources = data.frame(resource = 'mock'),",
    "       execution_ready = TRUE, config = cfg, analyses = 1L, contrasts = 0L)",
    "}",
    "plan_lisa_outputs <- function(path) { .mock_state$events <- c(.mock_state$events, 'plan'); list(summary = data.frame(plan = 'mock')) }",
    "run_lisa <- function(path) {",
    "  .mock_state$events <- c(.mock_state$events, 'run')",
    "  output_dir <- file.path(dirname(dirname(path)), 'runs', 'cptac-ccrcc-global-proteomic-standard')",
    "  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)",
    "  list(output_dir = output_dir)",
    "}",
    "verify_lisa_run <- function(output_dir) { .mock_state$events <- c(.mock_state$events, 'verify'); list(gate = 'PASS', findings = character()) }"
  ), file.path(mock_src, "R", "mock.R"))
  mock_lib <- file.path(driver_root, "mock-library")
  dir.create(mock_lib)
  install_status <- system2(
    file.path(R.home("bin"), "R"),
    c("CMD", "INSTALL", "--no-multiarch", "--no-docs", paste0("--library=", mock_lib), mock_src),
    stdout = FALSE, stderr = FALSE
  )
  expect_identical(install_status, 0L)

  driver <- file.path(driver_root, "driver.R")
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "script <- args[[1]]; bundle <- args[[2]]",
    "state <- get('.mock_state', envir = asNamespace('lisaR'))",
    "run_scenario <- function(project, prepare_only) {",
    "  Sys.setenv(LISAR_CPTAC_CCRCC_BUNDLE = bundle, LISAR_CPTAC_CCRCC_PROJECT_DIR = project,",
    "             LISAR_CPTAC_CCRCC_RESOURCE_CACHE = file.path(dirname(project), 'cache'),",
    "             LISAR_CPTAC_CCRCC_PREPARE_ONLY = if (prepare_only) '1' else '')",
    "  state$events <- character()",
    "  result <- tryCatch(source(script, local = new.env(parent = globalenv())), error = identity)",
    "  list(events = state$events, ok = !inherits(result, 'error'))",
    "}",
    "root <- tempfile('cptac-prepare-only-mock-run-')",
    "dir.create(root, recursive = TRUE)",
    "default_result <- run_scenario(file.path(root, 'default-project'), prepare_only = FALSE)",
    "prepare_result <- run_scenario(file.path(root, 'prepare-only-project'), prepare_only = TRUE)",
    "cat('SCENARIO=default EVENTS=', paste(default_result$events, collapse=','), ' OK=', default_result$ok, '\\n', sep='')",
    "cat('SCENARIO=prepare_only EVENTS=', paste(prepare_result$events, collapse=','), ' OK=', prepare_result$ok, '\\n', sep='')"
  ), driver)

  rscript <- file.path(R.home("bin"), "Rscript")
  run_env <- c(paste0("R_LIBS=", paste(c(mock_lib, .libPaths()), collapse = .Platform$path.sep)))
  out <- suppressWarnings(system2(
    rscript, c("--vanilla", driver, script, paths$bundle), env = run_env, stdout = TRUE, stderr = TRUE
  ))
  status <- attr(out, "status")
  expect_true((is.null(status) || identical(status, 0L)))

  default_line <- grep("^SCENARIO=default ", out, value = TRUE)
  prepare_line <- grep("^SCENARIO=prepare_only ", out, value = TRUE)
  expect_length(default_line, 1L)
  expect_length(prepare_line, 1L)
  expect_match(default_line, "EVENTS=validate,plan,run,verify", fixed = TRUE)
  expect_match(default_line, "OK=TRUE", fixed = TRUE)
  expect_match(prepare_line, "EVENTS=validate,plan ", fixed = TRUE)
  expect_match(prepare_line, "OK=TRUE", fixed = TRUE)
})
