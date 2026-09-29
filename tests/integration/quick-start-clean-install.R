#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
work_root <- if (length(args)) args[[1L]] else tempfile("lisaR-clean-install-")
expected_library <- Sys.getenv("LISAR_EXPECTED_LIBRARY", unset = "")
if (dir.exists(work_root) || file.exists(work_root)) {
  stop("Quick Start integration destination must not already exist: ", work_root,
       call. = FALSE)
}
dir.create(work_root, recursive = TRUE)
on.exit(unlink(work_root, recursive = TRUE, force = TRUE), add = TRUE)

Sys.setenv(
  http_proxy = "http://127.0.0.1:9",
  https_proxy = "http://127.0.0.1:9",
  HTTP_PROXY = "http://127.0.0.1:9",
  HTTPS_PROXY = "http://127.0.0.1:9"
)
library(lisaR)
installed_root <- normalizePath(
  system.file(package = "lisaR"), winslash = "/", mustWork = TRUE
)
if (nzchar(expected_library)) {
  expected_library <- normalizePath(
    expected_library, winslash = "/", mustWork = TRUE
  )
  if (!identical(dirname(installed_root), expected_library)) {
    stop(
      "Quick Start loaded lisaR outside the clean CI library: ",
      installed_root, call. = FALSE
    )
  }
}
project <- lisa_init_project(file.path(work_root, "project"))
config <- file.path(project, "study.yml")
validation <- validate_lisa_config(config)
stopifnot(
  isTRUE(validation$valid),
  identical(validation$analyses, 2L),
  identical(validation$contrasts, 1L),
  identical(as.integer(validation$config$pipeline$workers), 4L),
  !isTRUE(validation$config$pipeline$dry_run)
)
estimate <- plan_lisa_outputs(config)
stopifnot(identical(estimate$summary$mode, "full"))
result <- suppressWarnings(run_lisa(config))
check <- verify_lisa_run(result$output_dir)
parallel_path <- file.path(result$output_dir, "parallel_execution.tsv")
if (!file.exists(parallel_path)) {
  stop("Quick Start did not record parallel_execution.tsv.", call. = FALSE)
}
parallel_execution <- utils::read.delim(
  parallel_path, check.names = FALSE, stringsAsFactors = FALSE
)
required_parallel_columns <- c(
  "platform", "backend_effective", "workers_requested", "workers_effective"
)
if (!all(required_parallel_columns %in% names(parallel_execution)) ||
    !nrow(parallel_execution)) {
  stop("Quick Start parallel execution receipt is incomplete.", call. = FALSE)
}
runtime_platform <- lisaR:::lisa_runtime_platform()
workers_requested <- as.integer(parallel_execution$workers_requested)
workers_effective <- as.integer(parallel_execution$workers_effective)
backend_effective <- as.character(parallel_execution$backend_effective)
stopifnot(
  all(parallel_execution$platform == runtime_platform),
  all(workers_requested == 4L),
  all(workers_effective >= 1L),
  all(workers_effective <= workers_requested)
)
if (identical(runtime_platform, "linux")) {
  expected_backend <- ifelse(workers_effective > 1L, "multicore", "sequential")
  stopifnot(identical(backend_effective, unname(expected_backend)))
} else {
  stopifnot(
    all(workers_effective == 1L),
    all(backend_effective == "sequential")
  )
}
stopifnot(
  identical(check$gate, "PASS"),
  file.exists(file.path(result$output_dir, "report_index.html")),
  file.exists(file.path(result$output_dir, "report_output_policy.tsv")),
  dir.exists(result$extension_output_dir),
  file.exists(file.path(result$extension_output_dir, "index.html"))
)
cat("LISA_QUICK_START_GATE=PASS\n")
cat("package_version=", as.character(utils::packageVersion("lisaR")), "\n", sep = "")
cat("package_root=", installed_root, "\n", sep = "")
cat(
  "parallel_contract=platform:", runtime_platform,
  ";requested:", paste(sort(unique(workers_requested)), collapse = ","),
  ";effective:", paste(sort(unique(workers_effective)), collapse = ","),
  ";backend:", paste(sort(unique(backend_effective)), collapse = ","),
  "\n", sep = ""
)
cat("project=", project, "\n", sep = "")
cat("output=", result$output_dir, "\n", sep = "")
