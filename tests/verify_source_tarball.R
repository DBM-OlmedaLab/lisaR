#!/usr/bin/env Rscript
# Initial source-tarball boundary test.
# Usage: Rscript --vanilla tests/verify_source_tarball.R [package-source-root]

args <- commandArgs(trailingOnly = TRUE)
source_root <- if (length(args) == 0L) "." else args[[1L]]
if (length(args) == 0L && !file.exists(file.path(source_root, "DESCRIPTION"))) {
  message("SKIP: R CMD check is already running from the built source tarball")
  quit(save = "no", status = 0L)
}
source_root <- normalizePath(source_root, mustWork = TRUE)
work_dir <- tempfile("lisar-source-tarball-")
dir.create(work_dir, recursive = TRUE)
on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

r_command <- file.path(R.home("bin"), "R")
previous_working_directory <- getwd()
on.exit(setwd(previous_working_directory), add = TRUE)
setwd(work_dir)
build_output <- system2(
  r_command,
  c("--no-environ", "--no-site-file", "--no-init-file", "CMD", "build", "--no-build-vignettes", source_root),
  stdout = TRUE,
  stderr = TRUE,
  env = c("R_PROFILE_USER=/dev/null", "R_ENVIRON_USER=/dev/null", "R_TESTS=")
)
build_status <- attr(build_output, "status")
if (is.null(build_status)) {
  build_status <- 0L
}
if (!identical(as.integer(build_status), 0L)) {
  stop(
    paste(c("R CMD build failed in source-tarball boundary test:", build_output), collapse = "\n"),
    call. = FALSE
  )
}

tarballs <- list.files(work_dir, pattern = "^lisaR_.*\\.tar\\.gz$", full.names = TRUE)
if (length(tarballs) != 1L) {
  stop("source-tarball boundary test did not produce exactly one tarball", call. = FALSE)
}

contents <- utils::untar(tarballs[[1L]], list = TRUE)
forbidden <- c(
  "runs", "staging", "input_staging", "exports", "backups", "audits",
  "hpc_runs", "smoke-output",
  "state_last_[^/]*"
)
pattern <- paste0("/(?:", paste(forbidden, collapse = "|"), ")(?:/|$)")
if (any(grepl(pattern, contents, perl = TRUE))) {
  stop("source tarball contains an operational runtime path", call. = FALSE)
}

extract_dir <- file.path(work_dir, "extracted")
dir.create(extract_dir)
utils::untar(tarballs[[1L]], exdir = extract_dir)
text_extensions <- c(
  "R", "Rd", "Rmd", "md", "txt", "tsv", "csv", "json", "yaml", "yml",
  "DESCRIPTION", "NAMESPACE", "LICENSE"
)
files <- list.files(extract_dir, recursive = TRUE, full.names = TRUE,
                    all.files = TRUE, no.. = TRUE)
files <- files[file.info(files)$isdir %in% FALSE]
is_text <- basename(files) %in% text_extensions |
  tools::file_ext(files) %in% text_extensions
files <- files[is_text]

private_path_patterns <- c(
  "(?:/home|/Users)/[[:alnum:]_.-]+/",
  "[A-Za-z]:[/\\\\]Users[/\\\\][[:alnum:]_. -]+[/\\\\]"
)

hits <- character()
for (path in files) {
  text <- paste(suppressWarnings(readLines(path, warn = FALSE)), collapse = "\n")
  text <- iconv(text, from = "", to = "UTF-8", sub = "byte")
  # Portability tests deliberately contain one synthetic home directory.
  # The exemption is restricted to test files and this explicit fixture name.
  if (grepl("/lisaR/tests/", path, fixed = TRUE)) {
    text <- gsub("(?:/home|/Users)/fixture-user(?:/|\"|$)", "<fixture-home>/", text, perl = TRUE)
    text <- gsub("[A-Za-z]:[/\\\\]Users[/\\\\]+fixture-user", "<fixture-home>", text, perl = TRUE)
  }
  matched <- private_path_patterns[vapply(
    private_path_patterns, grepl, logical(1), x = text,
    ignore.case = TRUE, perl = TRUE
  )]
  if (length(matched)) {
    hits <- c(hits, sprintf("%s: %s", sub(paste0("^", extract_dir, "/"), "", path),
                            paste(matched, collapse = ", ")))
  }
}
if (length(hits)) {
  stop(
    paste(c("source tarball contains private/internal references:", hits),
          collapse = "\n"),
    call. = FALSE
  )
}

source(file.path(source_root, "tests", "verify_source_tarball_helpers.R"))
retired_layout_hits <- lisa_retired_layout_hits(files)
if (length(retired_layout_hits)) {
  stop(
    paste(
      c(
        "source tarball contains retired layout dependency references:",
        sub(paste0("^", extract_dir, "/"), "", retired_layout_hits)
      ),
      collapse = "\n"
    ),
    call. = FALSE
  )
}

message("PASS: source tarball contains no operational paths or private/internal references")
