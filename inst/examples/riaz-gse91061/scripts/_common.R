# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Shared helpers for the Riaz GSE91061 worked example.
#
# The scripts in this directory are designed to run both from a source checkout
# and from a copied example directory. All mutable files are written below the
# directory declared by LISAR_RIAZ_PROJECT_DIR. Original files are read from
# LISAR_RIAZ_INPUT_CACHE. Keeping these two roots separate makes it impossible for an
# analysis run to modify the immutable input cache by accident.

riaz_script_path <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  flag <- args[grepl("^--file=", args)]
  if (length(flag) != 1L) {
    stop("Unable to resolve the current script path from --file=.")
  }
  normalizePath(sub("^--file=", "", flag[[1]]), mustWork = TRUE)
}

riaz_example_dir <- function() {
  normalizePath(dirname(riaz_script_path()), mustWork = TRUE)
}

riaz_project_dir <- function() {
  value <- Sys.getenv("LISAR_RIAZ_PROJECT_DIR", unset = "")
  if (!nzchar(value)) {
    stop(
      "LISAR_RIAZ_PROJECT_DIR is required and must point to a writable project root."
    )
  }
  dir.create(value, recursive = TRUE, showWarnings = FALSE)
  normalizePath(value, mustWork = TRUE)
}

riaz_input_cache <- function() {
  value <- Sys.getenv("LISAR_RIAZ_INPUT_CACHE", unset = "")
  if (!nzchar(value)) {
    value <- file.path(riaz_project_dir(), "inputs", "original")
  }
  dir.create(value, recursive = TRUE, showWarnings = FALSE)
  normalizePath(value, mustWork = TRUE)
}

# Prepared Riaz projects use the ordinary C1 resource resolver: installed
# lisaR supplies the core/category-map bytes and TERM2GENE belongs in its
# regular external cache, never under the study project.
riaz_normal_resource_cache <- function() {
  riaz_require_packages("lisaR")
  value <- Sys.getenv("LISAR_RIAZ_RESOURCE_CACHE", unset = "")
  if (!nzchar(value)) {
    value <- get(
      "lisa_dictionary_cache_root", envir = asNamespace("lisaR"),
      inherits = FALSE
    )()
  }
  value <- path.expand(value)
  project <- riaz_project_dir()
  if (dir.exists(value)) {
    riaz_assert(!riaz_path_is_link(value),
                "LISAR_RIAZ_RESOURCE_CACHE cannot be a symbolic link.")
    value <- normalizePath(value, winslash = "/", mustWork = TRUE)
  }
  riaz_assert(
    !riaz_path_is_within(value, project),
    paste(
      "LISAR_RIAZ_RESOURCE_CACHE must be outside LISAR_RIAZ_PROJECT_DIR;",
      "prepared projects must not contain a scientific-resource cache."
    )
  )
  value
}

riaz_configure_normal_resources <- function() {
  cache_root <- riaz_normal_resource_cache()
  options(
    lisaR.dictionary_cache_root = cache_root,
    lisaR.dictionary_registry = NULL,
    lisaR.shared_dictionary_root = NULL
  )
  cache_root
}

riaz_output_dir <- function(...) {
  path <- file.path(riaz_project_dir(), "outputs", ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

riaz_write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data,
    path,
    sep = "\t",
    quote = FALSE,
    na = "NA"
  )
  invisible(path)
}

riaz_read_tsv <- function(path) {
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    connection <- gzfile(path, open = "rt")
    on.exit(close(connection), add = TRUE)
    return(utils::read.delim(
      connection,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      na.strings = c("", "NA")
    ))
  }
  data.table::fread(path, sep = "\t", data.table = FALSE, na.strings = c("", "NA"))
}

riaz_sha256 <- function(path) {
  if (!requireNamespace("lisaR", quietly = TRUE)) {
    stop(
      "LISA-SHA256-003 the installed lisaR namespace is required for SHA-256.",
      call. = FALSE
    )
  }
  get("lisa_sha256_file", envir = asNamespace("lisaR"), inherits = FALSE)(path)
}

riaz_require_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "Required R packages are unavailable: ",
      paste(missing, collapse = ", ")
    )
  }
}

riaz_assert <- function(condition, message) {
  if (length(condition) != 1L || is.na(condition) || !condition) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

riaz_path_is_link <- function(path) {
  target <- Sys.readlink(path)
  length(target) == 1L && !is.na(target) && nzchar(target)
}

riaz_path_entry_exists <- function(path) {
  file.exists(path) || dir.exists(path) || riaz_path_is_link(path)
}

# `normalizePath(mustWork = FALSE)` returns a path that does not exist yet
# unchanged, so a planned output keeps the caller spelling while an existing
# root is canonicalized. macOS reaches one directory as both `/var/...` and
# `/private/var/...`, and Windows as both an 8.3 short name and the long name,
# so the guard would read two spellings of the same directory as an escape.
# Resolve the deepest existing ancestor of each side and re-append the rest.
riaz_path_resolve <- function(value) {
  value <- sub("/+$", "", gsub("\\\\", "/", path.expand(as.character(value))))
  parent <- value
  remainder <- character()
  while (!dir.exists(parent) && !identical(parent, dirname(parent))) {
    remainder <- c(basename(parent), remainder)
    parent <- dirname(parent)
  }
  parent <- sub("/+$", "", gsub("\\\\", "/",
    normalizePath(parent, winslash = "/", mustWork = FALSE)))
  if (!nzchar(parent)) parent <- "/"
  if (!length(remainder)) return(parent)
  paste0(parent, if (identical(parent, "/")) "" else "/",
         paste(remainder, collapse = "/"))
}

riaz_path_is_within <- function(path, root) {
  path <- riaz_path_resolve(path)
  root <- riaz_path_resolve(root)
  identical(path, root) || startsWith(path, paste0(root, "/"))
}

riaz_new_output_path <- function(path, label = "output path") {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || grepl("[\r\n]", path)) {
    stop(label, " must be one non-empty path.", call. = FALSE)
  }
  expanded <- path.expand(path)
  parent <- dirname(expanded)
  reject_missing_parent <- function() {
    stop(
      label, " must have an existing parent directory: ", parent,
      call. = FALSE
    )
  }
  if (!dir.exists(parent)) reject_missing_parent()
  # Windows can report `dir.exists("missing/..")` as true even though the
  # missing component was never traversable. For a lexical traversal, inspect
  # ancestors back through the component immediately before the final `..`.
  # Ordinary paths (including UNC parents) retain the prior single check.
  pieces <- strsplit(gsub("\\\\", "/", parent), "/", fixed = TRUE)[[1L]]
  if (any(pieces == "..")) {
    ancestor <- parent
    repeat {
      if (!dir.exists(ancestor)) reject_missing_parent()
      ancestor_pieces <- strsplit(
        gsub("\\\\", "/", ancestor), "/", fixed = TRUE
      )[[1L]]
      if (!any(ancestor_pieces == "..")) break
      next_ancestor <- dirname(ancestor)
      if (identical(next_ancestor, ancestor)) break
      ancestor <- next_ancestor
    }
  }
  leaf <- basename(expanded)
  if (!nzchar(leaf) || leaf %in% c(".", "..")) {
    stop(label, " must name a new child directory.", call. = FALSE)
  }
  candidate <- file.path(
    normalizePath(parent, winslash = "/", mustWork = TRUE),
    leaf
  )
  if (riaz_path_entry_exists(candidate)) {
    stop(label, " must not already exist, including as a symbolic link.",
         call. = FALSE)
  }
  candidate
}

riaz_safe_relative_path <- function(path, label = "relative path") {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || grepl("[\r\n]", path)) {
    stop(label, " must be one non-empty relative path.", call. = FALSE)
  }
  value <- gsub("\\\\", "/", path)
  components <- strsplit(value, "/", fixed = TRUE)[[1]]
  if (startsWith(value, "/") || grepl("^[A-Za-z]:/", value) ||
      any(!nzchar(components) | components %in% c(".", ".."))) {
    stop(label, " contains an unsafe or non-relative component.", call. = FALSE)
  }
  value
}

riaz_effect_concordance <- function(current, benchmark) {
  required <- c("symbol", "log2FoldChange", "stat")
  riaz_assert(
    all(required %in% names(current)) && all(required %in% names(benchmark)),
    "Concordance inputs must contain symbol, log2FoldChange, and stat."
  )
  joined <- merge(
    current[, required, drop = FALSE],
    benchmark[, required, drop = FALSE],
    by = "symbol",
    suffixes = c("_current", "_benchmark")
  )
  metric <- function(current_values, benchmark_values) {
    keep <- is.finite(current_values) & is.finite(benchmark_values)
    riaz_assert(sum(keep) >= 3L, "Too few finite values for concordance.")
    current_values <- current_values[keep]
    benchmark_values <- benchmark_values[keep]
    data.frame(
      common_finite_values = length(current_values),
      spearman = unname(stats::cor(
        current_values, benchmark_values, method = "spearman"
      )),
      directional_concordance = mean(
        sign(current_values) == sign(benchmark_values)
      ),
      stringsAsFactors = FALSE
    )
  }
  lfc <- metric(
    joined$log2FoldChange_current,
    joined$log2FoldChange_benchmark
  )
  stat <- metric(joined$stat_current, joined$stat_benchmark)
  data.frame(
    common_symbols = nrow(joined),
    lfc_common_finite = lfc$common_finite_values,
    lfc_spearman = lfc$spearman,
    lfc_directional_concordance = lfc$directional_concordance,
    stat_common_finite = stat$common_finite_values,
    stat_spearman = stat$spearman,
    stat_directional_concordance = stat$directional_concordance,
    stringsAsFactors = FALSE
  )
}

riaz_write_package_manifest <- function(path) {
  installed <- as.data.frame(
    utils::installed.packages(noCache = TRUE),
    stringsAsFactors = FALSE
  )
  keep <- c("Package", "Version", "LibPath", "Priority")
  installed <- installed[, keep, drop = FALSE]
  names(installed) <- tolower(names(installed))
  installed <- installed[order(installed$package, installed$libpath), ]
  riaz_write_tsv(installed, path)
}

riaz_rscript_executable <- function() {
  suffix <- if (.Platform$OS.type == "windows") ".exe" else ""
  path <- file.path(R.home("bin"), paste0("Rscript", suffix))
  if (!file.exists(path) || dir.exists(path)) {
    stop("The Rscript executable is unavailable.", call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

riaz_run_script <- function(name) {
  script <- file.path(riaz_example_dir(), name)
  status <- system2(
    riaz_rscript_executable(),
    shQuote(script),
    stdout = "",
    stderr = ""
  )
  if (!identical(status, 0L)) {
    stop("Worked-example stage failed: ", name, " (exit ", status, ").")
  }
  invisible(TRUE)
}
