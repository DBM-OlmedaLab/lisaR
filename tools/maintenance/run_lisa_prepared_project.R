#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) if (is.null(x)) y else x

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  run_lisa_prepared_project.R --package-dir DIR --project-dir DIR",
      "    --dictionary FILE --category-map FILE --term2gene FILE",
      "    --lisa-project-root DIR",
      "    --report-title TITLE [--collections A,B,C]",
      "",
      "Development-only mutable runner for a reviewed lisaR source checkout.",
      "It writes into an existing prepared project and is not the transactional",
      "production entry point. Use run_lisa() for immutable production runs.",
      "The prepared project directory must contain:",
      "  config/de_index.tsv",
      "  config/contrast_index.tsv",
      "  input/*.tsv",
      "  input_matrices_vst/*.tsv",
      "",
      "Absolute de_index paths are relocated to the project directory when the",
      "referenced basename exists under input/ or input_matrices_vst/.",
      "--package-dir is an explicit trust boundary: it must be the reviewed local",
      "lisaR source tree. The runner rejects missing or unexpected R source files",
      "before sourcing any package code.",
      sep = "\n"
    ),
    "\n",
    file = stderr()
  )
}

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    name <- sub("^--", "", key)
    if (name %in% c("help", "dry-run")) {
      out[[name]] <- TRUE
      i <- i + 1
    } else {
      if (i == length(args)) stop("Missing value for ", key, call. = FALSE)
      out[[name]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

opt <- parse_args(args)
if (isTRUE(opt$help)) {
  usage()
  quit(status = 0)
}

required <- c(
  "package-dir", "project-dir", "dictionary", "category-map", "term2gene",
  "lisa-project-root", "report-title"
)
missing <- required[!vapply(required, function(x) nzchar(opt[[x]] %||% ""), logical(1))]
if (length(missing)) {
  usage()
  stop("Missing required argument(s): ", paste(missing, collapse = ", "), call. = FALSE)
}

raw_package_dir <- path.expand(opt[["package-dir"]])
raw_dictionary_path <- path.expand(opt[["dictionary"]])
raw_category_map_path <- path.expand(opt[["category-map"]])
raw_term2gene <- path.expand(opt[["term2gene"]])

runner_path_key <- function(path) {
  key <- gsub("\\", "/", as.character(path), fixed = TRUE)
  key <- sub("/+$", "", key)
  if (identical(.Platform$OS.type, "windows")) tolower(key) else key
}

runner_path_has_indirection <- function(path) {
  linked <- tryCatch(fs::link_exists(path), error = function(error) NA)
  if (length(linked) != 1L || is.na(linked) || isTRUE(linked)) return(TRUE)
  if (!file.exists(path) && !dir.exists(path)) return(TRUE)
  parent <- normalizePath(dirname(path), winslash = "/", mustWork = TRUE)
  leaf <- basename(path)
  lexical <- if (!nzchar(leaf) || identical(leaf, ".")) {
    parent
  } else if (identical(leaf, "..")) {
    dirname(parent)
  } else {
    file.path(parent, leaf)
  }
  resolved <- normalizePath(path, winslash = "/", mustWork = TRUE)
  !identical(runner_path_key(lexical), runner_path_key(resolved))
}

for (raw_path in c(
    raw_package_dir, raw_dictionary_path, raw_category_map_path,
    raw_term2gene
)) {
  if (runner_path_has_indirection(raw_path)) {
    stop(
      "The development runner rejects symbolic-link code and resource arguments: ",
      raw_path, call. = FALSE
    )
  }
}

package_dir <- normalizePath(raw_package_dir, winslash = "/", mustWork = TRUE)
project_dir <- normalizePath(opt[["project-dir"]], winslash = "/", mustWork = TRUE)
dictionary_path <- normalizePath(raw_dictionary_path, winslash = "/", mustWork = TRUE)
category_map_path <- normalizePath(raw_category_map_path, winslash = "/", mustWork = TRUE)
dictionary_dir <- dirname(dictionary_path)
term2gene <- normalizePath(raw_term2gene, winslash = "/", mustWork = TRUE)
lisa_project_root <- normalizePath(
  opt[["lisa-project-root"]], winslash = "/", mustWork = TRUE
)
report_title <- opt[["report-title"]]
collections <- strsplit(opt[["collections"]] %||% "GOBP-C2,GOMF,GOCC,PATHWAYS,HALLMARKS", ",", fixed = TRUE)[[1]]
collections <- trimws(collections[nzchar(trimws(collections))])

for (subdir in c("config", "input", "input_matrices_vst", "manifests", "logs")) {
  dir.create(file.path(project_dir, subdir), recursive = TRUE, showWarnings = FALSE)
}

description_path <- file.path(package_dir, "DESCRIPTION")
description <- tryCatch(read.dcf(description_path), error = function(error) NULL)
if (is.null(description) || nrow(description) != 1L ||
    !identical(unname(description[[1L, "Package"]]), "lisaR")) {
  stop("--package-dir must identify a reviewed lisaR source tree.", call. = FALSE)
}

# This inventory must equal R/ exactly: the check below refuses to run if any
# file is missing OR unexpected, so a new R/ file has to be reviewed and listed
# here. This explicit allowlist includes the current input adapters, portable
# path checks and explorer. Do not replace it with automatic directory discovery.
source_files <- c(
  "category_inference.R", "category_inference_pipeline.R", "support_grades.R", "hommel_support.R", "support_grades_pipeline.R",
  "calibration.R", "category_member_evidence.R", "category_evidence.R", "category_figure_links.R", "category_navigation.R", "category_plot_variants.R", "contrast_evidence.R", "contrast_navigation.R",
    "report_shell.R", "presentation_integration.R", "gene_evidence.R", "content_cache.R", "de_input_preparation.R", "collections.R", "config_pipeline.R", "contracts.R",
  "dependencies.R", "example_distribution.R", "example_project.R", "example_paths.R",
  "explore_app.R", "explore_contrast.R", "explore_engine.R", "explore_export.R",
  "explore_export_jobs.R", "explore_kegg.R", "explore_presentation.R",
  "explore_request.R", "explore_worker.R", "explore_workspace.R",
  "portable_paths.R", "figure_source_io.R",
  "full_product_completion.R", "full_product_inventory.R", "gene_level.R", "globals.R",
  "html_report.R", "input_contracts.R", "lisa_pipeline.R", "lollipop_reproduction.R", "manifest.R",
  "pathway_contracts.R", "post_lisa_products.R", "post_script_interfaces.R",
  "public_api.R", "report_architecture.R", "report_foundation.R",
  "report_packaging.R", "resource_install.R", "resources.R", "rng_contract.R", "run_LISA_DE.R",
  "run_defaults.R",
  "run_management.R", "scientific_report.R", "security.R", "trusted_rds.R",
  "utils_io.R", "validate_inputs.R"
)
source_dir <- file.path(package_dir, "R")
if (runner_path_has_indirection(source_dir)) {
  stop("The reviewed lisaR R/ directory must not be a symbolic link.",
       call. = FALSE)
}
source_dir <- normalizePath(source_dir, winslash = "/", mustWork = TRUE)
if (!startsWith(source_dir, paste0(package_dir, "/"))) {
  stop("The reviewed lisaR R/ directory resolves outside --package-dir.",
       call. = FALSE)
}
observed_source_files <- sort(list.files(
  source_dir, pattern = "[.]R$", full.names = FALSE
))
missing_source_files <- setdiff(source_files, observed_source_files)
unexpected_source_files <- setdiff(observed_source_files, source_files)
if (length(missing_source_files) || length(unexpected_source_files)) {
  stop(
    "The reviewed lisaR source inventory does not match --package-dir. Missing: ",
    paste(missing_source_files, collapse = ", "),
    "; unexpected: ", paste(unexpected_source_files, collapse = ", "),
    call. = FALSE
  )
}
source_paths <- file.path(source_dir, source_files)
source_indirection <- vapply(
  source_paths, runner_path_has_indirection, logical(1)
)
if (any(source_indirection)) {
  stop(
    "The reviewed lisaR source inventory contains a symbolic link: ",
    paste(basename(source_paths[source_indirection]),
          collapse = ", "),
    call. = FALSE
  )
}
resolved_source_paths <- normalizePath(
  source_paths, winslash = "/", mustWork = TRUE
)
if (any(dirname(resolved_source_paths) != source_dir) ||
    any(!vapply(
      resolved_source_paths,
      function(path) isTRUE(utils::file_test("-f", path)),
      logical(1)
    ))) {
  stop("The reviewed lisaR source inventory contains a non-regular or escaping file.",
       call. = FALSE)
}
runner_function_names <- ls(.GlobalEnv, all.names = TRUE)
runner_function_names <- runner_function_names[vapply(
  runner_function_names,
  function(name) is.function(get(name, envir = .GlobalEnv, inherits = FALSE)),
  logical(1)
)]
runner_function_snapshot <- setNames(lapply(
  runner_function_names,
  function(name) get(name, envir = .GlobalEnv, inherits = FALSE)
), runner_function_names)
for (source_file in file.path(source_dir, source_files)) {
  sys.source(source_file, envir = .GlobalEnv, keep.source = FALSE)
}
# Child Rscript processes resolve lisaR::: from the installed namespace. Keep
# an internal inventory of every function introduced or replaced by the
# reviewed source tree so the subprocess boundary can reject a different
# installed implementation before executing it.
source_runtime_names <- ls(.GlobalEnv, all.names = TRUE)
source_runtime_names <- source_runtime_names[vapply(
  source_runtime_names,
  function(name) is.function(get(name, envir = .GlobalEnv, inherits = FALSE)),
  logical(1)
)]
source_runtime_names <- source_runtime_names[vapply(
  source_runtime_names,
  function(name) {
    if (!name %in% names(runner_function_snapshot)) return(TRUE)
    before <- runner_function_snapshot[[name]]
    after <- get(name, envir = .GlobalEnv, inherits = FALSE)
    !identical(formals(before), formals(after)) ||
      !identical(body(before), body(after), ignore.srcref = TRUE)
  },
  logical(1)
)]
options(lisaR.source_runtime_functions = sort(source_runtime_names))

for (resource in c(dictionary_path, category_map_path, term2gene)) {
  if (runner_path_has_indirection(resource) ||
      !file.exists(resource) || dir.exists(resource) ||
      !isTRUE(utils::file_test("-f", resource))) {
    stop(
      "The development runner requires regular, non-symlink resource files: ",
      resource, call. = FALSE
    )
  }
}

relocate_project_path <- function(path, target_subdir) {
  if (is.na(path) || !nzchar(path)) return(path)
  if (file.exists(path)) return(normalizePath(path, mustWork = TRUE))
  candidate <- file.path(project_dir, target_subdir, basename(path))
  if (file.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))
  path
}

de_index_path <- file.path(project_dir, "config", "de_index.tsv")
contrast_index_path <- file.path(project_dir, "config", "contrast_index.tsv")
if (!file.exists(de_index_path)) stop("Missing de_index: ", de_index_path, call. = FALSE)
if (!file.exists(contrast_index_path)) stop("Missing contrast_index: ", contrast_index_path, call. = FALSE)

de_index <- read_lisa_tsv(de_index_path)
if (!all(c("analysis_id", "de_path", "counts_matrix_path") %in% names(de_index))) {
  stop("de_index lacks required columns: analysis_id, de_path, counts_matrix_path", call. = FALSE)
}
de_index$de_path <- vapply(de_index$de_path, relocate_project_path, character(1), target_subdir = "input")
de_index$counts_matrix_path <- vapply(
  de_index$counts_matrix_path,
  relocate_project_path,
  character(1),
  target_subdir = "input_matrices_vst"
)

unresolved_de <- de_index$de_path[!file.exists(de_index$de_path)]
unresolved_mtx <- de_index$counts_matrix_path[!file.exists(de_index$counts_matrix_path)]
if (length(unresolved_de) || length(unresolved_mtx)) {
  stop(
    "Unresolved input paths after relocation.\nDE: ",
    paste(unresolved_de, collapse = ", "),
    "\nMatrices: ",
    paste(unresolved_mtx, collapse = ", "),
    call. = FALSE
  )
}

resolved_de_index <- file.path(project_dir, "config", "de_index.resolved.tsv")
write_lisa_tsv(de_index, resolved_de_index)

write_lisa_tsv(
  data.frame(
    key = c(
      "package_dir", "project_dir", "dictionary", "dictionary_sha256",
      "category_map", "category_map_sha256", "term2gene",
      "term2gene_sha256", "lisa_project_root", "report_title", "collections"
    ),
    value = c(
      package_dir, project_dir, dictionary_path,
      lisa_sha256_file(dictionary_path), category_map_path,
      lisa_sha256_file(category_map_path), term2gene,
      lisa_sha256_file(term2gene), lisa_project_root, report_title,
      paste(collections, collapse = ",")
    ),
    stringsAsFactors = FALSE
  ),
  file.path(project_dir, "manifests", "development_runner_inputs.tsv")
)

run_lisa_pipeline(
  de_index = resolved_de_index,
  contrast_index = contrast_index_path,
  dictionary_dir = dictionary_dir,
  term2gene = term2gene,
  dictionary_path = dictionary_path,
  category_map_path = category_map_path,
  output_dir = project_dir,
  collections = collections,
  base_dir = project_dir,
  lisa_project_root = lisa_project_root,
  lisa_dictionary = "core",
  file_label_prefix = "semantic",
  plot_formats = c("png"),
  export_formats = c("tsv", "xlsx"),
  # This mutable development runner executes the same bounded canonical layer
  # as run_lisa(). Gene-level products belong to the post-run extension APIs.
  run_gene_level = FALSE,
  run_reports = TRUE,
  run_kegg_maps = TRUE,
  kegg_cache_root = Sys.getenv("LISAR_KEGG_CACHE", ""),
  kegg_snapshot_id = Sys.getenv("LISAR_KEGG_SNAPSHOT", ""),
  report_title = report_title,
  dry_run = isTRUE(opt[["dry-run"]]),
  .test_package_dir = package_dir
)
