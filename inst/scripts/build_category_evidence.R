#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Rebuild descriptive evidence from an existing completed collection without
# rerunning enrichment or requiring legacy gene-prioritisation products.
options(stringsAsFactors = FALSE)
parse_args <- function(args) {
  out <- list(project_dir = "", analysis_id = "", universe = "GOBP-C2", tier = "not recorded",
    gsea_padj_cutoff = .25, de_padj_cutoff = .05, positive_contrast = "not recorded",
    max_sets = 25L, max_genes = 40L, formats = "png", output_dir = "",
    category_dictionary = "", categories = "", source_data = "true", recipes = "true",
    overlap_export = "displayed",
    lisa_internal_renderer_sha256 = "")
  if (length(args) %% 2L) stop("Arguments must be --flag value pairs.", call. = FALSE)
  used <- character()
  while (length(args)) {
    if (!startsWith(args[[1L]], "--")) stop("Arguments must be --flag value pairs.", call. = FALSE)
    key <- gsub("-", "_", sub("^--", "", args[[1L]]))
    if (!key %in% names(out)) stop("Unknown argument: ", args[[1L]], call. = FALSE)
    if (key %in% used) stop("Duplicate argument: ", args[[1L]], call. = FALSE)
    used <- c(used, key); out[[key]] <- args[[2L]]; args <- args[-c(1L, 2L)]
  }
  for (key in c("project_dir", "analysis_id")) if (!nzchar(out[[key]])) stop("Required argument: --", gsub("_", "-", key), call. = FALSE)
  for (key in c("max_sets", "max_genes", "gsea_padj_cutoff", "de_padj_cutoff")) out[[key]] <- suppressWarnings(as.numeric(out[[key]]))
  for (key in c("source_data", "recipes")) {
    if (!tolower(out[[key]]) %in% c("true", "false")) stop("--", gsub("_", "-", key), " must be true or false.", call. = FALSE)
    out[[key]] <- identical(tolower(out[[key]]), "true")
  }
  out$formats <- unique(trimws(strsplit(tolower(out$formats), ",", fixed = TRUE)[[1L]]))
  out$formats <- out$formats[nzchar(out$formats)]
  if (any(!out$formats %in% c("png", "pdf", "svg"))) stop("--formats must contain png, pdf and/or svg (or be empty).", call. = FALSE)
  if (!out$overlap_export %in% c("displayed", "all", "none")) stop("--overlap-export must be displayed, all or none.", call. = FALSE)
  out
}

cfg <- parse_args(commandArgs(TRUE))
# The orchestrator supplies its candidate library to the child. No source-tree
# fallback may silently mix versions of the renderer and installed namespace.
ns <- asNamespace("lisaR")
if (!exists("build_lisa_category_evidence", envir = ns, inherits = FALSE))
  stop("Installed lisaR does not provide category evidence; install the candidate package first.", call. = FALSE)
collection_dir <- file.path(cfg$project_dir, "outputs", "single_de", cfg$analysis_id, paste0("collection_", cfg$universe))
gsea_paths <- list.files(file.path(collection_dir, "enrichment"), pattern = "_GSEA_.*_annotated[.]tsv$", full.names = TRUE)
gsea_paths <- gsea_paths[startsWith(basename(gsea_paths), paste0(cfg$analysis_id, "_GSEA_"))]
if (length(gsea_paths) != 1L) stop("Expected exactly one full annotated GSEA table in the completed collection; found ", length(gsea_paths), ".", call. = FALSE)
de_path <- file.path(collection_dir, "inputs", paste0(cfg$analysis_id, "_standardized_DE.tsv"))
ledger_path <- file.path(collection_dir, "qc", paste0(cfg$analysis_id, "_GSEA_universe_ledger.tsv"))
if (!nzchar(cfg$output_dir)) cfg$output_dir <- file.path(cfg$project_dir, "outputs", "gene_level", "single_de", cfg$analysis_id, paste0("collection_", cfg$universe), "category_evidence")
if (identical(cfg$positive_contrast, "not recorded")) {
  path <- file.path(cfg$project_dir, "config", "de_index.tsv")
  if (file.exists(path)) {
    index <- get("lisa_evidence_read", ns)(path)
    if (all(c("analysis_id", "positive_direction") %in% names(index))) {
      value <- as.character(index$positive_direction[index$analysis_id == cfg$analysis_id])
      if (length(value) == 1L && !is.na(value) && nzchar(value)) cfg$positive_contrast <- value
    }
  }
}
dictionary <- if (nzchar(cfg$category_dictionary)) get("lisa_evidence_read", ns)(cfg$category_dictionary) else NULL
if (is.null(dictionary)) {
  summaries <- list.files(file.path(collection_dir, "lisa_tables"), pattern = "_GSEA_category_summary[.]tsv$", full.names = TRUE)
  summaries <- summaries[startsWith(basename(summaries), paste0(cfg$analysis_id, "_"))]
  if (length(summaries) == 1L) dictionary <- get("lisa_evidence_read", ns)(summaries[[1L]])
}
if (!is.null(dictionary) && "universe" %in% names(dictionary)) dictionary <- dictionary[dictionary$universe %in% cfg$universe, , drop = FALSE]
evidence <- get("build_lisa_category_evidence", ns)(gsea_table = gsea_paths[[1L]],
  de_table = if (file.exists(de_path)) de_path else NULL,
  universe_ledger = if (file.exists(ledger_path)) ledger_path else NULL,
  category_dictionary = dictionary, analysis_id = cfg$analysis_id, collection = cfg$universe,
  tier = cfg$tier, positive_contrast = cfg$positive_contrast, gsea_padj_cutoff = cfg$gsea_padj_cutoff,
  de_padj_cutoff = cfg$de_padj_cutoff, max_sets = cfg$max_sets, max_genes = cfg$max_genes)
support <- get("lisa_hommel_read_results", ns)(collection_dir, cfg$analysis_id, cfg$universe)
if (!is.null(support)) evidence <- get("lisa_hommel_attach", ns)(evidence, support)
get("render_lisa_category_evidence", ns)(evidence, cfg$output_dir, formats = cfg$formats,
  source_data = cfg$source_data, recipes = cfg$recipes, overlap_export = cfg$overlap_export,
  categories = if (nzchar(cfg$categories)) trimws(strsplit(cfg$categories, ",", fixed = TRUE)[[1L]]) else NULL)
cat("Category evidence completed: ", nrow(evidence$categories), " categories; ", nrow(evidence$sets),
  " mapped set assignments; ", nrow(evidence$leading_edges), " leading-edge assignments.\n", sep = "")
