#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Connect existing A/B category evidence; this script performs no GSEA or DE.
options(stringsAsFactors = FALSE)
parse_args <- function(args) {
  out <- list(project_dir = "", contrast_id = "", analysis_a = "", analysis_b = "", label = "",
    universe = "GOBP-C2", tier = "core", evidence_a = "", evidence_b = "", contrast_summary = "", output_dir = "",
    max_sets = 25L, max_genes = 40L, formats = "png", variants = "clean,percentages,direction,dispersion",
    report_href = "", gene_explorer_href = "", lisa_internal_renderer_sha256 = "")
  if (length(args) %% 2L) stop("Arguments must be --flag value pairs.", call. = FALSE)
  used <- character()
  while (length(args)) {
    if (!startsWith(args[[1L]], "--")) stop("Arguments must be --flag value pairs.", call. = FALSE)
    key <- gsub("-", "_", sub("^--", "", args[[1L]]))
    if (!key %in% names(out)) stop("Unknown contrast evidence argument: ", args[[1L]], call. = FALSE)
    if (key %in% used) stop("Duplicate contrast evidence argument: ", args[[1L]], call. = FALSE)
    used <- c(used, key); out[[key]] <- args[[2L]]; args <- args[-c(1L, 2L)]
  }
  if (!nzchar(out$project_dir) || !nzchar(out$contrast_id)) stop("Required: --project-dir and --contrast-id.", call. = FALSE)
  for (key in c("max_sets", "max_genes")) out[[key]] <- suppressWarnings(as.numeric(out[[key]]))
  for (key in c("formats", "variants")) out[[key]] <- Filter(nzchar, trimws(strsplit(out[[key]], ",", fixed = TRUE)[[1L]]))
  out
}
cfg <- parse_args(commandArgs(TRUE)); ns <- asNamespace("lisaR")
if (!exists("build_lisa_contrast_evidence", envir = ns, inherits = FALSE)) stop("Installed lisaR lacks contrast evidence.", call. = FALSE)
index_path <- file.path(cfg$project_dir, "config", "contrast_index.tsv")
row <- NULL
if (file.exists(index_path)) {
  index <- get("lisa_evidence_read", ns)(index_path)
  if ("contrast_id" %in% names(index)) {
    match_id <- index$contrast_id == cfg$contrast_id
    if ("output_id" %in% names(index)) match_id <- match_id | paste(index$contrast_id, index$output_id, sep = "_") == cfg$contrast_id
    found <- index[match_id, , drop = FALSE]
    if (nrow(found) > 1L) stop("Contrast ID is ambiguous in contrast_index.tsv.", call. = FALSE)
    if (nrow(found)) row <- found
  }
}
for (side in c("a", "b")) {
  key <- paste0("analysis_", side)
  recorded <- if (is.null(row)) "" else as.character(get("lisa_evidence_column", ns)(row, c(key, paste0("contrast_", side))))[[1L]]
  if (!nzchar(cfg[[key]])) cfg[[key]] <- recorded
  if (nzchar(recorded) && !identical(cfg[[key]], recorded)) stop("Explicit analysis side conflicts with contrast index.", call. = FALSE)
  if (!nzchar(cfg[[key]]) || grepl("[/\\\\]", cfg[[key]])) stop("Exact valid analysis ID required for side ", side, ".", call. = FALSE)
  evidence_key <- paste0("evidence_", side)
  if (!nzchar(cfg[[evidence_key]])) cfg[[evidence_key]] <- file.path(cfg$project_dir, "outputs", "gene_level", "single_de", cfg[[key]], paste0("collection_", cfg$universe), "category_evidence")
}
if (!nzchar(cfg$label)) cfg$label <- if (is.null(row)) cfg$contrast_id else as.character(get("lisa_evidence_column", ns)(row, c("contrast_label", "label"), cfg$contrast_id))[[1L]]
if (!nzchar(cfg$contrast_summary)) {
  name <- if (!is.null(row) && "output_id" %in% names(row)) paste(row$contrast_id, row$output_id, sep = "_") else cfg$contrast_id
  directory <- file.path(cfg$project_dir, "outputs", "category_contrasts", name, paste0("collection_", cfg$universe), "lisa_tables")
  candidates <- list.files(directory, pattern = "_GSEA_category_contrast[.]tsv$", full.names = TRUE)
  if (length(candidates) != 1L) stop("Expected one canonical GSEA category contrast table.", call. = FALSE)
  cfg$contrast_summary <- candidates[[1L]]
}
evidence <- get("build_lisa_contrast_evidence", ns)(cfg$evidence_a, cfg$evidence_b,
  contrast_id = cfg$contrast_id, label = cfg$label, max_sets = cfg$max_sets, max_genes = cfg$max_genes,
  contrast_summary = cfg$contrast_summary)
for (key in c("analysis_a", "analysis_b", "tier")) if (!identical(as.character(evidence$metadata[[key]]), cfg[[key]])) stop("Rendered evidence conflicts with requested ", key, ".", call. = FALSE)
if (!identical(evidence$metadata$collection, cfg$universe)) stop("Rendered evidence collection does not match --universe.", call. = FALSE)
if (nzchar(cfg$gene_explorer_href)) evidence$metadata$gene_explorer_href <- cfg$gene_explorer_href
if (!nzchar(cfg$output_dir)) cfg$output_dir <- file.path(cfg$project_dir, "outputs", "gene_level", "category_contrasts", cfg$contrast_id, paste0("collection_", cfg$universe), "contrast_evidence")
get("render_lisa_contrast_evidence", ns)(evidence, cfg$output_dir, formats = cfg$formats, variants = cfg$variants, report_href = cfg$report_href)
cat("Descriptive contrast evidence completed: ", nrow(evidence$categories), " categories; ", nrow(evidence$sets), " paired set assignments. No enrichment or DE recalculated.\n", sep = "")
