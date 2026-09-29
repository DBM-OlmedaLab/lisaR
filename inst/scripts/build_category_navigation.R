#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Single-analysis category navigator. Reads a completed summary; no GSEA/DE work.
options(stringsAsFactors = FALSE)
parse_args <- function(args) {
  out <- list(project_dir = "", analysis_id = "", universe = "GOBP-C2", tier = "core",
    gsea_padj_cutoff = .25, positive_contrast = "not recorded", output_dir = "",
    lisa_internal_renderer_sha256 = "")
  if (length(args) %% 2L) stop("Arguments must be --flag value pairs.", call. = FALSE)
  seen <- character()
  while (length(args)) {
    if (!startsWith(args[[1L]], "--")) stop("Arguments must be --flag value pairs.", call. = FALSE)
    key <- gsub("-", "_", sub("^--", "", args[[1L]]))
    if (!key %in% names(out) || key %in% seen) stop("Unknown or duplicate argument: ", args[[1L]], call. = FALSE)
    out[[key]] <- args[[2L]]; seen <- c(seen, key); args <- args[-c(1L, 2L)]
  }
  for (key in c("project_dir", "analysis_id")) if (!nzchar(out[[key]])) stop("Required argument: --", gsub("_", "-", key), call. = FALSE)
  out$gsea_padj_cutoff <- suppressWarnings(as.numeric(out$gsea_padj_cutoff))
  out
}
cfg <- parse_args(commandArgs(TRUE))
ns <- asNamespace("lisaR")
get("lisa_navigation_scope_component", ns)(cfg$analysis_id, "analysis_id")
get("lisa_navigation_scope_component", ns)(cfg$universe, "universe")
collection_dir <- file.path(cfg$project_dir, "outputs", "single_de", cfg$analysis_id, paste0("collection_", cfg$universe))
summaries <- list.files(file.path(collection_dir, "lisa_tables"), pattern = "_GSEA_category_summary[.]tsv$", full.names = TRUE)
summaries <- summaries[startsWith(basename(summaries), paste0(cfg$analysis_id, "_"))]
if (length(summaries) != 1L) stop("Expected exactly one existing GSEA category summary; found ", length(summaries), ".", call. = FALSE)
if (!nzchar(cfg$output_dir)) cfg$output_dir <- file.path(cfg$project_dir, "report_pages", "category_navigation", cfg$analysis_id, cfg$universe)
navigation <- get("build_lisa_category_navigation", ns)(summaries[[1L]], analysis_id = cfg$analysis_id,
  collection = cfg$universe, tier = cfg$tier, gsea_padj_cutoff = cfg$gsea_padj_cutoff,
  positive_contrast = cfg$positive_contrast)
support <- get("lisa_hommel_read_results", ns)(collection_dir, cfg$analysis_id, cfg$universe)
if (!is.null(support)) navigation <- get("lisa_hommel_attach", ns)(navigation, support)
get("render_lisa_category_navigation", ns)(navigation, cfg$output_dir)
cat("Category navigator completed: ", nrow(navigation$categories), " exact category links.\n", sep = "")
