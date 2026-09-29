#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# A/B contrast category navigator. Reads a completed category-contrast summary;
# no GSEA/DE work, no new delta or FDR is computed here.
options(stringsAsFactors = FALSE)
parse_args <- function(args) {
  out <- list(project_dir = "", contrast_id = "", scope = "", analysis_a = "", analysis_b = "",
    universe = "GOBP-C2", tier = "core", gsea_padj_cutoff = .25, label = "",
    contrast_summary = "", output_dir = "", lisa_internal_renderer_sha256 = "")
  if (length(args) %% 2L) stop("Arguments must be --flag value pairs.", call. = FALSE)
  seen <- character()
  while (length(args)) {
    if (!startsWith(args[[1L]], "--")) stop("Arguments must be --flag value pairs.", call. = FALSE)
    key <- gsub("-", "_", sub("^--", "", args[[1L]]))
    if (!key %in% names(out) || key %in% seen) stop("Unknown or duplicate argument: ", args[[1L]], call. = FALSE)
    out[[key]] <- args[[2L]]; seen <- c(seen, key); args <- args[-c(1L, 2L)]
  }
  for (key in c("project_dir", "contrast_id", "scope", "analysis_a", "analysis_b"))
    if (!nzchar(out[[key]])) stop("Required argument: --", gsub("_", "-", key), call. = FALSE)
  out$gsea_padj_cutoff <- suppressWarnings(as.numeric(out$gsea_padj_cutoff))
  if (!nzchar(out$label)) out$label <- out$contrast_id
  out
}
cfg <- parse_args(commandArgs(TRUE))
ns <- asNamespace("lisaR")
get("lisa_navigation_scope_component", ns)(cfg$contrast_id, "contrast_id")
get("lisa_navigation_scope_component", ns)(cfg$scope, "scope")
get("lisa_navigation_scope_component", ns)(cfg$universe, "universe")
if (!nzchar(cfg$contrast_summary)) {
  directory <- file.path(cfg$project_dir, "outputs", "category_contrasts", cfg$scope,
    paste0("collection_", cfg$universe), "lisa_tables")
  candidates <- list.files(directory, pattern = "_GSEA_category_contrast[.]tsv$", full.names = TRUE)
  if (length(candidates) != 1L) stop("Expected one canonical GSEA category contrast table.", call. = FALSE)
  cfg$contrast_summary <- candidates[[1L]]
}
if (!nzchar(cfg$output_dir)) cfg$output_dir <- file.path(cfg$project_dir, "report_pages",
  "category_navigation", cfg$scope, cfg$universe)
navigation <- get("build_lisa_contrast_navigation", ns)(cfg$contrast_summary, contrast_id = cfg$contrast_id,
  scope = cfg$scope, analysis_a = cfg$analysis_a, analysis_b = cfg$analysis_b, collection = cfg$universe,
  tier = cfg$tier, gsea_padj_cutoff = cfg$gsea_padj_cutoff, label = cfg$label)
get("render_lisa_contrast_navigation", ns)(navigation, cfg$output_dir)
cat("Contrast navigator completed: ", nrow(navigation$categories), " exact category links.\n", sep = "")
