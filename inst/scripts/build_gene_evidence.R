#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Build the offline global gene explorer from completed individual analyses.
# No enrichment, DE, or contrast computation is performed here.
options(stringsAsFactors = FALSE)
parse_args <- function(args) {
  out <- list(project_dir = "", output_dir = "", term2gene = "", tier = "core",
    analysis_ids = "", gsea_padj_cutoff = .25, de_padj_cutoff = .05,
    lisa_internal_renderer_sha256 = "")
  if (length(args) %% 2L) stop("Arguments must be --flag value pairs.", call. = FALSE)
  used <- character()
  while (length(args)) {
    if (!startsWith(args[[1L]], "--")) stop("Arguments must be --flag value pairs.", call. = FALSE)
    key <- gsub("-", "_", sub("^--", "", args[[1L]]))
    if (!key %in% names(out)) stop("Unknown gene evidence argument: ", args[[1L]], call. = FALSE)
    if (key %in% used) stop("Duplicate gene evidence argument: ", args[[1L]], call. = FALSE)
    used <- c(used, key); out[[key]] <- args[[2L]]; args <- args[-c(1L, 2L)]
  }
  for (key in c("project_dir", "term2gene")) if (!nzchar(out[[key]])) stop("Required argument: --", gsub("_", "-", key), call. = FALSE)
  for (key in c("gsea_padj_cutoff", "de_padj_cutoff")) out[[key]] <- suppressWarnings(as.numeric(out[[key]]))
  if (!nzchar(out$output_dir)) out$output_dir <- file.path(out$project_dir, "report_pages", "gene_evidence")
  out
}
cfg <- parse_args(commandArgs(TRUE))
ns <- asNamespace("lisaR")
if (!exists("build_lisa_gene_evidence", ns, inherits = FALSE)) stop("Installed lisaR lacks the global gene evidence builder.", call. = FALSE)
read_table <- get("lisa_gene_evidence_read", ns)
single_root <- file.path(cfg$project_dir, "outputs", "single_de")
analysis_ids <- if (nzchar(cfg$analysis_ids)) unique(trimws(strsplit(cfg$analysis_ids, ",", fixed = TRUE)[[1L]])) else basename(list.dirs(single_root, full.names = TRUE, recursive = FALSE))
if (!length(analysis_ids) || any(!nzchar(analysis_ids)) || any(grepl("[/\\\\]", analysis_ids))) stop("Gene evidence needs valid individual-analysis IDs.", call. = FALSE)
index_path <- file.path(cfg$project_dir, "config", "de_index.tsv")
index <- if (file.exists(index_path)) read_table(index_path) else data.frame()
scopes <- list()
for (analysis_id in analysis_ids) {
  root <- file.path(single_root, analysis_id)
  if (!dir.exists(root)) stop("Completed individual analysis missing: ", analysis_id, call. = FALSE)
  dirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[startsWith(basename(dirs), "collection_") & basename(dirs) != "collection_HALLMARKS"]
  if (!length(dirs)) stop("No completed semantic collections for individual analysis: ", analysis_id, call. = FALSE)
  positive <- "not recorded"
  if (all(c("analysis_id", "positive_direction") %in% names(index))) {
    value <- unique(as.character(index$positive_direction[index$analysis_id == analysis_id]))
    value <- value[!is.na(value) & nzchar(value)]
    if (length(value) > 1L) stop("Conflicting positive direction in DE index.", call. = FALSE)
    if (length(value)) positive <- value[[1L]]
  }
  for (dir in sort(dirs)) {
    collection <- sub("^collection_", "", basename(dir))
    gsea <- list.files(file.path(dir, "enrichment"), pattern = "_GSEA_.*_annotated[.]tsv$", full.names = TRUE)
    gsea <- gsea[startsWith(basename(gsea), paste0(analysis_id, "_GSEA_"))]
    ledger <- file.path(dir, "qc", paste0(analysis_id, "_GSEA_universe_ledger.tsv"))
    de <- file.path(dir, "inputs", paste0(analysis_id, "_standardized_DE.tsv"))
    mapping <- file.path(dir, "qc", paste0(analysis_id, "_input_column_mapping.tsv"))
    if (length(gsea) != 1L || !file.exists(ledger) || !file.exists(de)) stop("Incomplete GSEA/DE evidence for ", analysis_id, " / ", collection, call. = FALSE)
    scopes[[length(scopes) + 1L]] <- list(analysis_id = analysis_id, collection = collection,
      positive_contrast = positive, analysis_type = "single_de", gsea_table = gsea[[1L]],
      universe_ledger = ledger, de_table = de,
      column_mapping = if (file.exists(mapping)) mapping else NULL)
  }
}
evidence <- get("build_lisa_gene_evidence", ns)(scopes, cfg$term2gene, tier = cfg$tier,
  gsea_padj_cutoff = cfg$gsea_padj_cutoff, de_padj_cutoff = cfg$de_padj_cutoff)
get("render_lisa_gene_evidence", ns)(evidence, cfg$output_dir)
cat("Global gene evidence completed: ", nrow(evidence$analyses), " individual analyses; ", nrow(evidence$scopes),
  " analysis/collection scopes; ", length(unique(evidence$identifiers$symbol)), " exact symbols; ", nrow(evidence$memberships),
  " full memberships stored once.\n", sep = "")
