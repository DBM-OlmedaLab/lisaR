#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Assemble the five model outputs into lisaR-ready DE tables and supporting
# expression data.
#
# lisaR receives the complete tested universe for each analysis. Significance
# filtering is deliberately deferred to lisaR. Supplying only significant
# genes would invalidate gene-set enrichment and background calculations.

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))
source(file.path(riaz_example_dir(), "_model_helpers.R"))

riaz_require_packages(c(
  "AnnotationDbi", "data.table", "DESeq2", "org.Hs.eg.db"
))

longitudinal_dir <- riaz_output_dir("models", "longitudinal")
time_dir <- riaz_output_dir("models", "time_specific")
output_dir <- riaz_output_dir("lisa_inputs")
analysis_ids <- c(
  "responders_on_vs_pre",
  "pd_on_vs_pre",
  "responders_vs_pd_pre",
  "responders_vs_pd_on",
  "differential_longitudinal_response"
)

source_dir <- setNames(
  c(longitudinal_dir, longitudinal_dir, time_dir, time_dir, longitudinal_dir),
  analysis_ids
)
summary_rows <- list()

for (analysis_id in analysis_ids) {
  source_path <- file.path(
    source_dir[[analysis_id]],
    paste0(analysis_id, "_symbol_full.tsv")
  )
  table <- riaz_read_tsv(source_path)
  export <- table[
    ,
    c(
      "symbol", "entrez_id", "baseMean", "log2FoldChange",
      "lfcSE", "stat", "pvalue", "padj", "beta_converged"
    ),
    drop = FALSE
  ]

  # Preserve the complete tested gene universe, but never let a coefficient
  # that failed DESeq2's convergence diagnostic contribute to a LISA ranking
  # or significance call. The unmodified values remain available in the
  # *_symbol_full.tsv model table, making this quality-control decision fully
  # auditable. This is a convergence filter, not a significance filter.
  nonconverged <- !is.na(export$beta_converged) & !export$beta_converged
  export$lisa_eligible <- !nonconverged & is.finite(export$stat)
  export$lisa_exclusion_reason <- ifelse(
    nonconverged,
    "DESeq2 beta coefficient did not converge",
    ifelse(!is.finite(export$stat), "No finite ranking statistic", NA_character_)
  )
  export[nonconverged, c(
    "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
  )] <- NA

  riaz_write_tsv(export, file.path(output_dir, paste0(analysis_id, ".tsv")))
  raw_fdr <- sum(table$padj < 0.05, na.rm = TRUE)
  raw_nonconverged_fdr <- sum(
    !table$beta_converged & table$padj < 0.05,
    na.rm = TRUE
  )
  summary_rows[[analysis_id]] <- data.frame(
    analysis_id = analysis_id,
    tested_entrez_rows = nrow(riaz_read_tsv(file.path(
      source_dir[[analysis_id]],
      paste0(analysis_id, "_entrez_full.tsv")
    ))),
    unique_symbols = nrow(export),
    finite_rank_values = sum(is.finite(export$stat)),
    fdr_lt_0_05 = sum(export$padj < 0.05, na.rm = TRUE),
    raw_model_fdr_lt_0_05 = raw_fdr,
    excluded_nonconverged_genes = sum(nonconverged),
    excluded_nonconverged_raw_fdr_lt_0_05 = raw_nonconverged_fdr,
    nonconverged_lisa_fdr_lt_0_05 = sum(
      !export$beta_converged & export$padj < 0.05, na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

# Use normalized counts from the longitudinal model for optional expression-
# supported lisaR products. One representative Entrez row per current symbol
# is selected by highest mean normalized abundance.
dds <- readRDS(file.path(longitudinal_dir, "longitudinal_deseq2.rds"))
normalized <- DESeq2::counts(dds, normalized = TRUE)
expression <- log2(normalized + 1)
meta <- data.frame(
  entrez_id = rownames(expression),
  symbol = AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = rownames(expression),
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  ),
  mean_normalized_count = rowMeans(normalized),
  stringsAsFactors = FALSE
)
meta <- meta[!is.na(meta$symbol) & nzchar(meta$symbol), , drop = FALSE]
meta <- meta[
  order(
    meta$symbol,
    -meta$mean_normalized_count,
    suppressWarnings(as.numeric(meta$entrez_id)),
    na.last = TRUE
  ),
  ,
  drop = FALSE
]
meta <- meta[!duplicated(meta$symbol), , drop = FALSE]
expression_export <- data.frame(
  symbol = meta$symbol,
  expression[meta$entrez_id, , drop = FALSE],
  check.names = FALSE
)
riaz_write_tsv(
  expression_export,
  file.path(output_dir, "paired_log2_normalized_counts_by_symbol.tsv")
)
riaz_write_tsv(
  riaz_read_tsv(file.path(
    riaz_output_dir("design"),
    "analysis_sample_design.tsv"
  )),
  file.path(output_dir, "paired_sample_annotation.tsv")
)
riaz_write_tsv(
  do.call(rbind, summary_rows),
  file.path(output_dir, "lisa_input_summary.tsv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

print(do.call(rbind, summary_rows), row.names = FALSE)
