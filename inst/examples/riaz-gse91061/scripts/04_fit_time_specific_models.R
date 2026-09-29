#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Fit response-group comparisons separately at PRE and ON:
#
#   responders minus PD at PRE
#   responders minus PD at ON
#
# These are intentionally separate cross-sectional models. A patient fixed
# effect cannot estimate a stable between-patient response-group difference,
# because each patient belongs permanently to one response group. At each visit
# there is one sample per patient, so the identifiable model is:
#
#   ~ prior_ipilimumab + response_group
#
# Positive exported effects always mean higher expression in responders.

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))
source(file.path(riaz_example_dir(), "_model_helpers.R"))

riaz_require_packages(c(
  "AnnotationDbi", "data.table", "DESeq2", "org.Hs.eg.db",
  "SummarizedExperiment"
))

all_design <- riaz_read_tsv(file.path(
  riaz_output_dir("design"),
  "analysis_sample_design.tsv"
))
all_counts <- riaz_read_counts(all_design$sample_id)
output_dir <- riaz_output_dir("models", "time_specific")
summary_rows <- list()

for (visit_name in c("Pre", "On")) {
  design <- all_design[all_design$visit == visit_name, , drop = FALSE]
  design <- design[match(
    colnames(all_counts)[colnames(all_counts) %in% design$sample_id],
    design$sample_id
  ), , drop = FALSE]
  counts <- all_counts[, design$sample_id, drop = FALSE]
  design$ipi_progressed <- as.integer(design$ipi_cohort == "NIV3-PROG")
  design$response_pd <- as.integer(design$analysis_response == "PD")

  model_matrix <- model.matrix(
    ~ ipi_progressed + response_pd,
    data = design
  )
  fit <- riaz_fit_deseq2_matrix(counts, design, model_matrix)

  # The fitted response_pd coefficient is PD minus responders. Negating it
  # provides the stated responders-minus-PD orientation.
  result <- riaz_numeric_contrast(fit$dds, c(response_pd = -1))
  analysis_id <- if (visit_name == "Pre") {
    "responders_vs_pd_pre"
  } else {
    "responders_vs_pd_on"
  }
  exported <- riaz_export_model_result(
    result,
    fit$dds,
    analysis_id,
    output_dir
  )
  saveRDS(
    fit$dds,
    file.path(output_dir, paste0(tolower(visit_name), "_deseq2.rds"))
  )
  summary_rows[[visit_name]] <- data.frame(
    analysis_id = analysis_id,
    visit = visit_name,
    samples = nrow(design),
    responders = sum(design$response_pd == 0L),
    progressive_disease = sum(design$response_pd == 1L),
    model_columns = ncol(model_matrix),
    model_rank = fit$rank,
    fdr_lt_0_05 = sum(exported$padj < 0.05, na.rm = TRUE),
    nonconverged_genes = sum(
      !SummarizedExperiment::mcols(fit$dds)$betaConv,
      na.rm = TRUE
    )
  )
}

riaz_write_tsv(
  do.call(rbind, summary_rows),
  file.path(output_dir, "time_specific_model_summary.tsv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

print(do.call(rbind, summary_rows), row.names = FALSE)

