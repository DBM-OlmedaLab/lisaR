#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Fit the paired longitudinal model and export three contrasts:
#
#   1. Responders: ON minus PRE.
#   2. Progressive disease: ON minus PRE.
#   3. Differential longitudinal response:
#        (responders ON−PRE) minus (PD ON−PRE).
#
# Patient fixed effects absorb every stable between-patient property, including
# response group and prior ipilimumab status. Two time interactions then allow
# treatment-associated change to differ by response and prior ipilimumab.
#
# The two prior-ipilimumab strata are balanced within each response group
# (4/8 responders and 9/18 PD patients are progressed). The displayed
# within-group effects are therefore marginal effects averaged with weight 0.5
# over the two prior-ipilimumab strata.

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))
source(file.path(riaz_example_dir(), "_model_helpers.R"))

riaz_require_packages(c(
  "AnnotationDbi", "data.table", "DESeq2", "org.Hs.eg.db",
  "SummarizedExperiment"
))

design <- riaz_read_tsv(file.path(
  riaz_output_dir("design"),
  "analysis_sample_design.tsv"
))

# Preserve the original GEO/count-matrix sample order. The statistical model
# is invariant to row order in exact arithmetic, but DESeq2's iterative
# optimisation can follow slightly different floating-point paths when rows
# are permuted. Keeping the public source order makes convergence diagnostics
# reproducible and matches the independently validated reference run.
design$patient_id <- factor(design$patient_id)
design$on_indicator <- as.integer(design$visit == "On")
design$response_pd <- as.integer(design$analysis_response == "PD")
design$ipi_progressed <- as.integer(design$ipi_cohort == "NIV3-PROG")
design$on_by_pd <- design$on_indicator * design$response_pd
design$on_by_ipi <- design$on_indicator * design$ipi_progressed

counts <- riaz_read_counts(design$sample_id)
design <- design[match(colnames(counts), design$sample_id), , drop = FALSE]

# Explicit numeric columns avoid adding response or prior-ipilimumab main
# effects, which are perfectly collinear with patient fixed effects.
model_matrix <- model.matrix(
  ~ patient_id + on_indicator + on_by_ipi + on_by_pd,
  data = design
)
fit <- riaz_fit_deseq2_matrix(counts, design, model_matrix)
dds <- fit$dds
output_dir <- riaz_output_dir("models", "longitudinal")

# Positive values mean a greater on-treatment increase (or smaller decrease)
# in the named group. The 0.5 weight is the marginal average over the balanced
# prior-ipilimumab strata.
responders <- riaz_numeric_contrast(
  dds,
  c(on_indicator = 1, on_by_ipi = 0.5)
)
pd <- riaz_numeric_contrast(
  dds,
  c(on_indicator = 1, on_by_ipi = 0.5, on_by_pd = 1)
)

# on_by_pd is PD minus responders for longitudinal change. Negating it gives
# the project-wide orientation: responders minus PD.
interaction <- riaz_numeric_contrast(dds, c(on_by_pd = -1))

responders_table <- riaz_export_model_result(
  responders,
  dds,
  "responders_on_vs_pre",
  output_dir
)
pd_table <- riaz_export_model_result(
  pd,
  dds,
  "pd_on_vs_pre",
  output_dir
)
interaction_table <- riaz_export_model_result(
  interaction,
  dds,
  "differential_longitudinal_response",
  output_dir
)

model_summary <- data.frame(
  model = "patient_blocked_longitudinal",
  samples = nrow(design),
  patients = nlevels(design$patient_id),
  responders = length(unique(design$patient_id[design$response_pd == 0L])),
  progressive_disease = length(unique(
    design$patient_id[design$response_pd == 1L]
  )),
  model_rows = nrow(model_matrix),
  model_columns = ncol(model_matrix),
  model_rank = fit$rank,
  nonconverged_genes = sum(
    !SummarizedExperiment::mcols(dds)$betaConv,
    na.rm = TRUE
  )
)
riaz_write_tsv(model_summary, file.path(output_dir, "model_summary.tsv"))
riaz_write_tsv(design, file.path(output_dir, "model_sample_design.tsv"))
riaz_write_tsv(
  as.data.frame(model_matrix),
  file.path(output_dir, "model_matrix.tsv")
)

saveRDS(dds, file.path(output_dir, "longitudinal_deseq2.rds"))
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

cat("model_rank=", fit$rank, "/", ncol(model_matrix), "\n", sep = "")
cat(
  "responders_fdr_lt_0.05=",
  sum(responders_table$padj < 0.05, na.rm = TRUE),
  "\n",
  sep = ""
)
cat(
  "pd_fdr_lt_0.05=",
  sum(pd_table$padj < 0.05, na.rm = TRUE),
  "\n",
  sep = ""
)
cat(
  "interaction_fdr_lt_0.05=",
  sum(interaction_table$padj < 0.05, na.rm = TRUE),
  "\n",
  sep = ""
)
