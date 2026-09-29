#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Define the exact worked-example cohort after reconstructing the public design.
#
# Selection is deliberately separated from metadata parsing. This makes every
# inclusion/exclusion rule visible in an audit table and prevents a later model
# script from silently changing the scientific population.

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))

riaz_require_packages("data.table")

design_dir <- riaz_output_dir("design")
sample <- riaz_read_tsv(file.path(
  design_dir,
  "reconstructed_sample_metadata.tsv"
))
patient <- riaz_read_tsv(file.path(
  design_dir,
  "reconstructed_patient_metadata.tsv"
))

patient$exclusion_reason <- NA_character_
patient$exclusion_reason[!patient$complete_pre_on_pair] <- "not_strict_pair"
patient$exclusion_reason[
  is.na(patient$exclusion_reason) &
    !(patient$analysis_response %in% c("CRPR", "PD"))
] <- "response_not_CRPR_or_PD"
patient$exclusion_reason[
  is.na(patient$exclusion_reason) & is.na(patient$ipi_cohort)
] <- "missing_prior_ipilimumab_cohort"
patient$exclusion_reason[
  is.na(patient$exclusion_reason) & patient$patient_id == "Pt3"
] <- "published_PCA_outlier_Pt3"
patient$included <- is.na(patient$exclusion_reason)

included_patients <- patient$patient_id[patient$included]
analysis <- sample[sample$patient_id %in% included_patients, , drop = FALSE]
analysis$response_pd <- as.integer(analysis$analysis_response == "PD")
analysis$on_indicator <- as.integer(analysis$visit == "On")
analysis$ipi_progressed <- as.integer(analysis$ipi_cohort == "NIV3-PROG")

# Every included patient must contribute exactly one PRE and one ON sample.
pair_table <- table(analysis$patient_id, analysis$visit)
if (
  !all(c("Pre", "On") %in% colnames(pair_table)) ||
    any(pair_table[, "Pre"] != 1L) ||
    any(pair_table[, "On"] != 1L)
) {
  stop("The selected analysis cohort is not strictly paired.")
}

group_counts <- unique(analysis[
  ,
  c("patient_id", "analysis_response", "ipi_cohort")
])
group_counts <- as.data.frame(
  xtabs(~ analysis_response + ipi_cohort, data = group_counts)
)
names(group_counts) <- c("response", "prior_ipilimumab", "patients")

expected <- c(CRPR = 8L, PD = 18L)
observed <- table(unique(analysis[, c("patient_id", "analysis_response")])$analysis_response)
if (!identical(as.integer(observed[names(expected)]), as.integer(expected))) {
  stop(
    "Unexpected final group sizes: ",
    paste(names(observed), observed, collapse = ", ")
  )
}

riaz_write_tsv(
  patient,
  file.path(design_dir, "patient_inclusion_audit.tsv")
)
riaz_write_tsv(
  analysis,
  file.path(design_dir, "analysis_sample_design.tsv")
)
riaz_write_tsv(
  group_counts,
  file.path(design_dir, "analysis_group_balance.tsv")
)

cat("included_patients=", length(included_patients), "\n", sep = "")
cat("included_samples=", nrow(analysis), "\n", sep = "")
cat("responders=8\n")
cat("progressive_disease=18\n")

