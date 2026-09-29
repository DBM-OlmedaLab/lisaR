#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# End-to-end driver for stages 00--09 of the Riaz GSE91061 worked example.
# Stage 10 is deliberately separate: it requires the verified run produced by
# stage 09 and a second, new output directory for the curated showcase.
#
# Required environment variables:
#   LISAR_RIAZ_PROJECT_DIR
#   LISAR_RIAZ_INPUT_CACHE
#   LISAR_RIAZ_DICTIONARY_SOURCE / LISAR_RIAZ_DICTIONARY_SHA256
#   LISAR_RIAZ_TERM2GENE_SOURCE / LISAR_RIAZ_TERM2GENE_SHA256
#   LISAR_RIAZ_CATEGORY_MAP_SOURCE / LISAR_RIAZ_CATEGORY_MAP_SHA256
#
# Large public inputs and licensed/versioned semantic resources remain outside
# the package. Every copied input is checksum-verified by its owning stage.

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))

stages <- c(
  "00_download_and_verify_inputs.R",
  "01_reconstruct_sample_metadata.R",
  "02_build_analysis_cohort.R",
  "03_fit_longitudinal_model.R",
  "04_fit_time_specific_models.R",
  "05_export_lisa_inputs.R",
  "06_validate_exported_inputs.R",
  "07_prepare_lisa_config.R",
  "08_prepare_lisa_resources.R",
  "09_run_lisa_example.R"
)

for (stage in stages) {
  cat("\n=== ", stage, " ===\n", sep = "")
  riaz_run_script(stage)
}
cat("\nRIAZ_GSE91061_WORKED_EXAMPLE=PASS\n")
cat(paste(
  "Stage 10 is not run automatically. Invoke 10_build_curated_showcase.R",
  "with <verified-lisa-run> and <new-output-directory>.\n"
))
