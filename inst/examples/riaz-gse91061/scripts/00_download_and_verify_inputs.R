#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Acquire and verify the original inputs used by the worked example.
#
# GEO files have stable public download URLs and are downloaded automatically
# when absent. Riaz supplementary spreadsheets are publisher-hosted materials:
# the script never redistributes them, but it verifies user-provided originals
# against the exact hashes used for this example.
#
# Environment:
#   LISAR_RIAZ_PROJECT_DIR  writable analysis root
#   LISAR_RIAZ_INPUT_CACHE  immutable/download cache (defaults inside project root)

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))

riaz_require_packages("data.table")

cache <- riaz_input_cache()
output_dir <- riaz_output_dir("input_audit")

manifest <- data.frame(
  file_name = c(
    "GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz",
    "GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz",
    "GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt.gz",
    "GSE91061_series_matrix.txt.gz",
    "Riaz_Cell_2017_Table_S1.xlsx",
    "Riaz_Cell_2017_Table_S2.xlsx",
    "Riaz_Cell_2017_Table_S3.xlsx",
    "Riaz_Cell_2017_Table_S4.xlsx",
    "Riaz_Cell_2017_Table_S5.xlsx",
    "Riaz_Cell_2017_Table_S6.xlsx"
  ),
  sha256 = c(
    "11f69d1cee7771304732c88a5bfa79fdb4fc364d840e75e77e3baa833d08a07a",
    "d05f37f3b0bba0b7ad5ea825af9b4a2eda1ee209668618f9f228213e70bc03fe",
    "1d2a65245f3f6aeeefbdf8dbc17c968598ac0d8e28b16df9c07bf03304567582",
    "f3ebd04d019afab9f12af0b92cfcafa516bb2b170f5b45d22d5e330e2eebdca0",
    "400ed0c2ef7b25763e68c4e2ff141e7c1e88eee703fe42e5cc3bca5c57b87fda",
    "0d84ca0a83c98b0dc1a56528eab0ff74c22e6f0851a470fe28788c6141f59619",
    "d4354f7a27f64bfc0fa0ba98a4f2056de2b4efa33edcf0be0e78919d8b26b82e",
    "08e6ded09d843dba57d9c3d4bb53ed2b8fba5ee1a788aa7c569a1f1ee03d0643",
    "5a315062214ce4c5b84d558c81b36d9d4f45d739dc5d86987537de3bca072dff",
    "364cf0aa5d4fd9612c8b51e78243e2dac3f21d79fef6fb1688cf01c948941717"
  ),
  source = c(
    rep("NCBI GEO GSE91061 supplementary files", 3),
    "NCBI GEO GSE91061 series matrix",
    rep("Riaz et al. Cell 2017 supplementary tables", 6)
  ),
  url = c(
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/suppl/GSE91061_BMS038109Sample_Cytolytic_Score_20161026.txt.gz",
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE91nnn/GSE91061/matrix/GSE91061_series_matrix.txt.gz",
    rep("https://doi.org/10.1016/j.cell.2017.09.028", 6)
  ),
  automatic_download = c(rep(TRUE, 4), rep(FALSE, 6)),
  stringsAsFactors = FALSE
)

audit <- manifest
audit$path <- file.path(cache, audit$file_name)
audit$downloaded <- FALSE
audit$observed_sha256 <- NA_character_
audit$status <- "pending"

for (i in seq_len(nrow(audit))) {
  path <- audit$path[[i]]
  if (!file.exists(path) && isTRUE(audit$automatic_download[[i]])) {
    message("Downloading ", audit$file_name[[i]])
    tmp <- paste0(path, ".partial")
    utils::download.file(audit$url[[i]], tmp, mode = "wb", quiet = FALSE)
    file.rename(tmp, path)
    audit$downloaded[[i]] <- TRUE
  }
  if (!file.exists(path)) {
    stop(
      "Required original input is absent: ", path, "\n",
      "Obtain it from ", audit$url[[i]], " and rerun this verification stage."
    )
  }
  observed <- riaz_sha256(path)
  audit$observed_sha256[[i]] <- observed
  audit$status[[i]] <- if (identical(observed, audit$sha256[[i]])) {
    "verified"
  } else {
    "hash_mismatch"
  }
}

riaz_write_tsv(audit, file.path(output_dir, "input_manifest.tsv"))
if (any(audit$status != "verified")) {
  stop(
    "One or more inputs failed SHA-256 verification: ",
    paste(audit$file_name[audit$status != "verified"], collapse = ", ")
  )
}

cat("verified_inputs=", nrow(audit), "\n", sep = "")
