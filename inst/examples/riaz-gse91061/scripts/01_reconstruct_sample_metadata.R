#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Reconstruct sample-, patient- and clinical-design tables from the original
# GEO series matrix plus Riaz supplementary Tables S2 and S4.
#
# This is intentionally an R implementation: no intermediate metadata table
# from the historical analysis is trusted. Every field used by the models is
# reconstructed from the public source files and cross-checked before export.

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))

riaz_require_packages(c("data.table", "readxl"))

cache <- riaz_input_cache()
output_dir <- riaz_output_dir("design")
series_path <- file.path(cache, "GSE91061_series_matrix.txt.gz")
count_path <- file.path(
  cache,
  "GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz"
)

parse_geo_vector <- function(line) {
  # GEO vectors are quoted, tab-delimited fields. scan() handles both details
  # while preserving spaces and punctuation in characteristic values.
  scan(
    text = line,
    what = character(),
    sep = "\t",
    quote = "\"",
    quiet = TRUE,
    comment.char = ""
  )
}

series_connection <- gzfile(series_path, open = "rt")
geo_lines <- readLines(series_connection, warn = FALSE)
close(series_connection)
title_line <- geo_lines[startsWith(geo_lines, "!Sample_title")]
accession_line <- geo_lines[startsWith(
  geo_lines,
  "!Sample_geo_accession"
)]
characteristic_lines <- geo_lines[startsWith(
  geo_lines,
  "!Sample_characteristics_ch1"
)]
if (length(title_line) != 1L || length(accession_line) != 1L) {
  stop("Required GEO title/accession vectors are missing or duplicated.")
}

titles <- parse_geo_vector(title_line)[-1]
accessions <- parse_geo_vector(accession_line)[-1]
if (length(titles) != length(accessions)) {
  stop("GEO title and accession vectors have different lengths.")
}

characteristics <- lapply(characteristic_lines, function(line) {
  values <- parse_geo_vector(line)[-1]
  if (length(values) != length(titles)) {
    stop("A GEO characteristic vector has an unexpected length.")
  }
  split <- strsplit(values, ":", fixed = TRUE)
  keys <- unique(tolower(trimws(vapply(split, `[[`, character(1), 1L))))
  if (length(keys) != 1L) {
    stop("Mixed keys occur inside one GEO characteristic vector.")
  }
  value <- vapply(split, function(x) {
    trimws(paste(x[-1], collapse = ":"))
  }, character(1))
  list(key = keys[[1]], value = value)
})
characteristic_map <- setNames(
  lapply(characteristics, `[[`, "value"),
  vapply(characteristics, `[[`, character(1), "key")
)

pattern <- "^Pt([0-9]+)_(Pre|On)_(.+)$"
matches <- regexec(pattern, titles)
parts <- regmatches(titles, matches)
if (any(lengths(parts) != 4L)) {
  stop("Unexpected GEO sample-title format.")
}

sample <- data.frame(
  sample_id = titles,
  geo_accession = accessions,
  patient_id = paste0(
    "Pt",
    as.integer(vapply(parts, `[[`, character(1), 2L))
  ),
  visit = vapply(parts, `[[`, character(1), 3L),
  response_geo = characteristic_map[["response"]],
  aliquot = vapply(parts, `[[`, character(1), 4L),
  stringsAsFactors = FALSE
)
declared_visit <- characteristic_map[["visit (pre or on treatment)"]]
if (!identical(sample$visit, declared_visit)) {
  stop("Sample titles and GEO visit characteristics disagree.")
}

# Read only the count-matrix header. The sample order in this header becomes
# the canonical order for every subsequent model and expression matrix.
count_connection <- gzfile(count_path, open = "rt")
count_header <- names(read.csv(
  count_connection,
  nrows = 1L,
  check.names = FALSE
))[-1]
close(count_connection)
if (!setequal(count_header, sample$sample_id)) {
  stop("GEO metadata and raw-count matrix contain different sample IDs.")
}
sample <- sample[match(count_header, sample$sample_id), , drop = FALSE]

patient_split <- split(sample, sample$patient_id)
patient <- do.call(rbind, lapply(patient_split, function(rows) {
  response <- unique(rows$response_geo)
  data.frame(
    patient_id = rows$patient_id[[1]],
    response_geo = if (length(response) == 1L) response else "CONFLICT",
    n_samples = nrow(rows),
    n_pre = sum(rows$visit == "Pre"),
    n_on = sum(rows$visit == "On"),
    complete_pre_on_pair = sum(rows$visit == "Pre") == 1L &&
      sum(rows$visit == "On") == 1L,
    same_visit_duplicate = sum(rows$visit == "Pre") > 1L ||
      sum(rows$visit == "On") > 1L,
    sample_ids = paste(rows$sample_id, collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
patient <- patient[
  order(as.integer(sub("^Pt", "", patient$patient_id))),
  ,
  drop = FALSE
]
rownames(patient) <- NULL

# Clinical response and prior-ipilimumab cohort are sourced from Table S2.
clinical <- as.data.frame(
  readxl::read_excel(
    file.path(cache, "Riaz_Cell_2017_Table_S2.xlsx"),
    sheet = "Table S2",
    skip = 2
  ),
  stringsAsFactors = FALSE
)
clinical <- clinical[, c("Patient", "Cohort", "Response"), drop = FALSE]
names(clinical) <- c("patient_id", "ipi_cohort", "published_response")

# Table S4 provides an independent assay-availability cross-check.
availability <- as.data.frame(
  readxl::read_excel(
    file.path(cache, "Riaz_Cell_2017_Table_S4.xlsx"),
    sheet = "Table S4",
    skip = 2
  ),
  stringsAsFactors = FALSE
)
availability <- availability[
  ,
  c("Patient", "Pre-treatment RNA-Seq", "On-treatment RNA-Seq"),
  drop = FALSE
]
names(availability) <- c(
  "patient_id",
  "published_pre_rnaseq",
  "published_on_rnaseq"
)

patient <- merge(patient, clinical, by = "patient_id", all.x = TRUE, sort = FALSE)
patient <- merge(
  patient,
  availability,
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)
patient <- patient[
  order(as.integer(sub("^Pt", "", patient$patient_id))),
  ,
  drop = FALSE
]

patient$analysis_response <- ifelse(
  patient$published_response %in% c("CR", "PR"),
  "CRPR",
  ifelse(
    patient$published_response == "SD",
    "SD",
    ifelse(patient$published_response == "PD", "PD", "UNEVALUABLE")
  )
)
patient$response_concordant <- with(patient, {
  geo_equivalent <- ifelse(
    published_response %in% c("CR", "PR"),
    "PRCR",
    ifelse(published_response == "NE", "UNK", published_response)
  )
  is.na(geo_equivalent) | response_geo == geo_equivalent
})
patient$availability_concordant <- with(patient, {
  (n_pre > 0) == (as.integer(published_pre_rnaseq) == 1L) &
    (n_on > 0) == (as.integer(published_on_rnaseq) == 1L)
})

sample <- merge(
  sample,
  patient[
    ,
    c(
      "patient_id", "ipi_cohort", "published_response",
      "analysis_response", "complete_pre_on_pair",
      "same_visit_duplicate"
    )
  ],
  by = "patient_id",
  all.x = TRUE,
  sort = FALSE
)
sample <- sample[match(count_header, sample$sample_id), , drop = FALSE]

summary <- data.frame(
  metric = c(
    "public_samples", "public_patients", "patients_with_pre",
    "patients_with_on", "strict_pairs",
    "strict_pairs_with_clinical_and_ipi",
    "same_visit_duplicate_patients", "response_disagreements",
    "availability_disagreements"
  ),
  value = c(
    nrow(sample),
    nrow(patient),
    sum(patient$n_pre > 0),
    sum(patient$n_on > 0),
    sum(patient$complete_pre_on_pair),
    sum(
      patient$complete_pre_on_pair &
        patient$analysis_response %in% c("CRPR", "SD", "PD") &
        !is.na(patient$ipi_cohort)
    ),
    sum(patient$same_visit_duplicate),
    sum(!patient$response_concordant, na.rm = TRUE),
    sum(!patient$availability_concordant, na.rm = TRUE)
  )
)

riaz_write_tsv(sample, file.path(output_dir, "reconstructed_sample_metadata.tsv"))
riaz_write_tsv(patient, file.path(output_dir, "reconstructed_patient_metadata.tsv"))
riaz_write_tsv(summary, file.path(output_dir, "reconstruction_summary.tsv"))

cat("samples=", nrow(sample), "\n", sep = "")
cat("patients=", nrow(patient), "\n", sep = "")
cat("strict_pairs=", sum(patient$complete_pre_on_pair), "\n", sep = "")
