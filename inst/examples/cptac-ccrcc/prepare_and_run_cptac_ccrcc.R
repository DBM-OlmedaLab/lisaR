#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 or later (GPL-3.0-or-later). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Prepare and execute one bounded local CPTAC ccRCC global-proteomic LISA run.
#
# This script deliberately accepts only an already-reviewed *local* C3 bundle.
# It neither downloads CPTAC material nor distributes it.  The bundle contains
# the frozen primary DE table and the minimum paired-design/QC receipts needed
# to audit that input; no clinical extras are requested.  It copies those
# verified files into a new project through a staging directory, writes a
# relative-path configuration, validates normal lisaR resources, then performs
# exactly one standard run.  It never refits limma or creates sensitivity,
# HALLMARKS, KEGG-map, or full-report products.
#
# Set LISAR_CPTAC_CCRCC_PREPARE_ONLY=1 to stop right after the verified
# project is promoted: every preparation step above still runs (bundle/
# manifest verification, staged copy, resource-cache resolution, config
# synthesis, `validate_lisa_config()`, `plan_lisa_outputs()`, atomic
# promotion), but `run_lisa()` is never called. The script instead prints the
# exact command to perform that run later. Leaving the variable unset keeps
# the default full run unchanged.

c3_fail <- function(...) stop(..., call. = FALSE)
c3_assert <- function(ok, ...) if (!isTRUE(ok)) c3_fail(...)
c3_path_is_link <- function(path) {
  value <- Sys.readlink(path)
  length(value) == 1L && !is.na(value) && nzchar(value)
}
c3_regular_file <- function(path) {
  file.exists(path) && !dir.exists(path) &&
    isTRUE(utils::file_test("-f", path)) && !c3_path_is_link(path)
}
c3_sha256 <- function(path) {
  get("lisa_sha256_file", envir = asNamespace("lisaR"), inherits = FALSE)(path)
}
c3_relative <- function(path) {
  valid <- is.character(path) && length(path) == 1L && !is.na(path) &&
    nzchar(path) && !grepl("[[:cntrl:]]", path) && !grepl("\\\\", path) &&
    !startsWith(path, "/")
  pieces <- if (valid) strsplit(path, "/", fixed = TRUE)[[1L]] else character()
  c3_assert(valid && length(pieces) && !any(pieces %in% c("", ".", "..")),
            "C3 bundle manifest contains an unsafe relative path.")
  path
}
c3_within <- function(path, root) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  identical(path, root) || startsWith(path, paste0(root, "/"))
}
c3_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE,
                     na = "NA")
}

required_packages <- c("jsonlite", "yaml", "lisaR")
missing_packages <- required_packages[!vapply(
  required_packages, requireNamespace, logical(1), quietly = TRUE
)]
c3_assert(!length(missing_packages), "Required installed R packages are unavailable: ",
          paste(missing_packages, collapse = ", "))

bundle_value <- Sys.getenv("LISAR_CPTAC_CCRCC_BUNDLE", unset = "")
project_value <- Sys.getenv("LISAR_CPTAC_CCRCC_PROJECT_DIR", unset = "")
c3_assert(nzchar(bundle_value),
          "LISAR_CPTAC_CCRCC_BUNDLE is required and must name a reviewed local C3 bundle.")
c3_assert(nzchar(project_value),
          "LISAR_CPTAC_CCRCC_PROJECT_DIR is required and must name a new project directory.")
c3_assert(dir.exists(bundle_value) && !c3_path_is_link(bundle_value),
          "LISAR_CPTAC_CCRCC_BUNDLE must be an existing non-symlink directory.")
bundle <- normalizePath(bundle_value, winslash = "/", mustWork = TRUE)
manifest_path <- file.path(bundle, "cptac_ccrcc_bundle_manifest.tsv")
c3_assert(c3_regular_file(manifest_path),
          "The local C3 bundle needs one regular cptac_ccrcc_bundle_manifest.tsv.")
manifest <- utils::read.delim(manifest_path, sep = "\t", quote = "",
                              comment.char = "", stringsAsFactors = FALSE,
                              check.names = FALSE, colClasses = "character")
expected <- data.frame(
  path = c(
    "de/de_primary_conservative_80pairs.tsv",
    "audit/paired_design_samples.tsv",
    "audit/checks_primary_conservative_80pairs.tsv",
    "audit/data_quality_summary.tsv",
    "audit/input_manifest.tsv",
    "audit/pairs_80_source_validated.tsv"
  ),
  role = c("primary_de", "paired_design", "primary_qc", "data_quality",
           "input_provenance", "cohort_pairs"), stringsAsFactors = FALSE
)
c3_assert(identical(names(manifest), c("path", "bytes", "sha256", "role")) &&
            nrow(manifest) == nrow(expected),
          "C3 bundle manifest must contain exactly path, bytes, sha256 and role for the six reviewed inputs.")
manifest$path <- vapply(manifest$path, c3_relative, character(1))
c3_assert(!anyDuplicated(manifest$path) &&
            identical(manifest$path, expected$path) &&
            identical(manifest$role, expected$role) &&
            all(grepl("^[0-9a-f]{64}$", manifest$sha256)) &&
            all(is.finite(suppressWarnings(as.numeric(manifest$bytes))) &
                suppressWarnings(as.numeric(manifest$bytes)) > 0),
          "C3 bundle manifest differs from the approved primary-DE/audit allowlist.")
bundle_files <- file.path(bundle, manifest$path)
c3_assert(all(vapply(bundle_files, c3_regular_file, logical(1))) &&
            all(vapply(bundle_files, c3_within, logical(1), root = bundle)),
          "A C3 bundle input is missing, unsafe, or escapes the bundle root.")
c3_assert(identical(as.numeric(file.info(bundle_files)$size), as.numeric(manifest$bytes)) &&
            identical(unname(vapply(bundle_files, c3_sha256, character(1))), manifest$sha256),
          "C3 bundle file sizes or SHA-256 values do not match its reviewed manifest.")
primary_sha256 <- "7c55e895cfddec99cae8856b503fcbd7b8714a9916f3295db4ee17ad61767dba"
c3_assert(identical(manifest$sha256[[match("de/de_primary_conservative_80pairs.tsv", manifest$path)]], primary_sha256),
          "The bundle does not contain the accepted C3 primary conservative 80-pair DE table.")

primary <- utils::read.delim(bundle_files[[match("de/de_primary_conservative_80pairs.tsv", manifest$path)]],
                             sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
required_de <- c("entrez_id", "gene_symbol", "tumor_minus_nat", "t",
                 "p_value", "fdr_bh_all_genes")
c3_assert(nrow(primary) == 6482L && all(required_de %in% names(primary)) &&
            !anyNA(primary$gene_symbol) && !anyDuplicated(toupper(primary$gene_symbol)) &&
            all(is.finite(primary$t)) && all(is.finite(primary$tumor_minus_nat)) &&
            all(is.finite(primary$p_value)) && all(is.finite(primary$fdr_bh_all_genes)),
          "The accepted primary DE table is not the complete 6,482-gene, unique-symbol C3 input.")

project_parent <- dirname(path.expand(project_value))
c3_assert(dir.exists(project_parent) && !c3_path_is_link(project_parent),
          "LISAR_CPTAC_CCRCC_PROJECT_DIR must have an existing non-symlink parent.")
project <- file.path(normalizePath(project_parent, winslash = "/", mustWork = TRUE),
                     basename(path.expand(project_value)))
c3_assert(!basename(project) %in% c("", ".", ".."),
          "LISAR_CPTAC_CCRCC_PROJECT_DIR must name a new child directory.")
c3_assert(!file.exists(project) && !dir.exists(project) && !c3_path_is_link(project),
          "The requested CPTAC project already exists; preserve it and use a new destination.")
stage <- file.path(dirname(project), paste0(".", basename(project), ".staging-", Sys.getpid()))
c3_assert(!file.exists(stage) && !dir.exists(stage) && !c3_path_is_link(stage),
          "The deterministic CPTAC project staging directory already exists; inspect it and use a new destination.")
dir.create(stage, recursive = FALSE)
# A failure before the promotion deliberately leaves this directory in place
# for inspection.  Avoid a top-level `on.exit()` handler: when sourced by a
# caller it binds to the caller frame and can falsely report retention after a
# successful promotion.

for (i in seq_len(nrow(manifest))) {
  destination <- file.path(stage, "inputs", manifest$path[[i]])
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  c3_assert(file.copy(bundle_files[[i]], destination, overwrite = FALSE, copy.mode = TRUE),
            "Could not copy verified C3 bundle input into project staging: ", manifest$path[[i]])
  c3_assert(identical(c3_sha256(destination), manifest$sha256[[i]]),
            "Copied project input failed the bundle SHA-256 verification: ", manifest$path[[i]])
}
dir.create(file.path(stage, "provenance"), recursive = TRUE, showWarnings = FALSE)
file.copy(manifest_path, file.path(stage, "provenance", "cptac_ccrcc_bundle_manifest.tsv"),
          overwrite = FALSE, copy.mode = TRUE)
bundle_manifest_sha256 <- c3_sha256(manifest_path)

resource_cache <- Sys.getenv("LISAR_CPTAC_CCRCC_RESOURCE_CACHE", unset = "")
if (!nzchar(resource_cache)) {
  resource_cache <- get("lisa_dictionary_cache_root", envir = asNamespace("lisaR"),
                        inherits = FALSE)()
}
resource_cache <- path.expand(resource_cache)
c3_assert(!c3_within(resource_cache, stage),
          "LISAR_CPTAC_CCRCC_RESOURCE_CACHE must be outside the prepared project.")
options(lisaR.dictionary_cache_root = resource_cache,
        lisaR.dictionary_registry = NULL,
        lisaR.shared_dictionary_root = NULL)

config <- list(
  pipeline = list(
    schema_version = "1.0.0", profile = "global proteomic", evidence_mode = "full_de",
    output_dir = "../runs/cptac-ccrcc-global-proteomic-standard", workers = 1L,
    gsea_padj_cutoff = 0.25, dictionary_resource = "lisa_dictionary_core@1.0.0",
    term2gene_resource = "msigdb_term2gene@2026.1",
    category_map_resource = "lisa_category_map@1.0.0", lisa_dictionary = "core",
    file_label_prefix = "CPTAC_CCRCC_C3", report_title = "CPTAC ccRCC — tumor versus adjacent non-tumoral tissue",
    run_ora = FALSE, run_kegg_maps = FALSE, run_hallmarks = FALSE, dry_run = FALSE,
    resume = FALSE, sufficiency_action = "warn",
    duplicate_policies = list(de_table_duplicate_policy = "error",
      matrix_duplicate_policy = "error", mapped_id_collision_policy = "error")
  ),
  report = list(mode = "standard", category_evidence = TRUE,
    category_nes_variants = "clean", legacy_gene_products = FALSE,
    evidence_max_sets = 15L, evidence_max_genes = 40L,
    formats = list(png = TRUE, svg = FALSE, pdf = FALSE), source_data = TRUE,
    recipes = TRUE),
  project = list(title = "CPTAC ccRCC — paired global proteomics"),
  collections = c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS"),
  single_de = list(list(
    analysis_id = "cptac_ccrcc_tumor_minus_nat", label = "ccRCC tumor vs adjacent non-tumoral tissue",
    comparison = "Paired tumor versus adjacent non-tumoral renal tissue abundance", positive_direction = "Positive values indicate higher tumour abundance.",
    model_note = "Frozen C3 primary conservative 80-pair limma output; this route does not refit differential abundance.",
    de_path = "../inputs/de/de_primary_conservative_80pairs.tsv", species = "Homo sapiens",
    symbol_col = "gene_symbol", rank_col = "t", logfc_col = "tumor_minus_nat",
    pvalue_col = "p_value", padj_col = "fdr_bh_all_genes"
  )),
  allowlisted_source_paths = c(
    "../inputs/de/de_primary_conservative_80pairs.tsv",
    "../inputs/audit/paired_design_samples.tsv",
    "../inputs/audit/checks_primary_conservative_80pairs.tsv",
    "../inputs/audit/data_quality_summary.tsv",
    "../inputs/audit/input_manifest.tsv",
    "../inputs/audit/pairs_80_source_validated.tsv"
  ),
  source_data = list(
    list(path = "../inputs/de/de_primary_conservative_80pairs.tsv", role = "frozen_primary_differential_abundance", target_subdir = "de_inputs"),
    list(path = "../inputs/audit/paired_design_samples.tsv", role = "paired_design_audit", target_subdir = "audit"),
    list(path = "../inputs/audit/checks_primary_conservative_80pairs.tsv", role = "primary_qc", target_subdir = "audit"),
    list(path = "../inputs/audit/data_quality_summary.tsv", role = "data_quality", target_subdir = "audit"),
    list(path = "../inputs/audit/input_manifest.tsv", role = "input_provenance", target_subdir = "audit"),
    list(path = "../inputs/audit/pairs_80_source_validated.tsv", role = "cohort_pairs", target_subdir = "audit")
  )
)
dir.create(file.path(stage, "config"), recursive = TRUE, showWarnings = FALSE)
config_path <- file.path(stage, "config", "cptac-ccrcc-global-proteomic.yml")
json_path <- file.path(stage, "config", "cptac-ccrcc-global-proteomic.json")
yaml::write_yaml(config, config_path)
jsonlite::write_json(config, json_path, auto_unbox = TRUE, pretty = TRUE,
                     null = "null", digits = NA)
c3_write_tsv(data.frame(format = c("YAML", "JSON"), path = file.path("config", basename(c(config_path, json_path))),
                        sha256 = c(c3_sha256(config_path), c3_sha256(json_path)), stringsAsFactors = FALSE),
             file.path(stage, "provenance", "config_manifest.tsv"))

validation <- lisaR::validate_lisa_config(config_path, check_files = TRUE, strict = TRUE)
c3_write_tsv(validation$readiness, file.path(stage, "logs", "validation_readiness.tsv"))
c3_write_tsv(validation$resources, file.path(stage, "logs", "resource_resolution.tsv"))
c3_assert(isTRUE(validation$execution_ready),
          "CPTAC global-proteomic configuration is not execution-ready; inspect staging logs before retrying.")
c3_assert(identical(validation$config$pipeline$profile, "global proteomic") &&
            validation$analyses == 1L && validation$contrasts == 0L &&
            isFALSE(validation$config$pipeline$run_kegg_maps) &&
            isFALSE(validation$config$pipeline$run_hallmarks),
          "CPTAC execution configuration drifted from the bounded primary-only contract.")

plan <- lisaR::plan_lisa_outputs(config_path)
c3_write_tsv(plan$summary, file.path(stage, "logs", "output_plan_summary.tsv"))
get("lisa_compact_example_project", asNamespace("lisaR"))(stage, "cptac", config_path)
# Promote the verified preparation before `run_lisa()` starts. Configuration
# paths are relative, so canonical run receipts use final project locations.
# If the scientific run fails, the prepared project remains available for audit.
c3_assert(file.rename(stage, project),
          "Could not promote the verified CPTAC preparation directory.")
config_path <- file.path(project, "study.yml")

prepare_only <- toupper(Sys.getenv("LISAR_CPTAC_CCRCC_PREPARE_ONLY", unset = "FALSE")) %in%
  c("1", "TRUE", "YES")
if (prepare_only) {
  # The preparation above (bundle/manifest verification, staged copy,
  # resource-cache resolution, config synthesis, validate_lisa_config(),
  # plan_lisa_outputs(), atomic promotion) is already complete and identical
  # to the full-run path. Stop here: no run_lisa(), no verify_lisa_run(), no
  # run receipt, and therefore no runs/ directory under `project`.
  next_run_command <- paste0(
    "Rscript --vanilla -e ", shQuote(paste0(
      "result <- lisaR::run_lisa(", deparse(config_path), "); ",
      "verification <- lisaR::verify_lisa_run(result$output_dir); ",
      "stopifnot(identical(verification$gate, 'PASS')); ",
      "cat('RUN_DIR=', result$output_dir, '\\n', sep = '')"
    ))
  )
  cat("CPTAC_CCRCC_LOCAL_ROUTE=PREPARE_ONLY_PASS\n")
  cat("PROJECT_DIR=", project, "\n", sep = "")
  cat("CONFIG_PATH=", config_path, "\n", sep = "")
  cat("NEXT_RUN_COMMAND=", next_run_command, "\n", sep = "")
} else {
  result <- lisaR::run_lisa(config_path)
  verification <- lisaR::verify_lisa_run(result$output_dir)
  c3_write_tsv(data.frame(gate = verification$gate,
                           findings = paste(verification$findings, collapse = " | "),
                           stringsAsFactors = FALSE),
               file.path(project, "logs", "run_verification.tsv"))
  c3_assert(identical(verification$gate, "PASS"),
            "The completed CPTAC LISA run did not pass verification; prepared project was retained for audit.")
  c3_write_tsv(data.frame(
    bundle_manifest_sha256 = bundle_manifest_sha256, primary_de_sha256 = primary_sha256,
    input_rows = nrow(primary), profile = "global proteomic", rank_column = "t",
    effect_column = "tumor_minus_nat", pvalue_column = "p_value",
    padj_column = "fdr_bh_all_genes", workers = 1L, run_gate = verification$gate,
    output_dir = "out", stringsAsFactors = FALSE
  ), file.path(project, "logs", "cptac_ccrcc_prepared_receipt.tsv"))
  cat("CPTAC_CCRCC_LOCAL_ROUTE=PASS\n")
  cat("PROJECT_DIR=", project, "\n", sep = "")
  cat("RUN_DIR=", file.path(project, "out"), "\n", sep = "")
}
