#!/usr/bin/env Rscript
# Build the reviewed local-only C2 Riaz prepared-input bundle.
#
# This maintainer utility copies no executable files and never creates an
# archive. It is intentionally not the public installer: recipients use
# inst/examples/riaz-gse91061/scripts/11_prepare_prepared_project.R.

source_root <- Sys.getenv("LISAR_RIAZ_CANONICAL_PROJECT", unset = "")
bundle_root <- Sys.getenv("LISAR_RIAZ_PREPARED_BUNDLE_OUT", unset = "")
if (!nzchar(source_root) || !nzchar(bundle_root)) {
  stop(
    "Set LISAR_RIAZ_CANONICAL_PROJECT and LISAR_RIAZ_PREPARED_BUNDLE_OUT.",
    call. = FALSE
  )
}
if (!dir.exists(source_root) || file.exists(bundle_root) || dir.exists(bundle_root)) {
  stop("The canonical source must exist and the bundle output must be new.", call. = FALSE)
}
source_root <- normalizePath(source_root, winslash = "/", mustWork = TRUE)
bundle_root <- normalizePath(bundle_root, winslash = "/", mustWork = FALSE)
required <- c(
  "outputs/design/analysis_sample_design.tsv",
  file.path("outputs/models/longitudinal", c(
    "responders_on_vs_pre_symbol_full.tsv", "pd_on_vs_pre_symbol_full.tsv",
    "differential_longitudinal_response_symbol_full.tsv"
  )),
  file.path("outputs/models/time_specific", c(
    "responders_vs_pd_pre_symbol_full.tsv", "responders_vs_pd_on_symbol_full.tsv"
  )),
  file.path("outputs/lisa_inputs", c(
    "responders_on_vs_pre.tsv", "pd_on_vs_pre.tsv", "responders_vs_pd_pre.tsv",
    "responders_vs_pd_on.tsv", "differential_longitudinal_response.tsv",
    "paired_log2_normalized_counts_by_symbol.tsv", "paired_sample_annotation.tsv"
  ))
)
all_inputs <- list.files(
  file.path(source_root, "outputs"), pattern = "[.]tsv$", recursive = TRUE,
  full.names = FALSE
)
all_inputs <- file.path("outputs", all_inputs)
if (!all(required %in% all_inputs) || !all(vapply(all_inputs, function(path) {
  any(startsWith(path, c(
    "outputs/design/", "outputs/models/", "outputs/lisa_inputs/"
  )))
}, logical(1)))) {
  stop("Canonical project lacks the reviewed Riaz prepared input layout.", call. = FALSE)
}
paths <- sort(all_inputs)
sources <- file.path(source_root, paths)
if (!all(file.exists(sources)) || any(dir.exists(sources))) {
  stop("Canonical prepared inputs contain a missing or nonregular path.", call. = FALSE)
}
dir.create(bundle_root, recursive = TRUE, mode = "0700")
for (i in seq_along(sources)) {
  destination <- file.path(bundle_root, paths[[i]])
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE, mode = "0700")
  if (!file.copy(sources[[i]], destination, overwrite = FALSE, copy.mode = TRUE)) {
    stop("Failed to copy canonical prepared input: ", paths[[i]], call. = FALSE)
  }
}
bundle_files <- file.path(bundle_root, paths)
# The publication-facing installer requires SHA-256. Create it through the
# installed lisaR primitive rather than accepting user-provided hashes.
if (!requireNamespace("lisaR", quietly = TRUE)) {
  stop("lisaR is required to create the reviewed SHA-256 manifest.", call. = FALSE)
}
sha256 <- vapply(bundle_files, function(path) {
  get("lisa_sha256_file", asNamespace("lisaR"), inherits = FALSE)(path)
}, character(1))
manifest <- data.frame(
  path = paths,
  bytes = as.numeric(file.info(bundle_files)$size),
  sha256 = sha256,
  role = ifelse(startsWith(paths, "outputs/design/"), "design",
                ifelse(startsWith(paths, "outputs/models/"), "model_output",
                       "lisa_input")),
  stringsAsFactors = FALSE
)
utils::write.table(manifest, file.path(bundle_root, "prepared_bundle_manifest.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
provenance <- data.frame(
  key = c(
    "bundle_schema", "study", "source_artifact", "source_evidence",
    "de_model_action", "longitudinal_design", "longitudinal_contrasts",
    "time_specific_design", "time_specific_contrasts", "prepared_de_columns",
    "prepared_environment_evidence", "resource_route", "public_endpoint"
  ),
  value = c(
    "riaz-prepared-bundle@2", "Riaz et al. GSE91061",
    "s5-riaz-canonical saved prepared project (2026-09-07)",
    "IDENTITY.txt and STEPS.log record stages 06-09; outputs are manifest-verified",
    "reused_saved_DE_outputs_no_DESeq2_refit",
    "DESeq2 NB Wald: ~ patient_id + on_indicator + on_by_ipi + on_by_pd; maxit=10000",
    "responders ON-PRE and PD ON-PRE marginal over prior-ipilimumab strata; interaction=(responders change)-(PD change)",
    "DESeq2 NB Wald at PRE and ON separately: ~ ipi_progressed + response_pd",
    "responders-PD at PRE and ON (negative response_pd coefficient)",
    "symbol, entrez_id, log2FoldChange, stat, pvalue, padj, beta_converged, lisa_eligible, lisa_exclusion_reason",
    "validation_package_versions.tsv; canonical saved sessionInfo_lisa_standard.txt records R 4.5.3, DESeq2 1.50.2, lisaR 0.99.0",
    "ordinary_C1_builtin_core_category_map_and_pinned_TERM2GENE",
    "pending_reviewed_stable_distribution"
  ),
  stringsAsFactors = FALSE
)
utils::write.table(provenance, file.path(bundle_root, "prepared_bundle_provenance.tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
cat("RIAZ_PREPARED_BUNDLE=PASS\n")
