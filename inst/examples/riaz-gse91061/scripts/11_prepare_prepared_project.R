#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Install a reviewed, prepared Riaz input bundle into one new project.
#
# This script accepts a local directory bundle only.  It never extracts an
# archive, fetches arbitrary code, or refits DESeq2.  A future public endpoint
# must be reviewed and pinned separately; do not substitute an unreviewed URL
# for LISAR_RIAZ_PREPARED_BUNDLE.

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))

riaz_require_packages(c("data.table", "lisaR"))

bundle_root <- Sys.getenv("LISAR_RIAZ_PREPARED_BUNDLE", unset = "")
riaz_assert(nzchar(bundle_root), paste(
  "LISAR_RIAZ_PREPARED_BUNDLE is required and must name a reviewed local",
  "directory bundle. No public prepared endpoint is configured yet."
))
riaz_assert(dir.exists(bundle_root) && !riaz_path_is_link(bundle_root),
            "LISAR_RIAZ_PREPARED_BUNDLE must be one existing non-symlink directory.")
bundle_root <- normalizePath(bundle_root, winslash = "/", mustWork = TRUE)
manifest_path <- file.path(bundle_root, "prepared_bundle_manifest.tsv")
provenance_path <- file.path(bundle_root, "prepared_bundle_provenance.tsv")
riaz_assert(
  file.exists(manifest_path) && !dir.exists(manifest_path) &&
    isTRUE(utils::file_test("-f", manifest_path)) && !riaz_path_is_link(manifest_path),
  "Prepared bundle lacks one regular prepared_bundle_manifest.tsv."
)
riaz_assert(
  file.exists(provenance_path) && !dir.exists(provenance_path) &&
    isTRUE(utils::file_test("-f", provenance_path)) && !riaz_path_is_link(provenance_path),
  "Prepared bundle lacks one regular prepared_bundle_provenance.tsv."
)

manifest <- riaz_read_tsv(manifest_path)
required_columns <- c("path", "bytes", "sha256", "role")
riaz_assert(
  all(required_columns %in% names(manifest)) && nrow(manifest) >= 13L,
  "Prepared bundle manifest is malformed or implausibly small."
)
manifest <- manifest[, required_columns, drop = FALSE]
manifest$path <- vapply(manifest$path, riaz_safe_relative_path, character(1))
riaz_assert(!anyDuplicated(manifest$path),
            "Prepared bundle manifest has duplicate paths.")
riaz_assert(
  all(grepl("^[0-9a-f]{64}$", manifest$sha256)) &&
    all(is.finite(manifest$bytes) & manifest$bytes > 0),
  "Prepared bundle manifest must contain reviewed lower-case SHA-256 and bytes."
)
allowed_prefixes <- c("outputs/design/", "outputs/models/", "outputs/lisa_inputs/")
riaz_assert(
  all(vapply(manifest$path, function(path) {
    any(startsWith(path, allowed_prefixes)) && grepl("\\.tsv$", path)
  }, logical(1))),
  "Prepared bundle manifest may contain only reviewed Riaz TSV inputs, never executable code."
)
required_paths <- c(
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
riaz_assert(all(required_paths %in% manifest$path), paste(
  "Prepared bundle is missing required design, model-provenance, or LISA-input",
  "TSVs. Rebuild it from the reviewed canonical project."
))

bundle_paths <- file.path(bundle_root, manifest$path)
riaz_assert(
  all(file.exists(bundle_paths)) && !any(dir.exists(bundle_paths)) &&
    all(vapply(bundle_paths, function(path) {
      isTRUE(utils::file_test("-f", path))
    }, logical(1))) &&
    !any(vapply(bundle_paths, riaz_path_is_link, logical(1))),
  "Prepared bundle has a missing, nonregular, or symbolic-link input."
)
observed_bytes <- as.numeric(file.info(bundle_paths)$size)
observed_sha256 <- vapply(bundle_paths, riaz_sha256, character(1))
riaz_assert(
  identical(observed_bytes, as.numeric(manifest$bytes)) &&
    identical(unname(observed_sha256), as.character(manifest$sha256)),
  "Prepared bundle inputs do not match their reviewed manifest; no project was created."
)
manifest_sha256 <- riaz_sha256(manifest_path)

riaz_prepared_next_command <- function(project, resource_cache) {
  cat(
    "Next, run: LISAR_RIAZ_PROJECT_DIR=", shQuote(project),
    " LISAR_RIAZ_RESOURCE_CACHE=", shQuote(resource_cache),
    " Rscript ", shQuote(file.path(
      riaz_example_dir(), "09_run_lisa_example.R"
    )), "\n", sep = ""
  )
}

riaz_validate_promoted_project <- function(project, manifest, manifest_sha256) {
  compact <- file.exists(file.path(project, "provenance", "path_map.tsv"))
  receipt_path <- file.path(project, "logs", "riaz_prepared_receipt.tsv")
  riaz_assert(
    dir.exists(project) && !riaz_path_is_link(project) &&
      file.exists(receipt_path) && !dir.exists(receipt_path) &&
      isTRUE(utils::file_test("-f", receipt_path)) && !riaz_path_is_link(receipt_path),
    "LISAR_RIAZ_PROJECT_DIR already exists and is not this completed prepared project; it will not be overwritten."
  )
  receipt <- riaz_read_tsv(receipt_path)
  required_receipt <- c(
    "bundle_manifest_sha256", "prepared_files", "resource_cache", "next_stage",
    "gsea_or_report_run", "public_prepared_endpoint"
  )
  riaz_assert(
    nrow(receipt) == 1L && all(required_receipt %in% names(receipt)) &&
      identical(receipt$bundle_manifest_sha256[[1]], manifest_sha256) &&
      identical(as.integer(receipt$prepared_files[[1]]), nrow(manifest)) &&
      identical(receipt$next_stage[[1]], "09_run_lisa_example.R") &&
      identical(as.logical(receipt$gsea_or_report_run[[1]]), FALSE) &&
      receipt$public_prepared_endpoint[[1]] %in% c(
        "pending_reviewed_stable_distribution",
        get("lisa_example_distribution_contract",
            envir = asNamespace("lisaR"))("riaz-gse91061")$source_url) &&
      nzchar(receipt$resource_cache[[1]], keepNA = FALSE) &&
      dir.exists(receipt$resource_cache[[1]]) &&
      !riaz_path_is_link(receipt$resource_cache[[1]]),
    "The existing prepared-project receipt is incomplete, changed, or has an unavailable normal resource cache; it will not be overwritten."
  )

  saved_manifest <- file.path(project, "provenance", "prepared_bundle_manifest.tsv")
  riaz_assert(
    file.exists(saved_manifest) && !dir.exists(saved_manifest) &&
      isTRUE(utils::file_test("-f", saved_manifest)) && !riaz_path_is_link(saved_manifest) &&
      identical(riaz_sha256(saved_manifest), manifest_sha256),
    "The promoted project no longer contains the verified prepared-bundle manifest; it will not be overwritten."
  )
  # Stage 06 deliberately regenerates its four validation/manifest receipts
  # from the copied scientific inputs.  The design, five model outputs and
  # seven raw LISA inputs remain the immutable promoted bundle content; the
  # generated receipts are checked separately below for their own contract.
  immutable_manifest <- manifest[manifest$path %in% required_paths, , drop = FALSE]
  riaz_assert(
    nrow(immutable_manifest) == length(required_paths),
    "The verified prepared bundle no longer identifies all immutable Riaz inputs; it will not be overwritten."
  )
  promoted_rel <- immutable_manifest$path
  if (compact) {
    get("lisa_verify_compact_example_project", asNamespace("lisaR"))(project)
    map <- riaz_read_tsv(file.path(project, "provenance", "path_map.tsv"))
    indices <- match(promoted_rel, map$source)
    riaz_assert(!anyNA(indices), "Incomplete mapping of immutable Riaz inputs.")
    promoted_rel <- map$path[indices]
  }
  promoted_paths <- file.path(project, promoted_rel)
  riaz_assert(
    all(file.exists(promoted_paths)) && !any(dir.exists(promoted_paths)) &&
      all(vapply(promoted_paths, function(path) isTRUE(utils::file_test("-f", path)), logical(1))) &&
      !any(vapply(promoted_paths, riaz_path_is_link, logical(1))) &&
      identical(as.numeric(file.info(promoted_paths)$size), as.numeric(immutable_manifest$bytes)) &&
      identical(unname(vapply(promoted_paths, riaz_sha256, character(1))),
                as.character(immutable_manifest$sha256)),
    "The promoted project differs from its verified prepared bundle; it will not be overwritten."
  )

  if (compact) {
    resolution <- riaz_read_tsv(file.path(project, "logs", "lisa_resource_resolution.tsv"))
    context <- riaz_read_tsv(file.path(project, "provenance", "resource_cache.tsv"))
    riaz_assert(nrow(resolution) == 3L &&
      all(resolution$resource_cache == receipt$resource_cache[[1L]]) &&
      identical(context$resource_cache[[1L]], receipt$resource_cache[[1L]]),
      "The resource receipts disagree; the existing project will not be overwritten.")
    return(receipt)
  }

  config_manifest_path <- file.path(project, "config", "config_manifest.tsv")
  input_manifest_path <- file.path(project, "outputs", "lisa_inputs", "lisa_input_manifest.tsv")
  resource_receipt_path <- file.path(project, "logs", "lisa_resource_resolution.tsv")
  riaz_assert(
    all(file.exists(c(config_manifest_path, input_manifest_path, resource_receipt_path))) &&
      !any(vapply(c(config_manifest_path, input_manifest_path, resource_receipt_path),
                  riaz_path_is_link, logical(1))),
    "The promoted project is missing a stage 06--08 receipt; it will not be overwritten."
  )
  config_manifest <- riaz_read_tsv(config_manifest_path)
  input_manifest <- riaz_read_tsv(input_manifest_path)
  resource_receipt <- riaz_read_tsv(resource_receipt_path)
  riaz_assert(
    all(c("format", "path", "sha256") %in% names(config_manifest)) &&
      identical(config_manifest$path, c("config/riaz-gse91061.yml", "config/riaz-gse91061.json")) &&
      all(vapply(file.path(project, config_manifest$path), riaz_sha256, character(1)) ==
          config_manifest$sha256) &&
      all(c("path", "bytes", "md5") %in% names(input_manifest)) &&
      all(startsWith(input_manifest$path, "outputs/lisa_inputs/")) &&
      all(file.exists(file.path(project, input_manifest$path))) &&
      nrow(resource_receipt) == 3L && "resource_cache" %in% names(resource_receipt) &&
      length(unique(resource_receipt$resource_cache)) == 1L &&
      identical(resource_receipt$resource_cache[[1]], receipt$resource_cache[[1]]),
    "The promoted project has stale or changed stage 06--08 receipts; it will not be overwritten."
  )
  receipt
}

requested_project <- Sys.getenv("LISAR_RIAZ_PROJECT_DIR", unset = "")
riaz_assert(nzchar(requested_project), "LISAR_RIAZ_PROJECT_DIR is required.")
project_parent <- dirname(path.expand(requested_project))
riaz_assert(dir.exists(project_parent) && !riaz_path_is_link(project_parent),
            "LISAR_RIAZ_PROJECT_DIR must have an existing non-symlink parent.")
project_target <- file.path(
  normalizePath(project_parent, winslash = "/", mustWork = TRUE),
  basename(path.expand(requested_project))
)
riaz_assert(!(basename(project_target) %in% c("", ".", "..")),
            "LISAR_RIAZ_PROJECT_DIR must name a new project directory.")

if (riaz_path_entry_exists(project_target)) {
  receipt <- riaz_validate_promoted_project(
    project_target, manifest, manifest_sha256
  )
  cat("RIAZ_PREPARED_PROJECT=ALREADY_PREPARED\n")
  riaz_prepared_next_command(project_target, receipt$resource_cache[[1]])
  quit(status = 0L)
}

stage <- file.path(
  dirname(project_target),
  paste0(".", basename(project_target), ".prepared-staging-", Sys.getpid())
)
riaz_assert(!riaz_path_entry_exists(stage),
            "Prepared-project staging path already exists; inspect it and retry with a new target.")
riaz_assert(dir.create(stage, mode = "0700"),
            "Could not create the private prepared-project staging directory.")
on.exit(if (dir.exists(stage)) {
  message("Prepared-project staging was retained for inspection: ", stage)
}, add = TRUE)

for (i in seq_along(bundle_paths)) {
  target <- file.path(stage, manifest$path[[i]])
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE, mode = "0700")
  riaz_assert(file.copy(bundle_paths[[i]], target, overwrite = FALSE, copy.mode = TRUE),
              paste("Could not stage prepared input:", manifest$path[[i]]))
}
staged_paths <- file.path(stage, manifest$path)
riaz_assert(
  identical(as.numeric(file.info(staged_paths)$size), as.numeric(manifest$bytes)) &&
    identical(unname(vapply(staged_paths, riaz_sha256, character(1))),
              as.character(manifest$sha256)),
  "Prepared inputs changed during staging; project promotion was refused."
)
dir.create(file.path(stage, "provenance"), mode = "0700")
riaz_assert(file.copy(manifest_path, file.path(stage, "provenance"), overwrite = FALSE),
            "Could not save the verified prepared-bundle manifest.")
riaz_assert(file.copy(provenance_path, file.path(stage, "provenance"), overwrite = FALSE),
            "Could not save prepared-bundle provenance.")

Sys.setenv(LISAR_RIAZ_PROJECT_DIR = stage)
riaz_run_script("06_validate_exported_inputs.R")
riaz_run_script("07_prepare_lisa_config.R")
riaz_run_script("08_prepare_lisa_resources.R")
resource_receipt <- riaz_read_tsv(file.path(
  stage, "logs", "lisa_resource_resolution.tsv"
))
riaz_assert(
  nrow(resource_receipt) == 3L && "resource_cache" %in% names(resource_receipt) &&
    length(unique(resource_receipt$resource_cache)) == 1L,
  "The ordinary C1 resource receipt is incomplete."
)
riaz_write_tsv(data.frame(
  bundle_manifest_sha256 = manifest_sha256,
  prepared_bundle = bundle_root,
  prepared_files = nrow(manifest),
  resource_cache = resource_receipt$resource_cache[[1]],
  next_stage = "09_run_lisa_example.R",
  gsea_or_report_run = FALSE,
  public_prepared_endpoint = get("lisa_example_distribution_contract",
    envir = asNamespace("lisaR"))("riaz-gse91061")$source_url,
  stringsAsFactors = FALSE
), file.path(stage, "logs", "riaz_prepared_receipt.tsv"))

# Validate the complete staged object before the atomic promotion.  All
# stage-06/07 paths are intentionally project-relative, so this same contract
# remains true after `file.rename()`.
riaz_validate_promoted_project(stage, manifest, manifest_sha256)
get("lisa_compact_example_project", asNamespace("lisaR"))(
  stage, "riaz", file.path(stage, "config", "riaz-gse91061.yml"),
  resource_cache = resource_receipt$resource_cache[[1L]])

riaz_assert(file.rename(stage, project_target), paste(
  "Prepared project passed validation but could not be promoted atomically;",
  "the intact staging directory was retained:", stage
))
riaz_validate_promoted_project(project_target, manifest, manifest_sha256)
cat("RIAZ_PREPARED_PROJECT=PASS\n")
riaz_prepared_next_command(project_target, resource_receipt$resource_cache[[1]])
