#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Validate the five lisaR inputs before configuration or enrichment.
#
# This script is strict about scientific and structural invariants. Exact
# annotation-universe and FDR counts can change when R or Bioconductor packages
# change, so those values are reported against a versioned reference instead
# of being treated as universal truths. Effect and rank concordance with the
# canonical run are checked quantitatively when an independently authorised
# external benchmark is supplied. Benchmark rows are not distributed by
# lisaR.

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))

analysis_ids <- c(
  "responders_on_vs_pre",
  "pd_on_vs_pre",
  "responders_vs_pd_pre",
  "responders_vs_pd_on",
  "differential_longitudinal_response"
)
input_dir <- riaz_output_dir("lisa_inputs")
design <- riaz_read_tsv(file.path(
  riaz_output_dir("design"), "analysis_sample_design.tsv"
))

riaz_assert(nrow(design) == 52L, "Expected 52 paired samples.")
riaz_assert(length(unique(design$patient_id)) == 26L, "Expected 26 patients.")
riaz_assert(
  length(unique(design$patient_id[design$analysis_response == "CRPR"])) == 8L,
  "Expected eight responder patients."
)
riaz_assert(
  length(unique(design$patient_id[design$analysis_response == "PD"])) == 18L,
  "Expected 18 progressive-disease patients."
)
riaz_assert(
  all(table(design$patient_id) == 2L),
  "Every patient must contribute exactly PRE and ON samples."
)
riaz_assert(
  all(vapply(split(design$visit, design$patient_id), function(x) {
    identical(sort(x), c("On", "Pre"))
  }, logical(1))),
  "Every patient must contribute one PRE and one ON sample."
)

benchmark_dir <- Sys.getenv("LISAR_RIAZ_BENCHMARK_DIR", unset = "")
benchmark_available <- nzchar(benchmark_dir)
benchmark_manifest <- NULL
benchmark_effects <- NULL
if (benchmark_available) {
  expected_hashes <- c(
    manifest = Sys.getenv("LISAR_RIAZ_BENCHMARK_MANIFEST_SHA256", unset = ""),
    effects = Sys.getenv("LISAR_RIAZ_BENCHMARK_EFFECTS_SHA256", unset = "")
  )
  riaz_assert(
    all(grepl("^[0-9a-f]{64}$", expected_hashes)),
    paste(
      "LISAR_RIAZ_BENCHMARK_MANIFEST_SHA256 and",
      "LISAR_RIAZ_BENCHMARK_EFFECTS_SHA256 must contain reviewed lower-case hashes."
    )
  )
  riaz_assert(
    dir.exists(benchmark_dir) && !nzchar(Sys.readlink(benchmark_dir)),
    "LISAR_RIAZ_BENCHMARK_DIR must be one existing non-symlink directory."
  )
  benchmark_dir <- normalizePath(benchmark_dir, mustWork = TRUE)
  benchmark_paths <- c(
    manifest = file.path(benchmark_dir, "canonical-model-benchmark-v1.tsv"),
    effects = file.path(benchmark_dir, "canonical-model-effects-v1.tsv.gz")
  )
  riaz_assert(
    all(file.exists(benchmark_paths)) &&
      !any(dir.exists(benchmark_paths)) &&
      !any(nzchar(Sys.readlink(benchmark_paths))),
    "The external benchmark directory lacks the two regular benchmark files."
  )
  observed_hashes <- vapply(benchmark_paths, riaz_sha256, character(1L))
  riaz_assert(
    identical(unname(observed_hashes), unname(expected_hashes)),
    "An external Riaz benchmark failed its reviewed SHA-256 check."
  )
  benchmark_manifest <- riaz_read_tsv(benchmark_paths[["manifest"]])
  benchmark_effects <- riaz_read_tsv(benchmark_paths[["effects"]])
  riaz_assert(
    setequal(benchmark_manifest$analysis_id, analysis_ids),
    "The external benchmark does not contain the five expected analyses."
  )
  riaz_assert(
    setequal(unique(benchmark_effects$analysis_id), analysis_ids),
    "The external effect benchmark does not contain the five expected analyses."
  )
}

validation <- list()
concordance <- list()

for (analysis_id in analysis_ids) {
  path <- file.path(input_dir, paste0(analysis_id, ".tsv"))
  riaz_assert(file.exists(path), paste("Missing lisaR input:", path))
  x <- riaz_read_tsv(path)
  required <- c(
    "symbol", "entrez_id", "log2FoldChange", "stat", "pvalue", "padj",
    "beta_converged", "lisa_eligible", "lisa_exclusion_reason"
  )
  riaz_assert(
    all(required %in% names(x)),
    paste("Required columns are absent from", analysis_id)
  )
  riaz_assert(
    nrow(x) > 10000L,
    paste("Implausibly small universe:", analysis_id)
  )
  riaz_assert(!anyDuplicated(x$symbol), paste("Duplicated symbols:", analysis_id))
  riaz_assert(
    all(!is.na(x$symbol) & nzchar(x$symbol)),
    paste("Missing gene symbols:", analysis_id)
  )
  riaz_assert(
    sum(is.finite(x$stat)) > 10000L,
    paste("Too few finite ranking values:", analysis_id)
  )

  model_dir <- if (analysis_id %in% c(
    "responders_vs_pd_pre", "responders_vs_pd_on"
  )) "time_specific" else "longitudinal"
  raw <- riaz_read_tsv(file.path(
    riaz_output_dir("models", model_dir),
    paste0(analysis_id, "_symbol_full.tsv")
  ))
  raw_fdr <- sum(raw$padj < 0.05, na.rm = TRUE)

  unstable <- !is.na(raw$beta_converged) & !raw$beta_converged
  unstable_symbols <- raw$symbol[unstable]
  if (length(unstable_symbols)) {
    exported_unstable <- x[x$symbol %in% unstable_symbols, , drop = FALSE]
    riaz_assert(
      all(is.na(exported_unstable$stat)) &&
        all(is.na(exported_unstable$pvalue)) &&
        all(is.na(exported_unstable$padj)) &&
        all(!exported_unstable$lisa_eligible),
      paste("A nonconverged coefficient entered a LISA rank:", analysis_id)
    )
  }

  expected <- if (benchmark_available) {
    benchmark_manifest[
      benchmark_manifest$analysis_id == analysis_id, , drop = FALSE
    ]
  } else {
    data.frame(
      universe_rows = NA_integer_, raw_model_fdr_lt_0_05 = NA_integer_
    )
  }
  validation[[analysis_id]] <- data.frame(
    analysis_id = analysis_id,
    universe_rows = nrow(x),
    benchmark_universe_rows = expected$universe_rows,
    universe_delta = nrow(x) - expected$universe_rows,
    finite_lisa_rank_values = sum(is.finite(x$stat)),
    raw_model_fdr_lt_0_05 = raw_fdr,
    benchmark_raw_model_fdr_lt_0_05 = expected$raw_model_fdr_lt_0_05,
    raw_model_fdr_delta = raw_fdr - expected$raw_model_fdr_lt_0_05,
    lisa_input_fdr_lt_0_05 = sum(x$padj < 0.05, na.rm = TRUE),
    excluded_nonconverged_coefficients = sum(unstable),
    nonconverged_coefficients_in_lisa_rank = sum(
      !x$lisa_eligible & is.finite(x$stat), na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )

  if (benchmark_available) {
    reference <- benchmark_effects[
      benchmark_effects$analysis_id == analysis_id, , drop = FALSE
    ]
    comparison <- riaz_effect_concordance(raw, reference)
    comparison$status <- "pass"
    comparison$analysis_id <- analysis_id
    comparison$minimum_common_symbol_fraction <- 0.98
    comparison$minimum_spearman <- 0.995
    comparison$minimum_directional_concordance <- 0.995
    comparison$common_symbol_fraction <- comparison$common_symbols /
      min(nrow(raw), nrow(reference))
    riaz_assert(
      comparison$common_symbol_fraction >=
        comparison$minimum_common_symbol_fraction,
      paste("Too few genes overlap the canonical benchmark:", analysis_id)
    )
    riaz_assert(
      comparison$lfc_spearman >= comparison$minimum_spearman &&
        comparison$stat_spearman >= comparison$minimum_spearman,
      paste("Effect or ranking concordance is below threshold:", analysis_id)
    )
    riaz_assert(
      comparison$lfc_directional_concordance >=
        comparison$minimum_directional_concordance &&
        comparison$stat_directional_concordance >=
          comparison$minimum_directional_concordance,
      paste("Directional concordance is below threshold:", analysis_id)
    )
  } else {
    comparison <- data.frame(
      status = "not_run_external_benchmark_not_supplied",
      analysis_id = analysis_id,
      common_symbols = NA_integer_,
      lfc_common_finite = NA_integer_, lfc_spearman = NA_real_,
      lfc_directional_concordance = NA_real_,
      stat_common_finite = NA_integer_, stat_spearman = NA_real_,
      stat_directional_concordance = NA_real_,
      minimum_common_symbol_fraction = 0.98,
      minimum_spearman = 0.995,
      minimum_directional_concordance = 0.995,
      common_symbol_fraction = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  concordance[[analysis_id]] <- comparison
}

validation <- do.call(rbind, validation)
riaz_assert(
  length(unique(validation$universe_rows)) == 1L,
  "The five exported analyses do not share one gene universe."
)
riaz_assert(
  all(validation$nonconverged_coefficients_in_lisa_rank == 0L),
  "At least one nonconverged coefficient remains in a LISA ranking."
)

# The longitudinal interaction must be Responders minus PD. This algebraic
# check is independent of labels, plot subtitles, and the reference benchmark.
longitudinal <- lapply(
  c(
    responders = "responders_on_vs_pre",
    pd = "pd_on_vs_pre",
    differential = "differential_longitudinal_response"
  ),
  function(id) {
    riaz_read_tsv(file.path(
      riaz_output_dir("models", "longitudinal"),
      paste0(id, "_symbol_full.tsv")
    ))[, c("symbol", "log2FoldChange"), drop = FALSE]
  }
)
orientation <- Reduce(
  function(left, right) merge(left, right, by = "symbol"),
  longitudinal
)
names(orientation) <- c(
  "symbol", "responders_log2fc", "pd_log2fc", "differential_log2fc"
)
orientation_error <- orientation$differential_log2fc -
  (orientation$responders_log2fc - orientation$pd_log2fc)
riaz_assert(
  max(abs(orientation_error), na.rm = TRUE) < 1e-10,
  "The longitudinal interaction is not Responders minus PD."
)

concordance <- do.call(rbind, concordance)
concordance <- concordance[, c(
  "analysis_id", setdiff(names(concordance), "analysis_id")
)]
riaz_write_tsv(
  validation,
  file.path(input_dir, "lisa_input_validation.tsv")
)
riaz_write_tsv(
  concordance,
  file.path(input_dir, "canonical_effect_concordance.tsv")
)
riaz_write_package_manifest(file.path(
  input_dir, "validation_package_versions.tsv"
))

manifest_paths <- c(
  file.path(input_dir, paste0(analysis_ids, ".tsv")),
  file.path(input_dir, "paired_log2_normalized_counts_by_symbol.tsv"),
  file.path(input_dir, "paired_sample_annotation.tsv"),
  file.path(input_dir, "lisa_input_validation.tsv"),
  file.path(input_dir, "canonical_effect_concordance.tsv"),
  file.path(input_dir, "validation_package_versions.tsv")
)
manifest <- data.frame(
  # Store project-relative paths.  Stage 11 runs this script in an atomic
  # sibling staging directory and later promotes that directory, so absolute
  # paths would become stale even though the verified bytes did not change.
  path = file.path("outputs", "lisa_inputs", basename(manifest_paths)),
  bytes = file.info(manifest_paths)$size,
  sha256 = unname(tools::md5sum(manifest_paths)),
  stringsAsFactors = FALSE
)
# The portable base-R checksum above is MD5 and is labelled explicitly. The
# project-level Chainrunner manifest adds SHA-256 using the operating system.
names(manifest)[names(manifest) == "sha256"] <- "md5"
riaz_write_tsv(manifest, file.path(input_dir, "lisa_input_manifest.tsv"))

print(validation, row.names = FALSE)
print(concordance, row.names = FALSE)
cat(sprintf(
  "RIAZ_CANONICAL_CONCORDANCE=%s\n",
  if (benchmark_available) "PASS" else "NOT_RUN_EXTERNAL_BENCHMARK_NOT_SUPPLIED"
))
cat("LISA_INPUT_GATE=PASS\n")
