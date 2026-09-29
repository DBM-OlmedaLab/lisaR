#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Validate, plan, dry-run, execute and verify the bounded Riaz lisaR workflow.
#
# The canonical configuration is never edited in place. A separate planning
# configuration is written with dry_run: true. Scientific computation starts
# only after validation, the output estimate, the fixed resource budget, and a
# planning artifact have all passed their gates.

source(file.path(dirname(normalizePath(sub(
  "^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
))), "_common.R"))

project_dir <- riaz_project_dir()
local_library <- file.path(project_dir, "rlib")
if (dir.exists(local_library)) {
  .libPaths(c(local_library, .libPaths()))
}
riaz_require_packages(c("data.table", "jsonlite", "lisaR", "yaml"))

active_version <- utils::packageVersion("lisaR")
riaz_assert(
  active_version >= "1.0.0" && active_version < "2.0.0",
  paste("This example requires lisaR 1.x (>= 1.0.0 and < 2.0.0); found", active_version)
)
expected_api <- c(
  "lisa_init_project", "validate_lisa_config", "plan_lisa_outputs",
  "run_lisa", "verify_lisa_run", "plan_lisa_extension",
  "render_lisa_categories", "run_lisa_de", "run_lisa_contrast"
)
riaz_assert(
  all(expected_api %in% getNamespaceExports("lisaR")),
  "The active lisaR installation does not expose the required public API."
)
active_package_path <- normalizePath(
  find.package("lisaR"), winslash = "/", mustWork = TRUE
)
active_description_path <- file.path(active_package_path, "DESCRIPTION")
riaz_assert(
  file.exists(active_description_path),
  "The active lisaR installation has no DESCRIPTION file."
)
active_description_sha256 <- riaz_sha256(active_description_path)

resource_cache <- riaz_configure_normal_resources()

source_config_path <- file.path(project_dir, "config", "riaz-gse91061.yml")
if (file.exists(file.path(project_dir, "study.yml")))
  source_config_path <- file.path(project_dir, "study.yml")
riaz_assert(
  file.exists(source_config_path) &&
    !dir.exists(source_config_path) &&
    isTRUE(utils::file_test("-f", source_config_path)) &&
    !riaz_path_is_link(source_config_path),
  paste(
    "Run 07_prepare_lisa_config.R first; its YAML must be one regular",
    "non-symlink file."
  )
)
source_config_path <- normalizePath(
  source_config_path, winslash = "/", mustWork = TRUE
)
riaz_assert(
  riaz_path_is_within(source_config_path, project_dir),
  "The runtime configuration must resolve inside LISAR_RIAZ_PROJECT_DIR."
)
source_config_sha256 <- riaz_sha256(source_config_path)
config_path <- file.path(
  dirname(source_config_path), "exec.yml"
)
riaz_assert(
  !riaz_path_entry_exists(config_path),
  paste(
    "The execution snapshot already exists; preserve this interrupted project",
    "and start from a new project directory:", config_path
  )
)
riaz_assert(
  file.copy(source_config_path, config_path, overwrite = FALSE),
  "Could not create the execution-configuration snapshot."
)
config_path <- normalizePath(config_path, winslash = "/", mustWork = TRUE)
execution_config_sha256 <- riaz_sha256(config_path)
riaz_assert(
  identical(source_config_sha256, execution_config_sha256),
  "The execution snapshot differs from the reviewed source configuration."
)
configuration <- yaml::read_yaml(config_path)
planned_output <- normalizePath(
  file.path(dirname(config_path), configuration$pipeline$output_dir),
  winslash = "/", mustWork = FALSE
)
riaz_assert(
  !riaz_path_entry_exists(planned_output),
  paste("Choose a new output directory; the configured destination exists:", planned_output)
)

validation <- lisaR::validate_lisa_config(
  config_path, check_files = TRUE, strict = TRUE
)
riaz_assert(isTRUE(validation$execution_ready), "The Riaz configuration is not execution-ready.")
riaz_assert(identical(validation$schema_version, "1.0.0"), "Expected configuration schema 1.0.0.")
prepared_tutorial <- file.exists(file.path(project_dir, "provenance", "path_map.tsv"))
riaz_assert(validation$analyses == if (prepared_tutorial) 2L else 5L,
            "Unexpected number of analyses for the selected Riaz route.")
riaz_assert(validation$contrasts == if (prepared_tutorial) 1L else 2L,
            "Unexpected number of contrasts for the selected Riaz route.")
if (prepared_tutorial) {
  riaz_assert(identical(vapply(validation$config$single_de, `[[`, character(1), "analysis_id"),
    c("responders_vs_pd_pre", "responders_vs_pd_on")), "The tutorial requires PRE and ON only.")
}
riaz_assert(
  identical(
    unlist(validation$config$collections, use.names = FALSE),
    c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS")
  ),
  "The Riaz example must use the four semantic collections in canonical order."
)
riaz_assert(
  identical(as.integer(validation$config$pipeline$workers), 4L),
  "The canonical Riaz example must exercise the 4-worker parallel path."
)
riaz_assert(
  identical(as.numeric(validation$config$pipeline$gsea_padj_cutoff), 0.25),
  "The canonical Riaz example must fix gsea_padj_cutoff at 0.25."
)
riaz_assert(
  isFALSE(validation$config$pipeline$run_hallmarks) &&
    isFALSE(validation$config$pipeline$run_kegg_maps),
  "HALLMARKS and KEGG must remain excluded from this bounded example."
)
riaz_assert(
  identical(validation$config$report$mode, "standard"),
  "The canonical Riaz example must use report.mode = standard."
)
riaz_formats <- validation$config$report$formats[c("png", "svg", "pdf")]
riaz_assert(
  identical(riaz_formats, list(png = TRUE, svg = FALSE, pdf = FALSE)),
  "The standard Riaz report must request PNG only."
)
riaz_assert(
  isTRUE(validation$config$report$source_data) &&
    isTRUE(validation$config$report$recipes),
  "The standard Riaz report must retain source tables and recipes."
)

normalize_resource_identity <- function(resources) {
  columns <- c("config_key", "resource_id", "path", "sha256")
  riaz_assert(
    all(columns %in% names(resources)) && nrow(resources) == 3L &&
      all(resources$checked %in% TRUE) && all(resources$ready %in% TRUE),
    "Strict validation did not resolve three checked, ready resources."
  )
  identity <- resources[, columns, drop = FALSE]
  identity$path <- vapply(
    identity$path, normalizePath, character(1),
    winslash = "/", mustWork = TRUE
  )
  identity <- identity[order(identity$config_key), , drop = FALSE]
  rownames(identity) <- NULL
  identity
}
initial_resource_identity <- normalize_resource_identity(validation$resources)

logs <- file.path(project_dir, "logs")
dir.create(logs, recursive = TRUE, showWarnings = FALSE)
riaz_write_tsv(
  data.frame(
    lisaR_version = as.character(active_version),
    package_path = active_package_path,
    description_sha256 = active_description_sha256,
    stringsAsFactors = FALSE
  ),
  file.path(logs, "lisa_standard_package_identity.tsv")
)
riaz_write_tsv(
  validation$readiness,
  file.path(logs, "lisa_standard_validation_readiness.tsv")
)
riaz_write_tsv(
  validation$inputs,
  file.path(logs, "lisa_standard_validation_inputs.tsv")
)
riaz_write_tsv(
  validation$dependencies,
  file.path(logs, "lisa_standard_validation_dependencies.tsv")
)
riaz_write_tsv(
  validation$resources,
  file.path(logs, "lisa_standard_validation_resources.tsv")
)
riaz_write_tsv(
  initial_resource_identity,
  file.path(logs, "lisa_standard_resource_identity.tsv")
)

resource_audit_path <- file.path(
  project_dir, "logs", "lisa_resource_resolution.tsv"
)
riaz_assert(
  file.exists(resource_audit_path),
  "Run 08_prepare_lisa_resources.R before planning the Riaz output."
)
resource_audit <- riaz_read_tsv(resource_audit_path)
resolved_dictionary <- validation$resources[
  validation$resources$config_key == "dictionary_resource",
  , drop = FALSE
]
riaz_assert(
  nrow(resolved_dictionary) == 1L &&
    isTRUE(resolved_dictionary$checked[[1]]) &&
    isTRUE(resolved_dictionary$ready[[1]]) &&
    file.exists(resolved_dictionary$path[[1]]),
  "Strict validation did not resolve one ready dictionary resource."
)
dictionary_row <- resource_audit[
  resource_audit$resource_id == resolved_dictionary$resource_id[[1]],
  , drop = FALSE
]
riaz_assert(
  nrow(dictionary_row) == 1L &&
    identical(
      normalizePath(
        dictionary_row$path[[1]], winslash = "/", mustWork = TRUE
      ),
      normalizePath(
        resolved_dictionary$path[[1]], winslash = "/", mustWork = TRUE
      )
    ) &&
    identical(dictionary_row$sha256[[1]], resolved_dictionary$sha256[[1]]),
  paste(
    "The resource-installation audit does not match the dictionary resolved",
    "by strict configuration validation."
  )
)
riaz_assert(
  identical(
    riaz_sha256(resolved_dictionary$path[[1]]),
    resolved_dictionary$sha256[[1]]
  ),
  "The resolved installed dictionary changed after strict validation."
)
planning_dictionary <- riaz_read_tsv(resolved_dictionary$path[[1]])
riaz_assert(
  all(c("universe", "category_id") %in% names(planning_dictionary)),
  "The installed dictionary lacks universe/category_id planning fields."
)
category_counts <- vapply(
  unlist(validation$config$collections, use.names = FALSE),
  function(collection) length(unique(
    planning_dictionary$category_id[
      planning_dictionary$universe == collection
    ]
  )),
  integer(1)
)
riaz_assert(
  all(category_counts > 0L),
  "The installed dictionary has no categories for a configured collection."
)
riaz_write_tsv(
  data.frame(
    collection = names(category_counts),
    expected_categories = unname(category_counts),
    dictionary_resource = resolved_dictionary$resource_id[[1]],
    dictionary_sha256 = resolved_dictionary$sha256[[1]],
    stringsAsFactors = FALSE
  ),
  file.path(logs, "lisa_standard_planning_category_counts.tsv")
)

output_plan <- lisaR::plan_lisa_outputs(
  config_path, category_counts = category_counts
)
riaz_write_tsv(
  output_plan$summary,
  file.path(logs, "lisa_standard_output_plan_summary.tsv")
)
riaz_write_tsv(
  output_plan$collections,
  file.path(logs, "lisa_standard_output_plan_collections.tsv")
)
riaz_write_tsv(
  output_plan$products,
  file.path(logs, "lisa_standard_output_plan_products.tsv")
)
riaz_write_tsv(
  output_plan$policy,
  file.path(logs, "lisa_standard_output_plan_policy.tsv")
)

max_planned_files <- 20000L
# The standard presentation includes member plots, four NES variants and
# paired contrast evidence with all/same/opposite subsets. Use the recorded
# estimate and leave bounded room for the shared gene index.
# Apply the same ceiling to both the plan and the completed directory below.
max_planned_bytes <- 4 * 1024^3
planned_figure_files <- as.numeric(output_plan$summary$expected_figure_files[[1]])
format_multiplier <- as.numeric(output_plan$summary$format_multiplier[[1]])
riaz_assert(
  is.finite(format_multiplier) && format_multiplier >= 1,
  "The output plan has no valid figure-format multiplier."
)
planned_units <- planned_figure_files / format_multiplier
planned_report_artifact_files <- planned_units * (
  format_multiplier +
    as.integer(validation$config$report$source_data) +
    as.integer(validation$config$report$recipes)
)
# Reserve room for canonical tables, manifests, HTML, ledgers and copied
# allowlisted inputs that are not represented as figure units by the planner.
planned_nonfigure_reserve <- 2000L
planned_files_with_reserve <- planned_report_artifact_files +
  planned_nonfigure_reserve
planned_bytes <- as.numeric(output_plan$summary$conservative_bytes[[1]])
riaz_assert(
  is.finite(planned_files_with_reserve) &&
    planned_files_with_reserve <= max_planned_files,
  paste(
    "Planned report artifacts plus the fixed non-figure reserve exceed the",
    "20,000-file Riaz safety ceiling:", planned_files_with_reserve
  )
)
riaz_assert(
  is.finite(planned_bytes) && planned_bytes <= max_planned_bytes,
  paste("Planned output exceeds the 4 GiB Riaz safety ceiling:", planned_bytes)
)

preflight_receipt <- data.frame(
  lisaR_version = as.character(active_version),
  schema_version = validation$schema_version,
  workers = 4L,
  report_mode = "standard",
  collections = paste(
    unlist(validation$config$collections, use.names = FALSE), collapse = ";"
  ),
  hallmarks_requested = FALSE,
  kegg_requested = FALSE,
  plot_formats = "png",
  figure_source_data = TRUE,
  figure_recipes = TRUE,
  canonical_table_format = "tsv",
  planned_figure_files = planned_figure_files,
  planned_report_artifact_files = planned_report_artifact_files,
  planned_nonfigure_reserve = planned_nonfigure_reserve,
  planned_files_with_reserve = planned_files_with_reserve,
  max_planned_files = max_planned_files,
  planned_conservative_bytes = planned_bytes,
  max_planned_bytes = max_planned_bytes,
  config_sha256 = execution_config_sha256,
  status = "PASS",
  stringsAsFactors = FALSE
)
riaz_write_tsv(
  preflight_receipt,
  file.path(logs, "lisa_standard_preflight_receipt.tsv")
)

# Preserve the execution configuration and create a separate planning file.
plan_configuration <- configuration
plan_configuration$pipeline$dry_run <- TRUE
plan_config_path <- file.path(
  dirname(config_path), "plan.yml"
)
riaz_assert(
  !riaz_path_entry_exists(plan_config_path),
  paste(
    "The planning configuration already exists; preserve this interrupted",
    "project and start from a new project directory:", plan_config_path
  )
)
yaml::write_yaml(plan_configuration, plan_config_path)
plan_config_path <- normalizePath(plan_config_path, winslash = "/", mustWork = TRUE)
riaz_write_tsv(
  data.frame(
    role = c("source", "execution_snapshot", "planning"),
    path = c(source_config_path, config_path, plan_config_path),
    dry_run = c(FALSE, FALSE, TRUE),
    sha256 = c(
      source_config_sha256, execution_config_sha256,
      riaz_sha256(plan_config_path)
    ),
    stringsAsFactors = FALSE
  ),
  file.path(logs, "lisa_standard_config_manifest.tsv")
)

plan_validation <- lisaR::validate_lisa_config(
  plan_config_path, check_files = TRUE, strict = TRUE
)
riaz_assert(
  isTRUE(plan_validation$execution_ready),
  "The derived planning configuration is not execution-ready."
)
riaz_assert(
  isTRUE(plan_validation$config$pipeline$dry_run),
  "The derived planning configuration must set dry_run: true."
)
plan_result <- lisaR::run_lisa(plan_config_path)
saveRDS(plan_result, file.path(logs, "lisa_standard_plan_result.rds"))
writeLines(
  capture.output(str(plan_result, max.level = 5)),
  file.path(logs, "lisa_standard_plan_result.txt")
)
riaz_assert(identical(plan_result$gate, "PLAN_PASS"), "The dry run did not return PLAN_PASS.")
riaz_assert(
  !riaz_path_entry_exists(plan_result$output_dir),
  "The dry run created the intended scientific destination."
)
riaz_assert(dir.exists(plan_result$plan_dir), "The dry-run planning directory is absent.")
plan_verification <- lisaR::verify_lisa_run(plan_result$plan_dir)
saveRDS(
  plan_verification,
  file.path(logs, "lisa_standard_plan_verification.rds")
)
writeLines(
  capture.output(str(plan_verification, max.level = 5)),
  file.path(logs, "lisa_standard_plan_verification.txt")
)
riaz_assert(
  identical(plan_verification$gate, "FAIL") &&
    identical(plan_verification$artifact_type, "plan") &&
    "planning_artifact_not_scientific_run" %in% plan_verification$findings,
  "The dry-run directory was not recognised as a non-scientific plan."
)
cat("LISA_PLAN_GATE=PASS\n")

# Revalidate the preserved execution configuration after planning, then run.
riaz_assert(
  identical(riaz_sha256(config_path), execution_config_sha256),
  paste(
    "The execution snapshot changed after planning; preserve the evidence",
    "and restart in a new project directory."
  )
)
riaz_assert(
  identical(
    normalizePath(find.package("lisaR"), winslash = "/", mustWork = TRUE),
    active_package_path
  ) &&
    identical(utils::packageVersion("lisaR"), active_version) &&
    identical(riaz_sha256(active_description_path), active_description_sha256),
  "The active lisaR installation changed after planning."
)
execution_validation <- lisaR::validate_lisa_config(
  config_path, check_files = TRUE, strict = TRUE
)
riaz_assert(
  isTRUE(execution_validation$execution_ready) &&
    isFALSE(execution_validation$config$pipeline$dry_run),
  "The preserved execution configuration is not ready for scientific computation."
)
execution_resource_identity <- normalize_resource_identity(
  execution_validation$resources
)
riaz_assert(
  identical(execution_resource_identity, initial_resource_identity),
  paste(
    "The resolved resource identity changed after planning; preserve the",
    "evidence and restart in a new project directory."
  )
)
riaz_assert(
  all(vapply(seq_len(nrow(execution_resource_identity)), function(i) {
    identical(
      riaz_sha256(execution_resource_identity$path[[i]]),
      execution_resource_identity$sha256[[i]]
    )
  }, logical(1))),
  "A resolved scientific resource changed after execution validation."
)
riaz_assert(
  !riaz_path_entry_exists(planned_output),
  "The scientific destination appeared after planning; choose a new destination."
)

result <- lisaR::run_lisa(config_path)
saveRDS(result, file.path(logs, "lisa_standard_result.rds"))
writeLines(
  capture.output(str(result, max.level = 5)),
  file.path(logs, "lisa_standard_result.txt")
)
riaz_assert(
  identical(result$gate, "PASS"),
  "The transactional lisaR run did not finish with gate PASS."
)

verification <- lisaR::verify_lisa_run(result$output_dir)
saveRDS(
  verification,
  file.path(logs, "lisa_standard_verification.rds")
)
writeLines(
  capture.output(str(verification, max.level = 5)),
  file.path(logs, "lisa_standard_verification.txt")
)
riaz_assert(
  identical(verification$gate, "PASS") &&
    identical(verification$artifact_type, "scientific_run"),
  "The completed result did not verify as a PASS scientific_run."
)

output_files <- list.files(
  result$output_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE,
  no.. = TRUE
)
output_files <- output_files[file.info(output_files)$isdir %in% FALSE]
actual_files <- length(output_files)
actual_bytes <- sum(as.numeric(file.info(output_files)$size))
riaz_assert(
  actual_files <= max_planned_files,
  paste("Completed output exceeds the 20,000-file safety ceiling:", actual_files)
)
riaz_assert(
  actual_bytes <= max_planned_bytes,
  paste("Completed output exceeds the 4 GiB safety ceiling:", actual_bytes)
)
riaz_write_tsv(
  data.frame(
    output_dir = normalizePath(result$output_dir, winslash = "/"),
    gate = verification$gate,
    artifact_type = verification$artifact_type,
    files = actual_files,
    bytes = actual_bytes,
    max_files = max_planned_files,
    max_bytes = max_planned_bytes,
    stringsAsFactors = FALSE
  ),
  file.path(logs, "lisa_standard_completion_receipt.tsv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(logs, "sessionInfo_lisa_standard.txt")
)
cat("LISA_STANDARD_GATE=PASS\n")
