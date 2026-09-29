# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Input contract configuration, evidence, duplicate-resolution, and preflight contracts.

lisa_pipeline_schema_version <- function() "1.0.0"

lisa_config_allowed_keys <- function(section) {
  switch(section,
    top_level = c(
      "schema_version", "pipeline", "project", "collections", "single_de",
      "de_index", "de_index_path", "contrasts", "contrast_index",
      "contrast_index_path", "source_data", "allowlisted_source_paths", "report"
    ),
    pipeline = c(
      "schema_version", "profile", "evidence_mode", "evidence",
      "duplicate_policies", "output_dir", "workers", "base_dir", "de_index_path",
      "contrast_index_path", "dictionary_resource", "term2gene_resource",
      "category_map_resource", "lisa_dictionary", "file_label_prefix",
      "gsea_padj_cutoff",
      # These legacy names remain recognizable only so that the migration
      # validator can issue its specific repair message.
      "plot_formats", "export_formats", "run_gene_level", "run_reports",
      "run_kegg_maps", "run_hallmarks", "run_ora", "pathways_category_display",
      "kegg_map_selection", "kegg_access_mode", "kegg_cache_root",
      "kegg_snapshot_id", "report_title", "dry_run",
      "resume", "run_id", "sufficiency_action", "sufficiency_ignore_reason",
      "cache_dir", "cache_mode", "cache_max_bytes",
      "lisa_project_root", "dictionary_dir", "term2gene",
      "dictionary_cache_root", "shared_dictionary_root", "dictionary_registry"
    ),
    project = c("title", "output_dir"),
    report = c("mode", "formats", "source_data", "recipes",
      "category_evidence", "legacy_gene_products", "evidence_max_sets", "evidence_max_genes",
      "category_nes_variants"),
    report_formats = c("png", "svg", "pdf"),
    source_data = c(
      "path", "source_path", "role", "source_role", "target_subdir",
      "allowlisted_paths"
    ),
    evidence = c("required_columns", "products"),
    duplicate_policies = c(
      "de_table_duplicate_policy", "matrix_duplicate_policy",
      "mapped_id_collision_policy"
    ),
    duplicate_policy = c(
      "type", "criterion", "tie_breaker", "column", "numeric_method",
      "non_numeric_method", "mapping_file", "source_row_col", "output_id_col"
    ),
    tie_breaker = "column",
    analysis = c(
      "analysis_id", "de_path", "species", "label", "comparison",
      "positive_direction", "model_note", "symbol_col", "rank_col",
      "logfc_col", "pvalue_col", "padj_col", "gene_id_col", "msigdb_mode",
      "term2gene_target_species", "matrix_path", "expression_matrix_path",
      "counts_matrix_path", "vst_matrix_path", "normalized_matrix_path",
      "tpm_matrix_path", "expression_matrix", "counts_matrix", "vst_matrix",
      "normalized_matrix", "tpm_matrix", "sample_counts_path", "sample_counts",
      "matrix_id_col", "matrix_feature_col", "sample_include_regex",
      "heatmap_sample_include_regex", "sample_exclude_regex",
      "heatmap_sample_exclude_regex", "mapped_id_path", "mapped_id_col",
      "source_id_col"
    ),
    contrast = c(
      "contrast_id", "output_id", "analysis_a", "analysis_b", "contrast_label",
      "comparison", "positive_direction"
    ),
    stop("Unknown configuration section: ", section, call. = FALSE)
  )
}

lisa_config_assert_object <- function(x, field) {
  if (!is.list(x) || is.data.frame(x) ||
      (length(x) && (is.null(names(x)) || any(!nzchar(names(x)))))) {
    lisa_structural_error(
      "LISA-CONFIG-TYPE-002", field, paste(class(x), collapse = ","),
      "configuration section must be a named object", "one JSON/YAML object",
      "Declare named fields inside the object."
    )
  }
  invisible(TRUE)
}

lisa_config_assert_known_keys <- function(x, allowed, field) {
  lisa_config_assert_object(x, field)
  unknown <- setdiff(names(x), allowed)
  if (length(unknown)) {
    bad_field <- paste0(field, ".", unknown[[1L]])
    lisa_structural_error(
      "LISA-CONFIG-UNKNOWN-001", bad_field, unknown[[1L]],
      "unknown configuration key", paste(allowed, collapse = ", "),
      "Correct the field name or remove it; undeclared metadata is not interpreted."
    )
  }
  invisible(TRUE)
}

lisa_config_assert_known_rows <- function(rows, allowed, field) {
  if (is.data.frame(rows)) {
    unknown <- setdiff(names(rows), allowed)
    if (length(unknown)) {
      lisa_structural_error(
        "LISA-CONFIG-UNKNOWN-001", paste0(field, ".", unknown[[1L]]), unknown[[1L]],
        "unknown configuration key", paste(allowed, collapse = ", "),
        "Correct the column name or remove it; undeclared metadata is not interpreted."
      )
    }
    return(invisible(TRUE))
  }
  if (!is.list(rows)) {
    lisa_structural_error(
      "LISA-CONFIG-TYPE-002", field, paste(class(rows), collapse = ","),
      "configuration rows must be an array of named objects",
      "one JSON/YAML array", "Declare each row as a named object."
    )
  }
  for (i in seq_along(rows)) {
    lisa_config_assert_known_keys(rows[[i]], allowed, sprintf("%s[%d]", field, i))
  }
  invisible(TRUE)
}

lisa_config_validate_duplicate_policy_keys <- function(policies) {
  lisa_config_assert_known_keys(
    policies, lisa_config_allowed_keys("duplicate_policies"),
    "pipeline.duplicate_policies"
  )
  for (key in names(policies)) {
    policy <- policies[[key]]
    if (!is.list(policy)) next
    path <- paste0("pipeline.duplicate_policies.", key)
    lisa_config_assert_known_keys(
      policy, lisa_config_allowed_keys("duplicate_policy"), path
    )
    if ("tie_breaker" %in% names(policy) && is.list(policy$tie_breaker)) {
      lisa_config_assert_known_keys(
        policy$tie_breaker, lisa_config_allowed_keys("tie_breaker"),
        paste0(path, ".tie_breaker")
      )
    }
  }
  invisible(TRUE)
}

lisa_validate_raw_config <- function(cfg, strict = TRUE) {
  lisa_config_assert_object(cfg, "config")
  if (isTRUE(strict)) {
    lisa_config_assert_known_keys(
      cfg, lisa_config_allowed_keys("top_level"), "config"
    )
  }

  pipeline <- if ("pipeline" %in% names(cfg)) cfg$pipeline else list()
  lisa_config_assert_object(pipeline, "pipeline")
  # Executable code is part of the installed lisaR package, never study data.
  # Keep this check outside the optional unknown-key gate so strict = FALSE
  # cannot turn a portable configuration into a code-selection capability.
  if ("package_dir" %in% names(pipeline)) {
    code <- "LISA-CONFIG-CODE-001"
    condition <- structure(
      list(
        message = paste(
          code,
          "field=pipeline.package_dir; reason=a study configuration cannot select the executable lisaR tree;",
          "expected=post-processing scripts from the installed lisaR package;",
          "repair=remove pipeline.package_dir and install the intended lisaR version before running."
        ),
        call = NULL,
        code = code,
        field = "pipeline.package_dir"
      ),
      class = c("lisa_code_selection_error", "lisa_config_error", "lisa_error", "error", "condition")
    )
    stop(condition)
  }
  if (isTRUE(strict)) {
    lisa_config_assert_known_keys(
      pipeline, lisa_config_allowed_keys("pipeline"), "pipeline"
    )
  }
  for (key in c("run_kegg_maps", "run_hallmarks", "run_ora", "dry_run", "resume")) {
    if (key %in% names(pipeline)) {
      lisa_config_bool(pipeline[[key]], paste0("pipeline.", key))
    }
  }
  if ("evidence" %in% names(pipeline) && !is.null(pipeline$evidence) && isTRUE(strict)) {
    lisa_config_assert_known_keys(
      pipeline$evidence, lisa_config_allowed_keys("evidence"), "pipeline.evidence"
    )
  }
  if ("duplicate_policies" %in% names(pipeline) &&
      !is.null(pipeline$duplicate_policies) && isTRUE(strict)) {
    lisa_config_validate_duplicate_policy_keys(pipeline$duplicate_policies)
  }

  if ("project" %in% names(cfg) && !is.null(cfg$project) && isTRUE(strict)) {
    lisa_config_assert_known_keys(
      cfg$project, lisa_config_allowed_keys("project"), "project"
    )
  }
  if ("report" %in% names(cfg) && !is.null(cfg$report)) {
    report <- cfg$report
    lisa_config_assert_object(report, "report")
    if (isTRUE(strict)) {
      lisa_config_assert_known_keys(
        report, lisa_config_allowed_keys("report"), "report"
      )
    }
    if ("formats" %in% names(report) && !is.null(report$formats)) {
      lisa_config_assert_object(report$formats, "report.formats")
      if (isTRUE(strict)) {
        lisa_config_assert_known_keys(
          report$formats, lisa_config_allowed_keys("report_formats"),
          "report.formats"
        )
      }
      for (key in intersect(names(report$formats), lisa_config_allowed_keys("report_formats"))) {
        lisa_config_bool(report$formats[[key]], paste0("report.formats.", key))
      }
    }
    for (key in intersect(names(report), c("source_data", "recipes", "category_evidence", "legacy_gene_products"))) {
      lisa_config_bool(report[[key]], paste0("report.", key))
    }
  }

  if (isTRUE(strict)) {
    if ("source_data" %in% names(cfg) && !is.null(cfg$source_data)) {
      lisa_config_assert_known_rows(
        cfg$source_data, lisa_config_allowed_keys("source_data"), "source_data"
      )
    }
    if ("single_de" %in% names(cfg) && !is.null(cfg$single_de)) {
      lisa_config_assert_known_rows(
        cfg$single_de, lisa_config_allowed_keys("analysis"), "single_de"
      )
    }
    if ("de_index" %in% names(cfg) && !is.null(cfg$de_index)) {
      lisa_config_assert_known_rows(
        cfg$de_index, lisa_config_allowed_keys("analysis"), "de_index"
      )
    }
    if ("contrasts" %in% names(cfg) && !is.null(cfg$contrasts)) {
      lisa_config_assert_known_rows(
        cfg$contrasts, lisa_config_allowed_keys("contrast"), "contrasts"
      )
    }
    if ("contrast_index" %in% names(cfg) && !is.null(cfg$contrast_index)) {
      lisa_config_assert_known_rows(
        cfg$contrast_index, lisa_config_allowed_keys("contrast"), "contrast_index"
      )
    }
  }
  invisible(TRUE)
}

lisa_validate_gsea_padj_cutoff <- function(value, field = "pipeline.gsea_padj_cutoff") {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0 || value > 1) {
    lisa_structural_error(
      "LISA-CONFIG-GSEA-001", field, paste(value, collapse = ","),
      "the GSEA adjusted-P cutoff must be one finite numeric value",
      "one number in the inclusive interval [0, 1]",
      "Set gsea_padj_cutoff to the single FDR threshold used by summaries, contrasts and figures."
    )
  }
  as.numeric(value)
}

lisa_validate_config_collections <- function(cfg, run_hallmarks = TRUE) {
  configured <- lisa_config_get(cfg, "collections", lisa_default_collections())
  collections <- as.character(unlist(configured, use.names = FALSE))
  canonical <- lisa_default_collections()

  if (!length(collections)) {
    stop(
      paste0(
        "LISA-COLLECTION-001 collections must contain at least one canonical ",
        "collection. Repair: select one or more of ",
        paste(canonical, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  if (anyNA(collections) || any(!nzchar(collections))) {
    stop(
      paste0(
        "LISA-COLLECTION-002 collections contains an empty value. Repair: ",
        "use only canonical collection IDs: ",
        paste(canonical, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }
  duplicates <- unique(collections[duplicated(collections)])
  if (length(duplicates)) {
    stop(
      paste0(
        "LISA-COLLECTION-003 duplicate collection ID(s): ",
        paste(duplicates, collapse = ", "),
        ". Repair: list each canonical collection at most once."
      ),
      call. = FALSE
    )
  }
  unknown <- setdiff(collections, canonical)
  if (length(unknown)) {
    stop(
      paste0(
        "LISA-COLLECTION-004 unknown collection ID(s): ",
        paste(unknown, collapse = ", "), ". Repair: use only ",
        paste(canonical, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }

  if (!isTRUE(run_hallmarks)) {
    collections <- collections[collections != "HALLMARKS"]
  }
  if (!length(collections)) {
    stop(
      paste0(
        "LISA-COLLECTION-005 run_hallmarks = false removed the only selected ",
        "collection. Repair: select at least one non-HALLMARKS collection or ",
        "set run_hallmarks = true."
      ),
      call. = FALSE
    )
  }
  collections
}

lisa_legacy_config_keys <- function() {
  c(
    plot_formats = "report.formats.{png,svg,pdf}",
    export_formats = "report.source_data and the mandatory canonical TSV tables",
    run_reports = "report.mode",
    run_gene_level = "render_lisa_categories() products"
  )
}

lisa_assert_no_legacy_config <- function(pipeline) {
  legacy <- intersect(names(lisa_legacy_config_keys()), names(pipeline))
  if (!length(legacy)) return(invisible(TRUE))
  replacements <- lisa_legacy_config_keys()[legacy]
  detail <- paste(sprintf("pipeline.%s -> %s", legacy, replacements), collapse = "; ")
  lisa_structural_error(
    "LISA-CONFIG-MIGRATION-001", "pipeline", paste(legacy, collapse = ","),
    "0.5.0 output keys are not accepted implicitly by schema 1.0.0",
    "the versioned 1.0.0 report contract",
    paste0("Migrate explicitly: ", detail, ".")
  )
}

lisa_config_bool <- function(value, field = "configuration boolean") {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    printable <- if (is.list(value)) {
      paste(unlist(value, recursive = TRUE, use.names = FALSE), collapse = ",")
    } else {
      paste(value, collapse = ",")
    }
    lisa_structural_error(
      "LISA-CONFIG-BOOL-001", field, printable,
      "configuration switches must be one JSON/YAML boolean",
      "one scalar true or false", "Set the switch explicitly to true or false."
    )
  }
  isTRUE(value)
}

lisa_report_bool <- function(value, field, default) {
  if (is.null(value)) return(default)
  lisa_config_bool(value, field)
}

lisa_validate_report_config <- function(cfg, allow_selected = FALSE,
                                        strict = TRUE) {
  report <- lisa_config_get(cfg, "report", list())
  if (!is.list(report)) {
    lisa_structural_error("LISA-CONFIG-REPORT-001", "report", class(report)[[1]],
      "report must be an object", "a named report object",
      "Declare report.mode, report.formats, report.source_data and report.recipes.")
  }
  allowed <- lisa_config_allowed_keys("report")
  if (isTRUE(strict)) {
    lisa_config_assert_known_keys(report, allowed, "report")
  }
  mode <- tolower(trimws(as.character(lisa_config_get(report, "mode", "standard")[[1]])))
  allowed_modes <- if (isTRUE(allow_selected)) c("selected", "full") else c("standard", "full")
  if (!mode %in% allowed_modes) {
    lisa_structural_error("LISA-CONFIG-REPORT-004", "report.mode", mode,
      "unsupported report mode", paste(allowed_modes, collapse = " or "),
      if (isTRUE(allow_selected))
        "Use selected for chosen additions or full for every applicable addition."
      else
        "Selected is not a study report mode. Use it with the extension API after a standard run has completed.")
  }
  formats <- lisa_config_get(report, "formats", list())
  if (!is.list(formats)) {
    lisa_structural_error("LISA-CONFIG-REPORT-005", "report.formats", class(formats)[[1]],
      "formats must be an object", "png, svg and pdf boolean switches",
      "Declare report.formats as a named object.")
  }
  if (isTRUE(strict)) {
    lisa_config_assert_known_keys(
      formats, c("png", "svg", "pdf"), "report.formats"
    )
  }
  format_flags <- c(
    png = lisa_report_bool(formats$png, "report.formats.png", TRUE),
    svg = lisa_report_bool(formats$svg, "report.formats.svg", FALSE),
    pdf = lisa_report_bool(formats$pdf, "report.formats.pdf", FALSE)
  )
  list(
    mode = mode,
    formats = format_flags,
    source_data = lisa_report_bool(report$source_data, "report.source_data", TRUE),
    recipes = lisa_report_bool(report$recipes, "report.recipes", FALSE),
    category_evidence = lisa_report_bool(report$category_evidence, "report.category_evidence", TRUE),
    legacy_gene_products = lisa_report_bool(report$legacy_gene_products, "report.legacy_gene_products", FALSE),
    evidence_max_sets = lisa_config_positive_integer(report$evidence_max_sets, "report.evidence_max_sets", 25L),
    evidence_max_genes = lisa_config_positive_integer(report$evidence_max_genes, "report.evidence_max_genes", 40L),
    category_nes_variants = lisa_config_category_nes_variants(
      if ("category_nes_variants" %in% names(report)) report$category_nes_variants else
        c("clean", "percentages", "direction", "dispersion"))
  )
}

# YAML/JSON arrays become unnamed lists; R callers may supply character vectors.
# Reject coercions, duplicates and empty selections instead of silently changing
# which figures a report requested. The order is preserved for the HTML selector.
lisa_config_category_nes_variants <- function(value) {
  if (is.list(value) && !is.data.frame(value) && is.null(names(value)) &&
      length(value) && all(vapply(value, function(x)
        is.character(x) && length(x) == 1L && !is.na(x), logical(1)))) {
    value <- unlist(value, use.names = FALSE)
  }
  if (!is.character(value) || !length(value) || anyNA(value) ||
      any(!value %in% c("clean", "percentages", "direction", "dispersion")) || anyDuplicated(value)) {
    lisa_structural_error("LISA-CONFIG-NES-VARIANTS-001", "report.category_nes_variants",
      paste(value, collapse = ","), "expected a nonempty unique selection of NES plot variants",
      "one to four of clean, percentages, direction, dispersion",
      "Select at least one variant, for example [clean], or all four: [clean, percentages, direction, dispersion].")
  }
  unname(value)
}

# Configuration controls have one numeric type in every supported encoding.
lisa_config_positive_integer <- function(value, field, default) {
  if (is.null(value)) return(default)
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != floor(value) ||
      value > .Machine$integer.max) {
    lisa_structural_error("LISA-CONFIG-INTEGER-001", field, paste(value, collapse = ","),
      "expected one positive integer", paste0("1 through ", .Machine$integer.max),
      "Use an integer-valued number, not a string or a vector.")
  }
  as.integer(value)
}

lisa_validate_cache_config <- function(pipeline) {
  mode <- lisa_config_get(pipeline, "cache_mode", "off")
  if (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
      !mode %in% c("off", "readwrite", "readonly", "refresh")) {
    stop("LISA-CONFIG-CACHE-001 pipeline.cache_mode must be off, readwrite, readonly or refresh.", call. = FALSE)
  }
  path <- lisa_config_get(pipeline, "cache_dir", NULL)
  if (!is.null(path) && (!is.character(path) || length(path) != 1L ||
      is.na(path) || !nzchar(trimws(path)))) {
    stop("LISA-CONFIG-CACHE-002 pipeline.cache_dir must be one non-empty path.", call. = FALSE)
  }
  if (mode != "off" && is.null(path)) {
    stop("LISA-CONFIG-CACHE-002 active caching requires an explicit pipeline.cache_dir.", call. = FALSE)
  }
  bytes <- lisa_config_get(pipeline, "cache_max_bytes", 536870912)
  if (!is.numeric(bytes) || length(bytes) != 1L || is.na(bytes) ||
      !is.finite(bytes) || bytes < 1 || bytes != floor(bytes) || bytes > 2^53 - 1) {
    stop("LISA-CONFIG-CACHE-003 pipeline.cache_max_bytes must be one positive exactly representable integer (at most 2^53 - 1).", call. = FALSE)
  }
  list(mode = mode, dir = path, max_bytes = as.numeric(bytes))
}

lisa_structural_error <- function(code, field, value = "", reason, expected, repair, rows = NULL, file = NULL) {
  location <- paste0("field=", field, if (!is.null(file)) paste0("; file=", file) else "")
  if (!is.null(rows)) location <- paste0(location, "; rows=", paste(rows, collapse = ","))
  stop(sprintf("%s: %s; value=%s; reason=%s; expected=%s; repair=%s", code, location,
    as.character(value), reason, expected, repair), call. = FALSE)
}

lisa_evidence_contract <- function(evidence_mode = "full_de", custom = NULL) {
  mode <- tolower(as.character(evidence_mode[[1]]))
  built_in <- list(full_de = c("effect", "significance"), rank_only = "rank", effect_only = "effect")
  if (identical(mode, "custom")) {
    if (!is.list(custom) || !length(custom$required_columns %||% character())) {
      lisa_structural_error("LISA-EVIDENCE-004", "pipeline.evidence.required_columns", mode,
        "custom evidence has no declared required_columns", "pipeline.evidence.required_columns as a non-empty list",
        "Declare the source columns and compatible products in pipeline.evidence.")
    }
    required <- as.character(unlist(custom$required_columns))
    allowed <- c("effect", "significance", "rank")
    if (any(!required %in% allowed)) {
      lisa_structural_error(
        "LISA-EVIDENCE-005", "pipeline.evidence.required_columns",
        paste(setdiff(required, allowed), collapse = ","),
        "custom evidence requirements must name semantic capabilities",
        paste(allowed, collapse = ", "),
        "Map source columns with logfc_col, padj_col and rank_col, then declare the matching capabilities."
      )
    }
    return(list(mode = mode, required = unique(required),
      products = unique(as.character(unlist(custom$products %||% character())))))
  }
  if (!mode %in% names(built_in)) lisa_structural_error("LISA-EVIDENCE-001", "evidence_mode", mode,
    "unsupported evidence mode", "full_de, rank_only, effect_only, or custom",
    "Set pipeline.evidence_mode to one approved value.")
  products <- switch(mode, full_de = c("gene_level", "volcano", "gene_cards", "kegg_maps"),
    rank_only = "ranked_enrichment", effect_only = "effect_summary")
  list(mode = mode, required = built_in[[mode]], products = products)
}

lisa_validate_pipeline_config <- function(cfg, strict = TRUE) {
  strict <- lisa_config_bool(strict, "strict")
  lisa_validate_raw_config(cfg, strict = strict)
  pipeline <- lisa_config_get(cfg, "pipeline", list())
  # Validate path-bearing identifiers before output destinations, locks or
  # staging names can be derived from them. Human-facing titles and labels are
  # deliberately outside this ASCII technical-ID contract.
  for (key in c("run_id", "file_label_prefix")) {
    value <- pipeline[[key]]
    if (!is.null(value)) {
      lisa_safe_id(value, paste0("pipeline.", key))
    }
  }
  lisa_assert_resume_disabled(
    lisa_config_bool(lisa_config_get(pipeline, "resume", FALSE))
  )
  version <- as.character(lisa_config_get(pipeline, "schema_version", lisa_config_get(cfg, "schema_version", "")))
  if (identical(version, "0.5.0")) {
    lisa_structural_error("LISA-CONFIG-MIGRATION-000", "pipeline.schema_version", version,
      "schema 0.5.0 requires explicit migration", lisa_pipeline_schema_version(),
      "Set schema_version to 1.0.0, replace plot_formats/export_formats/run_reports/run_gene_level with the report object, and review the resulting output plan.")
  }
  if (identical(version, "0.6.0")) {
    lisa_structural_error("LISA-CONFIG-MIGRATION-000", "pipeline.schema_version", version,
      "schema 0.6.0 requires explicit migration", lisa_pipeline_schema_version(),
      paste0("Review report.mode and the evidence options, select one to four report.category_nes_variants ",
        "(clean, percentages, direction, dispersion; all four by default), then set pipeline.schema_version to 1.0.0 ",
        "in a copy of the configuration and review plan_lisa_outputs(). No input files or scientific parameters need to change."))
  }
  if (!identical(version, lisa_pipeline_schema_version())) lisa_structural_error("LISA-CONFIG-001", "pipeline.schema_version", version,
    "configuration schema version is missing or unsupported", lisa_pipeline_schema_version(),
    "Migrate the configuration explicitly and set pipeline.schema_version to 1.0.0.")
  lisa_assert_no_legacy_config(pipeline)
  profile <- tolower(as.character(lisa_config_get(pipeline, "profile", "")))
  allowed_profiles <- c("transcriptomic/genomic", "global proteomic", "targeted", "custom")
  if (!profile %in% allowed_profiles) lisa_structural_error("LISA-CONFIG-002", "pipeline.profile", profile,
    "unsupported or missing scientific profile", paste(allowed_profiles, collapse = ", "),
    "Declare the profile appropriate to the input modality.")
  evidence <- lisa_evidence_contract(lisa_config_get(pipeline, "evidence_mode", ""), lisa_config_get(pipeline, "evidence", NULL))
  policies <- lisa_config_get(pipeline, "duplicate_policies", NULL)
  keys <- c("de_table_duplicate_policy", "matrix_duplicate_policy", "mapped_id_collision_policy")
  if (!is.list(policies) || !all(keys %in% names(policies))) lisa_structural_error("LISA-DUP-001", "pipeline.duplicate_policies", "",
    "each input class requires an explicit duplicate policy", paste(keys, collapse = ", "),
    "Declare all three policies; use error, select, aggregate, or mapping_file.")
  types <- stats::setNames(vapply(keys, function(key) lisa_validate_duplicate_policy(policies[[key]], key), character(1)), keys)
  if (identical(types[["de_table_duplicate_policy"]], "aggregate")) lisa_structural_error("LISA-DUP-010", "de_table_duplicate_policy", "aggregate",
    "aggregate is not an approved DE-table operation because it can create synthetic effect or significance statistics",
    "error, select, or mapping_file", "Select an observed row, provide a row-level mapping file, or resolve the source table upstream.")
  pathways_category_display <- lisa_config_choice(pipeline, "pathways_category_display", "all", "all")
  if (!identical(pathways_category_display$value, "all")) lisa_structural_error("LISA-PATHWAYS-006", "pipeline.pathways_category_display", pathways_category_display$value,
    "pathway testing requires the complete PATHWAYS universe", "all", "Set pathways_category_display to all; filtering may be performed only downstream of the complete status table.")
  kegg_map_selection <- lisa_config_choice(pipeline, "kegg_map_selection", "significant_only", "significant_only")
  if (!identical(kegg_map_selection$value, "significant_only")) lisa_structural_error("LISA-KEGG-015", "pipeline.kegg_map_selection", kegg_map_selection$value,
    "only existing significant KEGG results may be selected", "significant_only", "Set kegg_map_selection to significant_only.")
  hallmarks <- lisa_config_option(pipeline, "run_hallmarks", TRUE)
  kegg_maps <- lisa_config_option(pipeline, "run_kegg_maps", FALSE)
  kegg_access_mode <- lisa_config_choice(
    pipeline, "kegg_access_mode", "external", c("external", "cache_only")
  )
  kegg_cache_root <- lisa_config_get(pipeline, "kegg_cache_root", "")
  kegg_snapshot_id <- lisa_config_get(pipeline, "kegg_snapshot_id", "")
  if (isTRUE(kegg_maps$value)) {
    if (!identical(kegg_access_mode$value, "cache_only")) {
      lisa_structural_error("LISA-KEGG-016", "pipeline.kegg_access_mode",
        kegg_access_mode$value,
        "configured KEGG painting never retrieves pathway diagrams automatically",
        "cache_only",
        paste0("Use cache_only with a recipient-authorized immutable local KEGG snapshot. ",
          "Do not use this route to accept terms, fetch diagrams, or distribute map bytes."))
    }
    if (!is.character(kegg_cache_root) || length(kegg_cache_root) != 1L ||
        is.na(kegg_cache_root) || !nzchar(trimws(kegg_cache_root)) ||
        !is.character(kegg_snapshot_id) || length(kegg_snapshot_id) != 1L ||
        is.na(kegg_snapshot_id) || !nzchar(trimws(kegg_snapshot_id))) {
      lisa_structural_error("LISA-KEGG-017", "pipeline.kegg_cache_root/kegg_snapshot_id",
        "missing", "cache_only KEGG painting requires both an existing cache root and immutable snapshot ID",
        "non-empty kegg_cache_root and kegg_snapshot_id",
        "Point to the recipient-authorized local snapshot; lisaR will verify its per-resource SHA-256 metadata before rendering.")
    }
    lisa_safe_id(kegg_snapshot_id, "pipeline.kegg_snapshot_id")
  }
  collections <- lisa_validate_config_collections(
    cfg, run_hallmarks = hallmarks$value
  )
  cache <- lisa_validate_cache_config(pipeline)
  list(schema_version = version, profile = profile, evidence = evidence, duplicate_policies = policies,
    cache = cache,
    collections = collections,
    workers = lisa_validate_workers(lisa_config_get(pipeline, "workers", 4L)),
    gsea_padj_cutoff = lisa_validate_gsea_padj_cutoff(
      lisa_config_get(pipeline, "gsea_padj_cutoff", 0.25)
    ),
    report = lisa_validate_report_config(cfg, strict = strict),
    ora = lisa_config_option(pipeline, "run_ora", FALSE),
    hallmarks = hallmarks,
    kegg_maps = kegg_maps, kegg_access_mode = kegg_access_mode,
    kegg_cache_root = kegg_cache_root, kegg_snapshot_id = kegg_snapshot_id,
    pathways_category_display = pathways_category_display, kegg_map_selection = kegg_map_selection)
}

lisa_config_option <- function(x, key, default) {
  if (is.null(x[[key]])) return(list(value = default, source = "default"))
  list(value = lisa_config_bool(x[[key]], paste0("pipeline.", key)), source = "explicit")
}

lisa_config_choice <- function(x, key, default, allowed) {
  explicit <- !is.null(x[[key]])
  value <- if (explicit) tolower(trimws(as.character(x[[key]][[1]]))) else default
  if (!value %in% allowed) lisa_structural_error("LISA-CONFIG-014", paste0("pipeline.", key), value,
    "unsupported option", paste(allowed, collapse = ", "), "Set an approved explicit value or omit the option to use its default.")
  list(value = value, source = if (explicit) "explicit" else "default")
}

lisa_policy_type <- function(policy) tolower(as.character(if (is.character(policy)) policy[[1]] else policy$type %||% ""))

lisa_validate_duplicate_policy <- function(policy, field) {
  value <- lisa_policy_type(policy)
  if (!value %in% c("error", "select", "aggregate", "mapping_file")) lisa_structural_error("LISA-DUP-002", field, value,
    "unsupported duplicate policy", "error, select, aggregate, or mapping_file",
    "Declare one approved policy; implicit first is prohibited.")
  if (identical(value, "select")) {
    criterion <- tolower(as.character(policy$criterion %||% ""))
    tie <- policy$tie_breaker %||% ""
    if (!is.list(policy) || !criterion %in% c("min_padj", "max_abs_rank", "max_quality", "min", "max", "min_abs", "max_abs") || !nzchar(as.character(if (is.list(tie)) tie$column %||% "" else tie))) lisa_structural_error("LISA-DUP-003", field, value,
      "select requires a supported criterion and deterministic tie_breaker", "criterion plus tie_breaker=row_number or an existing column",
      "Provide both fields explicitly.")
    if (criterion %in% c("min", "max", "min_abs", "max_abs") && !nzchar(as.character(policy$column %||% ""))) lisa_structural_error("LISA-DUP-003", field, value,
      "generic select criteria require policy.column", "a source column name", "Set policy.column to the numeric ranking column.")
  }
  if (identical(value, "aggregate")) {
    numeric_method <- tolower(as.character(policy$numeric_method %||% ""))
    non_numeric_method <- tolower(as.character(policy$non_numeric_method %||% ""))
    if (!is.list(policy) || !numeric_method %in% c("mean", "sum", "median", "min", "max") || !non_numeric_method %in% c("constant", "collapse")) lisa_structural_error("LISA-DUP-007", field, value,
      "aggregate requires explicit numeric_method and non_numeric_method", "numeric_method=mean|sum|median|min|max and non_numeric_method=constant|collapse",
      "Declare deterministic aggregation for each column class.")
  }
  if (identical(value, "mapping_file") && (!is.list(policy) || !nzchar(as.character(policy$mapping_file %||% "")))) lisa_structural_error("LISA-DUP-008", field, value,
    "mapping_file requires a declared mapping_file", "mapping_file plus source_row/output_id columns when needed",
    "Provide a tab-separated mapping file and its column names.")
  value
}

lisa_preflight_de_table <- function(de_table, profile, evidence_mode, columns = list(), duplicate_policy, file = NULL, warning_action = "warn", ignore_reason = "") {
  df <- if (is.character(de_table) && length(de_table) == 1) read_lisa_tsv(de_table) else de_table
  id_col <- columns$id %||% columns$symbol %||% "symbol"
  effect_col <- columns$effect %||% columns$logfc %||% "log2FoldChange"
  significance_col <- columns$significance %||% columns$padj %||% "padj"
  rank_col <- columns$rank %||% "rank_value"
  contract <- lisa_evidence_contract(evidence_mode, columns$custom)
  need <- c(id = id_col)
  if ("effect" %in% contract$required) need <- c(need, effect = effect_col)
  if ("significance" %in% contract$required) need <- c(need, significance = significance_col)
  if ("rank" %in% contract$required) need <- c(need, rank = rank_col)
  missing <- need[!unname(need) %in% names(df)]
  if (length(missing)) lisa_structural_error("LISA-EVIDENCE-002", names(missing)[1], missing[[1]], "required evidence column is absent",
    paste(unname(need), collapse = ", "), "Add the source column or select an evidence_mode compatible with available evidence.", file = file)
  if ("significance" %in% names(need)) {
    vals <- suppressWarnings(as.numeric(df[[significance_col]]))
    if (any(!is.na(vals) & (vals < 0 | vals > 1))) lisa_structural_error("LISA-PREFLIGHT-003", significance_col, "outside [0,1]",
      "significance values are invalid", "numeric values in [0, 1]", "Correct p-value/FDR values in the source table.", file = file)
  }
  if ("rank" %in% names(need) && !any(is.finite(suppressWarnings(as.numeric(df[[rank_col]]))))) lisa_structural_error("LISA-PREFLIGHT-004", rank_col, "no finite values",
    "ranking cannot be interpreted numerically", "at least one finite numeric rank", "Provide a numeric ranking column.", file = file)
  resolution <- lisa_resolve_duplicates(df, id_col, duplicate_policy,
    input_class = "de_table", file = file, columns = list(padj = significance_col, rank = rank_col))
  warnings <- lisa_sufficiency_messages(profile, length(unique(resolution$data[[id_col]])))
  lisa_handle_preflight_warnings(warnings, warning_action, ignore_reason)
  list(data = resolution$data, duplicate_resolution = resolution$evidence, warnings = warnings, evidence = contract)
}

lisa_preflight_matrix <- function(matrix, profile, feature_col = "feature_id", duplicate_policy, file = NULL, warning_action = "warn", ignore_reason = "") {
  df <- if (is.character(matrix) && length(matrix) == 1) read_lisa_tsv(matrix) else as.data.frame(matrix, stringsAsFactors = FALSE)
  if (!feature_col %in% names(df)) lisa_structural_error("LISA-MATRIX-001", "feature_col", feature_col, "matrix feature column is absent",
    "a declared feature identifier column", "Set matrix_feature_col to a column present in the matrix.", file = file)
  value_cols <- setdiff(names(df), feature_col)
  numeric_columns <- vapply(df[value_cols], function(x) {
    all(is.na(x) | is.finite(suppressWarnings(as.numeric(x))))
  }, logical(1))
  if (!length(value_cols) || any(!numeric_columns)) lisa_structural_error("LISA-MATRIX-002", "matrix", "non-numeric values",
    "matrix measurements must be numeric", "at least one numeric sample column", "Provide a numeric expression or abundance matrix.", file = file)
  resolution <- lisa_resolve_duplicates(df, feature_col, duplicate_policy, input_class = "matrix", file = file)
  warnings <- lisa_sufficiency_messages(profile, length(unique(resolution$data[[feature_col]])))
  lisa_handle_preflight_warnings(warnings, warning_action, ignore_reason)
  list(data = resolution$data, duplicate_resolution = resolution$evidence, warnings = warnings)
}

lisa_resolve_mapped_id_collisions <- function(mapping, source_id_col = "source_id", mapped_id_col = "mapped_id", duplicate_policy, file = NULL) {
  if (!all(c(source_id_col, mapped_id_col) %in% names(mapping))) lisa_structural_error("LISA-MAP-001", "mapping", "missing columns",
    "mapped-ID collision check requires source and mapped identifier columns", paste(c(source_id_col, mapped_id_col), collapse = ", "),
    "Provide the declared mapping columns.", file = file)
  lisa_resolve_duplicates(mapping, mapped_id_col, duplicate_policy, input_class = "mapped_id", file = file)
}

lisa_sufficiency_messages <- function(profile, n_features) {
  limits <- switch(profile, "transcriptomic/genomic" = c(strong = 2000L, warning = 5000L), "global proteomic" = c(strong = 500L, warning = 1500L), c(strong = 0L, warning = 0L))
  if (identical(profile, "targeted")) return("LISA-SUFFICIENCY-003: targeted profile has a limited feature universe; no universal feature-count minimum is applied.")
  if (limits[["strong"]] > 0 && n_features < limits[["strong"]]) return(sprintf("LISA-SUFFICIENCY-002: strong warning: %s input has %d usable unique features, below the %d strong-warning threshold.", profile, n_features, limits[["strong"]]))
  if (limits[["warning"]] > 0 && n_features < limits[["warning"]]) return(sprintf("LISA-SUFFICIENCY-001: warning: %s input has %d usable unique features, below the %d recommended threshold.", profile, n_features, limits[["warning"]]))
  character()
}

lisa_handle_preflight_warnings <- function(warnings, action = "warn", ignore_reason = "") {
  if (!length(warnings)) return(invisible(NULL))
  action <- tolower(as.character(action[[1]]))
  if (!action %in% c("warn", "error", "ignore_with_reason")) lisa_structural_error("LISA-PREFLIGHT-005", "warning_action", action,
    "unsupported sufficiency action", "warn, error, or ignore_with_reason", "Choose an approved action.")
  if (identical(action, "error")) stop(paste(warnings, collapse = "\n"), call. = FALSE)
  if (identical(action, "ignore_with_reason") && !nzchar(trimws(as.character(ignore_reason[[1]])))) lisa_structural_error("LISA-PREFLIGHT-005", "warning_action", action,
    "ignored warnings need a justification", "ignore_with_reason plus non-empty reason", "Record the investigator justification.")
  if (identical(action, "warn")) warning(paste(warnings, collapse = "\n"), call. = FALSE)
  invisible(NULL)
}

lisa_duplicate_audit_empty <- function() data.frame(original_row = integer(), original_id = character(), output_id = character(), action = character(), policy = character(), parameters = character(), input_values = character(), output_values = character(), stringsAsFactors = FALSE)

lisa_canonical_duplicate_id <- function(ids) {
  # Match the package's established ranking normalization: symbols, including
  # mouse-native identifiers such as Stat3, are compared in uppercase only.
  # This is deliberately not species-specific identifier conversion.
  toupper(as.character(ids))
}

lisa_resolve_duplicates <- function(df, id_col, policy, input_class = c("de_table", "matrix", "mapped_id"), file = NULL, columns = list()) {
  input_class <- match.arg(input_class)
  if (!id_col %in% names(df)) lisa_structural_error("LISA-DUP-009", "id_col", id_col, "identifier column is absent", "a column in the input table", "Set the identifier column explicitly.", file = file)
  ids <- as.character(df[[id_col]])
  canonical_ids <- lisa_canonical_duplicate_id(ids)
  duplicate_ids <- unique(canonical_ids[duplicated(canonical_ids)])
  if (!length(duplicate_ids)) return(list(data = df, evidence = lisa_duplicate_audit_empty()))
  type <- lisa_validate_duplicate_policy(policy, input_class)
  rows <- which(canonical_ids %in% duplicate_ids)
  if (identical(type, "error")) lisa_structural_error("LISA-DUP-004", id_col, paste(duplicate_ids, collapse = ","), "duplicate identifiers require resolution",
    "unique identifiers or an explicit resolving policy", "Use an approved policy with deterministic parameters.", rows = rows, file = file)
  if (identical(type, "select")) return(lisa_duplicate_select(df, id_col, policy, input_class, file, columns = columns))
  if (identical(type, "aggregate")) return(lisa_duplicate_aggregate(df, id_col, policy, input_class, file))
  lisa_duplicate_mapping_file(df, id_col, policy, input_class, file)
}

lisa_duplicate_select <- function(df, id_col, policy, input_class, file, columns = list()) {
  criterion <- tolower(as.character(policy$criterion))
  score_col <- switch(criterion,
    min_padj = columns$padj %||% "padj",
    max_abs_rank = columns$rank %||% "rank_value",
    max_quality = "quality",
    policy$column
  )
  if (!score_col %in% names(df)) lisa_structural_error("LISA-DUP-006", id_col, score_col, "select criterion column is absent",
    score_col, "Add the criterion column or choose an applicable explicit policy.", file = file)
  tie <- policy$tie_breaker
  tie_col <- if (identical(tie, "row_number")) NULL else if (is.list(tie)) as.character(tie$column %||% "") else as.character(tie)
  if (!is.null(tie_col) && !tie_col %in% names(df)) lisa_structural_error("LISA-DUP-011", id_col, tie_col, "tie_breaker column is absent",
    "row_number or an existing source column", "Correct tie_breaker in the duplicate policy.", file = file)
  canonical_ids <- lisa_canonical_duplicate_id(df[[id_col]])
  groups <- split(seq_len(nrow(df)), canonical_ids)
  kept <- vapply(groups, function(group_rows) {
    values <- suppressWarnings(as.numeric(df[[score_col]][group_rows]))
    score <- switch(criterion, min_padj = values, min = values, min_abs = abs(values), -abs(values))
    ties <- if (is.null(tie_col)) group_rows else df[[tie_col]][group_rows]
    group_rows[order(score, ties, group_rows, na.last = TRUE)][1]
  }, integer(1))
  out <- df[sort(kept), , drop = FALSE]
  list(data = out, evidence = lisa_duplicate_audit(df, id_col, out, type = "select", parameters = paste0("input_class=", input_class, ";criterion=", criterion, ";column=", score_col, ";tie_breaker=", if (is.null(tie_col)) "row_number" else tie_col), selected_rows = kept), selected_rows = kept)
}

lisa_duplicate_aggregate <- function(df, id_col, policy, input_class, file) {
  if (identical(input_class, "de_table")) lisa_structural_error("LISA-DUP-010", id_col, "aggregate", "DE-table aggregation can synthesize statistics",
    "error, select, or mapping_file", "Resolve DE duplicates upstream or select an observed row.", file = file)
  numeric_method <- tolower(as.character(policy$numeric_method)); non_numeric_method <- tolower(as.character(policy$non_numeric_method))
  groups <- split(seq_len(nrow(df)), lisa_canonical_duplicate_id(df[[id_col]]))
  rows <- lapply(groups, function(group_rows) {
    out <- lapply(names(df), function(col) {
      values <- df[[col]][group_rows]
      if (identical(col, id_col)) return(as.character(values[[1]]))
      numeric_values <- suppressWarnings(as.numeric(values))
      if (all(is.na(values) | !is.na(numeric_values))) return(switch(numeric_method, mean = mean(numeric_values, na.rm = TRUE), sum = sum(numeric_values, na.rm = TRUE), median = stats::median(numeric_values, na.rm = TRUE), min = min(numeric_values, na.rm = TRUE), max = max(numeric_values, na.rm = TRUE)))
      unique_values <- unique(as.character(values))
      if (identical(non_numeric_method, "constant") && length(unique_values) > 1) lisa_structural_error("LISA-DUP-012", col, paste(unique_values, collapse = ","), "non-numeric duplicate values differ under constant aggregation",
        "identical values or non_numeric_method=collapse", "Choose an explicit non-numeric aggregation method.", rows = group_rows, file = file)
      if (identical(non_numeric_method, "collapse")) return(paste(unique_values, collapse = ";"))
      unique_values[[1]]
    })
    stats::setNames(as.list(out), names(df))
  })
  out <- as.data.frame(do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE)), stringsAsFactors = FALSE, check.names = FALSE)
  for (col in names(df)) if (is.numeric(df[[col]])) out[[col]] <- as.numeric(out[[col]])
  list(data = out, evidence = lisa_duplicate_audit(df, id_col, out, type = "aggregate", parameters = paste0("input_class=", input_class, ";numeric_method=", numeric_method, ";non_numeric_method=", non_numeric_method)))
}

lisa_duplicate_mapping_file <- function(df, id_col, policy, input_class, file) {
  mapping_path <- as.character(policy$mapping_file)
  if (!file.exists(mapping_path)) lisa_structural_error("LISA-DUP-013", "mapping_file", mapping_path, "mapping file does not exist", "an existing TSV mapping file", "Correct duplicate_policies.*.mapping_file.", file = file)
  map <- read_lisa_tsv(mapping_path)
  row_col <- as.character(policy$source_row_col %||% "source_row")
  output_col <- as.character(policy$output_id_col %||% "output_id")
  if (!all(c(row_col, output_col) %in% names(map))) lisa_structural_error("LISA-DUP-014", "mapping_file", mapping_path, "mapping file lacks required columns",
    paste(c(row_col, output_col), collapse = ", "), "Add one row per source row with its resolved output identifier.", file = file)
  source_rows <- suppressWarnings(as.integer(map[[row_col]]))
  if (any(!is.finite(source_rows)) || anyDuplicated(source_rows)) lisa_structural_error("LISA-DUP-015", row_col, "invalid", "mapping source rows must be unique integers",
    "unique positive source-row numbers", "Correct the mapping file source-row column.", file = file)
  canonical_ids <- lisa_canonical_duplicate_id(df[[id_col]])
  target_rows <- which(duplicated(canonical_ids) | duplicated(canonical_ids, fromLast = TRUE))
  if (!all(target_rows %in% source_rows)) lisa_structural_error("LISA-DUP-016", row_col, paste(setdiff(target_rows, source_rows), collapse = ","), "mapping file does not resolve every duplicate input row",
    "one mapping row for every duplicate input row", "Add the missing duplicate rows to the mapping file.", file = file)
  out <- df
  out[[id_col]][target_rows] <- as.character(map[[output_col]][match(target_rows, source_rows)])
  if (anyDuplicated(lisa_canonical_duplicate_id(out[[id_col]]))) lisa_structural_error("LISA-DUP-017", output_col, "duplicate output identifiers", "mapping_file did not produce unique identifiers",
    "unique resolved output identifiers", "Revise the mapping file or use aggregate/select where appropriate.", file = file)
  list(data = out, evidence = lisa_duplicate_audit(df, id_col, out, type = "mapping_file", parameters = paste0("input_class=", input_class, ";mapping_file=", normalizePath(mapping_path, mustWork = TRUE), ";source_row_col=", row_col, ";output_id_col=", output_col)))
}

lisa_duplicate_audit <- function(input, id_col, output, type, parameters, selected_rows = NULL) {
  original_rows <- seq_len(nrow(input))
  # A dropped duplicate remains auditable: it points to the retained/aggregated
  # output identifier rather than disappearing from the resolution evidence.
  output_ids <- as.character(input[[id_col]])
  action <- rep("retained", nrow(input))
  if (!is.null(selected_rows)) {
    selected_rows <- as.integer(selected_rows)
    canonical_ids <- lisa_canonical_duplicate_id(input[[id_col]])
    selected_ids <- as.character(input[[id_col]][selected_rows])
    output_ids <- selected_ids[match(canonical_ids, lisa_canonical_duplicate_id(selected_ids))]
    action <- ifelse(original_rows %in% selected_rows, "retained", "dropped")
  } else if (nrow(output) == nrow(input)) {
    output_ids <- as.character(output[[id_col]])
  }
  data.frame(original_row = original_rows, original_id = as.character(input[[id_col]]), output_id = output_ids, action = action, policy = type,
    parameters = parameters, input_values = vapply(original_rows, function(i) paste(input[i, ], collapse = ";"), character(1)),
    output_values = vapply(original_rows, function(i) { hit <- which(as.character(output[[id_col]]) == output_ids[[i]]); if (!length(hit)) "" else paste(output[hit[[1]], ], collapse = ";") }, character(1)), stringsAsFactors = FALSE)
}

lisa_product_applicability <- function(evidence, run_gene_level = TRUE, run_kegg_maps = FALSE) {
  rows <- data.frame(product = c("gene_level", "volcano", "gene_cards", "kegg_maps", "lisa_gps_provisional", "ranked_enrichment", "effect_summary"),
    enabled = c(isTRUE(run_gene_level), isTRUE(run_gene_level), isTRUE(run_gene_level), isTRUE(run_kegg_maps), TRUE, TRUE, TRUE), stringsAsFactors = FALSE)
  rows$required_evidence <- c("full_de", "full_de", "full_de", "full_de", "full_de", "rank_only", "effect_only")
  rows$status <- ifelse(!rows$enabled, "disabled", ifelse(rows$required_evidence == evidence$mode, "enabled", "not_applicable"))
  rows$reason <- ifelse(rows$status == "enabled", paste0("evidence_mode=", evidence$mode, " satisfies product evidence requirement"),
    ifelse(rows$status == "disabled", "disabled by configuration", paste0("requires evidence_mode=", rows$required_evidence, "; received ", evidence$mode)))
  rows
}
