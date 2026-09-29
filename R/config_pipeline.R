# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

read_lisa_pipeline_config <- function(config_path) {
  if (!file.exists(config_path)) {
    stop("Pipeline config does not exist: ", config_path, call. = FALSE)
  }
  ext <- tolower(tools::file_ext(config_path))
  if (ext == "json") {
    lisa_require_optional("jsonlite", "reading JSON pipeline configs")
    return(jsonlite::fromJSON(config_path, simplifyVector = FALSE))
  }
  if (ext %in% c("yaml", "yml")) {
    if (requireNamespace("yaml", quietly = TRUE)) {
      return(yaml::read_yaml(config_path))
    }
    return(read_lisa_yaml_via_python(config_path))
  }
  stop("Unsupported config extension: .", ext, ". Use .json, .yaml or .yml.", call. = FALSE)
}

lisa_tag_duplicate_audit <- function(audit, input_class) {
  audit <- as.data.frame(audit, stringsAsFactors = FALSE)
  cbind(
    input_class = rep(as.character(input_class), nrow(audit)),
    audit,
    stringsAsFactors = FALSE
  )
}

lisa_materialize_effective_input <- function(data, analysis_id, input_class,
                                             source_path, output_dir, policy,
                                             duplicate_resolution_path,
                                             source_rows) {
  effective_path <- file.path(
    output_dir, "effective_inputs",
    paste0(analysis_id, "_", input_class, ".tsv")
  )
  write_lisa_tsv(data, effective_path)
  list(
    path = effective_path,
    receipt = data.frame(
      analysis_id = analysis_id,
      input_class = input_class,
      source_path = as.character(source_path),
      effective_path = effective_path,
      sha256 = lisa_sha256_file(effective_path),
      duplicate_policy = lisa_policy_type(policy),
      source_rows = as.integer(source_rows),
      effective_rows = nrow(data),
      duplicate_resolution = duplicate_resolution_path,
      stringsAsFactors = FALSE
    )
  )
}

# Resolve expression scale before materialization changes filenames. The same
# helper serves the installed heatmap script and effective-input preparation.
# *_input_scale fields are generated provenance, not new configuration knobs.
lisa_expression_matrix_source <- function(path, index_row = data.frame(), index_dir = ".") {
  columns <- c("matrix_path", "expression_matrix", "expression_matrix_path",
    "counts_matrix", "counts_matrix_path", "vst_matrix", "vst_matrix_path",
    "normalized_matrix", "normalized_matrix_path", "tpm_matrix", "tpm_matrix_path",
    "sample_counts_path", "sample_counts")
  typed_scales <- c(vst_matrix = "transformed", vst_matrix_path = "transformed",
    counts_matrix = "count_like", counts_matrix_path = "count_like",
    tpm_matrix = "count_like", tpm_matrix_path = "count_like",
    sample_counts_path = "count_like", sample_counts = "count_like")
  scales <- origins <- character()
  for (column in if (nrow(index_row)) intersect(columns, names(index_row)) else character()) {
    declared <- as.character(index_row[[column]][[1L]])
    if (is.na(declared) || !nzchar(trimws(declared))) next
    declared <- lisa_norm_path(declared, index_dir)
    if (!file.exists(declared) || !file.exists(path) ||
        !identical(normalizePath(declared), normalizePath(path))) next
    field <- paste0(column, "_input_scale")
    if (field %in% names(index_row) && !is.na(index_row[[field]][[1L]]) &&
        nzchar(as.character(index_row[[field]][[1L]]))) {
      scale <- as.character(index_row[[field]][[1L]])
      origin_field <- paste0(column, "_input_scale_source")
      origin <- if (origin_field %in% names(index_row)) as.character(index_row[[origin_field]][[1L]]) else NA_character_
      if (!scale %in% c("transformed", "count_like") || is.na(origin) || !nzchar(origin)) {
        stop("Invalid materialized expression-scale provenance for ", column, ".", call. = FALSE)
      }
      scales <- c(scales, scale)
      origins <- c(origins, origin)
    } else if (column %in% names(typed_scales)) {
      scales <- c(scales, unname(typed_scales[[column]]))
      origins <- c(origins, paste0("de_index:", column))
    }
  }
  if (length(unique(scales)) > 1L) {
    stop("Conflicting expression scales declared for the same matrix in de_index: ",
         paste(unique(origins), collapse = ", "),
         ". Retain only declarations consistent with the stored matrix scale.", call. = FALSE)
  }
  input_scale <- if (length(scales)) scales[[1L]] else if (
    grepl("vst|rlog|rld|log2", basename(path), ignore.case = TRUE)
  ) "transformed" else "count_like"
  list(path = path, input_scale = input_scale,
       input_scale_source = if (length(origins)) paste(unique(origins), collapse = ";") else "legacy_basename",
       expression_transform = if (input_scale == "transformed") "none" else "log2(x+1)")
}

read_lisa_yaml_via_python <- function(config_path) {
  py <- Sys.which("python3")
  if (!nzchar(py)) {
    stop("YAML config requires R package 'yaml' or python3 with PyYAML.", call. = FALSE)
  }
  code <- paste(
    "import json, sys, yaml",
    "with open(sys.argv[1], 'r', encoding='utf-8') as handle:",
    "    data = yaml.safe_load(handle) or {}",
    "print(json.dumps(data))",
    sep = "\n"
  )
  # The parser receives code with -c and the config as a separate argument: no shell
  # interpolation or temporary script is required.
  parser_root <- lisa_run_root(file.path(tempdir(), "lisaR-config-parser"))
  json <- lisa_run_subprocess(py, c("-c", code, config_path), required = FALSE,
                              stage = "yaml_config_parser", run_root = parser_root)
  if (json$exit_code != 0L) stop("Could not parse YAML config via python3/PyYAML: ", json$diagnostics, call. = FALSE)
  lisa_require_optional("jsonlite", "reading YAML pipeline configs")
  jsonlite::fromJSON(paste(json$stdout, collapse = "\n"), simplifyVector = FALSE)
}

run_lisa_pipeline_from_config <- function(config_path) {
  resolved_config <- lisa_resolve_config(config_path, .execution = TRUE)
  config_path <- resolved_config$config_path
  cfg <- resolved_config$cfg
  validated_config <- resolved_config$contract
  config_dir <- resolved_config$config_dir
  pipeline <- resolved_config$pipeline
  project <- lisa_config_get(cfg, "project", list())
  report_contract <- validated_config$report
  if (identical(report_contract$mode, "selected")) {
    stop("LISA-REPORT-MODE-002 selected is an extension of a completed standard run. Use plan_lisa_extension() and render_lisa_categories() with that run.", call. = FALSE)
  }

  output_dir <- lisa_config_get(pipeline, "output_dir", lisa_config_get(project, "output_dir", NULL))
  if (is.null(output_dir) || !nzchar(as.character(output_dir))) {
    stop("Config must define pipeline.output_dir or project.output_dir.", call. = FALSE)
  }
  output_dir <- lisa_config_managed_path(output_dir, config_dir)
  output_dir <- lisa_managed_destination(
    output_dir, create_parent = TRUE
  )$path
  dry_run <- lisa_config_bool(
    lisa_config_get(pipeline, "dry_run", TRUE), "pipeline.dry_run"
  )
  old_run_root <- options("lisaR.run_root")
  on.exit({
    if ("lisaR.run_root" %in% names(old_run_root)) {
      options(lisaR.run_root = old_run_root$lisaR.run_root)
    } else {
      options(lisaR.run_root = NULL)
    }
  }, add = TRUE)

  base_dir <- lisa_config_path(lisa_config_get(pipeline, "base_dir", output_dir), config_dir)
  report_title <- as.character(lisa_config_get(pipeline, "report_title", lisa_config_get(project, "title", "LISA report")))
  collections <- validated_config$collections
  # The resolver already applied policy and destination guards before reading
  # indexes. Recheck the destination at this boundary in case it appeared since.
  lisa_assert_new_final_output(output_dir)

  de_index <- resolved_config$de_index
  contrast_index <- resolved_config$contrast_index
  lisa_check_run_paths(output_dir, de_index, contrast_index, collections,
    as.character(lisa_config_get(pipeline, "file_label_prefix", "semantic")),
    report_contract)
  species <- unique(as.character(de_index$species))
  if (length(species) != 1L) stop("LISA-SPECIES-006 canonical execution requires one declared species per run. Repair: split mixed-species analyses into separate runs.", call. = FALSE)
  dependency_status <- lisa_assert_config_dependencies(
    pipeline, species[[1L]], report = report_contract
  )
  resource_resolution <- lisa_resolve_pipeline_resources(
    lisa_config_get(resolved_config$raw_config, "pipeline", list()), species[[1L]], stop_on_error = TRUE
  )
  dictionary_resource <- resource_resolution$resolved$dictionary_resource
  term2gene_resource <- resource_resolution$resolved$term2gene_resource
  category_map_resource <- resource_resolution$resolved$category_map_resource

  # All material writes occur in an exclusive staging directory. Neither the
  # final path nor an abandoned staging tree is reused. A retry starts with a
  # fresh run identity and leaves the failed staging tree available for audit.
  final_output_dir <- output_dir
  contract <- lisa_configured_run_contract(
    config_path = config_path,
    cfg = cfg,
    config_dir = config_dir,
    de_index = de_index,
    contrast_index = contrast_index,
    duplicate_policies = validated_config$duplicate_policies,
    resources = resource_resolution$resolved,
    dependency_status = dependency_status
  )
  configured_run_id <- lisa_config_get(pipeline, "run_id", NULL)
  plan_run_id <- if (isTRUE(dry_run)) {
    configured_run_id %||% lisa_new_run_id(contract[["contract"]])
  } else {
    NULL
  }
  transaction_output_dir <- if (isTRUE(dry_run)) {
    file.path(dirname(final_output_dir), paste0("plan-",
      substr(digest::digest(paste(final_output_dir, plan_run_id), algo = "sha256", serialize = FALSE), 1L, 12L)))
  } else {
    final_output_dir
  }
  tx <- lisa_transaction_begin(
    transaction_output_dir, contract,
    run_id = if (isTRUE(dry_run)) plan_run_id else configured_run_id,
    resume = lisa_config_bool(lisa_config_get(pipeline, "resume", FALSE), "pipeline.resume"),
    artifact_kind = if (isTRUE(dry_run)) "plan" else "scientific_run"
  )
  promoted <- FALSE
  on.exit(if (!promoted) lisa_transaction_abort(tx, "run interrupted or failed"), add = TRUE)
  output_dir <- tx$staging_dir
  options(lisaR.run_root = output_dir)
  lisa_transaction_event(tx, "planned")
  if (isTRUE(dry_run)) {
    write_lisa_tsv(data.frame(
      artifact_type = "planning_only",
      gate = "PLAN_PASS",
      scientific_complete = FALSE,
      intended_output_dir = final_output_dir,
      stringsAsFactors = FALSE
    ), file.path(output_dir, "lisa_plan_status.tsv"))
  }

  source_manifest <- lisa_copy_config_source_data(cfg, output_dir, config_dir)
  de_index <- lisa_remap_registered_source_paths(de_index, source_manifest)
  contrast_index <- lisa_remap_registered_source_paths(contrast_index, source_manifest)
  lisa_validate_run_identity(de_index, contrast_index)

  preflight <- lapply(seq_len(nrow(de_index)), function(i) {
    row <- de_index[i, , drop = FALSE]
    columns <- list(
      symbol = if ("symbol_col" %in% names(row)) row$symbol_col[[1]] else "symbol",
      logfc = if ("logfc_col" %in% names(row)) row$logfc_col[[1]] else "log2FoldChange",
      padj = if ("padj_col" %in% names(row)) row$padj_col[[1]] else "padj",
      rank = if ("rank_col" %in% names(row)) row$rank_col[[1]] else "rank_value",
      custom = if (identical(validated_config$evidence$mode, "custom")) {
        list(
          required_columns = validated_config$evidence$required,
          products = validated_config$evidence$products
        )
      } else {
        NULL
      }
    )
    result <- lisa_preflight_de_table(row$de_path[[1]], validated_config$profile, validated_config$evidence$mode,
      columns = columns, duplicate_policy = validated_config$duplicate_policies$de_table_duplicate_policy,
      file = row$de_path[[1]], warning_action = as.character(lisa_config_get(pipeline, "sufficiency_action", "warn")),
      ignore_reason = as.character(lisa_config_get(pipeline, "sufficiency_ignore_reason", "")))
    duplicate_path <- file.path(output_dir, "duplicate_resolution", paste0(row$analysis_id[[1]], "_duplicate_resolution.tsv"))
    lisa_guarded_dir_create(dirname(duplicate_path))
    duplicate_audit <- lisa_tag_duplicate_audit(result$duplicate_resolution, "de_table")
    receipts <- list()
    de_effective <- lisa_materialize_effective_input(
      result$data, row$analysis_id[[1]], "de_table", row$de_path[[1]], output_dir,
      validated_config$duplicate_policies$de_table_duplicate_policy, duplicate_path,
      source_rows = nrow(read_lisa_tsv(row$de_path[[1]]))
    )
    row$de_path[[1]] <- de_effective$path
    receipts[[length(receipts) + 1L]] <- de_effective$receipt
    matrix_path_cols <- intersect(c("matrix_path", "expression_matrix_path", "counts_matrix_path", "vst_matrix_path", "normalized_matrix_path"), names(row))
    matrix_source_row <- row
    for (matrix_path_col in matrix_path_cols) {
      # Keep a stable row schema even when an analysis omits this matrix.
      row[[paste0(matrix_path_col, "_input_scale")]] <- ""
      row[[paste0(matrix_path_col, "_input_scale_source")]] <- ""
      matrix_path <- as.character(row[[matrix_path_col]][[1]])
      if (is.na(matrix_path) || !nzchar(matrix_path)) next
      matrix_source <- lisa_expression_matrix_source(matrix_path, matrix_source_row, config_dir)
      matrix_feature_col <- if (identical(matrix_path_col, "matrix_path")) {
        if ("matrix_id_col" %in% names(row) && nzchar(as.character(row$matrix_id_col[[1]]))) row$matrix_id_col[[1]] else "symbol"
      } else if ("matrix_feature_col" %in% names(row) && nzchar(as.character(row$matrix_feature_col[[1]]))) {
        row$matrix_feature_col[[1]]
      } else {
        "feature_id"
      }
      matrix_result <- lisa_preflight_matrix(matrix_path, validated_config$profile,
        feature_col = matrix_feature_col,
        duplicate_policy = validated_config$duplicate_policies$matrix_duplicate_policy, file = matrix_path,
        warning_action = as.character(lisa_config_get(pipeline, "sufficiency_action", "warn")),
        ignore_reason = as.character(lisa_config_get(pipeline, "sufficiency_ignore_reason", "")))
      matrix_audit <- lisa_tag_duplicate_audit(matrix_result$duplicate_resolution, matrix_path_col)
      duplicate_audit <- rbind(duplicate_audit, matrix_audit)
      matrix_effective <- lisa_materialize_effective_input(
        matrix_result$data, row$analysis_id[[1]], matrix_path_col, matrix_path, output_dir,
        validated_config$duplicate_policies$matrix_duplicate_policy, duplicate_path,
        source_rows = nrow(read_lisa_tsv(matrix_path))
      )
      row[[matrix_path_col]][[1]] <- matrix_effective$path
      row[[paste0(matrix_path_col, "_input_scale")]] <- matrix_source$input_scale
      row[[paste0(matrix_path_col, "_input_scale_source")]] <- matrix_source$input_scale_source
      receipts[[length(receipts) + 1L]] <- matrix_effective$receipt
    }
    write_lisa_tsv(duplicate_audit, duplicate_path)
    if ("mapped_id_path" %in% names(row) && nzchar(as.character(row$mapped_id_path[[1]]))) {
      mapped_id_path <- as.character(row$mapped_id_path[[1]])
      mapped_id_col <- if ("mapped_id_col" %in% names(row) && nzchar(as.character(row$mapped_id_col[[1]]))) row$mapped_id_col[[1]] else "mapped_id"
      source_id_col <- if ("source_id_col" %in% names(row) && nzchar(as.character(row$source_id_col[[1]]))) row$source_id_col[[1]] else "source_id"
      mapped_result <- lisa_resolve_mapped_id_collisions(read_lisa_tsv(mapped_id_path), source_id_col = source_id_col,
        mapped_id_col = mapped_id_col, duplicate_policy = validated_config$duplicate_policies$mapped_id_collision_policy, file = mapped_id_path)
      mapped_audit <- lisa_tag_duplicate_audit(mapped_result$evidence, "mapped_id")
      duplicate_audit <- rbind(duplicate_audit, mapped_audit)
      write_lisa_tsv(duplicate_audit, duplicate_path)
      mapped_effective <- lisa_materialize_effective_input(
        mapped_result$data, row$analysis_id[[1]], "mapped_id", mapped_id_path, output_dir,
        validated_config$duplicate_policies$mapped_id_collision_policy, duplicate_path,
        source_rows = nrow(read_lisa_tsv(mapped_id_path))
      )
      row$mapped_id_path[[1]] <- mapped_effective$path
      receipts[[length(receipts) + 1L]] <- mapped_effective$receipt
    }
    list(
      row = row,
      preflight = data.frame(analysis_id = row$analysis_id[[1]], evidence_mode = result$evidence$mode,
        warnings = paste(result$warnings, collapse = " | "), duplicate_resolution = duplicate_path, stringsAsFactors = FALSE),
      receipts = do.call(rbind, receipts)
    )
  })
  de_index <- do.call(rbind, lapply(preflight, `[[`, "row"))
  write_lisa_tsv(do.call(rbind, lapply(preflight, `[[`, "preflight")), file.path(output_dir, "scientific_preflight.tsv"))
  write_lisa_tsv(do.call(rbind, lapply(preflight, `[[`, "receipts")), file.path(output_dir, "effective_inputs.tsv"))
  # The scientific computation is always the bounded standard layer. With
  # report mode "full" the FULL presentation products, and with run_kegg_maps
  # the native KEGG maps, are rendered from the saved DE/GSEA evidence inside
  # this same run before its one report and its one promotion.
  product_plan <- lisa_product_plan(validated_config$evidence,
    run_gene_level = FALSE,
    run_kegg_maps = isTRUE(validated_config$kegg_maps$value))
  write_lisa_tsv(product_plan, file.path(output_dir, "product_plan.tsv"))
  write_lisa_tsv(product_plan, file.path(output_dir, "product_applicability.tsv"))

  config_out <- file.path(output_dir, "config")
  lisa_guarded_dir_create(config_out)
  write_lisa_tsv(de_index, file.path(config_out, "de_index.tsv"))
  if (nrow(contrast_index) > 0) {
    write_lisa_tsv(contrast_index, file.path(config_out, "contrast_index.tsv"))
  } else {
    write_lisa_tsv(data.frame(), file.path(config_out, "contrast_index.tsv"))
  }
  lisa_guarded_copy(config_path, file.path(config_out, basename(config_path)), overwrite = TRUE)
  normalized_cfg <- cfg
  normalized_cfg$collections <- as.list(validated_config$collections)
  normalized_cfg$pipeline$workers <- validated_config$workers
  normalized_cfg$pipeline$gsea_padj_cutoff <- validated_config$gsea_padj_cutoff
  for (resource_key in resource_resolution$status$config_key) {
    row <- resource_resolution$status$config_key == resource_key
    normalized_cfg$pipeline[[resource_key]] <-
      resource_resolution$status$resource_id[row][[1L]]
  }
  if (identical(resource_resolution$dictionary_tier, "custom")) {
    normalized_cfg$pipeline$lisa_dictionary <- NULL
  } else {
    normalized_cfg$pipeline$lisa_dictionary <-
      resource_resolution$dictionary_tier
  }
  normalized_cfg$report <- report_contract
  normalized_cfg$report$formats <- as.list(report_contract$formats)
  normalized_cfg$report$category_nes_variants <- as.list(report_contract$category_nes_variants)
  write_lisa_tsv(resolved_config$provenance, file.path(config_out, "configuration_provenance.tsv"))
  lisa_write_normalized_config(normalized_cfg, file.path(config_out, "pipeline_config.normalized.json"))
  write_lisa_tsv(
    resource_resolution$status,
    file.path(config_out, "resource_resolution.tsv")
  )
  # An ordinary uniform FULL run has no owner-subset concept: every
  # analysis/contrast in this run's own de_index/contrast_index is implicitly
  # "selected" for full-scope extras. Declaring that full set here lets
  # report_full_scope_owner() distinguish a genuinely unselected owner
  # (empty_state "not_requested") from a selected owner whose expected
  # products never materialized (must fail closed as "missing"), instead of
  # inferring selection purely from directory existence. Contrast identities
  # must match the on-disk owner id the report generator actually uses
  # (lisa_contrast_name(): "<contrast_id>_<output_id-or-contrast_id>"), not
  # the bare contrast_id.
  full_scope_analyses <- if (identical(report_contract$mode, "full") && nrow(de_index) > 0) {
    paste(as.character(de_index$analysis_id), collapse = ";")
  } else ""
  full_scope_contrasts <- if (identical(report_contract$mode, "full") && nrow(contrast_index) > 0) {
    paste(vapply(seq_len(nrow(contrast_index)), function(i) {
      lisa_contrast_name(contrast_index[i, , drop = FALSE])
    }, character(1)), collapse = ";")
  } else ""
  write_lisa_tsv(data.frame(
    key = c("configuration_schema_version", "profile", "evidence_mode", "gsea_padj_cutoff", "report_mode", "report_png", "report_svg", "report_pdf", "report_source_data", "report_recipes", "report_category_evidence", "report_legacy_gene_products", "report_evidence_max_sets", "report_evidence_max_genes", "report_category_nes_variants", "run_ora", "run_ora_source", "run_hallmarks", "run_hallmarks_source", "run_kegg_maps", "run_kegg_maps_source", "kegg_access_mode", "kegg_access_mode_source", "kegg_snapshot_id", "lisa_gps_development_contract", "pathways_category_display", "pathways_category_display_source", "kegg_map_selection", "kegg_map_selection_source", "full_scope_analyses", "full_scope_contrasts"),
    value = c(validated_config$schema_version, validated_config$profile, validated_config$evidence$mode,
      lisa_config_text(validated_config$gsea_padj_cutoff),
      report_contract$mode, as.character(report_contract$formats),
      as.character(report_contract$source_data), as.character(report_contract$recipes),
      as.character(report_contract$category_evidence), as.character(report_contract$legacy_gene_products),
      as.character(report_contract$evidence_max_sets), as.character(report_contract$evidence_max_genes),
      paste(report_contract$category_nes_variants, collapse = ","),
      as.character(validated_config$ora$value), validated_config$ora$source,
      as.character(validated_config$hallmarks$value), validated_config$hallmarks$source,
      as.character(validated_config$kegg_maps$value), validated_config$kegg_maps$source,
      validated_config$kegg_access_mode$value, validated_config$kegg_access_mode$source,
      if (isTRUE(validated_config$kegg_maps$value)) validated_config$kegg_snapshot_id else "",
      lisa_gps_development_contract()$development_contract_version,
      validated_config$pathways_category_display$value, validated_config$pathways_category_display$source,
      validated_config$kegg_map_selection$value, validated_config$kegg_map_selection$source,
      full_scope_analyses, full_scope_contrasts),
    stringsAsFactors = FALSE), file.path(output_dir, "contract_manifest.tsv"))

  result <- run_lisa_pipeline(
    de_index = file.path(config_out, "de_index.tsv"),
    contrast_index = file.path(config_out, "contrast_index.tsv"),
    dictionary_dir = dirname(dictionary_resource$path),
    term2gene = term2gene_resource$path,
    dictionary_path = dictionary_resource$path,
    category_map_path = category_map_resource$path,
    output_dir = output_dir,
    collections = collections,
    base_dir = output_dir,
    lisa_project_root = lisa_config_path(lisa_config_get(pipeline, "lisa_project_root", NULL), config_dir),
    lisa_dictionary = resource_resolution$dictionary_tier,
    file_label_prefix = as.character(lisa_config_get(pipeline, "file_label_prefix", "semantic")),
    plot_formats = names(report_contract$formats)[report_contract$formats],
    export_formats = "tsv",
    run_gene_level = FALSE,
    run_reports = TRUE,
    run_kegg_maps = FALSE,
    kegg_access_mode = validated_config$kegg_access_mode$value,
    kegg_cache_root = if (isTRUE(validated_config$kegg_maps$value)) lisa_config_path(validated_config$kegg_cache_root, config_dir) else "",
    kegg_snapshot_id = if (isTRUE(validated_config$kegg_maps$value)) validated_config$kegg_snapshot_id else "",
    run_hallmarks = validated_config$hallmarks$value,
    run_ora = validated_config$ora$value,
    gsea_padj_cutoff = validated_config$gsea_padj_cutoff,
    report_title = report_title,
    report_mode = "standard",
    category_evidence = report_contract$category_evidence,
    category_nes_variants = report_contract$category_nes_variants,
    evidence_max_sets = report_contract$evidence_max_sets,
    evidence_max_genes = report_contract$evidence_max_genes,
    cache_dir = pipeline$cache_dir,
    cache_mode = pipeline$cache_mode,
    cache_max_bytes = pipeline$cache_max_bytes,
    source_data = report_contract$source_data,
    recipes = report_contract$recipes,
    registered_category_map = TRUE,
    workers = validated_config$workers,
    dry_run = dry_run,
    .presentation = list(
      report_mode = report_contract$mode,
      report = list(
        mode = report_contract$mode, formats = as.list(report_contract$formats),
        source_data = report_contract$source_data, recipes = report_contract$recipes,
        category_evidence = report_contract$category_evidence,
        category_nes_variants = as.list(report_contract$category_nes_variants),
        legacy_gene_products = report_contract$legacy_gene_products,
        evidence_max_sets = report_contract$evidence_max_sets,
        evidence_max_genes = report_contract$evidence_max_genes),
      kegg_maps = isTRUE(validated_config$kegg_maps$value),
      kegg_cache_root = if (isTRUE(validated_config$kegg_maps$value))
        lisa_config_path(validated_config$kegg_cache_root, config_dir) else "",
      kegg_snapshot_id = if (isTRUE(validated_config$kegg_maps$value))
        validated_config$kegg_snapshot_id else "")
  )
  lisa_transaction_event(tx, if (isTRUE(dry_run)) "plan_computed" else "computed")
  promoted_result <- lisa_transaction_promote(tx, workers = validated_config$workers)
  promoted <- TRUE
  promoted_source_manifest <- file.path(
    promoted_result$output_dir, "source_data", "source_data_manifest.tsv"
  )
  invisible(c(result, list(config_path = config_path, output_dir = final_output_dir,
    source_manifest = promoted_source_manifest,
    plan_dir = if (isTRUE(dry_run)) promoted_result$output_dir else "",
    plan_path = if (isTRUE(dry_run)) file.path(promoted_result$output_dir, "lisa_pipeline_plan.tsv") else "",
    run_id = promoted_result$run_id, gate = promoted_result$gate)))
}

# Render native KEGG diagrams strictly as an extension of a completed canonical
# run.  The extension stages only the saved DE and GSEA tables needed to build
# the painter's gene-level evidence; it never alters the source run or reruns
# DE/GSEA.  KEGG resources are resolved by the installed painter in cache_only
# mode, which validates snapshot metadata and per-resource hashes.
lisa_render_kegg_maps_extension <- function(source_run, output_dir,
                                             kegg_cache_root,
                                             kegg_snapshot_id,
                                             analysis_ids = NULL,
                                             collections = NULL,
                                             selection_policy = "significant_only",
                                             max_pathways_per_collection = 20L) {
  source_run <- lisa_assert_run_tree_safe(source_run)
  if (!is.character(kegg_cache_root) || length(kegg_cache_root) != 1L ||
      is.na(kegg_cache_root) || !nzchar(kegg_cache_root) ||
      !is.character(kegg_snapshot_id) || length(kegg_snapshot_id) != 1L ||
      is.na(kegg_snapshot_id) || !nzchar(kegg_snapshot_id)) {
    stop("LISA-KEGG-010 cache_only KEGG maps require an explicit cache root and immutable snapshot ID.", call. = FALSE)
  }
  if (!identical(selection_policy, "significant_only")) {
    stop("LISA-KEGG-015 saved-evidence map extension supports only significant_only selection.", call. = FALSE)
  }
  if (length(max_pathways_per_collection) != 1L || is.na(max_pathways_per_collection) ||
      !is.numeric(max_pathways_per_collection) ||
      max_pathways_per_collection != as.integer(max_pathways_per_collection) ||
      max_pathways_per_collection < 1L) {
    stop("LISA-KEGG-024 max_pathways_per_collection must be one positive integer.", call. = FALSE)
  }
  max_pathways_per_collection <- as.integer(max_pathways_per_collection)
  destination <- lisa_managed_destination(output_dir, create_parent = TRUE)
  output_dir <- destination$path
  if (lisa_path_within(output_dir, source_run)) {
    stop("LISA-KEGG-019 map extension output_dir must be outside the immutable source run.", call. = FALSE)
  }
  if (lisa_path_entry_exists(output_dir)) {
    stop("LISA-KEGG-020 map extension output_dir already exists and is immutable.", call. = FALSE)
  }
  extension_id <- paste0("lisa-kegg-extension-", substr(lisa_sha256_text(paste(source_run, kegg_snapshot_id, sep = "\n")), 1L, 16L))
  staging <- lisa_short_staging_path(output_dir, extension_id)
  if (lisa_path_entry_exists(staging)) {
    stop("LISA-KEGG-021 deterministic map-extension staging directory already exists.", call. = FALSE)
  }
  staging <- lisa_run_root(staging)
  old_run_root <- getOption("lisaR.run_root", NULL)
  options(lisaR.run_root = staging)
  on.exit(options(lisaR.run_root = old_run_root), add = TRUE)
  promoted <- FALSE
  on.exit(if (!promoted && lisa_path_entry_exists(staging)) {
    lisa_guarded_delete(staging, recursive = TRUE, run_root = NULL)
  }, add = TRUE)

  work_root <- file.path(staging, "work")
  lisa_guarded_dir_create(work_root, staging)
  de_index_path <- file.path(source_run, "config", "de_index.tsv")
  if (!file.exists(de_index_path)) {
    stop("LISA-KEGG-022 saved-evidence extension requires config/de_index.tsv in the completed source run.", call. = FALSE)
  }
  de_index <- read_lisa_tsv(de_index_path)
  if (!all(c("analysis_id", "species") %in% names(de_index))) {
    stop("LISA-KEGG-022 saved-evidence extension requires analysis_id and species in config/de_index.tsv.", call. = FALSE)
  }
  package_dir <- lisa_resolve_package_dir()
  code_ledger <- lisa_write_code_identity_ledger(
    package_dir, "build_single_de_kegg_pathway_painter.R",
    file.path(work_root, "code_identity.tsv")
  )
  collection_roots <- Sys.glob(file.path(source_run, "outputs", "single_de", "*", "collection_*"))
  rows <- list()
  for (collection_root in sort(collection_roots)) {
    analysis_id <- basename(dirname(collection_root))
    collection <- sub("^collection_", "", basename(collection_root))
    if (!is.null(analysis_ids) && !analysis_id %in% as.character(analysis_ids)) next
    if (!is.null(collections) && !collection %in% as.character(collections)) next
    de_path <- file.path(collection_root, "inputs", paste0(analysis_id, "_standardized_DE.tsv"))
    gsea_paths <- Sys.glob(file.path(collection_root, "enrichment", paste0(analysis_id, "_GSEA_*_annotated.tsv")))
    if (!file.exists(de_path) || length(gsea_paths) != 1L) next
    species <- unique(as.character(de_index$species[de_index$analysis_id == analysis_id]))
    if (length(species) != 1L || is.na(species) || !nzchar(species)) {
      stop("LISA-SPECIES-005 saved-evidence map extension cannot determine one species for analysis ", analysis_id, ".", call. = FALSE)
    }
    rel <- vapply(c(de_path, gsea_paths), lisa_extension_relative, character(1), root = source_run)
    lisa_extension_stage_files(source_run, work_root, rel)
    prefix <- lisa_extension_file_label_prefix(gsea_paths[[1L]], analysis_id)
    built <- lisa_build_single_gene_level_tables(analysis_id, collection,
      file.path(work_root, "outputs"), prefix)
    if (!identical(built$status, "completed")) {
      rows[[length(rows) + 1L]] <- lisa_post_status_row("single_de_kegg_pathway_painter", analysis_id, "", collection,
        built$status, built$output_path, built$message)
      next
    }
    painted <- lisa_run_post_script(package_dir, "build_single_de_kegg_pathway_painter.R",
      c("--project-dir", work_root, "--analysis-id", analysis_id,
        "--universe", collection, "--species", species,
        "--kegg-cache-root", kegg_cache_root,
        "--kegg-snapshot-id", kegg_snapshot_id,
        "--kegg-access-mode", "cache_only",
        # Preserve the painter's established bounded product policy. The
        # earlier one-map value was proof scaffolding, not product semantics.
        "--top-pathways", as.character(max_pathways_per_collection)),
      trusted_run_root = work_root, code_ledger = code_ledger)
    rows[[length(rows) + 1L]] <- lisa_post_status_row("single_de_kegg_pathway_painter", analysis_id, "", collection,
      painted$status, painted$output_path, painted$message)
  }
  status <- if (length(rows)) do.call(rbind, rows) else data.frame()
  failed <- if (nrow(status) && "status" %in% names(status)) status$status == "failed" else logical()
  if (any(failed)) {
    detail <- paste(paste(status$analysis_id[failed], status$collection[failed],
      status$message[failed], sep = ":"), collapse = " | ")
    stop("LISA-KEGG-025 required cache-only painter unit failed: ", detail, call. = FALSE)
  }
  if (!nrow(status) || !any(status$status == "completed")) {
    detail <- if (nrow(status)) paste(paste(status$analysis_id, status$collection,
      status$status, status$message, sep = ":"), collapse = " | ") else "no eligible saved evidence"
    stop("LISA-KEGG-023 saved-evidence map extension produced no completed native painter output: ", detail, call. = FALSE)
  }
  write_lisa_tsv(status, file.path(work_root, "kegg_extension_status.tsv"))
  write_lisa_tsv(data.frame(
    extension_id = extension_id, source_run = source_run,
    kegg_access_mode = "cache_only", kegg_snapshot_id = kegg_snapshot_id,
    selection_policy = selection_policy,
    max_pathways_per_collection = max_pathways_per_collection,
    completed_painters = sum(status$status == "completed"), stringsAsFactors = FALSE
  ), file.path(work_root, "kegg_extension_receipt.tsv"))
  lisa_rebase_kegg_extension_paths(work_root, output_dir)
  report_path <- lisa_write_kegg_extension_report(
    work_root, output_dir, source_run,
    selection_policy = selection_policy,
    max_pathways_per_collection = max_pathways_per_collection
  )
  lisa_assert_run_tree_safe(work_root)
  lisa_promote_managed_directory(work_root, output_dir)
  lisa_guarded_delete(staging, recursive = TRUE, run_root = NULL)
  promoted <- TRUE
  invisible(list(output_dir = output_dir, report_path = file.path(output_dir, basename(report_path)),
    status = status))
}

# Painter scripts write absolute paths while operating in the protected staging
# tree.  Rebase only their small path-bearing receipts before the atomic rename;
# large node/gene-map audit tables that do not contain file paths are never
# loaded merely for promotion.
lisa_rebase_kegg_extension_paths <- function(work_root, output_dir) {
  # Receipts store whichever spelling of the work root the writer used. On a
  # platform where the root is reachable under an alias (macOS `/private/var`,
  # a Windows 8.3 short name) or where `TMPDIR` introduces a doubled
  # separator, rebasing only the canonical spelling silently leaves the staged
  # paths pointing at the temporary run. Substitute every known spelling,
  # longest first so an alias can never partially consume another.
  requested_spelling <- sub("[/\\\\]+$", "",
                            path.expand(as.character(work_root)))
  requested_root <- gsub("\\\\", "/", requested_spelling)
  work_root <- normalizePath(work_root, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  source_roots <- unique(c(requested_spelling, work_root,
                           lisa_path_aliases(requested_root), requested_root))
  source_roots <- source_roots[nzchar(source_roots)]
  # R's Windows path constructors can write native backslashes into receipt
  # values even though canonical managed paths use forward slashes. Match both
  # separator spellings, then always write the canonical promoted root.
  if (.Platform$OS.type == "windows") {
    source_roots <- unique(c(
      source_roots,
      chartr("/", intToUtf8(92L), source_roots)
    ))
  }
  source_roots <- source_roots[order(nchar(source_roots), decreasing = TRUE)]
  paths <- c(
    file.path(work_root, "kegg_extension_status.tsv"),
    Sys.glob(file.path(work_root, "outputs", "gene_level", "single_de", "*",
      "collection_*", "kegg_painter", "*_kegg_pathway_painter_index.tsv")),
    Sys.glob(file.path(work_root, "outputs", "gene_level", "single_de", "*",
      "collection_*", "kegg_painter", "*_kegg_pathway_painter_nodes.tsv"))
  )
  for (path in paths[file.exists(paths)]) {
    table <- read_lisa_tsv(path)
    for (column in names(table)) {
      if (is.character(table[[column]])) {
        for (source_root in source_roots) {
          table[[column]] <- gsub(source_root, output_dir, table[[column]],
                                  fixed = TRUE)
        }
      }
    }
    write_lisa_tsv(table, path)
  }
  invisible(TRUE)
}

lisa_kegg_extension_inventory <- function(work_root, output_dir) {
  work_root <- normalizePath(work_root, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  index_files <- sort(Sys.glob(file.path(work_root, "outputs", "gene_level", "single_de", "*",
    "collection_*", "kegg_painter", "*_kegg_pathway_painter_index.tsv")))
  if (!length(index_files)) {
    stop("LISA-REPORT-KEGG-005 KEGG extension has no painter index.", call. = FALSE)
  }
  rows <- list()
  for (index_path in index_files) {
    index <- read_lisa_tsv(index_path)
    required <- c("rank", "analysis_id", "universe", "kegg_id", "kegg_title",
      "output_png", "output_pdf")
    if (!nrow(index) || !all(required %in% names(index))) {
      stop("LISA-REPORT-KEGG-006 KEGG painter index is empty or incomplete: ", index_path, call. = FALSE)
    }
    ids <- as.character(index$kegg_id)
    titles <- as.character(index$kegg_title)
    if (any(is.na(ids) | !grepl("^(hsa|mmu)[0-9]{5}$", ids)) ||
        any(is.na(titles) | !nzchar(trimws(titles))) ||
        anyDuplicated(paste(index$analysis_id, index$universe, ids, sep = "\r"))) {
      stop("LISA-REPORT-KEGG-007 KEGG painter index lacks unique genuine map identifiers/titles.", call. = FALSE)
    }
    painter_dir <- dirname(index_path)
    for (i in seq_len(nrow(index))) {
      output_names <- basename(c(as.character(index$output_png[[i]]),
        as.character(index$output_pdf[[i]])))
      local_outputs <- file.path(painter_dir, output_names)
      if (any(!nzchar(output_names)) || any(!file.exists(local_outputs))) {
        stop("LISA-REPORT-KEGG-008 selected map ", ids[[i]], " (", titles[[i]],
          ") has missing local output links. Verify the declared cache snapshot and painter receipt.",
          call. = FALSE)
      }
      stem <- sub("_painted[.]png$", "", output_names[[1L]])
      sidecars <- file.path(painter_dir, paste0(stem,
        c("_painted_source.tsv", "_painted_recipe.R")))
      if (!all(file.exists(sidecars))) {
        stop("LISA-REPORT-KEGG-009 selected map ", ids[[i]], " (", titles[[i]],
          ") lacks its figure-specific source or executable recipe.", call. = FALSE)
      }
      relative <- function(path) {
        path <- normalizePath(path, winslash = "/", mustWork = TRUE)
        prefix <- paste0(work_root, "/")
        if (!startsWith(path, prefix)) {
          stop("LISA-REPORT-KEGG-010 KEGG report asset escaped its extension root.", call. = FALSE)
        }
        substring(path, nchar(prefix) + 1L)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        rank = as.integer(index$rank[[i]]),
        analysis_id = as.character(index$analysis_id[[i]]),
        collection = as.character(index$universe[[i]]),
        kegg_id = ids[[i]],
        kegg_title = titles[[i]],
        png = relative(local_outputs[[1L]]),
        pdf = relative(local_outputs[[2L]]),
        source_data = relative(sidecars[[1L]]),
        recipe = relative(sidecars[[2L]]),
        painter_index = relative(index_path),
        stringsAsFactors = FALSE
      )
    }
  }
  inventory <- do.call(rbind, rows)
  inventory <- inventory[order(inventory$analysis_id, inventory$collection, inventory$rank,
    inventory$kegg_id), , drop = FALSE]
  rownames(inventory) <- NULL
  inventory
}

lisa_write_kegg_extension_report <- function(work_root, output_dir, source_run,
                                             selection_policy = "significant_only",
                                             max_pathways_per_collection = 20L) {
  inventory <- lisa_kegg_extension_inventory(work_root, output_dir)
  source_report <- file.path(source_run, "report_index.html")
  if (!file.exists(source_report)) {
    stop("LISA-REPORT-KEGG-011 canonical source report is absent; cannot create a connected KEGG report flow.", call. = FALSE)
  }
  report_path <- file.path(work_root, "report_index.html")
  canonical_href <- lisa_presentation_relative(source_report, output_dir)
  slash <- function(x) gsub("\\\\", "/", x)
  esc <- lisa_html_escape
  map_cards <- vapply(seq_len(nrow(inventory)), function(i) {
    row <- inventory[i, , drop = FALSE]
    paste0(
      '<article class="map-card" data-kegg-id="', esc(row$kegg_id), '">',
      '<h2><code>', esc(row$kegg_id), '</code> \u2014 ', esc(row$kegg_title), '</h2>',
      '<p>', esc(row$analysis_id), ' \u00b7 ', esc(row$collection), ' \u00b7 rank ', row$rank, '</p>',
      '<a href="', esc(slash(row$png)), '"><img loading="lazy" src="', esc(slash(row$png)),
      '" alt="KEGG ', esc(row$kegg_id), ' \u2014 ', esc(row$kegg_title), '"></a>',
      '<nav><a download href="', esc(slash(row$png)), '">PNG</a> \u00b7 ',
      '<a download href="', esc(slash(row$pdf)), '">PDF</a> \u00b7 ',
      '<a download href="', esc(slash(row$source_data)), '">source data</a> \u00b7 ',
      '<a download href="', esc(slash(row$recipe)), '">recipe</a> \u00b7 ',
      '<a download href="', esc(slash(row$painter_index)), '">painter index</a></nav>',
      '</article>'
    )
  }, character(1L))
  html <- c(
    '<!doctype html>', '<html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1">',
    '<title>LISA native KEGG pathway maps</title>',
    '<style>body{font-family:system-ui,sans-serif;max-width:1200px;margin:auto;padding:2rem;color:#172033}header{border-bottom:1px solid #ccd5e0;margin-bottom:2rem}.maps{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:1.25rem}.map-card{border:1px solid #ccd5e0;border-radius:10px;padding:1rem;background:#fff}.map-card img{width:100%;height:auto;border:1px solid #e3e8ef}nav{margin-top:.7rem}code{font-size:.9em}</style>',
    '</head><body><header><p><a href="', esc(slash(canonical_href)),
    '">\u2190 Back to the canonical LISA report</a></p>',
    '<h1>Native KEGG pathway maps</h1>',
    '<p>Postcanonical cache-only extension. Maps are identified from the validated painter index, not filenames. Selection policy: <code>',
    esc(selection_policy), '</code>; maximum ', as.integer(max_pathways_per_collection),
    ' maps per analysis/collection. Every selected map is required in the declared immutable snapshot.</p></header>',
    '<main class="maps">', map_cards, '</main></body></html>'
  )
  lisa_guarded_write(report_path, function(target) writeLines(html, target, useBytes = TRUE))
  write_lisa_tsv(inventory, file.path(work_root, "kegg_report_manifest.tsv"))
  hrefs <- c(canonical_href, unlist(inventory[c("png", "pdf", "source_data", "recipe", "painter_index")],
    use.names = FALSE))
  targets <- c(source_report, file.path(work_root, hrefs[-1L]))
  link_audit <- data.frame(
    href = slash(hrefs),
    target_kind = c("canonical_report", rep("extension_asset", length(hrefs) - 1L)),
    target_exists = file.exists(targets),
    stringsAsFactors = FALSE
  )
  write_lisa_tsv(link_audit, file.path(work_root, "kegg_report_link_audit.tsv"))
  if (any(!link_audit$target_exists)) {
    stop("LISA-REPORT-KEGG-012 KEGG report contains a broken local link.", call. = FALSE)
  }
  report_path
}

lisa_config_get <- function(x, key, default = NULL) {
  if (is.null(x) || is.null(x[[key]])) return(default)
  x[[key]]
}

lisa_config_path <- function(path, base_dir) {
  if (is.null(path) || length(path) == 0 || is.na(path) || !nzchar(as.character(path))) return(NULL)
  path <- as.character(path)
  if (lisa_is_absolute_path(path)) return(path)
  normalizePath(file.path(base_dir, path), winslash = "/", mustWork = FALSE)
}

lisa_config_managed_path <- function(path, base_dir) {
  if (is.null(path) || length(path) != 1L || is.na(path) ||
      !nzchar(as.character(path)) || grepl("[[:cntrl:]]", path)) {
    stop("Unsafe configured managed destination.", call. = FALSE)
  }
  if (length(base_dir) != 1L || is.na(base_dir) ||
      !dir.exists(base_dir)) {
    stop(
      "Configured managed destination base must be one existing directory.",
      call. = FALSE
    )
  }
  path <- as.character(path)
  if (lisa_is_absolute_path(path)) return(path)

  current <- normalizePath(
    base_dir, winslash = "/", mustWork = TRUE
  )
  pieces <- strsplit(
    gsub("\\\\", "/", path), "/", fixed = TRUE
  )[[1L]]
  unresolved <- FALSE

  for (piece in pieces) {
    if (!nzchar(piece) || identical(piece, ".")) next
    if (!isTRUE(unresolved) && lisa_path_entry_exists(current) &&
        !dir.exists(current)) {
      stop(
        paste(
          "Configured managed destinations cannot traverse through a",
          "filesystem entry that is not a directory."
        ),
        call. = FALSE
      )
    }
    if (identical(piece, "..")) {
      if (isTRUE(unresolved)) {
        stop(
          paste(
            "Configured managed destinations may not traverse through a",
            "filesystem component that does not exist."
          ),
          call. = FALSE
        )
      }
      current <- dirname(current)
      next
    }

    candidate <- file.path(current, piece)
    if (!isTRUE(unresolved) && lisa_path_entry_exists(candidate)) {
      if (lisa_path_is_link(candidate)) {
        if (!file.exists(candidate) && !dir.exists(candidate)) {
          stop(
            "Configured managed destination contains a broken symbolic link: ",
            candidate,
            call. = FALSE
          )
        }
        current <- normalizePath(
          candidate, winslash = "/", mustWork = TRUE
        )
      } else {
        current <- candidate
      }
    } else {
      unresolved <- TRUE
      current <- candidate
    }
  }

  gsub("\\\\", "/", current)
}

# Use round-trip double precision rather than display precision. This also
# keeps numeric list cells semantically identical to explicit TSV cells.
lisa_config_text <- function(value) {
  if (is.null(value) || !length(value)) return("")
  if (is.list(value)) value <- unlist(value, use.names = FALSE)
  if (length(value) == 1L && is.na(value)) return("")
  if (is.numeric(value)) return(paste(sprintf("%.17g", value), collapse = ";"))
  paste(as.character(value), collapse = ";")
}

lisa_config_table <- function(rows = NULL, path = NULL, config_dir = getwd(), required = TRUE) {
  if (!is.null(path) && nzchar(as.character(path))) {
    path <- lisa_config_path(path, config_dir)
    if (!file.exists(path) || file.info(path)$size == 0) {
      return(read_lisa_tsv(path, required = required))
    }
    # IDs, labels and column names are identifiers, not numbers. In particular,
    # automatic type conversion must not turn an analysis ID "001" into "1".
    rows <- utils::read.delim(path, sep = "\t", header = TRUE, quote = "",
      comment.char = "", check.names = FALSE, stringsAsFactors = FALSE,
      colClasses = "character", na.strings = character())
  }
  if (is.null(rows)) {
    if (isTRUE(required)) stop("Config table is missing.", call. = FALSE)
    return(data.frame())
  }
  if (is.data.frame(rows)) {
    rows[] <- lapply(rows, function(column) vapply(seq_along(column),
      function(i) if (is.na(column[[i]])) "" else lisa_config_text(column[[i]]), character(1)))
    row.names(rows) <- NULL
    return(rows)
  }
  if (!is.list(rows) || length(rows) == 0) return(data.frame())
  keys <- unique(unlist(lapply(rows, names), use.names = FALSE))
  out <- do.call(rbind, lapply(rows, function(row) {
    vals <- stats::setNames(as.list(rep("", length(keys))), keys)
    for (key in names(row)) vals[[key]] <- lisa_config_text(row[[key]])
    as.data.frame(vals, stringsAsFactors = FALSE, check.names = FALSE)
  }))
  row.names(out) <- NULL
  out
}

lisa_config_resolve_path_columns <- function(df, config_dir) {
  if (is.null(df) || nrow(df) == 0) return(df)
  path_cols <- unique(c(grep("(^path$|_path$|_dir$|_file$|_excel$|_rds$)", names(df), value = TRUE),
    intersect(c("expression_matrix", "counts_matrix", "vst_matrix", "normalized_matrix", "tpm_matrix", "sample_counts"), names(df))))
  for (col in path_cols) {
    df[[col]] <- vapply(df[[col]], function(x) {
      if (is.na(x) || !nzchar(as.character(x))) return("")
      lisa_config_path(x, config_dir)
    }, character(1))
  }
  df
}

# These declarative prohibitions must win over missing input/index errors.
# They are shared by validation, planning and execution, without performing I/O.
lisa_assert_config_resource_policy <- function(pipeline) {
  if (!is.null(lisa_config_get(pipeline, "dictionary_dir", NULL)) || !is.null(lisa_config_get(pipeline, "term2gene", NULL))) stop("LISA-RESOURCE-013 free-form dictionary_dir and term2gene paths are not permitted. Repair: declare registered logical resources.", call. = FALSE)
  forbidden_resource_roots <- c("dictionary_cache_root", "shared_dictionary_root", "dictionary_registry")
  configured_roots <- forbidden_resource_roots[vapply(forbidden_resource_roots, function(key) !is.null(lisa_config_get(pipeline, key, NULL)), logical(1))]
  if (length(configured_roots)) {
    stop(sprintf("LISA-RESOURCE-015 per-run resource root(s) are forbidden: %s. Repair: select registered logical IDs in the run config; administrators must configure registry and cache roots through lisaR options or LISAR_* environment variables.", paste(configured_roots, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

# Shared semantic boundary for validation, planning and configured execution.
# Public YAML/JSON/list + TSV index forms remain unchanged. Resolution is
# read-only: it neither creates output directories nor loads scientific data.
lisa_resolve_config <- function(config, strict = TRUE, config_dir = NULL,
                                .execution = FALSE) {
  config_path <- NULL
  if (is.character(config) && length(config) == 1L && !is.na(config)) {
    config_path <- normalizePath(config, winslash = "/", mustWork = TRUE)
    raw <- read_lisa_pipeline_config(config_path)
    config_dir <- dirname(config_path)
  } else {
    if (!is.list(config) || is.null(names(config))) {
      stop("`config` must be a YAML/JSON path or a named configuration list.", call. = FALSE)
    }
    raw <- config
    if (is.null(config_dir)) config_dir <- getwd()
    config_dir <- normalizePath(config_dir, winslash = "/", mustWork = TRUE)
  }
  contract <- lisa_validate_pipeline_config(raw, strict = strict)
  pipeline <- lisa_config_get(raw, "pipeline", list())
  lisa_assert_config_resource_policy(pipeline)
  project <- lisa_config_get(raw, "project", list())
  provenance <- list()
  add <- function(field, value, source, origin = field) {
    provenance[[length(provenance) + 1L]] <<- data.frame(field = field,
      value = lisa_config_text(value), source = source, origin = origin,
      stringsAsFactors = FALSE)
    value
  }
  effective <- raw
  effective$pipeline <- pipeline
  put <- function(key, value, origin = paste0("pipeline.", key), explicit = !is.null(pipeline[[key]])) {
    effective$pipeline[key] <<- list(add(paste0("pipeline.", key), value,
      if (explicit) "explicit" else "default", origin))
  }
  put("schema_version", contract$schema_version,
    if (is.null(pipeline$schema_version)) "schema_version" else "pipeline.schema_version",
    explicit = !is.null(pipeline$schema_version) || !is.null(raw$schema_version))
  effective$schema_version <- NULL
  put("profile", contract$profile)
  put("evidence_mode", contract$evidence$mode)
  put("workers", contract$workers)
  put("gsea_padj_cutoff", contract$gsea_padj_cutoff)
  put("run_ora", contract$ora$value)
  put("run_hallmarks", contract$hallmarks$value)
  put("run_kegg_maps", contract$kegg_maps$value)
  put("kegg_access_mode", contract$kegg_access_mode$value)
  if (isTRUE(contract$kegg_maps$value)) {
    put("kegg_cache_root", lisa_config_path(contract$kegg_cache_root, config_dir))
    put("kegg_snapshot_id", contract$kegg_snapshot_id)
  }
  put("pathways_category_display", contract$pathways_category_display$value)
  put("kegg_map_selection", contract$kegg_map_selection$value)
  put("dry_run", lisa_config_bool(lisa_config_get(pipeline, "dry_run", TRUE), "pipeline.dry_run"))
  put("resume", lisa_config_bool(lisa_config_get(pipeline, "resume", FALSE), "pipeline.resume"))
  put("cache_mode", contract$cache$mode)
  # Cache roots are future managed destinations, not existing scientific
  # inputs. Resolve safe parent traversal before the strict cache write guard.
  put("cache_dir", if (is.null(contract$cache$dir)) NULL else
    lisa_config_managed_path(contract$cache$dir, config_dir))
  put("cache_max_bytes", contract$cache$max_bytes)
  for (key in c("file_label_prefix", "sufficiency_action", "sufficiency_ignore_reason")) {
    defaults <- c(file_label_prefix = "semantic",
      sufficiency_action = "warn", sufficiency_ignore_reason = "")
    put(key, as.character(lisa_config_get(pipeline, key, defaults[[key]])))
  }
  output <- lisa_config_get(pipeline, "output_dir", project$output_dir)
  if (!is.null(output)) put("output_dir", lisa_config_managed_path(output, config_dir),
    if (is.null(pipeline$output_dir)) "project.output_dir" else "pipeline.output_dir", TRUE)
  put("base_dir", lisa_config_path(lisa_config_get(pipeline, "base_dir", effective$pipeline$output_dir), config_dir))
  put("report_title", as.character(lisa_config_get(pipeline, "report_title", lisa_config_get(project, "title", "LISA report"))),
    if (is.null(pipeline$report_title) && !is.null(project$title)) "project.title" else "pipeline.report_title",
    !is.null(pipeline$report_title) || !is.null(project$title))
  if (!is.null(pipeline$lisa_project_root)) put("lisa_project_root", lisa_config_path(pipeline$lisa_project_root, config_dir))

  policies <- contract$duplicate_policies
  for (key in names(policies)) {
    policy <- policies[[key]]
    if (identical(lisa_policy_type(policy), "mapping_file")) {
      policy$mapping_file <- lisa_config_path(policy$mapping_file, config_dir)
      policies[[key]] <- policy
    }
    add(paste0("pipeline.duplicate_policies.", key), lisa_policy_type(policy), "explicit")
  }
  effective$pipeline$duplicate_policies <- contract$duplicate_policies <- policies
  effective$collections <- as.list(contract$collections)
  add("collections", contract$collections, if (is.null(raw$collections)) "default" else "explicit")
  effective$report <- contract$report
  effective$report$formats <- as.list(contract$report$formats)
  effective$report$category_nes_variants <- as.list(contract$report$category_nes_variants)
  for (key in setdiff(names(contract$report), "formats")) {
    add(paste0("report.", key), contract$report[[key]],
      if (is.null(raw$report[[key]])) "default" else "explicit")
  }
  for (key in names(contract$report$formats)) {
    add(paste0("report.formats.", key), contract$report$formats[[key]],
      if (is.null(raw$report$formats[[key]])) "default" else "explicit")
  }

  resolve_index <- function(row_key, alias, path_key, kind, required) {
    rows <- lisa_config_get(raw, row_key, raw[[alias]])
    path <- lisa_config_get(pipeline, path_key, raw[[path_key]])
    out <- lisa_config_table(rows, path, config_dir, required)
    if (isTRUE(strict)) lisa_config_assert_known_rows(out, lisa_config_allowed_keys(kind), row_key)
    out <- lisa_config_resolve_path_columns(out, config_dir)
    add(row_key, nrow(out), "explicit", if (!is.null(path)) lisa_config_path(path, config_dir) else
      if (!is.null(raw[[row_key]])) row_key else alias)
    out
  }
  # Execution adds an early immutable-destination admission check, not a
  # separate interpretation of configuration values. Read-only validation and
  # planning may still inspect a config whose final output already exists.
  if (isTRUE(.execution)) {
    destination <- effective$pipeline$output_dir
    if (is.null(destination) || !nzchar(as.character(destination))) {
      stop("Config must define pipeline.output_dir or project.output_dir.", call. = FALSE)
    }
    lisa_assert_new_final_output(destination)
  }
  de_index <- resolve_index("single_de", "de_index", "de_index_path", "analysis", TRUE)
  de_index <- validate_lisa_de_index(de_index, base_dir = config_dir, require_files = FALSE)
  contrast_index <- resolve_index("contrasts", "contrast_index", "contrast_index_path", "contrast", FALSE)
  contrast_index <- validate_lisa_contrast_index(contrast_index, de_index)
  as_rows <- function(x) lapply(seq_len(nrow(x)), function(i) as.list(x[i, , drop = FALSE]))
  effective$single_de <- as_rows(de_index)
  effective$contrasts <- as_rows(contrast_index)
  for (key in c("de_index", "de_index_path", "contrast_index", "contrast_index_path")) effective[[key]] <- NULL
  effective$pipeline$de_index_path <- effective$pipeline$contrast_index_path <- NULL
  selection <- lisa_select_pipeline_resource_references(pipeline)
  for (key in names(selection$references)) {
    ref <- selection$references[[key]]
    effective$pipeline[[key]] <- ref$resource_id
    add(paste0("pipeline.", key), ref$resource_id,
      if (is.null(pipeline[[key]]) && is.null(pipeline$lisa_dictionary)) "default" else "explicit", ref$selection_source)
  }
  effective$pipeline$lisa_dictionary <- if (selection$dictionary_tier == "custom") NULL else selection$dictionary_tier
  effective <- lisa_config_canonical_object(effective)
  contract$cache$dir <- effective$pipeline$cache_dir
  structure(list(cfg = effective, raw_config = raw, contract = contract,
    pipeline = effective$pipeline, config_path = config_path, config_dir = config_dir,
    de_index = de_index, contrast_index = contrast_index,
    provenance = do.call(rbind, provenance)), class = "lisa_resolved_config")
}

# Named object keys are unordered; unnamed arrays (including analysis and
# collection order) are not. Never sort scientific vectors or task sequences.
lisa_config_canonical_object <- function(x) {
  if (!is.list(x) || is.data.frame(x)) return(x)
  if (!is.null(names(x)) && all(nzchar(names(x)))) x <- x[sort(names(x))]
  lapply(x, lisa_config_canonical_object)
}

lisa_write_normalized_config <- function(cfg, path) {
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    lisa_guarded_write(path, function(target) writeLines(jsonlite::toJSON(cfg, auto_unbox = TRUE, pretty = TRUE, digits = 17, null = "null", na = "null"), target, useBytes = TRUE))
  }
  invisible(path)
}

lisa_source_abort <- function(code, message, source_path = NULL,
                              allowlisted_paths = NULL) {
  condition <- structure(
    list(
      message = paste(code, message),
      call = NULL,
      code = code,
      source_path = source_path,
      allowlisted_paths = allowlisted_paths
    ),
    class = c("lisa_source_data_error", "lisa_error", "error", "condition")
  )
  stop(condition)
}

lisa_prepare_config_source_data <- function(cfg, config_dir) {
  rows <- lisa_config_get(cfg, "source_data", list())
  if (!is.list(rows)) {
    lisa_source_abort(
      "LISA-SOURCE-003",
      "source_data must be a list of explicitly allowlisted file entries."
    )
  }
  if (!length(rows)) return(list())

  lapply(seq_along(rows), function(index) {
    row <- rows[[index]]
    if (!is.list(row)) {
      lisa_source_abort(
        "LISA-SOURCE-003",
        sprintf("source_data entry %d must be a named object.", index)
      )
    }
    source_value <- lisa_config_get(
      row, "path", lisa_config_get(row, "source_path", NULL)
    )
    if (length(source_value) != 1L || is.na(source_value) ||
        !nzchar(as.character(source_value))) {
      lisa_source_abort(
        "LISA-SOURCE-003",
        sprintf("source_data entry %d must declare one non-empty path.", index)
      )
    }
    src <- lisa_config_path(source_value, config_dir)
    role <- as.character(lisa_config_get(row, "role", lisa_config_get(row, "source_role", "")))
    target_subdir <- as.character(lisa_config_get(row, "target_subdir", "registered_sources"))
    row_has_allowlist <- "allowlisted_paths" %in% names(row)
    allowed_value <- if (row_has_allowlist) {
      row[["allowlisted_paths"]]
    } else {
      lisa_config_get(cfg, "allowlisted_source_paths", NULL)
    }
    allowed <- unlist(allowed_value, use.names = FALSE)
    if (!is.character(allowed) || !length(allowed) || anyNA(allowed) ||
        any(!nzchar(trimws(allowed)))) {
      lisa_source_abort(
        "LISA-SOURCE-003",
        paste0(
          "source_data entry ", index,
          " requires an explicit non-empty allowlisted_paths entry or ",
          "allowlisted_source_paths list before any source file is copied."
        ),
        source_path = as.character(src)
      )
    }
    if (!file.exists(src) || dir.exists(src)) {
      lisa_source_abort(
        "LISA-SOURCE-004",
        paste0("configured source_data file does not exist or is not a regular file: ", src),
        source_path = as.character(src),
        allowlisted_paths = allowed
      )
    }
    src <- normalizePath(src, winslash = "/", mustWork = TRUE)
    allowed <- normalizePath(
      vapply(allowed, lisa_config_path, character(1), base_dir = config_dir),
      winslash = "/", mustWork = FALSE
    )
    if (!src %in% allowed) {
      lisa_source_abort(
        "LISA-SOURCE-001",
        paste0("source_data path is not an exact allowlist match: ", src),
        source_path = src,
        allowlisted_paths = allowed
      )
    }
    lisa_safe_id(target_subdir, "source_data target_subdir")
    list(
      source_path = src,
      source_role = role,
      target_subdir = target_subdir,
      allowlisted_paths = allowed
    )
  })
}

lisa_copy_config_source_data <- function(cfg, output_dir, config_dir) {
  prepared <- lisa_prepare_config_source_data(cfg, config_dir)
  manifest_path <- file.path(output_dir, "source_data", "source_data_manifest.tsv")
  lisa_guarded_dir_create(dirname(manifest_path))
  if (!length(prepared)) {
    write_lisa_tsv(data.frame(), manifest_path)
    return(manifest_path)
  }
  # Validation is deliberately complete before this copy phase. A bad later
  # entry therefore cannot leave an earlier source file partially imported.
  manifest <- lapply(prepared, function(entry) {
    src <- entry$source_path
    dest_dir <- file.path(output_dir, "source_data", entry$target_subdir)
    lisa_guarded_dir_create(dest_dir)
    dest <- file.path(dest_dir, basename(src))
    pre_sha <- lisa_sha256_file(src)
    if (!identical(normalizePath(src, mustWork = TRUE), normalizePath(dest, mustWork = FALSE))) {
      lisa_guarded_copy(src, dest, overwrite = TRUE)
    }
    post_sha <- lisa_sha256_file(dest)
    if (!identical(pre_sha, post_sha)) {
      lisa_source_abort(
        "LISA-SOURCE-002",
        paste0("source copy checksum mismatch: ", src),
        source_path = src,
        allowlisted_paths = entry$allowlisted_paths
      )
    }
    data.frame(
      source_path = src,
      source_role = entry$source_role,
      copied_path = dest,
      bytes = file.info(dest)$size,
      source_sha256 = pre_sha,
      sha256 = post_sha,
      stringsAsFactors = FALSE
    )
  })
  write_lisa_tsv(do.call(rbind, manifest), manifest_path)
  manifest_path
}

lisa_remap_registered_source_paths <- function(df, source_manifest_path) {
  if (is.null(df) || nrow(df) == 0 || !file.exists(source_manifest_path)) return(df)
  manifest <- read_lisa_tsv(source_manifest_path)
  required <- c("source_path", "copied_path")
  if (nrow(manifest) == 0 || !all(required %in% names(manifest))) return(df)
  path_cols <- unique(c(grep("(^path$|_path$|_dir$|_file$|_excel$|_rds$)", names(df), value = TRUE),
    intersect(c("expression_matrix", "counts_matrix", "vst_matrix", "normalized_matrix", "tpm_matrix", "sample_counts"), names(df))))
  if (length(path_cols) == 0) return(df)

  normalize_key <- function(x) {
    if (is.na(x) || !nzchar(as.character(x))) return("")
    normalizePath(as.character(x), mustWork = FALSE)
  }
  source_map <- stats::setNames(
    normalizePath(manifest$copied_path, mustWork = FALSE),
    vapply(manifest$source_path, normalize_key, character(1))
  )
  for (col in path_cols) {
    df[[col]] <- vapply(df[[col]], function(x) {
      key <- normalize_key(x)
      if (nzchar(key) && key %in% names(source_map)) return(unname(source_map[[key]]))
      as.character(x)
    }, character(1))
  }
  df
}

lisa_sha256_abort <- function(code, message, paths = character()) {
  condition <- structure(
    list(
      message = paste(code, message), call = NULL, code = code,
      paths = as.character(paths)
    ),
    class = c("lisa_sha256_error", "lisa_error", "error", "condition")
  )
  stop(condition)
}

lisa_sha256_is_valid <- function(x) {
  if (!is.character(x)) return(rep(FALSE, length(x)))
  valid <- !is.na(x) & grepl("^[0-9a-f]{64}$", x)
  valid[is.na(valid)] <- FALSE
  valid
}

lisa_sha256_all_valid <- function(x, expected_length = length(x)) {
  is.character(x) && length(x) == expected_length &&
    all(lisa_sha256_is_valid(x))
}

lisa_sha256_regular_files <- function(paths, run_root = NULL) {
  if (!is.character(paths) || anyNA(paths) || any(!nzchar(paths))) {
    lisa_sha256_abort(
      "LISA-SHA256-001",
      "paths must be non-empty character values naming regular files."
    )
  }
  if (!length(paths)) return(character())
  if (!is.null(run_root)) {
    root <- lisa_assert_run_tree_safe(run_root)
    paths <- vapply(paths, lisa_guarded_path, character(1), run_root = root)
  }
  regular <- file.exists(paths) & !dir.exists(paths) & utils::file_test("-f", paths)
  regular[is.na(regular)] <- FALSE
  if (!all(regular)) {
    lisa_sha256_abort(
      "LISA-SHA256-002",
      "every path must identify an existing regular file.",
      paths[!regular]
    )
  }
  unname(normalizePath(paths, winslash = "/", mustWork = TRUE))
}

# digest is an Imports dependency because it hashes file bytes in-process on
# every supported R platform.  No operating-system command or PATH lookup is
# part of the integrity boundary.
lisa_sha256_digest_backend <- function(paths) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    lisa_sha256_abort(
      "LISA-SHA256-003",
      "the required digest package is unavailable.", paths
    )
  }
  hashes <- vapply(
    paths,
    function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE),
    character(1), USE.NAMES = FALSE
  )
  stats::setNames(hashes, paths)
}

lisa_sha256_validate_backend_result <- function(hashes, paths) {
  named_correctly <- is.character(hashes) &&
    !is.null(names(hashes)) && identical(unname(names(hashes)), paths)
  if (!named_correctly || length(hashes) != length(paths)) {
    lisa_sha256_abort(
      "LISA-SHA256-004",
      "the SHA-256 backend returned the wrong cardinality or file names.",
      paths
    )
  }
  if (!lisa_sha256_all_valid(hashes, length(paths))) {
    lisa_sha256_abort(
      "LISA-SHA256-004",
      "the SHA-256 backend returned an empty, NA, or malformed digest.",
      paths
    )
  }
  hashes
}

lisa_sha256_files_primitive <- function(paths, run_root = NULL,
                                         .backend = lisa_sha256_digest_backend) {
  paths <- lisa_sha256_regular_files(paths, run_root = run_root)
  if (!length(paths)) return(stats::setNames(character(), character()))
  if (!is.function(.backend)) {
    lisa_sha256_abort("LISA-SHA256-003", "the SHA-256 backend is not callable.", paths)
  }
  hashes <- tryCatch(
    .backend(paths),
    error = function(error) {
      if (inherits(error, "lisa_sha256_error")) stop(error)
      lisa_sha256_abort(
        "LISA-SHA256-003",
        paste0("the SHA-256 backend failed: ", conditionMessage(error)),
        paths
      )
    }
  )
  lisa_sha256_validate_backend_result(hashes, paths)
}

lisa_sha256_file <- function(path, .backend = lisa_sha256_digest_backend) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    lisa_sha256_abort(
      "LISA-SHA256-001",
      "path must be one non-empty character value naming a regular file."
    )
  }
  unname(lisa_sha256_files_primitive(path, .backend = .backend)[[1L]])
}
