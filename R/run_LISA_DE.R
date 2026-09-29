# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# LISA: LLM-Inferred Semantic Annotation
# First standalone DE -> enrichment -> LISA annotation workflow.

run_LISA_DE <- function(
    input,
    input_type = c("auto", "de_table", "deseq2_dds", "deseq2_results", "edger"),
    output_dir,
    comparison_name = NULL,
    display_title = NULL,
    comparison_subtitle = NULL,
    positive_direction = NULL,
    species = "Homo sapiens",
    msigdb_mode = c("human", "mouse_hs_ortholog", "mouse_native"),
    lisa_dictionary = "core",
    lisa_dictionary_path = NULL,
    term2gene_path = NULL,
    category_map_path = NULL,
    universes = lisa_default_collections(),
    c2_sources = c("BIOCARTA", "KEGG_LEGACY", "KEGG_MEDICUS", "PID", "REACTOME", "WIKIPATHWAYS", "CP_OTHER"),
    symbol_col = NULL,
    rank_col = NULL,
    logfc_col = NULL,
    padj_col = NULL,
    pvalue_col = NULL,
    gene_id_col = NULL,
    deseq2_contrast = NULL,
    deseq2_name = NULL,
    trusted_rds = FALSE,
    rds_max_bytes = 512 * 1024^2,
    gsea_padj_cutoff = 0.25,
    ora_padj_cutoff = 0.05,
    ora_lfc_cutoff = 0.58,
    min_gs_size = 10,
    max_gs_size = 600,
    n_threads = 1,
    fgsea_nperm = NULL,
    fgsea_eps = 0,
    random_seed = 1729L,
    run_gsea = TRUE,
    run_ora = FALSE,
    outputs = c("tables", "qc", "plots", "integrated"),
    plots = c("lollipop", "direction_lollipop", "category_pathways", "barplot"),
    plot_order = c("supracategory", "fixed", "mean_NES", "n_genesets", "consistency"),
    include_empty_categories = TRUE,
    group_by_supracategory = TRUE,
    export_formats = c("tsv", "xlsx"),
    plot_formats = c("png"),
    figure_source_data = TRUE,
    figure_recipes = TRUE,
    plot_bg = "white",
    palette = "lisa_default",
    file_label_prefix = "LISA",
    make_heatmaps = FALSE,
    dds_for_heatmaps = NULL,
    heatmap_ctrl_regex = NULL,
    heatmap_trt_regex = NULL,
    lisa_project_root = NULL,
    allow_missing_dictionary = FALSE,
    registered_category_map = FALSE,
    overwrite = TRUE,
    verbose = TRUE,
    input_receipt = NULL,
    .receipt_object_verified = FALSE,
    cache_dir = NULL,
    cache_mode = c("off", "readwrite", "readonly", "refresh"),
    cache_max_bytes = 512 * 1024^2,
    category_nes_variants = c("clean", "percentages", "direction", "dispersion")
) {
  category_nes_variants <- lisa_nes_variants_validate(category_nes_variants)
  input_source <- input
  input_type <- match.arg(input_type)
  cache_mode <- match.arg(cache_mode)
  cache_options <- lisa_content_cache_options(cache_dir, cache_mode, cache_max_bytes)
  prepared_input <- lisa_prepare_de_input(
    input = input,
    input_type = input_type,
    trusted_rds = trusted_rds,
    rds_max_bytes = rds_max_bytes,
    input_receipt = input_receipt,
    .receipt_object_verified = .receipt_object_verified
  )
  input <- prepared_input$input
  input_type <- prepared_input$input_type
  input_receipt <- prepared_input$receipt
  if (is.null(comparison_name)) comparison_name <- infer_comparison_name(input_source)
  gsea_padj_cutoff <- lisa_validate_gsea_padj_cutoff(
    gsea_padj_cutoff,
    field = "gsea_padj_cutoff"
  )
  valid_lisa_dictionaries <- c("core", "expanded", "broad_sensitivity")
  if (is.null(lisa_dictionary_path)) {
    if (length(lisa_dictionary) == 1 && identical(lisa_dictionary, "all")) {
      lisa_dictionary <- valid_lisa_dictionaries
    }
    lisa_dictionary <- normalize_lisa_dictionary(lisa_dictionary)
    if (length(lisa_dictionary) > 1) {
      bad <- setdiff(lisa_dictionary, valid_lisa_dictionaries)
      if (length(bad) > 0) {
        stop(sprintf("Unknown LISA dictionary: %s", paste(bad, collapse = ", ")), call. = FALSE)
      }
      lisa_guarded_dir_create(output_dir, output_dir)
      call_args <- as.list(environment())
      formal_names <- names(formals(run_LISA_DE))
      formal_names <- formal_names[formal_names %in% names(call_args)]
      call_args <- call_args[formal_names]
      stopifnot(!anyNA(names(call_args)), all(nzchar(names(call_args))))
      call_args$.receipt_object_verified <- !is.null(input_receipt)
      results <- list()
      for (dict in lisa_dictionary) {
        dict_output_dir <- file.path(output_dir, paste0("dictionary_", dict))
        log_msg(verbose, "Running separate LISA output for dictionary: %s", dict)
        call_args$lisa_dictionary <- dict
        call_args$output_dir <- dict_output_dir
        results[[dict]] <- do.call(run_LISA_DE, call_args)
      }
      return(invisible(results))
    }
  } else if (length(lisa_dictionary) != 1) {
    stop("When lisa_dictionary_path is provided, lisa_dictionary must be a single label.", call. = FALSE)
  }

  msigdb_mode <- match.arg(msigdb_mode)
  if (is.null(lisa_dictionary_path)) {
    lisa_dictionary <- match.arg(lisa_dictionary, valid_lisa_dictionaries)
  } else {
    lisa_dictionary <- safe_file_label(as.character(lisa_dictionary)[1])
  }
  universes <- normalize_lisa_universes(universes)
  resource_paths <- lisa_resolve_low_level_resource_paths(
    species = species,
    lisa_dictionary = lisa_dictionary,
    universes = universes,
    lisa_dictionary_path = lisa_dictionary_path,
    term2gene_path = term2gene_path,
    category_map_path = category_map_path
  )
  lisa_dictionary_path <- resource_paths$lisa_dictionary_path
  term2gene_path <- resource_paths$term2gene_path
  category_map_path <- resource_paths$category_map_path
  registered_category_map <- isTRUE(registered_category_map) ||
    isTRUE(resource_paths$registered_category_map)
  if (length(universes) > 1) {
    lisa_guarded_dir_create(output_dir, output_dir)
    call_args <- as.list(environment())
    formal_names <- names(formals(run_LISA_DE))
    formal_names <- formal_names[formal_names %in% names(call_args)]
    call_args <- call_args[formal_names]
    stopifnot(!anyNA(names(call_args)), all(nzchar(names(call_args))))
    call_args$.receipt_object_verified <- !is.null(input_receipt)
    results <- list()
    for (universe in universes) {
      universe_output_dir <- file.path(output_dir, paste0("collection_", safe_file_label(universe)))
      log_msg(verbose, "Running separate LISA output for collection: %s", universe)
      call_args$universes <- universe
      call_args$output_dir <- universe_output_dir
      results[[universe]] <- do.call(run_LISA_DE, call_args)
    }
    return(invisible(results))
  }
  plot_order <- match.arg(plot_order)
  outputs <- unique(outputs)
  plots <- normalize_lisa_plots(plots)
  export_formats <- unique(export_formats)
  plot_formats <- unique(tolower(plot_formats))
  supported_plot_formats <- c("pdf", "svg", "png", "tiff", "tif", "jpeg", "jpg")
  unsupported_plot_formats <- setdiff(plot_formats, supported_plot_formats)
  if (length(unsupported_plot_formats) > 0) {
    stop(sprintf("Unsupported plot_formats: %s", paste(unsupported_plot_formats, collapse = ", ")), call. = FALSE)
  }
  if (!is.null(fgsea_nperm)) {
    if (length(fgsea_nperm) != 1 || is.na(fgsea_nperm) || fgsea_nperm <= 0) {
      stop("fgsea_nperm must be a positive integer, or NULL to use fgseaMultilevel.", call. = FALSE)
    }
    fgsea_nperm <- as.integer(fgsea_nperm)
  }

  project_root <- lisa_project_root %||% getwd()
  universe <- universes[1]
  is_hallmarks_direct <- identical(universe, "HALLMARKS")
  term2gene_resolved <- resolve_lisa_term2gene(
    project_root = project_root,
    term2gene_path = term2gene_path,
    msigdb_mode = msigdb_mode,
    species = species,
    verbose = verbose
  )
  term2gene_path <- term2gene_resolved$path
  lisa_guarded_dir_create(output_dir, output_dir)
  dirs <- make_output_dirs(output_dir)

  if (!is_hallmarks_direct && !file.exists(lisa_dictionary_path) && isTRUE(allow_missing_dictionary)) stop("LISA-RESOURCE-012 selected dictionary is missing. Repair: install the registered resource before analysis.", call. = FALSE)

  if (!is_hallmarks_direct) require_file(lisa_dictionary_path, "LISA dictionary")
  require_file(term2gene_path, "TERM2GENE")
  if (!is_hallmarks_direct) require_file(category_map_path, "category macrogroup map")

  term2gene_raw <- read_tsv_local(term2gene_path)
  if (is_hallmarks_direct) {
    log_msg(verbose, "Building direct MSigDB HALLMARKS collection")
    hallmarks <- build_hallmarks_lisa_inputs(term2gene_raw)
    lisa_dict <- hallmarks$lisa_dict
    category_map <- hallmarks$category_map
    run_ora <- FALSE
    lisa_dictionary_path <- "MSIGDB_HALLMARKS_DIRECT"
  } else {
    log_msg(verbose, "Loading LISA dictionary: %s", lisa_dictionary_path)
    lisa_dict <- read_tsv_local(lisa_dictionary_path)
    # Transitional reader (C6): a legacy `lisa_dictionary@1` resource still has
    # the historical `LISA_score`. Drop it here, once, so everything downstream
    # sees the canonical score-free contract and an old resource produces
    # results identical to the equivalent new one.
    lisa_dict <- lisa_dictionary_runtime_projection(lisa_dict)
    lisa_dict <- lisa_dict[lisa_dict$universe %in% universes, , drop = FALSE]
    lisa_dict$gene_set_id <- as.character(lisa_dict$gene_set_id)
    lisa_dict$category_id <- as.character(lisa_dict$category_id)
    lisa_dict$source_family <- infer_source_family(lisa_dict$gene_set_id)
    keep_c2 <- lisa_dict$universe != "GOBP-C2" | lisa_dict$source_family %in% c("GO:BP", c2_sources)
    lisa_dict <- lisa_dict[keep_c2, , drop = FALSE]

    category_map <- read_tsv_local(category_map_path)
    if (isTRUE(registered_category_map)) {
      category_map <- prepare_registered_category_map(category_map, lisa_dict)
    } else {
      category_map <- augment_category_map(category_map, lisa_dict)
    }
    if (!isTRUE(registered_category_map) && identical(universe, "PATHWAYS")) {
      if ("pathway_id" %in% colnames(category_map)) {
        category_map <- category_map[
          !is.na(category_map$pathway_id) & category_map$pathway_id != "",
          , drop = FALSE
        ]
      } else {
        category_map <- category_map[category_map$category_id %in% unique(lisa_dict$category_id), , drop = FALSE]
      }
      category_map <- category_map[
        order(as.numeric(category_map$macrogroup_order),
              as.numeric(category_map$category_order_within_macrogroup),
              category_map$category_id),
        , drop = FALSE
      ]
      rownames(category_map) <- NULL
    } else if (!isTRUE(registered_category_map)) {
      category_map <- alphabetize_category_order(category_map)
    }
    category_map$color <- build_lisa_palette(category_map, palette)
  }
  write_tsv_local(category_map, file.path(dirs$qc, "lisa_category_palette.tsv"))
  write_tsv_local(unique(lisa_dict[c("gene_set_id", "category_id")]),
    file.path(dirs$qc, "lisa_category_membership.tsv"))

  term2gene <- select_lisa_term2gene_universe(
    term2gene_raw = term2gene_raw,
    universe = universe,
    c2_sources = c2_sources,
    dictionary_gene_sets = unique(lisa_dict$gene_set_id)
  )
  term2gene_coverage <- build_lisa_term2gene_coverage(lisa_dict, term2gene)
  term2gene_lisa <- term2gene[
    term2gene$gs_name %in% unique(lisa_dict$gene_set_id),
    , drop = FALSE
  ]
  if (nrow(term2gene) == 0) {
    stop(sprintf("TERM2GENE has zero gene sets for LISA universe %s: %s", universe, term2gene_path), call. = FALSE)
  }

  log_msg(verbose, "Reading DE input")
  input_type_resolved <- resolve_input_type(input, input_type)
  de <- extract_de_table(input, input_type_resolved, deseq2_contrast, deseq2_name)
  std <- standardize_de_table(
    de,
    symbol_col = symbol_col,
    rank_col = rank_col,
    logfc_col = logfc_col,
    padj_col = padj_col,
    pvalue_col = pvalue_col,
    gene_id_col = gene_id_col
  )
  comparison_name <- safe_file_label(comparison_name)
  file_label_prefix <- safe_file_label(file_label_prefix %||% "LISA")
  input_receipt_path <- NULL
  if (!is.null(input_receipt)) {
    input_receipt_path <- file.path(
      dirs$qc, paste0(comparison_name, "_trusted_RDS_input_receipt.tsv")
    )
    write_tsv_local(lisa_rds_receipt_table(input_receipt), input_receipt_path)
  }
  gsea_task_id <- paste("gsea", comparison_name, universe, sep = ":")
  runtime_manifest <- lisa_runtime_manifest(
    root_seed = random_seed,
    task_ids = gsea_task_id,
    backend = "sequential",
    workers = n_threads
  )
  runtime_values <- c(
    root_seed = runtime_manifest$root_seed,
    effective_task_seed = runtime_manifest$task_seeds$effective_seed,
    task_id = runtime_manifest$task_seeds$task_id,
    rng_kind = paste(runtime_manifest$rng_kind, collapse = ";"),
    backend = runtime_manifest$backend,
    workers = runtime_manifest$workers,
    backend_requested = runtime_manifest$backend_requested,
    backend_effective = runtime_manifest$backend_effective,
    workers_requested = runtime_manifest$workers_requested,
    workers_effective = runtime_manifest$workers_effective,
    cap_reason = runtime_manifest$cap_reason,
    outer_task_id = runtime_manifest$outer_task_id,
    outer_root_seed = runtime_manifest$outer_root_seed,
    outer_effective_task_seed = runtime_manifest$outer_task_seed,
    inner_backend = runtime_manifest$inner_backend,
    inner_workers = runtime_manifest$inner_workers,
    inner_threads = runtime_manifest$inner_threads,
    platform = runtime_manifest$platform,
    locale = runtime_manifest$locale,
    timezone = runtime_manifest$timezone,
    runtime_manifest$thread_environment
  )
  thread_keys <- paste0("thread_", names(runtime_manifest$thread_environment))
  names(runtime_values)[seq.int(length(runtime_values) - length(thread_keys) + 1L, length(runtime_values))] <- thread_keys
  write_lisa_tsv(
    data.frame(
      key = names(runtime_values),
      value = unname(runtime_values),
      stringsAsFactors = FALSE
    ),
    file.path(dirs$qc, paste0(comparison_name, "_runtime_rng.tsv"))
  )
  write_lisa_tsv(runtime_manifest$versions, file.path(dirs$qc, paste0(comparison_name, "_runtime_versions.tsv")))
  write_lisa_tsv(runtime_manifest$lock_hashes, file.path(dirs$qc, paste0(comparison_name, "_runtime_lock_hashes.tsv")))
  write_tsv_local(term2gene_coverage$summary, file.path(dirs$qc, paste0(comparison_name, "_term2gene_source_family_coverage.tsv")))
  write_tsv_local(term2gene_coverage$missing, file.path(dirs$qc, paste0(comparison_name, "_term2gene_missing_gene_sets.tsv")))

  write_tsv_local(std$de, file.path(dirs$inputs, paste0(comparison_name, "_standardized_DE.tsv")))
  write_tsv_local(std$column_mapping, file.path(dirs$qc, paste0(comparison_name, "_input_column_mapping.tsv")))

  ranks <- make_rank_vector(std$de)
  write_tsv_local(
    data.frame(symbol = names(ranks), rank_value = as.numeric(ranks), stringsAsFactors = FALSE),
    file.path(dirs$inputs, paste0(comparison_name, "_ranked_genes.tsv"))
  )

  ranking_qc <- data.frame(
    comparison = comparison_name,
    input_type = input_type_resolved,
    species = species,
    rank_col = std$column_mapping$source_column[std$column_mapping$standard_column == "rank_value"][1],
    symbol_col = std$column_mapping$source_column[std$column_mapping$standard_column == "symbol"][1],
    n_input_rows = nrow(std$de_raw),
    n_standardized_rows = nrow(std$de),
    n_ranked_symbols = length(ranks),
    n_lisa_dictionary_rows = nrow(lisa_dict),
    n_lisa_gene_sets = length(unique(lisa_dict$gene_set_id)),
    n_available_term2gene_gene_sets = length(unique(term2gene$gs_name)),
    n_classified_term2gene_gene_sets = length(unique(term2gene_lisa$gs_name)),
    n_unclassified_term2gene_gene_sets = length(setdiff(unique(term2gene$gs_name), unique(lisa_dict$gene_set_id))),
    n_missing_term2gene_gene_sets = length(setdiff(unique(lisa_dict$gene_set_id), unique(term2gene$gs_name))),
    n_term2gene_rows = nrow(term2gene),
    lisa_dictionary = lisa_dictionary,
    lisa_dictionary_path = lisa_dictionary_path,
    lisa_universe = paste(universes, collapse = ","),
    msigdb_mode = term2gene_resolved$msigdb_mode,
    msigdb_db_species = term2gene_resolved$db_species,
    msigdb_target_species = term2gene_resolved$target_species,
    has_canonical_kegg_gene_sets = any(grepl("^KEGG_", unique(term2gene$gs_name)) & !grepl("^KEGG_MEDICUS", unique(term2gene$gs_name))),
    fgsea_nperm = ifelse(is.null(fgsea_nperm), NA_integer_, as.integer(fgsea_nperm)),
    fgsea_mode = ifelse(is.null(fgsea_nperm), "multilevel", "permutation"),
    gsea_padj_cutoff = gsea_padj_cutoff,
    stringsAsFactors = FALSE
  )
  write_tsv_local(ranking_qc, file.path(dirs$qc, paste0(comparison_name, "_ranking_qc.tsv")))

  gsea_all <- data.frame()
  gsea_base <- data.frame()
  gsea_ledger <- data.frame()
  gsea_unclassified <- data.frame()
  ora_all_tested <- data.frame()
  ora_all <- data.frame()
  cache_qc <- data.frame(
    stage = "fgsea", mode = cache_mode, status = "skipped", reason = "run_gsea_false",
    cache_key = "", artifact_sha256 = "", write_status = "not_requested",
    elapsed_seconds = NA_real_, stringsAsFactors = FALSE
  )
  for (component in c("ranks", "memberships", "universe", "parameters", "rng", "versions")) {
    cache_qc[[paste0(component, "_sha256")]] <- ""
  }

  if (run_gsea) {
    log_msg(verbose, "Running GSEA")
    lisa_require_optional("fgsea", "GSEA; alternatively call run_LISA_DE(..., run_gsea = FALSE) for ORA/QC only")
    pathways <- split(term2gene$gene_symbol, term2gene$gs_name)
    pathways <- lapply(pathways, unique)
    gsea_ledger <- build_lisa_gsea_ledger(
      term2gene = term2gene,
      lisa_dict = lisa_dict,
      ranks = ranks,
      min_gs_size = min_gs_size,
      max_gs_size = max_gs_size,
      universe = universe
    )
    cached_gsea <- lisa_cached_fgsea(
      pathways, ranks, min_gs_size, max_gs_size, n_threads, fgsea_nperm, fgsea_eps,
      random_seed, gsea_task_id,
      universe = list(collection = universe, c2_sources = c2_sources,
                      species = species, msigdb_mode = term2gene_resolved$msigdb_mode),
      cache_options = cache_options
    )
    gsea_result <- cached_gsea$result
    cache_qc <- cached_gsea$qc
    gsea_notes <- attr(gsea_result, "lisa_notes")
    gsea_notes <- unique(as.character(gsea_notes[!is.na(gsea_notes) & nzchar(gsea_notes)]))
    gsea_notes <- trimws(gsub("[\r\n]+", " ", gsea_notes))
    gsea_notes_path <- file.path(dirs$qc, paste0(comparison_name, "_GSEA_notes.tsv"))
    gsea_notes_table <- data.frame(
      analysis_id = rep(comparison_name, length(gsea_notes)),
      collection = rep(universe, length(gsea_notes)),
      source = rep("fgsea", length(gsea_notes)),
      note = gsea_notes,
      stringsAsFactors = FALSE
    )
    write_tsv_local(gsea_notes_table, gsea_notes_path)
    if (length(gsea_notes) > 0L) {
      log_msg(verbose, "GSEA recorded %d note(s) in %s", length(gsea_notes), gsea_notes_path)
    }
    if (nrow(gsea_result) > 0 && "leadingEdge" %in% names(gsea_result)) {
      gsea_result$leadingEdge <- vapply(gsea_result$leadingEdge, paste, character(1), collapse = "/")
    }
    gsea_base <- complete_lisa_gsea_results(gsea_result, gsea_ledger)
    gsea_ledger <- gsea_base
    gsea_unclassified <- lisa_unclassified_gsea_rows(gsea_base)
    gsea_all <- annotate_lisa(gsea_base, lisa_dict, by_col = "pathway")
    gsea_all <- add_lisa_category_metadata(gsea_all, category_map)
    gsea_all <- add_pathway_order_metadata(gsea_all)
    write_tsv_local(gsea_all, file.path(dirs$enrichment, paste0(comparison_name, "_GSEA_", file_label_prefix, "_annotated.tsv")))
    write_tsv_local(gsea_base, file.path(dirs$qc, paste0(comparison_name, "_GSEA_universe_ledger.tsv")))
  }
  cache_qc$analysis_id <- comparison_name
  cache_qc$collection <- universe
  cache_qc_path <- file.path(dirs$qc, paste0(comparison_name, "_content_cache.tsv"))
  write_tsv_local(cache_qc, cache_qc_path)

  if (run_ora) {
    log_msg(verbose, "Running ORA")
    ora_all_tested <- run_lisa_ora(
      std$de, term2gene_lisa, lisa_dict, ora_padj_cutoff, ora_lfc_cutoff,
      analysis_id = comparison_name, collection = universe
    )
    ora_all_tested <- add_lisa_category_metadata(ora_all_tested, category_map)
    if (nrow(ora_all_tested) > 0) ora_all_tested <- add_pathway_order_metadata(ora_all_tested)
    write_tsv_local(
      ora_all_tested,
      file.path(dirs$enrichment, paste0(comparison_name, "_ORA_all_tested_annotated.tsv"))
    )
    # Keep the established ORA artifact as the post-BH significant view used
    # by downstream LISA summaries and gene-level products.
    ora_all <- ora_all_tested[ora_all_tested$selected %in% TRUE, , drop = FALSE]
    write_tsv_local(ora_all, file.path(dirs$enrichment, paste0(comparison_name, "_ORA_", file_label_prefix, "_annotated.tsv")))
  }

  gsea_classified <- lisa_classified_gsea_rows(gsea_all)
  summaries <- build_lisa_summaries(
    gsea_all = gsea_classified,
    ora_all = ora_all,
    category_map = category_map,
    include_empty_categories = include_empty_categories,
    gsea_padj_cutoff = gsea_padj_cutoff,
    ora_padj_cutoff = ora_padj_cutoff,
    gsea_ledger = gsea_ledger
  )
  write_tsv_local(summaries$gsea, file.path(dirs$lisa_tables, paste0(comparison_name, "_", file_label_prefix, "_GSEA_category_summary.tsv")))
  write_tsv_local(summaries$ora, file.path(dirs$lisa_tables, paste0(comparison_name, "_", file_label_prefix, "_ORA_category_summary.tsv")))
  write_tsv_local(summaries$integrated, file.path(dirs$lisa_tables, paste0(comparison_name, "_", file_label_prefix, "_integrated_summary.tsv")))
  write_tsv_local(summaries$coverage, file.path(dirs$qc, paste0(comparison_name, "_lisa_mapping_coverage.tsv")))
  residual_path <- file.path(
    dirs$lisa_tables,
    paste0(comparison_name, "_", file_label_prefix, "_GSEA_OTHER_UNCLASSIFIED.tsv")
  )
  legacy_unmapped_path <- file.path(dirs$qc, paste0(comparison_name, "_unmapped_genesets.tsv"))
  exports_residual_table <- isTRUE(run_gsea) &&
    universe %in% c("GOBP-C2", "GOMF", "GOCC")
  if (isTRUE(exports_residual_table)) {
    write_tsv_local(gsea_unclassified, residual_path)
    write_tsv_local(gsea_unclassified, legacy_unmapped_path)
  } else if (isTRUE(overwrite)) {
    for (stale in c(residual_path, legacy_unmapped_path)) {
      if (file.exists(stale)) lisa_guarded_delete(stale, run_root = output_dir)
    }
  }

  if ("xlsx" %in% export_formats) {
    lisa_require_optional("openxlsx", "XLSX export")
    wb <- openxlsx::createWorkbook()
    add_wb_sheet(wb, "Ranking_QC", ranking_qc)
    add_wb_sheet(wb, "GSEA_LISA", gsea_all)
    add_wb_sheet(wb, "ORA_LISA", ora_all)
    add_wb_sheet(wb, "GSEA_Category", summaries$gsea)
    add_wb_sheet(wb, "ORA_Category", summaries$ora)
    add_wb_sheet(wb, "Integrated", summaries$integrated)
    add_wb_sheet(wb, "Coverage", summaries$coverage)
    if (isTRUE(exports_residual_table)) {
      add_wb_sheet(wb, "GSEA_Unclassified", gsea_unclassified)
    }
    openxlsx::saveWorkbook(wb, file.path(output_dir, paste0(comparison_name, "_", file_label_prefix, "_results.xlsx")), overwrite = overwrite)
  }

  nes_variant_manifest <- lisa_nes_runtime_status(category_nes_variants,
    if (!isTRUE(run_gsea)) "not_requested_gsea" else if (is_hallmarks_direct)
      "not_applicable_direct_gene_sets" else "not_requested_plot_output")
  if ("plots" %in% outputs) {
    lisa_require_optional("ggplot2", "plot output")
    log_msg(verbose, "Writing plots")
    if (isTRUE(run_gsea) && !is_hallmarks_direct) {
      prepared_nes <- lisa_prepare_category_nes_variants(
        prepare_category_plot_df(summaries$gsea, plot_order), gsea_classified,
        gsea_padj_cutoff = gsea_padj_cutoff, title = display_title %||% comparison_name,
        subtitle = paste(c(comparison_subtitle, positive_direction), collapse = "\n"))
      nes_variant_manifest <- lisa_write_runtime_nes_variants(prepared_nes,
        dirs, comparison_name, universe, file_label_prefix, category_nes_variants,
        plot_formats, scope = "single")
    }
    if ("lollipop" %in% plots && nrow(summaries$gsea) > 0) {
      p <- plot_lisa_gsea_lollipop(summaries$gsea, group_by_supracategory, plot_order)
      p <- lisa_apply_plot_context(p, display_title, comparison_subtitle, positive_direction)
      save_plot_multi(p, file.path(dirs$plots, paste0(comparison_name, "_", file_label_prefix, "_GSEA_lollipop")), plot_formats, width = 11, height = 12, bg = plot_bg)
      if (isTRUE(figure_source_data)) lisa_write_lollipop_figure_source(
        summaries$gsea, file.path(dirs$plots, paste0(comparison_name, "_", file_label_prefix, "_GSEA_lollipop_source.tsv")),
        group_by_supracategory, plot_order, display_title, comparison_subtitle,
        positive_direction, plot_bg = plot_bg)
    }
    if ("direction_lollipop" %in% plots && nrow(summaries$gsea) > 0) {
      p <- plot_lisa_gsea_direction_lollipop(summaries$gsea, group_by_supracategory, plot_order)
      p <- lisa_apply_plot_context(p, display_title, comparison_subtitle, positive_direction)
      save_plot_multi(p, file.path(dirs$plots, paste0(comparison_name, "_", file_label_prefix, "_GSEA_lollipop_direction_stats")), plot_formats, width = 11.5, height = 12, bg = plot_bg)
      if (isTRUE(figure_source_data)) lisa_write_lollipop_figure_source(
        summaries$gsea, file.path(dirs$plots, paste0(comparison_name, "_", file_label_prefix, "_GSEA_lollipop_direction_stats_source.tsv")),
        group_by_supracategory, plot_order, display_title, comparison_subtitle,
        positive_direction, annotation_variant = "direction", figure_width = 11.5,
        plot_bg = plot_bg)
    }
    if ("dumbbell" %in% plots && nrow(summaries$gsea) > 0) {
      p <- plot_lisa_gsea_lollipop(summaries$gsea, group_by_supracategory, plot_order)
      p <- lisa_apply_plot_context(p, display_title, comparison_subtitle, positive_direction)
      save_plot_multi(p, file.path(dirs$plots, paste0(comparison_name, "_", file_label_prefix, "_GSEA_dumbbell")), plot_formats, width = 11, height = 12, bg = plot_bg)
      if (isTRUE(figure_source_data)) lisa_write_lollipop_figure_source(
        summaries$gsea, file.path(dirs$plots, paste0(comparison_name, "_", file_label_prefix, "_GSEA_dumbbell_source.tsv")),
        group_by_supracategory, plot_order, display_title, comparison_subtitle,
        positive_direction, plot_bg = plot_bg)
    }
    if ("barplot" %in% plots && nrow(summaries$ora) > 0) {
      p <- plot_lisa_ora_barplot(summaries$ora, group_by_supracategory, plot_order)
      p <- lisa_apply_plot_context(p, display_title, comparison_subtitle, positive_direction)
      save_plot_multi(p, file.path(dirs$plots, paste0(comparison_name, "_", file_label_prefix, "_ORA_barplot")), plot_formats, width = 11, height = 12, bg = plot_bg)
    }
    if ("dotplot" %in% plots && nrow(gsea_classified) > 0) {
      p <- plot_lisa_pathway_dotplot(gsea_classified, gsea_padj_cutoff, group_by_supracategory, plot_order)
      p <- lisa_apply_plot_context(p, display_title, comparison_subtitle, positive_direction)
      save_plot_multi(p, file.path(dirs$plots, paste0(comparison_name, "_", file_label_prefix, "_GSEA_pathway_dotplot")), plot_formats, width = 12, height = 10, bg = plot_bg)
    }
    if ("category_pathways" %in% plots && nrow(gsea_classified) > 0 && nrow(summaries$gsea) > 0) {
      write_lisa_gsea_category_pathway_plots(
        gsea_all = gsea_classified,
        category_summary = summaries$gsea,
        gsea_padj_cutoff = gsea_padj_cutoff,
        plot_order = plot_order,
        output_dir = file.path(dirs$plots, paste0(comparison_name, "_", file_label_prefix, "_GSEA_category_pathways")),
        plot_formats = plot_formats,
        bg = plot_bg,
        comparison_subtitle = comparison_subtitle,
        positive_direction = positive_direction,
        source_data = figure_source_data,
        recipes = figure_recipes
      )
    }
  }

  write_tsv_local(nes_variant_manifest,
    file.path(dirs$qc, paste0(comparison_name, "_category_nes_variants.tsv")))

  write_run_report(
    path = file.path(output_dir, paste0(comparison_name, "_", file_label_prefix, "_RUN_REPORT.txt")),
    comparison_name = comparison_name,
    input_type = input_type_resolved,
    dictionary = lisa_dictionary_path,
    term2gene = term2gene_path,
    term2gene_resolved = term2gene_resolved,
    ranking_qc = ranking_qc,
    gsea_rows = nrow(gsea_all),
    gsea_universe_gene_sets = nrow(gsea_ledger),
    gsea_unclassified_gene_sets = nrow(gsea_unclassified),
    ora_rows = nrow(ora_all),
    outputs = dirs
  )

  if (make_heatmaps) {
    warning("make_heatmaps is reserved for the next slice. No heatmaps were generated in this first run_LISA_DE implementation.", call. = FALSE)
  }

  invisible(list(
    comparison_name = comparison_name,
    output_dir = output_dir,
    ranking_qc = ranking_qc,
    gsea = gsea_all,
    gsea_ledger = gsea_ledger,
    gsea_unclassified = gsea_unclassified,
    ora_all_tested = ora_all_tested,
    ora = ora_all,
    summaries = summaries,
    nes_variant_manifest = nes_variant_manifest,
    content_cache = cache_qc,
    content_cache_path = cache_qc_path,
    input_receipt = input_receipt,
    input_receipt_path = input_receipt_path
  ))
}

run_LISA_contrast <- function(
    contrast_a,
    contrast_b,
    output_dir,
    comparison_name = NULL,
    contrast_a_label = "Contrast A",
    contrast_b_label = "Contrast B",
    contrast_a_title = contrast_a_label,
    contrast_b_title = contrast_b_label,
    comparison_subtitle = NULL,
    positive_direction = NULL,
    universes = lisa_default_collections(),
    plot_sets = c("all", "same_direction", "opposite_direction"),
    plot_order = c("supracategory", "fixed", "delta_abs", "mean_abs"),
    include_missing_categories = TRUE,
    direction_epsilon = 0,
    group_by_supracategory = TRUE,
    export_formats = c("tsv", "xlsx"),
    plot_formats = c("png"),
    plot_bg = "white",
    annotate_gene_sets = TRUE,
    gsea_padj_cutoff = 0.25,
    figure_source_data = TRUE,
    figure_recipes = TRUE,
    file_label_prefix = "LISA",
    plot_title_prefix = "LISA category contrast dumbbell",
    overwrite = TRUE,
    verbose = TRUE,
    category_nes_variants = c("clean", "percentages", "direction", "dispersion")
) {
  category_nes_variants <- lisa_nes_variants_validate(category_nes_variants)
  gsea_padj_cutoff <- lisa_validate_gsea_padj_cutoff(
    gsea_padj_cutoff,
    field = "gsea_padj_cutoff"
  )
  universes <- normalize_lisa_universes(universes)
  plot_order <- match.arg(plot_order)
  plot_sets <- normalize_lisa_contrast_plot_sets(plot_sets)
  export_formats <- unique(tolower(export_formats))
  plot_formats <- unique(tolower(plot_formats))
  supported_plot_formats <- c("pdf", "svg", "png", "tiff", "tif", "jpeg", "jpg")
  unsupported_plot_formats <- setdiff(plot_formats, supported_plot_formats)
  if (length(unsupported_plot_formats) > 0) {
    stop(sprintf("Unsupported plot_formats: %s", paste(unsupported_plot_formats, collapse = ", ")), call. = FALSE)
  }
  if (is.null(comparison_name)) {
    comparison_name <- paste(safe_file_label(contrast_a_label), "vs", safe_file_label(contrast_b_label), sep = "_")
  }
  comparison_name <- safe_file_label(comparison_name)
  file_label_prefix <- safe_file_label(file_label_prefix %||% "LISA")
  plot_title_prefix <- plot_title_prefix %||% "LISA category contrast dumbbell"
  subtitle_parts <- c(
    if (!is.null(comparison_subtitle) && nzchar(as.character(comparison_subtitle[[1]]))) {
      as.character(comparison_subtitle[[1]])
    } else {
      sprintf("%s minus %s (A - B)", contrast_a_title, contrast_b_title)
    },
    if (!is.null(positive_direction) && nzchar(as.character(positive_direction[[1]]))) {
      as.character(positive_direction[[1]])
    }
  )
  contrast_subtitle <- paste(subtitle_parts, collapse = "\n")
  plot_contrast_subtitle <- lisa_plot_device_text(contrast_subtitle)

  if (length(universes) > 1) {
    lisa_guarded_dir_create(output_dir, output_dir)
    call_args <- as.list(environment())
    call_args <- call_args[names(formals(run_LISA_contrast))]
    results <- list()
    for (universe in universes) {
      universe_output_dir <- file.path(output_dir, paste0("collection_", safe_file_label(universe)))
      log_msg(verbose, "Running LISA contrast for collection: %s", universe)
      call_args$universes <- universe
      call_args$output_dir <- universe_output_dir
      results[[universe]] <- do.call(run_LISA_contrast, call_args)
    }
    return(invisible(results))
  }
  universe <- universes[1]

  lisa_guarded_dir_create(output_dir, output_dir)
  dirs <- make_output_dirs(output_dir)
  log_msg(verbose, "Loading LISA category summaries for contrast: %s", comparison_name)
  path_a <- find_lisa_category_summary(contrast_a, universe = universe, analysis = "GSEA")
  path_b <- find_lisa_category_summary(contrast_b, universe = universe, analysis = "GSEA")
  summary_a <- read_lisa_contrast_summary(path_a, contrast_a_label, gsea_padj_cutoff)
  summary_b <- read_lisa_contrast_summary(path_b, contrast_b_label, gsea_padj_cutoff)

  contrast_summary <- build_lisa_contrast_summary(
    summary_a = summary_a,
    summary_b = summary_b,
    contrast_a_label = contrast_a_label,
    contrast_b_label = contrast_b_label,
    universe = universe,
    include_missing_categories = include_missing_categories,
    direction_epsilon = direction_epsilon,
    gsea_padj_cutoff = gsea_padj_cutoff
  )
  # The canonical table and all-category plot retain category-map order. The
  # optional signal-only views may still use the requested effect-size order.
  contrast_summary <- order_lisa_contrast_summary(contrast_summary, "supracategory")

  table_stem <- file.path(dirs$lisa_tables, paste0(comparison_name, "_", file_label_prefix, "_GSEA_category_contrast"))
  write_tsv_local(contrast_summary, paste0(table_stem, ".tsv"))

  manifest <- data.frame()
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    for (plot_set in plot_sets) {
      plot_df <- filter_lisa_contrast_plot_set(contrast_summary, plot_set)
      if (!identical(plot_set, "all")) {
        plot_df <- order_lisa_contrast_summary(plot_df, plot_order)
      }
      plot_title <- lisa_contrast_plot_title(plot_set, universe, title_prefix = plot_title_prefix)
      stem <- file.path(dirs$plots,
        lisa_contrast_plot_stem(comparison_name, file_label_prefix, plot_set))
      height <- lisa_contrast_plot_height(nrow(plot_df))
      p_plain <- plot_lisa_contrast_dumbbell(
        contrast_summary = plot_df,
        contrast_a_label = contrast_a_label,
        contrast_b_label = contrast_b_label,
        contrast_a_title = contrast_a_title,
        contrast_b_title = contrast_b_title,
        title = plot_title,
        subtitle = plot_contrast_subtitle,
        group_by_supracategory = group_by_supracategory,
        annotate_gene_sets = FALSE
      )
      save_plot_multi(p_plain, stem, plot_formats, width = 12, height = height, bg = plot_bg)
      source_stem <- paste0(stem, "_source.tsv")
      recipe_stem <- paste0(stem, "_recipe.R")
      if (isTRUE(figure_source_data)) {
        lisa_write_contrast_figure_source(
          plot_df, source_stem, comparison_name, plot_title,
          plot_contrast_subtitle, plot_set, "plain", contrast_a_label,
          contrast_b_label, gsea_padj_cutoff,
          group_by_supracategory = group_by_supracategory
        )
      }
      if (isTRUE(figure_recipes)) {
        lisa_install_figure_recipe(recipe_stem)
      }
      manifest <- rbind(
        manifest,
        data.frame(
          comparison_name = comparison_name,
          universe = universe,
          plot_set = plot_set,
          annotation_variant = "plain",
          n_categories = nrow(plot_df),
          n_categories_with_signal = sum(plot_df$plot_has_any_significant_support, na.rm = TRUE),
          gsea_padj_cutoff = gsea_padj_cutoff,
          plot_stem = stem,
          source_tsv = if (isTRUE(figure_source_data)) source_stem else "",
          recipe_r = if (isTRUE(figure_recipes)) recipe_stem else "",
          stringsAsFactors = FALSE
        )
      )
      if (isTRUE(annotate_gene_sets)) {
        p_annotated <- plot_lisa_contrast_dumbbell(
          contrast_summary = plot_df,
          contrast_a_label = contrast_a_label,
          contrast_b_label = contrast_b_label,
          contrast_a_title = contrast_a_title,
          contrast_b_title = contrast_b_title,
          title = plot_title,
          subtitle = plot_contrast_subtitle,
          group_by_supracategory = group_by_supracategory,
          annotate_gene_sets = TRUE
        )
        annotated_stem <- paste0(stem, "_annotated")
        save_plot_multi(p_annotated, annotated_stem, plot_formats, width = 15.5, height = height * 1.08, bg = plot_bg)
        annotated_source <- paste0(annotated_stem, "_source.tsv")
        annotated_recipe <- paste0(annotated_stem, "_recipe.R")
        if (isTRUE(figure_source_data)) {
          lisa_write_contrast_figure_source(
            plot_df, annotated_source, comparison_name, plot_title,
            plot_contrast_subtitle, plot_set, "annotated", contrast_a_label,
            contrast_b_label, gsea_padj_cutoff,
            group_by_supracategory = group_by_supracategory
          )
        }
        if (isTRUE(figure_recipes)) {
          lisa_install_figure_recipe(annotated_recipe)
        }
        manifest <- rbind(
          manifest,
          data.frame(
            comparison_name = comparison_name,
            universe = universe,
            plot_set = plot_set,
            annotation_variant = "annotated",
            n_categories = nrow(plot_df),
            n_categories_with_signal = sum(plot_df$plot_has_any_significant_support, na.rm = TRUE),
            gsea_padj_cutoff = gsea_padj_cutoff,
            plot_stem = annotated_stem,
            source_tsv = if (isTRUE(figure_source_data)) annotated_source else "",
            recipe_r = if (isTRUE(figure_recipes)) annotated_recipe else "",
            stringsAsFactors = FALSE
          )
        )
      }
    }
  } else {
    warning("Package 'ggplot2' is required to write LISA contrast plots. Tables were still written.", call. = FALSE)
  }
  write_tsv_local(manifest, file.path(dirs$qc, paste0(comparison_name, "_plot_manifest.tsv")))

  # Supplementary descriptive variants reuse member NES values, not GSEA.
  # Summary-only legacy inputs remain supported, but cannot yield honest
  # medians/quantiles without their corresponding member tables.
  nes_variant_manifest <- lisa_nes_runtime_status(category_nes_variants,
    "unavailable_missing_member_table")
  if (identical(universe, "HALLMARKS")) {
    nes_variant_manifest$status <- "not_applicable_direct_gene_sets"
  } else if (!requireNamespace("ggplot2", quietly = TRUE)) {
    nes_variant_manifest$status <- "unavailable_plot_dependency"
  } else {
    members_a <- lisa_nes_members_for_summary(path_a, universe)
    members_b <- lisa_nes_members_for_summary(path_b, universe)
    if (!is.null(members_a) && !is.null(members_b)) {
      aligned_a <- lisa_nes_align_contrast_side(summary_a, contrast_summary)
      aligned_b <- lisa_nes_align_contrast_side(summary_b, contrast_summary)
      prepared_nes <- lisa_prepare_contrast_nes_variants(
        contrast_summary, aligned_a, aligned_b, members_a, members_b,
        gsea_padj_cutoff = gsea_padj_cutoff,
        label_a = contrast_a_title, label_b = contrast_b_title,
        title = paste0(plot_title_prefix, " - ", universe),
        subtitle = plot_contrast_subtitle)
      nes_variant_manifest <- lisa_write_runtime_nes_variants(prepared_nes,
        dirs, comparison_name, universe, file_label_prefix, category_nes_variants,
        plot_formats, scope = "contrast", plot_sets = plot_sets)
    } else {
      log_msg(verbose, "NES distribution variants unavailable: exact annotated GSEA member table missing for a summary-only contrast input.")
    }
  }
  write_tsv_local(nes_variant_manifest,
    file.path(dirs$qc, paste0(comparison_name, "_category_nes_variants.tsv")))

  input_manifest <- data.frame(
    comparison_name = comparison_name,
    universe = universe,
    contrast_a_label = contrast_a_label,
    contrast_b_label = contrast_b_label,
    contrast_a_title = contrast_a_title,
    contrast_b_title = contrast_b_title,
    comparison_subtitle = contrast_subtitle,
    positive_direction = as.character(positive_direction %||% ""),
    contrast_a_input = as.character(contrast_a),
    contrast_b_input = as.character(contrast_b),
    contrast_a_summary = path_a,
    contrast_b_summary = path_b,
    direction_epsilon = direction_epsilon,
    annotate_gene_sets = annotate_gene_sets,
    gsea_padj_cutoff = gsea_padj_cutoff,
    figure_source_data = figure_source_data,
    figure_recipes = figure_recipes,
    category_nes_variants = paste(category_nes_variants, collapse = ","),
    stringsAsFactors = FALSE
  )
  write_tsv_local(input_manifest, file.path(dirs$qc, paste0(comparison_name, "_input_manifest.tsv")))

  if ("xlsx" %in% export_formats) {
    lisa_require_optional("openxlsx", "XLSX contrast export")
    wb <- openxlsx::createWorkbook()
    add_wb_sheet(wb, "Category_Contrast", contrast_summary)
    add_wb_sheet(wb, "Input_Manifest", input_manifest)
    add_wb_sheet(wb, "Plot_Manifest", manifest)
    openxlsx::saveWorkbook(wb, file.path(output_dir, paste0(comparison_name, "_", file_label_prefix, "_contrast_results.xlsx")), overwrite = overwrite)
  }

  write_lisa_contrast_report(
    path = file.path(output_dir, paste0(comparison_name, "_", file_label_prefix, "_CONTRAST_REPORT.txt")),
    comparison_name = comparison_name,
    universe = universe,
    input_manifest = input_manifest,
    contrast_summary = contrast_summary,
    plot_sets = plot_sets,
    outputs = dirs
  )

  invisible(list(
    comparison_name = comparison_name,
    output_dir = output_dir,
    universe = universe,
    summary = contrast_summary,
    manifest = manifest,
    nes_variant_manifest = nes_variant_manifest,
    input_manifest = input_manifest
  ))
}

# Runtime wiring for the shared NES presentation module. All paths recorded in
# the supplemental manifest are relative to the collection output directory.
lisa_nes_runtime_stem <- function(comparison, prefix, contrast = FALSE) {
  # Exact identifiers stay in metadata. Bound presentation basenames so the
  # extra view/subset/settings suffixes remain portable in deep report trees.
  identity <- jsonlite::toJSON(list(comparison = comparison, prefix = prefix,
    scope = if (contrast) "contrast" else "single"), auto_unbox = TRUE)
  label <- gsub("[^A-Za-z0-9_-]", "_", as.character(comparison))
  paste0(substr(label, 1L, 8L), "_nes_",
    substr(digest::digest(as.character(identity), algo = "sha256", serialize = FALSE), 1L, 12L))
}

lisa_nes_runtime_status <- function(variants, status) {
  data.frame(variant = variants, status = status, format = "", path = "",
    source = "", settings = "", recipe = "", plot_set = "all", stringsAsFactors = FALSE)
}

lisa_write_runtime_nes_variants <- function(prepared, dirs, comparison_name,
    universe, file_label_prefix, variants, plot_formats, scope = "single", plot_sets = "all") {
  variants <- lisa_nes_variants_validate(variants)
  formats <- intersect(plot_formats, c("png", "pdf", "svg"))
  if (!length(formats)) return(lisa_nes_runtime_status(variants, "not_requested_supported_format"))
  prefix <- lisa_nes_runtime_stem(comparison_name, file_label_prefix, identical(scope, "contrast"))
  folder <- file.path(dirs$plots, "category_nes")
  lisa_guarded_dir_create(folder)
  prepared$metadata$analysis_id <- comparison_name
  prepared$metadata$collection <- universe
  prepared$metadata$category_nes_variants <- variants
  files <- render_lisa_category_nes_variants(prepared, folder,
    variants = variants, prefix = prefix, formats = formats, plot_sets = plot_sets)
  recipe_path <- file.path(folder, paste0(prefix, "_recipe.R"))
  renderer <- lisa_post_script_path(lisa_resolve_package_dir(),
    "reproduce_lisa_category_nes.R")
  lisa_copy_verified_code_file(renderer, recipe_path, lisa_sha256_file(renderer))
  files$status <- "rendered"
  files$recipe <- recipe_path
  for (field in c("path", "source", "settings", "recipe")) {
    files[[field]] <- file.path("plots", "category_nes", basename(files[[field]]))
  }
  files[, c("variant", "status", "format", "path", "source", "settings", "recipe", "plot_set"), drop = FALSE]
}

lisa_nes_members_for_summary <- function(summary_path, universe) {
  # Do not pair an arbitrary neighbouring enrichment table with a summary.
  # Only the exact canonical file naming relationship is accepted.
  root <- dirname(dirname(summary_path))
  folder <- file.path(root, "enrichment")
  if (!identical(basename(dirname(summary_path)), "lisa_tables") || !dir.exists(folder)) return(NULL)
  paths <- list.files(folder, pattern = "_GSEA_.+_annotated\\.tsv$", full.names = TRUE)
  matching_summary <- sub("_GSEA_(.+)_annotated\\.tsv$", "_\\1_GSEA_category_summary.tsv", basename(paths))
  paths <- paths[matching_summary == basename(summary_path)]
  if (!length(paths)) return(NULL)
  if (length(paths) != 1L) stop("Ambiguous annotated GSEA table for NES variants.", call. = FALSE)
  members <- utils::read.delim(paths[[1L]], sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character",
    na.strings = character())
  if (!all(c("category_id", "pathway", "NES", "padj") %in% names(members)))
    stop("Annotated GSEA members lack the NES-variant contract.", call. = FALSE)
  classified <- !is.na(members$category_id) & nzchar(members$category_id) &
    members$category_id != "OTHER_UNCLASSIFIED"
  if ("universe" %in% names(members) && any(classified &
      (is.na(members$universe) | members$universe != universe)))
    stop("Annotated GSEA members disagree with the contrast collection.", call. = FALSE)
  members
}

lisa_nes_align_contrast_side <- function(summary, contrast_summary) {
  ids <- as.character(contrast_summary$category_id)
  index <- match(ids, as.character(summary$category_id))
  side <- summary[index, , drop = FALSE]
  missing <- is.na(index)
  side$category_id <- ids
  # The existing full-outer contrast already distinguishes an absent category
  # from a contextual endpoint. Fill only absence counts; never invent NES=0.
  if (any(missing)) {
    for (field in intersect(c("n_genesets", "n_genesets_mapped", "n_genesets_evaluable",
      "n_genesets_significant", "n_pos_genesets", "n_neg_genesets"), names(side)))
      side[[field]][missing] <- 0L
    for (field in intersect(c("mean_NES", "mean_NES_contextual", "median_NES", "min_padj"), names(side)))
      side[[field]][missing] <- NA_real_
    for (field in intersect(c("category_display_name", "macrogroup_name"), names(side)))
      side[[field]][missing] <- as.character(contrast_summary[[field]][missing])
  }
  rownames(side) <- NULL
  side
}

`%||%` <- function(x, y) if (is.null(x)) y else x

normalize_lisa_dictionary <- function(x) {
  x <- as.character(x)
  x <- ifelse(tolower(x) %in% c("complete", "full", "completo"), "broad_sensitivity", x)
  x
}

normalize_lisa_universes <- function(x) {
  normalize_lisa_collections(x)
}

normalize_lisa_plots <- function(x) {
  x <- unique(tolower(as.character(x)))
  aliases <- c(
    directional_lollipop = "direction_lollipop",
    direction_stats = "direction_lollipop",
    lollipop_stats = "direction_lollipop",
    stats_lollipop = "direction_lollipop",
    pathway_dotplot = "dotplot",
    gsea_pathway_dotplot = "dotplot",
    category_subplots = "category_pathways",
    category_pathway_subplots = "category_pathways",
    category_pathway_plots = "category_pathways",
    gsea_category_pathways = "category_pathways",
    gsea_category_subplots = "category_pathways",
    ora_barplot = "barplot"
  )
  x[x %in% names(aliases)] <- unname(aliases[x[x %in% names(aliases)]])
  valid <- c("lollipop", "direction_lollipop", "dumbbell", "barplot", "dotplot", "category_pathways")
  bad <- setdiff(x, valid)
  if (length(bad) > 0) {
    stop(sprintf("Unsupported plots: %s", paste(bad, collapse = ", ")), call. = FALSE)
  }
  unique(x)
}

normalize_lisa_contrast_plot_sets <- function(x) {
  x <- unique(tolower(as.character(x)))
  aliases <- c(
    all_categories = "all",
    every_category = "all",
    same = "same_direction",
    concordant = "same_direction",
    same_sign = "same_direction",
    same_behavior = "same_direction",
    opposite = "opposite_direction",
    discordant = "opposite_direction",
    opposite_sign = "opposite_direction",
    opposite_behavior = "opposite_direction"
  )
  x[x %in% names(aliases)] <- unname(aliases[x[x %in% names(aliases)]])
  valid <- c("all", "same_direction", "opposite_direction")
  bad <- setdiff(x, valid)
  if (length(bad) > 0) {
    stop(sprintf("Unsupported LISA contrast plot_sets: %s", paste(bad, collapse = ", ")), call. = FALSE)
  }
  unique(x)
}

lisa_dictionary_file <- function(lisa_dictionary, universes) {
  suffix <- switch(
    lisa_dictionary,
    core = "core",
    expanded = "expanded",
    broad_sensitivity = "broad_sensitivity"
  )
  if (length(universes) == 1) {
    prefix <- switch(
      universes,
      `GOBP-C2` = "lisa_gobpc2_dictionary",
      GOMF = "lisa_gomf_dictionary",
      GOCC = "lisa_gocc_dictionary",
      PATHWAYS = "lisa_pathways_dictionary",
      NULL
    )
    if (!is.null(prefix)) return(sprintf("%s_%s_v0_1.tsv", prefix, suffix))
  }
  sprintf("lisa_dictionary_%s_v0_1.tsv", suffix)
}

require_file <- function(path, label) {
  if (!file.exists(path)) stop(sprintf("%s not found: %s", label, path), call. = FALSE)
}

write_lisa_pending_dictionary_result <- function(output_dir, dirs, comparison_name, file_label_prefix,
                                                 universe, lisa_dictionary, lisa_dictionary_path,
                                                 term2gene_path, project_root) {
  comparison_name <- safe_file_label(comparison_name)
  file_label_prefix <- safe_file_label(file_label_prefix %||% "LISA")
  pending <- data.frame(
    status = "pending_dictionary",
    universe = universe,
    lisa_dictionary = lisa_dictionary,
    lisa_dictionary_path = lisa_dictionary_path,
    term2gene_path = term2gene_path,
    project_root = project_root,
    stringsAsFactors = FALSE
  )
  write_tsv_local(pending, file.path(dirs$qc, paste0(comparison_name, "_pending_dictionary.tsv")))
  lisa_guarded_write(
    file.path(output_dir, paste0(comparison_name, "_", file_label_prefix, "_RUN_REPORT.txt")),
    function(path) writeLines(c(
      sprintf("LISA run status: pending_dictionary"),
      sprintf("collection: %s", universe),
      sprintf("dictionary tier: %s", lisa_dictionary),
      sprintf("expected dictionary path: %s", lisa_dictionary_path),
      sprintf("project_root: %s", project_root)
    ), con = path),
    run_root = output_dir
  )
  invisible(list(
    status = "pending_dictionary",
    comparison_name = comparison_name,
    output_dir = output_dir,
    pending = pending
  ))
}

build_hallmarks_lisa_inputs <- function(term2gene_raw) {
  required <- c("gs_collection", "gs_name", "gene_symbol")
  lisa_require_columns(term2gene_raw, required, "TERM2GENE")

  is_h <- term2gene_raw$gs_collection %in% c("H", "MH")
  if ("collection_requested" %in% colnames(term2gene_raw)) {
    is_h <- is_h | term2gene_raw$collection_requested %in% c("H", "MH")
  }
  term2gene_h <- term2gene_raw[is_h, , drop = FALSE]
  term2gene_h <- term2gene_h[!is.na(term2gene_h$gs_name) & term2gene_h$gs_name != "" &
                               !is.na(term2gene_h$gene_symbol) & term2gene_h$gene_symbol != "", , drop = FALSE]
  if (nrow(term2gene_h) == 0) {
    stop("No MSigDB H/HALLMARK rows found in TERM2GENE.", call. = FALSE)
  }

  gene_sets <- unique(term2gene_h$gs_name)
  display <- gsub("^HALLMARK_", "", gene_sets)
  display <- gsub("_", " ", display)
  display <- tools::toTitleCase(tolower(display))
  colors <- rep("#7A3E48", length(gene_sets))

  lisa_dict <- data.frame(
    universe = "HALLMARKS",
    gene_set_id = gene_sets,
    gene_set_name = display,
    source_id = gene_sets,
    category_id = gene_sets,
    category_display_name = display,
    # No LISA_score: the direct HALLMARKS collection only ever carried a
    # constant 4 to satisfy the old common format, never to classify anything.
    # Its real tier is `msigdb_direct`.
    tier = "msigdb_direct",
    source_family = "HALLMARKS",
    stringsAsFactors = FALSE
  )
  category_map <- data.frame(
    category_id = gene_sets,
    display_name = display,
    macrogroup_id = "HALLMARKS",
    macrogroup_name = "MSigDB Hallmarks",
    macrogroup_order = 1,
    category_order_within_macrogroup = seq_along(gene_sets),
    color = colors,
    stringsAsFactors = FALSE
  )
  list(lisa_dict = lisa_dict, category_map = category_map)
}

select_lisa_term2gene_universe <- function(term2gene_raw, universe,
                                           c2_sources,
                                           dictionary_gene_sets = character()) {
  required <- c(
    "gs_collection", "gs_subcollection", "gs_name", "gs_exact_source",
    "gene_symbol"
  )
  lisa_require_columns(term2gene_raw, required, "TERM2GENE")
  universe <- normalize_lisa_universes(universe)[[1L]]
  x <- term2gene_raw[, required, drop = FALSE]
  for (column in required) x[[column]] <- as.character(x[[column]])
  x <- x[
    !is.na(x$gs_name) & nzchar(x$gs_name) &
      !is.na(x$gene_symbol) & nzchar(x$gene_symbol),
    , drop = FALSE
  ]
  x$gene_symbol <- toupper(x$gene_symbol)

  collection <- toupper(trimws(x$gs_collection))
  subcollection <- toupper(trimws(x$gs_subcollection))
  selected <- collection == toupper(universe)
  rule <- "explicit_analysis_collection"

  if (!any(selected)) {
    source_family <- infer_source_family(x$gs_name)
    is_go_bp <- collection %in% c("C5", "M5") &
      subcollection %in% c("GO:BP", "BP")
    is_go_mf <- collection %in% c("C5", "M5") &
      subcollection %in% c("GO:MF", "MF")
    is_go_cc <- collection %in% c("C5", "M5") &
      subcollection %in% c("GO:CC", "CC")
    is_selected_cp <- collection %in% c("C2", "M2") &
      grepl("^CP($|:)", subcollection) & source_family %in% c2_sources
    selected <- switch(
      universe,
      `GOBP-C2` = is_go_bp | is_selected_cp,
      GOMF = is_go_mf,
      GOCC = is_go_cc,
      PATHWAYS = is_go_bp | is_selected_cp,
      HALLMARKS = collection %in% c("H", "MH"),
      rep(FALSE, nrow(x))
    )
    rule <- "msigdb_collection_metadata"
  }

  if (!any(selected) && length(dictionary_gene_sets)) {
    selected <- x$gs_name %in% as.character(dictionary_gene_sets)
    rule <- "dictionary_fallback_for_custom_metadata"
  }
  out <- unique(x[selected, required, drop = FALSE])
  out <- out[order(out$gs_name, out$gene_symbol), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "lisa_universe_selection_rule") <- rule
  out
}

build_lisa_gsea_ledger <- function(term2gene, lisa_dict, ranks,
                                   min_gs_size, max_gs_size, universe) {
  required <- c(
    "gs_collection", "gs_subcollection", "gs_name", "gs_exact_source",
    "gene_symbol"
  )
  lisa_require_columns(term2gene, required, "TERM2GENE universe")
  metadata <- unique(term2gene[c(
    "gs_collection", "gs_subcollection", "gs_name", "gs_exact_source"
  )])
  if (anyDuplicated(metadata$gs_name)) {
    stop(
      "LISA-GSEA-UNIVERSE-001 TERM2GENE has conflicting metadata for one gs_name.",
      call. = FALSE
    )
  }
  membership <- split(
    toupper(as.character(term2gene$gene_symbol)),
    as.character(term2gene$gs_name)
  )
  membership <- lapply(membership, function(value) {
    unique(value[!is.na(value) & nzchar(value)])
  })
  ids <- sort(names(membership))
  metadata <- metadata[match(ids, metadata$gs_name), , drop = FALSE]
  ranked <- names(ranks)
  total_size <- vapply(membership[ids], length, integer(1L))
  ranked_size <- vapply(
    membership[ids], function(value) length(intersect(value, ranked)),
    integer(1L)
  )
  eligible <- ranked_size >= min_gs_size & ranked_size <= max_gs_size
  eligibility_status <- ifelse(
    ranked_size == 0L,
    "no_rank_overlap",
    ifelse(
      ranked_size < min_gs_size,
      "below_min_size",
      ifelse(ranked_size > max_gs_size, "above_max_size", "eligible")
    )
  )

  assignments <- unique(lisa_dict[c("gene_set_id", "category_id")])
  assignment_counts <- table(as.character(assignments$gene_set_id))
  n_categories <- as.integer(assignment_counts[ids])
  n_categories[is.na(n_categories)] <- 0L
  classified <- n_categories > 0L

  data.frame(
    pathway = ids,
    analysis_collection = rep(as.character(universe), length(ids)),
    gs_collection = metadata$gs_collection,
    gs_subcollection = metadata$gs_subcollection,
    gs_exact_source = metadata$gs_exact_source,
    term2gene_source_family = infer_source_family(ids),
    source_gene_count = total_size,
    ranked_gene_count = ranked_size,
    min_gs_size = rep(as.integer(min_gs_size), length(ids)),
    max_gs_size = rep(as.integer(max_gs_size), length(ids)),
    included_in_gsea_universe = rep(TRUE, length(ids)),
    eligible_for_gsea = eligible,
    gsea_eligibility_status = eligibility_status,
    n_lisa_categories = n_categories,
    classification_status = ifelse(classified, "classified", "unclassified"),
    classification_bucket = ifelse(classified, "", "OTHER_UNCLASSIFIED"),
    stringsAsFactors = FALSE
  )
}

complete_lisa_gsea_results <- function(gsea_result, ledger) {
  if (!nrow(ledger)) return(ledger)
  expected <- as.character(ledger$pathway[ledger$eligible_for_gsea %in% TRUE])
  if (!nrow(gsea_result)) {
    actual <- character()
    gsea_result <- data.frame(pathway = character(), stringsAsFactors = FALSE)
  } else {
    if (!"pathway" %in% names(gsea_result)) {
      stop("LISA-GSEA-UNIVERSE-002 fgsea results lack pathway.", call. = FALSE)
    }
    actual <- as.character(gsea_result$pathway)
  }
  if (anyDuplicated(actual)) {
    stop("LISA-GSEA-UNIVERSE-003 fgsea returned duplicate pathway rows.", call. = FALSE)
  }
  missing <- setdiff(expected, actual)
  unexpected <- setdiff(actual, expected)
  if (length(missing) || length(unexpected)) {
    stop(
      paste0(
        "LISA-GSEA-UNIVERSE-004 fgsea result coverage differs from the ",
        "eligible universe; missing=", length(missing),
        ", unexpected=", length(unexpected), "."
      ),
      call. = FALSE
    )
  }
  out <- ledger
  result_index <- match(out$pathway, actual)
  for (column in setdiff(names(gsea_result), "pathway")) {
    out[[column]] <- gsea_result[[column]][result_index]
  }
  for (column in c("pval", "padj", "ES", "NES", "size")) {
    if (!column %in% names(out)) out[[column]] <- rep(NA_real_, nrow(out))
  }
  if (!"leadingEdge" %in% names(out)) {
    out$leadingEdge <- rep(NA_character_, nrow(out))
  }
  out$gsea_result_status <- ifelse(
    out$eligible_for_gsea,
    "tested",
    out$gsea_eligibility_status
  )
  out$analysis <- "GSEA"
  out
}

lisa_unclassified_gsea_rows <- function(df) {
  if (!nrow(df) || !"classification_status" %in% names(df)) {
    return(df[0, , drop = FALSE])
  }
  out <- df[
    as.character(df$classification_status) == "unclassified",
    , drop = FALSE
  ]
  out[order(out$pathway), , drop = FALSE]
}

lisa_classified_gsea_rows <- function(df) {
  if (!nrow(df) || !"category_id" %in% names(df)) {
    return(df[0, , drop = FALSE])
  }
  category <- as.character(df$category_id)
  keep <- !is.na(category) & nzchar(category) &
    category != "OTHER_UNCLASSIFIED"
  if ("classification_status" %in% names(df)) {
    keep <- keep & as.character(df$classification_status) == "classified"
  }
  df[keep, , drop = FALSE]
}

lisa_assert_classified_plot_rows <- function(df, context = "plot") {
  if (!nrow(df)) return(invisible(TRUE))
  category <- if ("category_id" %in% names(df)) {
    as.character(df$category_id)
  } else {
    rep(NA_character_, nrow(df))
  }
  invalid <- is.na(category) | !nzchar(category) |
    category == "OTHER_UNCLASSIFIED"
  if ("classification_status" %in% names(df)) {
    invalid <- invalid |
      as.character(df$classification_status) != "classified"
  }
  if (any(invalid)) {
    stop(
      "LISA-PLOT-UNCLASSIFIED-001 ", context,
      " received unclassified gene-set rows.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

log_msg <- function(verbose, ...) {
  if (isTRUE(verbose)) message(sprintf(...))
}

read_tsv_local <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE)
}

write_tsv_local <- function(x, path) {
  write_lisa_tsv(x, path)
}

resolve_lisa_term2gene <- function(project_root, term2gene_path, msigdb_mode, species, verbose = TRUE) {
  if (!is.null(term2gene_path)) {
    return(list(
      path = term2gene_path,
      msigdb_mode = paste0(msigdb_mode, "_custom_term2gene"),
      db_species = "custom",
      target_species = species,
      generated = FALSE
    ))
  }

  universe_dir <- file.path(project_root, "inputs", "msigdb_universe")
  cfg <- switch(
    msigdb_mode,
    human = list(
      path = file.path(universe_dir, "msigdb_hs_term2gene.tsv"),
      db_species = "HS",
      target_species = "Homo sapiens",
      collections = c("C2", "C5")
    ),
    mouse_hs_ortholog = list(
      path = file.path(universe_dir, "msigdb_mm_hs_ortholog_term2gene.tsv"),
      db_species = "HS",
      target_species = "Mus musculus",
      collections = c("C2", "C5")
    ),
    mouse_native = list(
      path = file.path(universe_dir, "msigdb_mm_native_term2gene.tsv"),
      db_species = "MM",
      target_species = "Mus musculus",
      collections = c("M2", "M5")
    )
  )

  generated <- FALSE
  if (!file.exists(cfg$path)) {
    log_msg(verbose, "TERM2GENE cache not found for msigdb_mode=%s; exporting with msigdbr: %s", msigdb_mode, cfg$path)
    export_lisa_msigdb_term2gene(
      out_path = cfg$path,
      db_species = cfg$db_species,
      target_species = cfg$target_species,
      collections = cfg$collections,
      msigdb_mode = msigdb_mode,
      verbose = verbose
    )
    generated <- TRUE
  }

  list(
    path = cfg$path,
    msigdb_mode = msigdb_mode,
    db_species = cfg$db_species,
    target_species = cfg$target_species,
    generated = generated
  )
}

export_lisa_msigdb_term2gene <- function(out_path, db_species, target_species, collections, msigdb_mode, verbose = TRUE) {
  lisa_require_optional("msigdbr", "exporting TERM2GENE caches; alternatively provide term2gene_path explicitly")

  fetch_one <- function(collection) {
    log_msg(verbose, "Fetching msigdbr db_species=%s species=%s collection=%s", db_species, target_species, collection)
    tryCatch(
      msigdbr::msigdbr(db_species = db_species, species = target_species, collection = collection),
      error = function(e) msigdbr::msigdbr(db_species = db_species, species = target_species, category = collection)
    )
  }

  raw <- do.call(rbind, lapply(collections, fetch_one))
  required <- c("gs_collection", "gs_subcollection", "gs_name", "gs_exact_source", "gene_symbol")
  for (cc in required) {
    if (!cc %in% colnames(raw)) raw[[cc]] <- NA_character_
  }
  term2gene <- unique(raw[, required, drop = FALSE])
  term2gene <- term2gene[!is.na(term2gene$gs_name) & term2gene$gs_name != "" &
                           !is.na(term2gene$gene_symbol) & term2gene$gene_symbol != "", , drop = FALSE]
  term2gene <- term2gene[order(term2gene$gs_name, term2gene$gene_symbol), , drop = FALSE]
  write_tsv_local(term2gene, out_path)

  metadata <- data.frame(
    field = c("msigdb_mode", "db_species", "target_species", "collections", "n_gene_sets", "n_term2gene_rows", "export_time", "msigdbr_version"),
    value = c(
      msigdb_mode,
      db_species,
      target_species,
      paste(collections, collapse = ","),
      as.character(length(unique(term2gene$gs_name))),
      as.character(nrow(term2gene)),
      as.character(Sys.time()),
      as.character(utils::packageVersion("msigdbr"))
    ),
    stringsAsFactors = FALSE
  )
  write_tsv_local(metadata, sub("\\.tsv$", "_metadata.tsv", out_path))
  invisible(out_path)
}

build_lisa_term2gene_coverage <- function(lisa_dict, term2gene_raw) {
  if (!"gs_name" %in% colnames(term2gene_raw)) {
    stop("TERM2GENE input must contain a gs_name column.", call. = FALSE)
  }
  dict_sets <- unique(as.character(lisa_dict$gene_set_id))
  available_sets <- unique(as.character(term2gene_raw$gs_name))
  source_family <- infer_source_family(dict_sets)
  detail <- data.frame(
    gene_set_id = dict_sets,
    source_family = source_family,
    available_in_term2gene = dict_sets %in% available_sets,
    stringsAsFactors = FALSE
  )
  if (nrow(detail) == 0L) {
    stop("The selected LISA dictionary contains zero rows for the requested collection.", call. = FALSE)
  }
  summary <- do.call(rbind, lapply(split(detail, detail$source_family), function(x) {
    data.frame(
      source_family = x$source_family[1],
      dictionary_gene_sets = nrow(x),
      available_term2gene_gene_sets = sum(x$available_in_term2gene),
      missing_term2gene_gene_sets = sum(!x$available_in_term2gene),
      pct_available = round(100 * sum(x$available_in_term2gene) / nrow(x), 2),
      stringsAsFactors = FALSE
    )
  }))
  summary <- summary[order(summary$source_family), , drop = FALSE]
  missing <- detail[!detail$available_in_term2gene, , drop = FALSE]
  list(summary = summary, missing = missing)
}

make_output_dirs <- function(output_dir) {
  dirs <- list(
    inputs = file.path(output_dir, "inputs"),
    enrichment = file.path(output_dir, "enrichment"),
    lisa_tables = file.path(output_dir, "lisa_tables"),
    qc = file.path(output_dir, "qc"),
    plots = file.path(output_dir, "plots")
  )
  invisible(lapply(dirs, lisa_guarded_dir_create, run_root = output_dir))
  dirs
}

resolve_input_type <- function(input, input_type) {
  if (input_type != "auto") return(input_type)
  if (is.character(input) && length(input) == 1 && grepl("\\.(tsv|txt|csv)$", input, ignore.case = TRUE)) return("de_table")
  if (is.character(input) && length(input) == 1 && grepl("\\.rds$", input, ignore.case = TRUE)) {
    lisa_rds_abort(
      "LISA-RDS-010",
      "an RDS path reached type resolution without passing the explicit trusted-RDS boundary."
    )
  }
  if (inherits(input, "DESeqDataSet")) return("deseq2_dds")
  if (inherits(input, "DESeqResults")) return("deseq2_results")
  if (inherits(input, "DGEList")) {
    stop(
      "LISA-INPUT-EDGER-001: `DGEList` contains counts, not a differential-expression result. Supply DGELRT, DGEExact or TopTags.",
      call. = FALSE
    )
  }
  if (inherits(input, c("DGELRT", "DGEExact", "TopTags"))) return("edger")
  if (is.data.frame(input) || is.matrix(input)) return("de_table")
  "de_table"
}

extract_de_table <- function(input, input_type, deseq2_contrast, deseq2_name) {
  obj <- input
  if (is.character(input) && length(input) == 1) {
    if (grepl("\\.rds$", input, ignore.case = TRUE)) {
      lisa_rds_abort(
        "LISA-RDS-010",
        "an RDS path reached extraction without passing the explicit trusted-RDS boundary."
      )
    }
    if (grepl("\\.csv$", input, ignore.case = TRUE)) return(utils::read.csv(input, stringsAsFactors = FALSE, check.names = FALSE))
    if (grepl("\\.(tsv|txt)$", input, ignore.case = TRUE)) return(read_tsv_local(input))
  }
  if (input_type == "deseq2_dds") {
    lisa_require_optional("DESeq2", "extracting results from a DESeq2 object")
    has_contrast <- !is.null(deseq2_contrast)
    has_name <- !is.null(deseq2_name)
    if (!has_contrast && !has_name) {
      stop(
        paste0(
          "LISA-INPUT-DESEQ2-001: a `DESeqDataSet` requires exactly one result selector. ",
          "Provide either `deseq2_name` or `deseq2_contrast`."
        ),
        call. = FALSE
      )
    }
    if (has_contrast && has_name) {
      stop(
        paste0(
          "LISA-INPUT-DESEQ2-002: `deseq2_name` and `deseq2_contrast` are mutually exclusive. ",
          "Provide exactly one result selector."
        ),
        call. = FALSE
      )
    }
    if (has_name) {
      if (!is.character(deseq2_name) || length(deseq2_name) != 1L ||
          is.na(deseq2_name) || !nzchar(trimws(deseq2_name))) {
        stop(
          "LISA-INPUT-DESEQ2-003: `deseq2_name` must be one non-empty character value.",
          call. = FALSE
        )
      }
      available_names <- DESeq2::resultsNames(obj)
      if (!deseq2_name %in% available_names) {
        stop(
          sprintf(
            "LISA-INPUT-DESEQ2-004: unknown `deseq2_name` '%s'. Available names: %s.",
            deseq2_name,
            paste(available_names, collapse = ", ")
          ),
          call. = FALSE
        )
      }
      return(as.data.frame(DESeq2::results(obj, name = deseq2_name)))
    }
    valid_character_contrast <- is.character(deseq2_contrast) &&
      length(deseq2_contrast) == 3L &&
      !anyNA(deseq2_contrast) &&
      all(nzchar(trimws(deseq2_contrast)))
    valid_list_contrast <- is.list(deseq2_contrast) &&
      length(deseq2_contrast) %in% c(1L, 2L) &&
      all(vapply(
        deseq2_contrast,
        function(value) {
          is.character(value) && length(value) > 0L &&
            !anyNA(value) && all(nzchar(trimws(value)))
        },
        logical(1)
      ))
    valid_numeric_contrast <- is.numeric(deseq2_contrast) &&
      length(deseq2_contrast) > 0L &&
      !anyNA(deseq2_contrast) &&
      all(is.finite(deseq2_contrast))
    if (!valid_character_contrast && !valid_list_contrast && !valid_numeric_contrast) {
      stop(
        paste0(
          "LISA-INPUT-DESEQ2-005: `deseq2_contrast` must be a non-empty DESeq2 contrast: ",
          "a three-value character vector, a one- or two-part character list, or a finite numeric vector."
        ),
        call. = FALSE
      )
    }
    return(as.data.frame(DESeq2::results(obj, contrast = deseq2_contrast)))
  }
  if (input_type == "deseq2_results") return(as.data.frame(obj))
  if (input_type == "edger") {
    if (inherits(obj, "DGEList")) {
      stop(
        "LISA-INPUT-EDGER-001: `DGEList` contains counts, not a differential-expression result. Supply DGELRT, DGEExact or TopTags.",
        call. = FALSE
      )
    }
    return(as.data.frame(lisa_edger_result_table(obj)))
  }
  as.data.frame(obj)
}

standardize_de_table <- function(de, symbol_col, rank_col, logfc_col, padj_col, pvalue_col, gene_id_col) {
  de_raw <- as.data.frame(de)
  de <- de_raw
  if (!is.null(rownames(de)) && !any(grepl("^rowname$", tolower(colnames(de))))) de$rowname <- rownames(de)

  pick <- function(user_col, candidates, required = TRUE) {
    if (!is.null(user_col)) {
      if (!user_col %in% colnames(de)) stop(sprintf("Column not found: %s", user_col), call. = FALSE)
      return(user_col)
    }
    idx <- match(tolower(candidates), tolower(colnames(de)))
    idx <- idx[!is.na(idx)]
    if (length(idx) > 0) return(colnames(de)[idx[1]])
    if (required) stop(sprintf("Required DE column not found. Candidates: %s", paste(candidates, collapse = ", ")), call. = FALSE)
    NA_character_
  }

  symbol_source <- pick(symbol_col, c("human_symbol", "symbol", "gene_symbol", "gene_name", "SYMBOL", "external_gene_name"))
  logfc_source <- pick(logfc_col, c("log2FoldChange", "logFC", "avg_log2FC", "log2FC"), required = FALSE)
  padj_source <- pick(padj_col, c("padj", "FDR", "adj.P.Val", "qvalue", "p.adjust"), required = FALSE)
  pval_source <- pick(pvalue_col, c("pvalue", "PValue", "P.Value", "p_val"), required = FALSE)
  gene_id_source <- pick(gene_id_col, c("gene_id", "ensembl", "ensembl_id", "id", "rowname"), required = FALSE)

  if (!is.null(rank_col)) {
    rank_source <- pick(rank_col, rank_col)
    rank_mode <- "explicit"
  } else {
    rank_source <- pick(NULL, c("wald_stat", "stat", "Wald", "score", "t", "z"), required = FALSE)
    rank_mode <- "stat"
    if (is.na(rank_source) && !is.na(logfc_source)) {
      rank_source <- logfc_source
      rank_mode <- "logfc"
    }
    if (is.na(rank_source) && !is.na(logfc_source) && !is.na(pval_source)) rank_mode <- "signed_logp"
    if (is.na(rank_source) && rank_mode != "signed_logp") {
      stop("Could not infer ranking column. Provide rank_col, or logfc_col plus pvalue_col.", call. = FALSE)
    }
  }

  out <- data.frame(
    gene_id = if (!is.na(gene_id_source)) as.character(de[[gene_id_source]]) else rownames(de),
    symbol = as.character(de[[symbol_source]]),
    log2FoldChange = if (!is.na(logfc_source)) suppressWarnings(as.numeric(de[[logfc_source]])) else NA_real_,
    padj = if (!is.na(padj_source)) suppressWarnings(as.numeric(de[[padj_source]])) else NA_real_,
    pvalue = if (!is.na(pval_source)) suppressWarnings(as.numeric(de[[pval_source]])) else NA_real_,
    stringsAsFactors = FALSE
  )
  if (rank_mode == "signed_logp") {
    out$rank_value <- sign(out$log2FoldChange) * -log10(pmax(out$pvalue, 1e-300))
  } else {
    out$rank_value <- suppressWarnings(as.numeric(de[[rank_source]]))
  }
  out <- out[!is.na(out$symbol) & out$symbol != "" & !is.na(out$rank_value), , drop = FALSE]
  out$symbol <- toupper(out$symbol)
  duplicate_symbols <- unique(out$symbol[duplicated(out$symbol)])
  if (length(duplicate_symbols)) {
    stop("LISA-DUP-018: residual canonical duplicate symbols reached standardize_de_table: ",
      paste(duplicate_symbols, collapse = ", "),
      ". Resolve duplicates before analysis; even identical rank values are not silently discarded.",
      call. = FALSE)
  }

  mapping <- data.frame(
    standard_column = c("gene_id", "symbol", "rank_value", "log2FoldChange", "padj", "pvalue"),
    source_column = c(gene_id_source, symbol_source, rank_source, logfc_source, padj_source, pval_source),
    stringsAsFactors = FALSE
  )
  list(de = out, de_raw = de_raw, column_mapping = mapping)
}

make_rank_vector <- function(de) {
  de$symbol_match <- toupper(as.character(de$symbol))
  duplicate_symbols <- unique(de$symbol_match[duplicated(de$symbol_match)])
  if (length(duplicate_symbols)) {
    stop("LISA-DUP-018: residual canonical duplicate symbols reached make_rank_vector: ",
      paste(duplicate_symbols, collapse = ", "),
      ". Resolve duplicates before ranking; rank values are never averaged.", call. = FALSE)
  }
  ranks <- suppressWarnings(as.numeric(de$rank_value))
  names(ranks) <- de$symbol_match
  sort(ranks[is.finite(ranks)], decreasing = TRUE)
}

run_fgsea_lisa <- function(pathways, ranks, min_gs_size, max_gs_size, n_threads, fgsea_nperm, fgsea_eps) {
  pathways <- lapply(pathways, function(x) intersect(toupper(x), names(ranks)))
  pathways <- pathways[lengths(pathways) >= min_gs_size & lengths(pathways) <= max_gs_size]
  if (length(pathways) == 0) return(data.frame())
  notes <- character()
  res <- withCallingHandlers({
    if (!is.null(fgsea_nperm)) {
      if (exists("fgseaSimple", envir = asNamespace("fgsea"), inherits = FALSE)) {
        fgsea::fgseaSimple(pathways = pathways, stats = ranks, nperm = fgsea_nperm,
                           minSize = min_gs_size, maxSize = max_gs_size, nproc = n_threads)
      } else {
        fgsea::fgsea(pathways = pathways, stats = ranks, nperm = fgsea_nperm,
                     minSize = min_gs_size, maxSize = max_gs_size, nproc = n_threads)
      }
    } else {
      fgsea::fgseaMultilevel(pathways = pathways, stats = ranks,
                             minSize = min_gs_size, maxSize = max_gs_size,
                             nproc = n_threads, eps = fgsea_eps)
    }
  }, warning = function(w) {
    notes <<- c(notes, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  out <- as.data.frame(res)
  attr(out, "lisa_notes") <- unique(notes)
  out
}

run_lisa_ora <- function(de, term2gene, lisa_dict, ora_padj_cutoff, ora_lfc_cutoff,
                         analysis_id = "unspecified", collection = "unspecified") {
  if (!all(c("log2FoldChange", "padj") %in% colnames(de))) return(data.frame())
  bg <- unique(de$symbol[!is.na(de$symbol) & de$symbol != ""])
  up <- unique(de$symbol[!is.na(de$padj) & de$padj <= ora_padj_cutoff & !is.na(de$log2FoldChange) & de$log2FoldChange >= ora_lfc_cutoff])
  down <- unique(de$symbol[!is.na(de$padj) & de$padj <= ora_padj_cutoff & !is.na(de$log2FoldChange) & de$log2FoldChange <= -ora_lfc_cutoff])
  rbind(
    run_ora_one(up, bg, term2gene, lisa_dict, "UP", ora_padj_cutoff,
                analysis_id = analysis_id, collection = collection),
    run_ora_one(down, bg, term2gene, lisa_dict, "DOWN", ora_padj_cutoff,
                analysis_id = analysis_id, collection = collection)
  )
}

empty_ora_family <- function() {
  data.frame(
    pathway = character(),
    analysis_id = character(),
    collection = character(),
    family_id = character(),
    direction = character(),
    overlap = integer(),
    set_size = integer(),
    query_size = integer(),
    pvalue = numeric(),
    padj = numeric(),
    genes = character(),
    n_eligible = integer(),
    n_tested = integer(),
    method = character(),
    selected = logical(),
    analysis = character(),
    stringsAsFactors = FALSE
  )
}

run_ora_one <- function(genes, bg, term2gene, lisa_dict, direction, padj_cutoff,
                        analysis_id = "unspecified", collection = "unspecified") {
  direction <- toupper(as.character(direction)[1])
  analysis_id <- as.character(analysis_id)[1]
  collection <- as.character(collection)[1]
  bg <- unique(toupper(as.character(bg)))
  bg <- bg[!is.na(bg) & nzchar(bg)]
  genes <- intersect(unique(toupper(as.character(genes))), bg)
  empty <- empty_ora_family()
  if (length(genes) < 5 || length(bg) < 10) {
    return(annotate_lisa(empty, lisa_dict, by_col = "pathway"))
  }
  t2g <- term2gene
  t2g <- t2g[
    !is.na(t2g$gs_name) & nzchar(as.character(t2g$gs_name)) &
      !is.na(t2g$gene_symbol) & nzchar(as.character(t2g$gene_symbol)),
    , drop = FALSE
  ]
  t2g$gs_name <- as.character(t2g$gs_name)
  t2g$gene_symbol <- toupper(as.character(t2g$gene_symbol))
  term_list <- split(t2g$gene_symbol, t2g$gs_name)
  term_list <- lapply(term_list, function(x) intersect(unique(x), bg))
  term_list <- term_list[lengths(term_list) > 0]
  if (length(term_list) == 0) {
    return(annotate_lisa(empty, lisa_dict, by_col = "pathway"))
  }
  # This unique, background-eligible term list is the hypothesis family.
  # Retain zero-overlap sets through BH; LISA annotation may fan one tested
  # gene set out to multiple category rows only after adjustment.
  hits <- lapply(term_list, function(x) intersect(x, genes))
  k <- lengths(hits)
  m <- lengths(term_list)
  n <- length(bg) - m
  q <- length(genes)
  pval <- stats::phyper(k - 1, m, n, q, lower.tail = FALSE)
  pval[k == 0L] <- 1
  padj <- stats::p.adjust(pval, method = "BH")
  n_eligible <- length(term_list)
  n_tested <- length(pval)
  out <- data.frame(
    pathway = names(term_list),
    analysis_id = analysis_id,
    collection = collection,
    family_id = paste("ORA", analysis_id, collection, direction, sep = ":"),
    direction = direction,
    overlap = as.integer(k),
    set_size = as.integer(m),
    query_size = q,
    pvalue = pval,
    padj = padj,
    genes = vapply(hits, paste, character(1), collapse = "/"),
    n_eligible = as.integer(n_eligible),
    n_tested = as.integer(n_tested),
    method = "BH",
    selected = !is.na(padj) & padj <= padj_cutoff,
    stringsAsFactors = FALSE
  )
  out$analysis <- paste0("ORA_", direction)
  annotate_lisa(out, lisa_dict, by_col = "pathway")
}

annotate_lisa <- function(df, lisa_dict, by_col) {
  # The merge is by gene set, never by score. C6 removed `LISA_score` from the
  # selection; every remaining column, the join key and the row multiplicity are
  # unchanged.
  dict_keep <- lisa_dict[, c("universe", "gene_set_id", "gene_set_name", "source_id", "category_id",
                             "category_display_name", "tier", "source_family"), drop = FALSE]
  merge(df, dict_keep, by.x = by_col, by.y = "gene_set_id", all.x = TRUE)
}

add_lisa_category_metadata <- function(df, category_map) {
  if (nrow(df) == 0 || !"category_id" %in% colnames(df)) return(df)
  map <- category_map[, c("category_id", "macrogroup_id", "macrogroup_name",
                          "macrogroup_order", "category_order_within_macrogroup", "color"), drop = FALSE]
  metadata_cols <- setdiff(colnames(map), "category_id")
  df <- df[, setdiff(colnames(df), metadata_cols), drop = FALSE]
  merge(df, map, by = "category_id", all.x = TRUE)
}

build_lisa_summaries <- function(gsea_all, ora_all, category_map,
                                 include_empty_categories,
                                 gsea_padj_cutoff, ora_padj_cutoff,
                                 gsea_ledger = NULL) {
  gsea_padj_cutoff <- lisa_validate_gsea_padj_cutoff(
    gsea_padj_cutoff,
    field = "gsea_padj_cutoff"
  )
  base <- category_map[, c("category_id", "display_name", "macrogroup_id", "macrogroup_name",
                           "macrogroup_order", "category_order_within_macrogroup", "color"), drop = FALSE]
  names(base)[names(base) == "display_name"] <- "category_display_name"

  # Category means are intentionally calculated independently for each
  # analysis.  `mean_NES` keeps its historical meaning (significant mapped
  # member sets), while the contextual mean preserves evaluable pre-threshold
  # information for an A/B endpoint without turning absent support into zero.
  # Evaluable means that both NES and adjusted p-value are finite; an incomplete
  # statistical result stays mapped, but cannot provide a contextual endpoint.
  gsea_mapped <- if (nrow(gsea_all) > 0) {
    gsea_all[!is.na(gsea_all$category_id) & nzchar(as.character(gsea_all$category_id)), , drop = FALSE]
  } else {
    data.frame()
  }
  if (nrow(gsea_mapped) > 0) {
    gsea <- do.call(rbind, lapply(split(gsea_mapped, gsea_mapped$category_id), function(x) {
      if ("pathway" %in% names(x)) {
        pathway <- as.character(x$pathway)
        keep <- !is.na(pathway) & nzchar(pathway) & !duplicated(pathway)
        x <- x[keep, , drop = FALSE]
      }
      nes <- suppressWarnings(as.numeric(x$NES))
      padj <- suppressWarnings(as.numeric(x$padj))
      evaluable <- is.finite(nes) & is.finite(padj)
      significant <- evaluable & padj <= gsea_padj_cutoff
      sig_nes <- nes[significant]
      contextual_nes <- nes[evaluable]
      pos_nes <- sig_nes[sig_nes > 0]
      neg_nes <- sig_nes[sig_nes < 0]
      n_mapped <- nrow(x)
      n_evaluable <- sum(evaluable)
      n_significant <- sum(significant)
      significant_mean <- if (n_significant > 0L) mean(sig_nes) else NA_real_
      data.frame(
        category_id = x$category_id[1],
        n_genesets = n_significant,
        n_genesets_mapped = n_mapped,
        n_genesets_evaluable = n_evaluable,
        n_genesets_significant = n_significant,
        n_pos_genesets = sum(sig_nes > 0),
        n_neg_genesets = sum(sig_nes < 0),
        mean_NES = significant_mean,
        mean_NES_contextual = if (n_evaluable > 0L) mean(contextual_nes) else NA_real_,
        mean_pos_NES = if (length(pos_nes) > 0) mean(pos_nes) else NA_real_,
        mean_neg_NES = if (length(neg_nes) > 0) mean(neg_nes) else NA_real_,
        median_NES = if (n_significant > 0L) stats::median(sig_nes) else NA_real_,
        consistency = if (n_significant > 0L) mean(sign(sig_nes) == sign(significant_mean)) else NA_real_,
        min_padj = if (n_significant > 0L) min(padj[significant]) else NA_real_,
        min_padj_mapped = if (n_evaluable > 0L) min(padj[evaluable]) else NA_real_,
        has_significant_support = n_significant > 0L,
        gsea_padj_cutoff = gsea_padj_cutoff,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    gsea <- data.frame(
      category_id = character(), n_genesets = integer(),
      n_genesets_mapped = integer(), n_genesets_evaluable = integer(),
      n_genesets_significant = integer(), n_pos_genesets = integer(),
      n_neg_genesets = integer(), mean_NES = numeric(),
      mean_NES_contextual = numeric(), mean_pos_NES = numeric(),
      mean_neg_NES = numeric(), median_NES = numeric(), consistency = numeric(),
      min_padj = numeric(), min_padj_mapped = numeric(),
      has_significant_support = logical(), gsea_padj_cutoff = numeric()
    )
  }
  gsea <- merge(base, gsea, by = "category_id", all.x = include_empty_categories, all.y = TRUE)
  gsea <- fill_summary_zeros(gsea, c(
    "n_genesets", "n_genesets_mapped", "n_genesets_evaluable",
    "n_genesets_significant", "n_pos_genesets", "n_neg_genesets"
  ))
  for (cc in c(
    "n_genesets", "n_genesets_mapped", "n_genesets_evaluable",
    "n_genesets_significant", "n_pos_genesets", "n_neg_genesets"
  )) {
    gsea[[cc]] <- as.integer(gsea[[cc]])
  }
  gsea$has_significant_support <- gsea$n_genesets_significant > 0L
  gsea$mean_NES[!gsea$has_significant_support] <- NA_real_
  gsea$gsea_padj_cutoff <- gsea_padj_cutoff
  gsea$same_direction_pct <- ifelse(gsea$n_genesets > 0, round(100 * gsea$consistency, 1), NA_real_)
  gsea$mean_NES_direction <- ifelse(
    gsea$n_genesets > 0,
    ifelse(gsea$mean_NES >= 0, "Positive NES", "Negative NES"),
    NA_character_
  )
  gsea$analysis <- "GSEA"
  gsea <- order_category_summary(gsea)

  ora_sig <- if (nrow(ora_all) > 0) ora_all[!is.na(ora_all$category_id) & !is.na(ora_all$padj) & ora_all$padj <= ora_padj_cutoff, , drop = FALSE] else data.frame()
  if (nrow(ora_sig) > 0) {
    ora <- do.call(rbind, lapply(split(ora_sig, paste(ora_sig$category_id, ora_sig$direction, sep = "|")), function(x) {
      data.frame(
        category_id = x$category_id[1],
        direction = x$direction[1],
        n_genesets = length(unique(x$pathway)),
        total_overlap = sum(x$overlap, na.rm = TRUE),
        min_padj = min(x$padj, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    ora <- data.frame(category_id = character(), direction = character(), n_genesets = integer(),
                      total_overlap = integer(), min_padj = numeric())
  }
  if (include_empty_categories) {
    expanded_base <- rbind(transform(base, direction = "UP"), transform(base, direction = "DOWN"))
    ora <- merge(expanded_base, ora, by = c("category_id", "direction"), all.x = TRUE, all.y = TRUE)
  } else {
    ora <- merge(base, ora, by = "category_id", all.y = TRUE)
  }
  ora <- fill_summary_zeros(ora, c("n_genesets", "total_overlap"))
  ora$analysis <- paste0("ORA_", ora$direction)
  ora <- order_category_summary(ora)

  integrated <- rbind_fill_local(list(
    data.frame(gsea[, intersect(colnames(gsea), c(
      "analysis", "category_id", "category_display_name", "macrogroup_id",
      "macrogroup_name", "n_genesets", "n_genesets_mapped",
      "n_genesets_evaluable", "n_genesets_significant", "n_pos_genesets",
      "n_neg_genesets", "mean_NES", "mean_NES_contextual",
      "mean_pos_NES", "mean_neg_NES", "median_NES", "consistency",
      "same_direction_pct", "mean_NES_direction", "min_padj",
      "min_padj_mapped", "has_significant_support", "gsea_padj_cutoff", "color"
    ))], stringsAsFactors = FALSE),
    data.frame(ora[, intersect(colnames(ora), c("analysis", "direction", "category_id", "category_display_name", "macrogroup_id", "macrogroup_name", "n_genesets", "total_overlap", "min_padj", "color"))], stringsAsFactors = FALSE)
  ))

  coverage <- coverage_table(gsea_all, ora_all, gsea_ledger = gsea_ledger)
  unmapped <- if (!is.null(gsea_ledger) && nrow(gsea_ledger)) {
    lisa_unclassified_gsea_rows(gsea_ledger)
  } else {
    rbind(
      unmapped_table(gsea_all, "GSEA"),
      unmapped_table(ora_all, "ORA")
    )
  }
  list(gsea = gsea, ora = ora, integrated = integrated, coverage = coverage, unmapped = unmapped)
}

find_lisa_category_summary <- function(input, universe, analysis = "GSEA") {
  input <- as.character(input)
  if (length(input) != 1 || is.na(input) || input == "") {
    stop("Each LISA contrast input must be a single file or directory.", call. = FALSE)
  }
  if (file.exists(input) && !dir.exists(input)) return(normalizePath(input, mustWork = TRUE))
  if (!dir.exists(input)) stop(sprintf("LISA contrast input not found: %s", input), call. = FALSE)

  # The label between the comparison name and the summary suffix is
  # user-configurable (`file_label_prefix`).  Match the stable suffix instead
  # of assuming the historical "LISA" or "semantic" prefixes.  Collection
  # selection below keeps recursive searches unambiguous.
  pattern <- sprintf("_%s_category_summary\\.tsv$", analysis)
  collection_dir <- file.path(input, paste0("collection_", safe_file_label(universe)), "lisa_tables")
  candidates <- character()
  if (dir.exists(collection_dir)) {
    candidates <- list.files(collection_dir, pattern = pattern, full.names = TRUE)
  }
  if (length(candidates) == 0) {
    direct_dir <- file.path(input, "lisa_tables")
    if (dir.exists(direct_dir)) candidates <- list.files(direct_dir, pattern = pattern, full.names = TRUE)
  }
  if (length(candidates) == 0) {
    all_candidates <- list.files(input, pattern = pattern, recursive = TRUE, full.names = TRUE)
    universe_token <- paste0("collection_", safe_file_label(universe))
    universe_candidates <- all_candidates[grepl(universe_token, all_candidates, fixed = TRUE)]
    candidates <- if (length(universe_candidates) > 0) universe_candidates else all_candidates
  }
  if (length(candidates) == 0) {
    stop(sprintf("No LISA %s category summary found under: %s", analysis, input), call. = FALSE)
  }
  if (length(candidates) > 1) {
    stop(sprintf(
      "Multiple LISA %s category summaries found for %s. Pass the exact summary file instead:\n%s",
      analysis, input, paste(candidates, collapse = "\n")
    ), call. = FALSE)
  }
  normalizePath(candidates[1], mustWork = TRUE)
}

read_lisa_contrast_summary <- function(path, contrast_label, gsea_padj_cutoff = 0.25) {
  require_file(path, sprintf("LISA category summary for %s", contrast_label))
  df <- read_tsv_local(path)
  if (!"category_id" %in% colnames(df)) {
    stop(sprintf("LISA category summary lacks category_id: %s", path), call. = FALSE)
  }
  normalize_lisa_gsea_summary_support(
    df,
    gsea_padj_cutoff = gsea_padj_cutoff,
    source = path
  )
}

normalize_lisa_gsea_summary_support <- function(df, gsea_padj_cutoff = 0.25,
                                                source = "in-memory summary") {
  gsea_padj_cutoff <- lisa_validate_gsea_padj_cutoff(
    gsea_padj_cutoff,
    field = "gsea_padj_cutoff"
  )
  if (!"category_id" %in% colnames(df)) {
    stop(sprintf("LISA category summary lacks category_id: %s", source), call. = FALSE)
  }
  if (!"category_display_name" %in% colnames(df)) df$category_display_name <- df$category_id
  if (!"macrogroup_id" %in% colnames(df)) df$macrogroup_id <- "UNCLASSIFIED"
  if (!"macrogroup_name" %in% colnames(df)) df$macrogroup_name <- "Other or unclassified"
  if (!"macrogroup_order" %in% colnames(df)) df$macrogroup_order <- 99
  if (!"category_order_within_macrogroup" %in% colnames(df)) df$category_order_within_macrogroup <- seq_len(nrow(df))
  if (!"color" %in% colnames(df)) df$color <- "#737373"
  if (!"n_genesets" %in% colnames(df)) df$n_genesets <- 0
  if (!"mean_NES" %in% colnames(df)) df$mean_NES <- NA_real_
  if (!"median_NES" %in% colnames(df)) df$median_NES <- NA_real_
  if (!"consistency" %in% colnames(df)) df$consistency <- NA_real_

  legacy_n <- suppressWarnings(as.numeric(df$n_genesets))
  legacy_n[!is.finite(legacy_n)] <- 0
  if (!"n_genesets_significant" %in% colnames(df)) df$n_genesets_significant <- legacy_n
  if (!"n_genesets_evaluable" %in% colnames(df)) df$n_genesets_evaluable <- legacy_n
  if (!"n_genesets_mapped" %in% colnames(df)) df$n_genesets_mapped <- df$n_genesets_evaluable
  count_cols <- c("n_genesets_mapped", "n_genesets_evaluable", "n_genesets_significant")
  for (cc in count_cols) {
    value <- suppressWarnings(as.numeric(df[[cc]]))
    value[!is.finite(value)] <- 0
    if (any(value < 0 | value != floor(value))) {
      stop(sprintf("LISA category summary has invalid %s in %s.", cc, source), call. = FALSE)
    }
    df[[cc]] <- as.integer(value)
  }
  invalid_counts <- df$n_genesets_significant > df$n_genesets_evaluable |
    df$n_genesets_evaluable > df$n_genesets_mapped
  if (any(invalid_counts)) {
    stop(sprintf(
      "LISA category summary support counts must satisfy significant <= evaluable <= mapped: %s",
      source
    ), call. = FALSE)
  }
  df$n_genesets <- df$n_genesets_significant
  df$mean_NES <- suppressWarnings(as.numeric(df$mean_NES))
  if (!"mean_NES_contextual" %in% colnames(df)) {
    df$mean_NES_contextual <- ifelse(df$n_genesets_evaluable > 0L, df$mean_NES, NA_real_)
  }
  df$mean_NES_contextual <- suppressWarnings(as.numeric(df$mean_NES_contextual))
  if (!"has_significant_support" %in% colnames(df)) {
    df$has_significant_support <- df$n_genesets_significant > 0L
  }
  df$has_significant_support <- df$n_genesets_significant > 0L
  if (any(df$has_significant_support & !is.finite(df$mean_NES))) {
    stop(sprintf(
      "LISA category summary has significant support but no finite mean_NES: %s",
      source
    ), call. = FALSE)
  }
  if (any(df$n_genesets_evaluable > 0L & !is.finite(df$mean_NES_contextual))) {
    stop(sprintf(
      "LISA category summary has evaluable member sets but no finite contextual mean_NES: %s",
      source
    ), call. = FALSE)
  }
  # Historical empty summaries encoded absence as mean_NES=0. The count is the
  # authoritative compatibility discriminator; unsupported aliases are NA.
  df$mean_NES[!df$has_significant_support] <- NA_real_

  cutoff_source <- "contrast_argument_legacy_summary"
  if ("gsea_padj_cutoff" %in% colnames(df)) {
    observed <- suppressWarnings(as.numeric(df$gsea_padj_cutoff))
    observed <- unique(observed[is.finite(observed)])
    if (length(observed) != 1L || !isTRUE(all.equal(observed[[1]], gsea_padj_cutoff,
      tolerance = .Machine$double.eps^0.5))) {
      stop(sprintf(
        "LISA category summary cutoff does not match gsea_padj_cutoff=%s: %s",
        format(gsea_padj_cutoff), source
      ), call. = FALSE)
    }
    cutoff_source <- "single_de_summary"
  }
  df$gsea_padj_cutoff <- gsea_padj_cutoff
  if (!"gsea_padj_cutoff_source" %in% colnames(df)) {
    df$gsea_padj_cutoff_source <- cutoff_source
  }
  if (!"same_direction_pct" %in% colnames(df)) {
    df$same_direction_pct <- ifelse(suppressWarnings(as.numeric(df$n_genesets)) > 0, 100 * suppressWarnings(as.numeric(df$consistency)), NA_real_)
  }
  if (!"mean_NES_direction" %in% colnames(df)) {
    df$mean_NES_direction <- ifelse(
      df$has_significant_support,
      ifelse(df$mean_NES >= 0, "Positive NES", "Negative NES"),
      NA_character_
    )
  }
  if (!"min_padj" %in% colnames(df)) df$min_padj <- NA_real_
  df
}

build_lisa_contrast_summary <- function(summary_a, summary_b, contrast_a_label, contrast_b_label,
                                        universe, include_missing_categories, direction_epsilon,
                                        gsea_padj_cutoff = 0.25) {
  summary_a <- normalize_lisa_gsea_summary_support(
    summary_a, gsea_padj_cutoff, source = "summary_a"
  )
  summary_b <- normalize_lisa_gsea_summary_support(
    summary_b, gsea_padj_cutoff, source = "summary_b"
  )
  a <- suffix_lisa_contrast_summary(summary_a, "A")
  b <- suffix_lisa_contrast_summary(summary_b, "B")
  merged <- merge(a, b, by = "category_id", all = include_missing_categories)
  meta_cols <- c("category_display_name", "macrogroup_id", "macrogroup_name",
                 "macrogroup_order", "category_order_within_macrogroup", "color")
  for (cc in meta_cols) {
    a_meta <- merged[[paste0(cc, "_A")]]
    b_meta <- merged[[paste0(cc, "_B")]]
    conflicts <- !is.na(a_meta) & !is.na(b_meta) &
      as.character(a_meta) != as.character(b_meta)
    if (any(conflicts)) {
      stop(
        "LISA-CONTRAST-SUPPORT-001: A and B contain incompatible category metadata for ",
        paste(merged$category_id[conflicts], collapse = ", "),
        " (field ", cc, "). Use the same versioned category map on both sides.",
        call. = FALSE
      )
    }
    merged[[cc]] <- coalesce_local(a_meta, b_meta)
  }
  merged$universe <- universe
  merged$contrast_a_label <- contrast_a_label
  merged$contrast_b_label <- contrast_b_label
  merged$gsea_padj_cutoff <- gsea_padj_cutoff

  numeric_cols <- c(
    "n_genesets", "n_genesets_mapped", "n_genesets_evaluable",
    "n_genesets_significant", "mean_NES", "mean_NES_contextual",
    "median_NES", "consistency", "same_direction_pct", "min_padj",
    "min_padj_mapped"
  )
  for (prefix in c("A", "B")) {
    for (cc in numeric_cols) {
      nm <- paste0(cc, "_", prefix)
      if (nm %in% colnames(merged)) merged[[nm]] <- suppressWarnings(as.numeric(merged[[nm]]))
    }
    for (n_col in paste0(c(
      "n_genesets_", "n_genesets_mapped_", "n_genesets_evaluable_",
      "n_genesets_significant_"
    ), prefix)) {
      merged[[n_col]][is.na(merged[[n_col]])] <- 0
      merged[[n_col]] <- as.integer(merged[[n_col]])
    }
    support_col <- paste0("has_significant_support_", prefix)
    merged[[support_col]] <- merged[[paste0("n_genesets_significant_", prefix)]] > 0L
    merged[[paste0("n_genesets_", prefix)]] <- merged[[paste0("n_genesets_significant_", prefix)]]
    mean_col <- paste0("mean_NES_", prefix)
    merged[[mean_col]][!merged[[support_col]]] <- NA_real_
  }

  merged$plot_has_any_significant_support <-
    merged$has_significant_support_A | merged$has_significant_support_B
  merged$display_mean_NES_A <- ifelse(
    merged$plot_has_any_significant_support,
    ifelse(merged$has_significant_support_A, merged$mean_NES_A, merged$mean_NES_contextual_A),
    NA_real_
  )
  merged$display_mean_NES_B <- ifelse(
    merged$plot_has_any_significant_support,
    ifelse(merged$has_significant_support_B, merged$mean_NES_B, merged$mean_NES_contextual_B),
    NA_real_
  )
  incomplete_display <- merged$plot_has_any_significant_support &
    (!is.finite(merged$display_mean_NES_A) | !is.finite(merged$display_mean_NES_B))
  if (any(incomplete_display)) {
    stop(
      "LISA-CONTRAST-SUPPORT-002: a category with significant support on one side lacks a finite contextual mean on the other: ",
      paste(merged$category_id[incomplete_display], collapse = ", "),
      ". Rebuild both single-analysis summaries from evaluable mapped gene sets; absence is never imputed as zero.",
      call. = FALSE
    )
  }
  merged$endpoint_source_A <- ifelse(
    !merged$plot_has_any_significant_support, "hidden_no_category_signal",
    ifelse(merged$has_significant_support_A, "significant_mean",
      ifelse(is.finite(merged$display_mean_NES_A), "contextual_mean_not_significant", "context_unavailable"))
  )
  merged$endpoint_source_B <- ifelse(
    !merged$plot_has_any_significant_support, "hidden_no_category_signal",
    ifelse(merged$has_significant_support_B, "significant_mean",
      ifelse(is.finite(merged$display_mean_NES_B), "contextual_mean_not_significant", "context_unavailable"))
  )
  visible_a <- merged$plot_has_any_significant_support & is.finite(merged$display_mean_NES_A)
  visible_b <- merged$plot_has_any_significant_support & is.finite(merged$display_mean_NES_B)
  merged$direction_sign_A <- lisa_direction_sign(merged$display_mean_NES_A, visible_a, direction_epsilon)
  merged$direction_sign_B <- lisa_direction_sign(merged$display_mean_NES_B, visible_b, direction_epsilon)
  merged$direction_A <- lisa_direction_label(merged$direction_sign_A, visible_a)
  merged$direction_B <- lisa_direction_label(merged$direction_sign_B, visible_b)
  merged$direction_class <- ifelse(
    !merged$plot_has_any_significant_support,
    "no_significant_support",
    ifelse(
      !visible_a | !visible_b,
      "context_unavailable",
      ifelse(
        merged$direction_sign_A == 0L | merged$direction_sign_B == 0L,
        "one_or_both_zero",
        ifelse(
          merged$direction_sign_A == merged$direction_sign_B & merged$direction_sign_A > 0,
          "same_positive",
          ifelse(
            merged$direction_sign_A == merged$direction_sign_B & merged$direction_sign_A < 0,
            "same_negative",
            "opposite"
          )
        )
      )
    )
  )
  merged$is_same_direction <- merged$direction_class %in% c("same_positive", "same_negative")
  merged$is_opposite_direction <- merged$direction_class == "opposite"
  # A contrast is oriented exactly as declared: contrast A minus contrast B.
  # This also matches the report-foundation label and makes a positive delta
  # mean "stronger in A" rather than silently reversing the configured order.
  comparable <- visible_a & visible_b
  merged$delta_mean_NES <- ifelse(
    comparable,
    merged$display_mean_NES_A - merged$display_mean_NES_B,
    NA_real_
  )
  merged$abs_delta_mean_NES <- abs(merged$delta_mean_NES)
  merged$mean_abs_NES <- ifelse(
    comparable,
    (abs(merged$display_mean_NES_A) + abs(merged$display_mean_NES_B)) / 2,
    NA_real_
  )

  keep <- c(
    "universe", "category_id", "category_display_name", "macrogroup_id", "macrogroup_name",
    "macrogroup_order", "category_order_within_macrogroup", "color",
    "contrast_a_label", "contrast_b_label", "gsea_padj_cutoff",
    "n_genesets_A", "n_genesets_mapped_A", "n_genesets_evaluable_A", "n_genesets_significant_A",
    "n_pos_genesets_A", "n_neg_genesets_A", "mean_NES_A", "mean_NES_contextual_A",
    "mean_pos_NES_A", "mean_neg_NES_A", "median_NES_A", "consistency_A",
    "same_direction_pct_A", "mean_NES_direction_A", "min_padj_A", "min_padj_mapped_A",
    "has_significant_support_A", "display_mean_NES_A", "endpoint_source_A",
    "n_genesets_B", "n_genesets_mapped_B", "n_genesets_evaluable_B", "n_genesets_significant_B",
    "n_pos_genesets_B", "n_neg_genesets_B", "mean_NES_B", "mean_NES_contextual_B",
    "mean_pos_NES_B", "mean_neg_NES_B", "median_NES_B", "consistency_B",
    "same_direction_pct_B", "mean_NES_direction_B", "min_padj_B", "min_padj_mapped_B",
    "has_significant_support_B", "display_mean_NES_B", "endpoint_source_B",
    "plot_has_any_significant_support",
    "direction_A", "direction_B", "direction_class", "is_same_direction", "is_opposite_direction",
    "delta_mean_NES", "abs_delta_mean_NES", "mean_abs_NES"
  )
  keep <- intersect(keep, colnames(merged))
  out <- merged[, keep, drop = FALSE]
  order_lisa_contrast_summary(out, "supracategory")
}

suffix_lisa_contrast_summary <- function(df, suffix) {
  meta_cols <- c("category_id", "category_display_name", "macrogroup_id", "macrogroup_name",
                 "macrogroup_order", "category_order_within_macrogroup", "color")
  metric_cols <- c("n_genesets", "mean_NES", "median_NES", "consistency", "same_direction_pct",
                   "mean_NES_direction", "n_pos_genesets", "n_neg_genesets",
                   "mean_pos_NES", "mean_neg_NES", "min_padj",
                   "n_genesets_mapped", "n_genesets_evaluable",
                   "n_genesets_significant", "mean_NES_contextual",
                   "min_padj_mapped", "has_significant_support",
                   "gsea_padj_cutoff_source")
  count_cols <- c("n_genesets", "n_pos_genesets", "n_neg_genesets",
                  "n_genesets_mapped", "n_genesets_evaluable",
                  "n_genesets_significant")
  for (col in setdiff(metric_cols, colnames(df))) {
    if (col %in% count_cols) {
      df[[col]] <- 0
    } else {
      df[[col]] <- NA
    }
  }
  keep <- intersect(c(meta_cols, metric_cols), colnames(df))
  out <- df[, keep, drop = FALSE]
  names(out)[names(out) != "category_id"] <- paste0(names(out)[names(out) != "category_id"], "_", suffix)
  out
}

coalesce_local <- function(a, b) {
  out <- a
  idx <- is.na(out) | out == ""
  out[idx] <- b[idx]
  out
}

lisa_direction_sign <- function(mean_nes, visible, epsilon = 0) {
  mean_nes <- suppressWarnings(as.numeric(mean_nes))
  visible <- as.logical(visible)
  visible[is.na(visible)] <- FALSE
  out <- rep(0L, length(mean_nes))
  out[visible & is.finite(mean_nes) & mean_nes > epsilon] <- 1L
  out[visible & is.finite(mean_nes) & mean_nes < -epsilon] <- -1L
  out
}

lisa_direction_label <- function(direction_sign, visible = rep(TRUE, length(direction_sign))) {
  ifelse(!visible, "not_displayed",
    ifelse(direction_sign > 0, "positive", ifelse(direction_sign < 0, "negative", "zero")))
}

order_lisa_contrast_summary <- function(df, plot_order) {
  if (nrow(df) == 0) return(df)
  if (!"macrogroup_order" %in% colnames(df)) df$macrogroup_order <- 99
  if (!"category_order_within_macrogroup" %in% colnames(df)) df$category_order_within_macrogroup <- 99
  if (!"category_display_name" %in% colnames(df)) df$category_display_name <- df$category_id
  if (!"macrogroup_name" %in% colnames(df)) df$macrogroup_name <- "Other or unclassified"
  df$macrogroup_order <- suppressWarnings(as.numeric(df$macrogroup_order))
  df$category_order_within_macrogroup <- suppressWarnings(as.numeric(df$category_order_within_macrogroup))
  df$macrogroup_order[is.na(df$macrogroup_order)] <- 99
  df$category_order_within_macrogroup[is.na(df$category_order_within_macrogroup)] <- 99
  if (plot_order %in% c("supracategory", "fixed")) {
    df <- df[order(df$macrogroup_order, df$category_order_within_macrogroup, df$category_id), , drop = FALSE]
  } else if (plot_order == "delta_abs") {
    df <- df[order(-df$abs_delta_mean_NES, df$category_id), , drop = FALSE]
  } else if (plot_order == "mean_abs") {
    df <- df[order(-df$mean_abs_NES, df$category_id), , drop = FALSE]
  }
  rownames(df) <- NULL
  df
}

filter_lisa_contrast_plot_set <- function(df, plot_set) {
  if (nrow(df) == 0) return(df)
  if (plot_set == "same_direction") return(df[df$is_same_direction, , drop = FALSE])
  if (plot_set == "opposite_direction") return(df[df$is_opposite_direction, , drop = FALSE])
  # The canonical all-category view is a dictionary map, not a significance
  # filter. Rows with no supported side remain as labelled blank rows.
  df
}

plot_lisa_contrast_dumbbell <- function(contrast_summary, contrast_a_label, contrast_b_label,
                                        contrast_a_title, contrast_b_title, title,
                                        subtitle = sprintf("%s minus %s (A - B)", contrast_a_title, contrast_b_title),
                                        group_by_supracategory,
                                        annotate_gene_sets = TRUE) {
  lisa_assert_classified_plot_rows(contrast_summary, "contrast dumbbell")
  if (nrow(contrast_summary) == 0) {
    return(
      ggplot2::ggplot(data.frame(x = 0, y = 1, label = "No categories match this filter"),
                      ggplot2::aes(x = x, y = y, label = label)) +
        ggplot2::geom_text(size = 4, color = "grey35") +
        ggplot2::labs(title = title, subtitle = subtitle, x = "Mean NES", y = NULL) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_blank(),
          plot.background = ggplot2::element_rect(fill = "white", color = NA),
          panel.background = ggplot2::element_rect(fill = "white", color = NA)
        )
    )
  }
  df <- contrast_summary
  if (!"n_genesets_A" %in% names(df)) {
    df$n_genesets_A <- as.integer(is.finite(suppressWarnings(as.numeric(df$mean_NES_A))))
  }
  if (!"n_genesets_B" %in% names(df)) {
    df$n_genesets_B <- as.integer(is.finite(suppressWarnings(as.numeric(df$mean_NES_B))))
  }
  if (!"same_direction_pct_A" %in% names(df)) df$same_direction_pct_A <- NA_real_
  if (!"same_direction_pct_B" %in% names(df)) df$same_direction_pct_B <- NA_real_
  if (!"has_significant_support_A" %in% names(df)) {
    df$has_significant_support_A <- suppressWarnings(as.numeric(df$n_genesets_A)) > 0
  }
  if (!"has_significant_support_B" %in% names(df)) {
    df$has_significant_support_B <- suppressWarnings(as.numeric(df$n_genesets_B)) > 0
  }
  if (!"plot_has_any_significant_support" %in% names(df)) {
    df$plot_has_any_significant_support <-
      df$has_significant_support_A | df$has_significant_support_B
  }
  if (!"display_mean_NES_A" %in% names(df)) {
    contextual <- if ("mean_NES_contextual_A" %in% names(df)) df$mean_NES_contextual_A else df$mean_NES_A
    df$display_mean_NES_A <- ifelse(
      df$plot_has_any_significant_support,
      ifelse(df$has_significant_support_A, df$mean_NES_A, contextual),
      NA_real_
    )
  }
  if (!"display_mean_NES_B" %in% names(df)) {
    contextual <- if ("mean_NES_contextual_B" %in% names(df)) df$mean_NES_contextual_B else df$mean_NES_B
    df$display_mean_NES_B <- ifelse(
      df$plot_has_any_significant_support,
      ifelse(df$has_significant_support_B, df$mean_NES_B, contextual),
      NA_real_
    )
  }
  if (!"n_genesets_evaluable_A" %in% names(df)) df$n_genesets_evaluable_A <- df$n_genesets_A
  if (!"n_genesets_evaluable_B" %in% names(df)) df$n_genesets_evaluable_B <- df$n_genesets_B
  if (!"gsea_padj_cutoff" %in% names(df)) df$gsea_padj_cutoff <- NA_real_
  if (!"category_display_name" %in% names(df)) df$category_display_name <- df$category_id
  if (!"macrogroup_name" %in% names(df)) df$macrogroup_name <- "Other or unclassified"
  if (!"color" %in% names(df)) df$color <- "#737373"
  df$category_display_name[is.na(df$category_display_name) | df$category_display_name == ""] <- df$category_id[is.na(df$category_display_name) | df$category_display_name == ""]
  df$macrogroup_name[is.na(df$macrogroup_name) | df$macrogroup_name == ""] <- "Other or unclassified"
  df$color[is.na(df$color) | df$color == ""] <- "#737373"
  df$category_plot_id <- make.unique(paste(df$category_id, df$category_display_name, sep = "__"))
  df$category_plot_id <- factor(df$category_plot_id, levels = rev(unique(df$category_plot_id)))
  df$macrogroup_name <- factor(df$macrogroup_name, levels = unique(df$macrogroup_name))
  label_values <- stats::setNames(df$category_display_name, as.character(df$category_plot_id))
  finite_nes <- c(df$display_mean_NES_A, df$display_mean_NES_B)
  finite_nes <- finite_nes[is.finite(finite_nes)]
  max_abs_nes <- if (length(finite_nes) > 0) max(abs(finite_nes), na.rm = TRUE) else 1
  if (!is.finite(max_abs_nes) || max_abs_nes == 0) max_abs_nes <- 1
  x_limit_multiplier <- if (isTRUE(annotate_gene_sets)) 1.52 else 1.16
  x_limits <- c(-max_abs_nes, max_abs_nes) * x_limit_multiplier
  label_pad <- max_abs_nes * 0.055
  a_is_left_endpoint <- df$display_mean_NES_A <= df$display_mean_NES_B
  a_is_left_endpoint[is.na(a_is_left_endpoint)] <- TRUE

  point_df <- rbind(
    data.frame(
      category_plot_id = df$category_plot_id,
      macrogroup_name = df$macrogroup_name,
      mean_NES = df$display_mean_NES_A,
      n_genesets = df$n_genesets_A,
      n_genesets_evaluable = df$n_genesets_evaluable_A,
      same_direction_pct = df$same_direction_pct_A,
      has_significant_support = df$has_significant_support_A,
      category_has_signal = df$plot_has_any_significant_support,
      contrast = contrast_a_label,
      label_side = ifelse(a_is_left_endpoint, -1, 1),
      point_fill = vapply(df$color, lighten_color, character(1), amount = 0.55),
      stringsAsFactors = FALSE
    ),
    data.frame(
      category_plot_id = df$category_plot_id,
      macrogroup_name = df$macrogroup_name,
      mean_NES = df$display_mean_NES_B,
      n_genesets = df$n_genesets_B,
      n_genesets_evaluable = df$n_genesets_evaluable_B,
      same_direction_pct = df$same_direction_pct_B,
      has_significant_support = df$has_significant_support_B,
      category_has_signal = df$plot_has_any_significant_support,
      contrast = contrast_b_label,
      label_side = ifelse(a_is_left_endpoint, 1, -1),
      point_fill = df$color,
      stringsAsFactors = FALSE
    )
  )
  point_df$contrast <- factor(point_df$contrast, levels = c(contrast_a_label, contrast_b_label))
  point_df$n_genesets <- suppressWarnings(as.numeric(point_df$n_genesets))
  point_df$n_genesets_evaluable <- suppressWarnings(as.numeric(point_df$n_genesets_evaluable))
  point_df$same_direction_pct <- suppressWarnings(as.numeric(point_df$same_direction_pct))
  point_df$endpoint_visible <- point_df$category_has_signal & is.finite(point_df$mean_NES)
  point_df$point_label <- ifelse(
    point_df$has_significant_support & is.finite(point_df$same_direction_pct),
    sprintf("n_sig=%s; %.0f%%", format(point_df$n_genesets, trim = TRUE, scientific = FALSE), point_df$same_direction_pct),
    ifelse(
      point_df$has_significant_support,
      sprintf("n_sig=%s", format(point_df$n_genesets, trim = TRUE, scientific = FALSE)),
      sprintf("NS; n_eval=%s", format(point_df$n_genesets_evaluable, trim = TRUE, scientific = FALSE))
    )
  )
  point_df$label_x <- point_df$mean_NES + point_df$label_side * label_pad
  point_df$label_hjust <- ifelse(point_df$label_side < 0, 1, 0)

  segment_df <- df[
    df$plot_has_any_significant_support &
      is.finite(df$display_mean_NES_A) & is.finite(df$display_mean_NES_B),
    , drop = FALSE
  ]
  visible_points <- point_df[point_df$endpoint_visible, , drop = FALSE]
  supported_points <- visible_points[visible_points$has_significant_support, , drop = FALSE]
  contextual_points <- visible_points[!visible_points$has_significant_support, , drop = FALSE]

  p <- ggplot2::ggplot(df, ggplot2::aes(y = category_plot_id)) +
    ggplot2::geom_blank(ggplot2::aes(x = 0, y = category_plot_id)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.35) +
    ggplot2::geom_segment(
      data = segment_df,
      ggplot2::aes(x = display_mean_NES_A, xend = display_mean_NES_B,
        y = category_plot_id, yend = category_plot_id, color = color),
      inherit.aes = FALSE,
      linewidth = 0.55,
      alpha = 0.70
    ) +
    ggplot2::geom_point(
      data = supported_points,
      ggplot2::aes(x = mean_NES, y = category_plot_id, shape = contrast, size = n_genesets, fill = point_fill),
      inherit.aes = FALSE,
      color = "grey18",
      stroke = 0.28,
      alpha = 0.96
    ) +
    ggplot2::geom_point(
      data = contextual_points,
      ggplot2::aes(x = mean_NES, y = category_plot_id, shape = contrast, size = n_genesets),
      inherit.aes = FALSE,
      fill = NA,
      color = "grey35",
      stroke = 0.75,
      alpha = 0.38
    )
  if (isTRUE(annotate_gene_sets)) {
    p <- p +
      ggplot2::geom_text(
        data = visible_points[visible_points$contrast == contrast_a_label, , drop = FALSE],
        ggplot2::aes(x = label_x, y = category_plot_id, label = point_label, hjust = label_hjust),
        inherit.aes = FALSE,
        nudge_y = 0.17,
        size = 3.25,
        color = "grey25",
        lineheight = 0.9,
        show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = visible_points[visible_points$contrast == contrast_b_label, , drop = FALSE],
        ggplot2::aes(x = label_x, y = category_plot_id, label = point_label, hjust = label_hjust),
        inherit.aes = FALSE,
        nudge_y = -0.17,
        size = 3.25,
        color = "grey25",
        lineheight = 0.9,
        show.legend = FALSE
    )
  }
  contrast_levels <- c(contrast_a_label, contrast_b_label)
  p <- p +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_shape_manual(
      values = stats::setNames(c(21, 22), contrast_levels),
      limits = contrast_levels,
      drop = FALSE
    ) +
    ggplot2::scale_size_continuous(range = c(1.7, 6.2), breaks = c(0, 1, 5, 10, 20), name = "n gene sets") +
    # With free-y facets, unused factor levels must be dropped per panel.
    # Keeping every global category level in every macrogroup duplicates the
    # y axis (N categories x M macrogroup panels) and makes the plot illegible.
    # geom_blank() above still trains and preserves genuine blank category
    # rows within the macrogroup to which they belong.
    ggplot2::scale_y_discrete(labels = label_values, drop = TRUE) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "Mean NES (significant-set mean; contextual pre-threshold mean for an unsupported side)",
      y = NULL,
      shape = "Contrast",
      caption = sprintf(
        paste0(
          "Filled: significant support at GSEA FDR <= %s.\n",
          "Hollow/translucent: contextual mean without significant support. Blank row: neither side significant."
        ),
        format(unique(df$gsea_padj_cutoff)[1], trim = TRUE)
      )
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(size = 5.2, fill = "grey75", color = "grey18", stroke = 0.45)),
      size = ggplot2::guide_legend(override.aes = list(shape = 21, fill = "grey75", color = "grey18"))
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey35"),
      plot.caption = ggplot2::element_text(
        hjust = 0, size = 7.5, lineheight = 0.95,
        margin = ggplot2::margin(t = 5)
      ),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = grid::unit(c(5.5, 18, 5.5, 5.5), "pt")
    )
  if (isTRUE(group_by_supracategory) && "macrogroup_name" %in% colnames(df)) {
    p <- p +
      ggplot2::facet_grid(macrogroup_name ~ ., scales = "free_y", space = "free_y", switch = "y") +
      ggplot2::theme(
        strip.placement = "outside",
        strip.background.y = ggplot2::element_rect(fill = "grey94", color = "grey80"),
        strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 8, face = "bold"),
        panel.spacing.y = grid::unit(0.25, "lines")
      )
  }
  p
}

lisa_contrast_plot_title <- function(plot_set, universe, title_prefix = "LISA category contrast dumbbell") {
  label <- switch(
    plot_set,
    all = "all categories",
    same_direction = "same-direction categories",
    opposite_direction = "opposite-direction categories",
    plot_set
  )
  sprintf("%s: %s (%s)", title_prefix, label, universe)
}

lisa_contrast_plot_height <- function(n_categories) {
  max(4.8, 3.2 + 0.23 * max(1, n_categories))
}

lisa_install_figure_recipe <- function(target, code_ledger = NULL) {
  package_dir <- lisa_resolve_package_dir()
  renderer <- lisa_post_script_path(package_dir, "reproduce_lisa_figure.R")
  if (is.null(code_ledger)) {
    expected_sha256 <- lisa_sha256_file(renderer)
  } else {
    recipe_identity <- lisa_verify_code_recipe_ledger(
      code_ledger, package_dir, "reproduce_lisa_figure.R"
    )
    renderer <- recipe_identity$path
    expected_sha256 <- recipe_identity$sha256
  }
  lisa_copy_verified_code_file(renderer, target, expected_sha256)
  invisible(target)
}

lisa_encode_figure_source_text <- function(x) {
  x <- enc2utf8(as.character(x))
  x <- gsub("%", "%25", x, fixed = TRUE)
  x <- gsub("\t", "%09", x, fixed = TRUE)
  x <- gsub("\r", "%0D", x, fixed = TRUE)
  gsub("\n", "%0A", x, fixed = TRUE)
}

lisa_write_contrast_figure_source <- function(df, path, comparison_name, title,
                                              subtitle, plot_set,
                                              annotation_variant,
                                              contrast_a_label,
                                              contrast_b_label,
                                              gsea_padj_cutoff,
                                              group_by_supracategory = TRUE,
                                              figure_width = NULL,
                                              figure_height = NULL,
                                              figure_dpi = 300) {
  if (length(plot_set) != 1L || is.na(plot_set) || !nzchar(as.character(plot_set)) ||
      length(annotation_variant) != 1L || is.na(annotation_variant) ||
      !nzchar(as.character(annotation_variant)) ||
      length(group_by_supracategory) != 1L || is.na(group_by_supracategory)) {
    stop(
      "LISA-FIGURE-SOURCE-001 dumbbell source metadata requires one plot set and annotation variant.",
      call. = FALSE
    )
  }
  plot_set <- as.character(plot_set)
  annotation_variant <- tolower(as.character(annotation_variant))
  if (!annotation_variant %in% c("plain", "annotated")) {
    stop(
      "LISA-FIGURE-SOURCE-001 dumbbell annotation variant must be plain or annotated.",
      call. = FALSE
    )
  }
  n_categories <- nrow(df)
  if (is.null(figure_width)) {
    figure_width <- if (identical(annotation_variant, "annotated")) 15.5 else 12
  }
  if (is.null(figure_height)) {
    figure_height <- lisa_contrast_plot_height(n_categories)
    if (identical(annotation_variant, "annotated")) figure_height <- figure_height * 1.08
  }
  dimensions <- c(
    figure_width = suppressWarnings(as.numeric(figure_width)),
    figure_height = suppressWarnings(as.numeric(figure_height)),
    figure_dpi = suppressWarnings(as.numeric(figure_dpi))
  )
  if (any(lengths(list(figure_width, figure_height, figure_dpi)) != 1L) ||
      any(!is.finite(dimensions)) || any(dimensions <= 0)) {
    stop(
      "LISA-FIGURE-SOURCE-001 dumbbell figure dimensions and DPI must be positive scalars.",
      call. = FALSE
    )
  }
  text_fields <- list(
    comparison_name = comparison_name,
    title = title,
    subtitle = subtitle,
    contrast_a_label = contrast_a_label,
    contrast_b_label = contrast_b_label
  )
  text_ok <- vapply(text_fields, function(value) {
    length(value) == 1L && !is.na(value) && nzchar(as.character(value))
  }, logical(1))
  if (!all(text_ok)) {
    stop(
      "LISA-FIGURE-SOURCE-001 dumbbell identity, title, subtitle and contrast labels must be non-empty scalars.",
      call. = FALSE
    )
  }
  title_encoded <- lisa_encode_figure_source_text(title)
  subtitle_encoded <- lisa_encode_figure_source_text(subtitle)
  contrast_a_label_encoded <- lisa_encode_figure_source_text(contrast_a_label)
  contrast_b_label_encoded <- lisa_encode_figure_source_text(contrast_b_label)
  source <- df
  empty_plot <- nrow(source) == 0L
  if (isTRUE(empty_plot)) {
    # A header-only TSV cannot retain the renderer, titles or A/B labels.
    # Preserve one explicitly non-scientific metadata row so an empty filtered
    # plot remains exactly reproducible without inventing a category or NES.
    source <- as.data.frame(
      lapply(source, function(column) column[NA_integer_]),
      stringsAsFactors = FALSE, optional = TRUE
    )
  }
  row_count <- nrow(source)
  source$figure_id <- rep(paste0(
    "lisa_dumbbell__", comparison_name, "__", plot_set, "__", annotation_variant
  ), row_count)
  source$figure_type <- rep("lisa_dumbbell", row_count)
  source$source_row_order <- if (isTRUE(empty_plot)) NA_integer_ else seq_len(row_count)
  source$selected_for_plot <- rep(!empty_plot, row_count)
  source$figure_record_kind <- rep(
    if (isTRUE(empty_plot)) "empty_plot_metadata" else "category",
    row_count
  )
  source$figure_title <- rep(title_encoded, row_count)
  source$figure_subtitle <- rep(subtitle_encoded, row_count)
  source$plot_set <- rep(plot_set, row_count)
  source$annotation_variant <- rep(annotation_variant, row_count)
  source$group_by_supracategory <- rep(isTRUE(group_by_supracategory), row_count)
  source$figure_width <- rep(unname(dimensions[["figure_width"]]), row_count)
  source$figure_height <- rep(unname(dimensions[["figure_height"]]), row_count)
  source$figure_dpi <- rep(unname(dimensions[["figure_dpi"]]), row_count)
  source$contrast_a_label <- rep(contrast_a_label_encoded, row_count)
  source$contrast_b_label <- rep(contrast_b_label_encoded, row_count)
  source$gsea_padj_cutoff <- rep(gsea_padj_cutoff, row_count)
  write_tsv_local(source, path)
  persisted <- read_tsv_local(path)
  # Decimal TSV serialization is not guaranteed to preserve a double's exact
  # bit pattern. Keep a tight finite numeric comparison while all textual and
  # logical contract fields below remain exact.
  numeric_metadata_matches <- function(observed, expected) {
    observed <- suppressWarnings(as.numeric(observed))
    expected <- suppressWarnings(as.numeric(expected))
    length(observed) == row_count &&
      length(expected) == 1L &&
      all(is.finite(observed)) &&
      is.finite(expected) &&
      isTRUE(all.equal(
        observed, rep(expected, row_count),
        tolerance = .Machine$double.eps^0.5,
        check.attributes = FALSE
      ))
  }
  metadata_ok <- nrow(persisted) == row_count &&
    all(persisted$plot_set == plot_set) &&
    all(persisted$annotation_variant == annotation_variant) &&
    all(persisted$group_by_supracategory == isTRUE(group_by_supracategory)) &&
    numeric_metadata_matches(
      persisted$figure_width, dimensions[["figure_width"]]
    ) &&
    numeric_metadata_matches(
      persisted$figure_height, dimensions[["figure_height"]]
    ) &&
    numeric_metadata_matches(
      persisted$figure_dpi, dimensions[["figure_dpi"]]
    ) &&
    all(persisted$figure_title == title_encoded) &&
    all(persisted$figure_subtitle == subtitle_encoded) &&
    all(persisted$contrast_a_label == contrast_a_label_encoded) &&
    all(persisted$contrast_b_label == contrast_b_label_encoded)
  if (!isTRUE(metadata_ok)) {
    stop(
      "LISA-FIGURE-SOURCE-002 persisted dumbbell source metadata differs from the rendered figure contract.",
      call. = FALSE
    )
  }
  invisible(path)
}

write_lisa_contrast_report <- function(path, comparison_name, universe, input_manifest,
                                       contrast_summary, plot_sets, outputs) {
  lines <- c(
    sprintf("LISA contrast report: %s", comparison_name),
    sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("Collection: %s", universe),
    "",
    sprintf("Contrast A: %s", input_manifest$contrast_a_title[1]),
    sprintf("Contrast A summary: %s", input_manifest$contrast_a_summary[1]),
    sprintf("Contrast B: %s", input_manifest$contrast_b_title[1]),
    sprintf("Contrast B summary: %s", input_manifest$contrast_b_summary[1]),
    sprintf("GSEA adjusted-P cutoff: %s", input_manifest$gsea_padj_cutoff[1]),
    "",
    sprintf("Categories: %s", nrow(contrast_summary)),
    sprintf("Categories with signal on at least one side: %s",
      sum(contrast_summary$plot_has_any_significant_support, na.rm = TRUE)),
    sprintf("Blank categories (neither side significant): %s",
      sum(!contrast_summary$plot_has_any_significant_support, na.rm = TRUE)),
    sprintf("Same direction: %s", sum(contrast_summary$is_same_direction, na.rm = TRUE)),
    sprintf("Opposite direction: %s", sum(contrast_summary$is_opposite_direction, na.rm = TRUE)),
    sprintf("Supported zero endpoint: %s", sum(
      contrast_summary$plot_has_any_significant_support &
        (contrast_summary$direction_A == "zero" | contrast_summary$direction_B == "zero"),
      na.rm = TRUE
    )),
    sprintf("Plot sets: %s", paste(plot_sets, collapse = ", ")),
    "",
    "Output directories:",
    sprintf("lisa_tables: %s", outputs$lisa_tables),
    sprintf("qc: %s", outputs$qc),
    sprintf("plots: %s", outputs$plots)
  )
  lisa_guarded_write(path, function(target) writeLines(lines, target))
}

fill_summary_zeros <- function(df, cols) {
  for (cc in cols) if (cc %in% colnames(df)) df[[cc]][is.na(df[[cc]])] <- 0
  if ("min_padj" %in% colnames(df)) df$min_padj[is.infinite(df$min_padj)] <- NA_real_
  df
}

rbind_fill_local <- function(xs) {
  xs <- xs[vapply(xs, nrow, integer(1)) > 0]
  if (length(xs) == 0) return(data.frame())
  cols <- unique(unlist(lapply(xs, colnames)))
  xs <- lapply(xs, function(x) {
    missing <- setdiff(cols, colnames(x))
    for (cc in missing) x[[cc]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, xs)
}

order_category_summary <- function(df) {
  df[order(as.numeric(df$macrogroup_order), as.numeric(df$category_order_within_macrogroup), df$category_id), , drop = FALSE]
}

alphabetize_category_order <- function(category_map) {
  if (!"display_name" %in% colnames(category_map)) {
    category_map$display_name <- category_map$category_id
  }
  category_map$macrogroup_name[is.na(category_map$macrogroup_name) | category_map$macrogroup_name == ""] <- "Other or unclassified"
  category_map$display_name[is.na(category_map$display_name) | category_map$display_name == ""] <- category_map$category_id[is.na(category_map$display_name) | category_map$display_name == ""]
  category_map$macrogroup_sort_key <- tolower(category_map$macrogroup_name)
  category_map$category_sort_key <- tolower(category_map$display_name)
  macro_levels <- sort(unique(category_map$macrogroup_sort_key))
  macro_order <- stats::setNames(seq_along(macro_levels), macro_levels)
  category_map$macrogroup_order <- unname(macro_order[category_map$macrogroup_sort_key])
  category_map <- category_map[order(category_map$macrogroup_order, category_map$category_sort_key, category_map$category_id), , drop = FALSE]
  category_map$category_order_within_macrogroup <- ave(
    seq_len(nrow(category_map)),
    category_map$macrogroup_order,
    FUN = seq_along
  )
  category_map$macrogroup_sort_key <- NULL
  category_map$category_sort_key <- NULL
  rownames(category_map) <- NULL
  category_map
}

coverage_table <- function(gsea_all, ora_all, gsea_ledger = NULL) {
  one <- function(df, label) {
    if (nrow(df) == 0) return(data.frame(
      analysis = label,
      total_rows = 0L,
      mapped_rows = 0L,
      unmapped_rows = 0L,
      pct_mapped = NA_real_,
      eligible_rows = NA_integer_,
      tested_rows = NA_integer_
    ))
    data.frame(
      analysis = label,
      total_rows = nrow(df),
      mapped_rows = sum(!is.na(df$category_id)),
      unmapped_rows = sum(is.na(df$category_id)),
      pct_mapped = round(100 * sum(!is.na(df$category_id)) / nrow(df), 2),
      eligible_rows = NA_integer_,
      tested_rows = NA_integer_,
      stringsAsFactors = FALSE
    )
  }
  gsea <- if (!is.null(gsea_ledger) && nrow(gsea_ledger)) {
    classified <- as.character(gsea_ledger$classification_status) ==
      "classified"
    data.frame(
      analysis = "GSEA",
      total_rows = nrow(gsea_ledger),
      mapped_rows = sum(classified),
      unmapped_rows = sum(!classified),
      pct_mapped = round(100 * sum(classified) / nrow(gsea_ledger), 2),
      eligible_rows = sum(gsea_ledger$eligible_for_gsea %in% TRUE),
      tested_rows = sum(
        as.character(gsea_ledger$gsea_result_status) == "tested"
      ),
      stringsAsFactors = FALSE
    )
  } else {
    one(gsea_all, "GSEA")
  }
  rbind(gsea, one(ora_all, "ORA"))
}

unmapped_table <- function(df, label) {
  if (nrow(df) == 0 || !"category_id" %in% colnames(df)) return(data.frame())
  out <- df[is.na(df$category_id), , drop = FALSE]
  if (nrow(out) == 0) return(data.frame())
  out$analysis_label <- label
  out
}

infer_source_family <- function(gene_set_id) {
  x <- as.character(gene_set_id)
  out <- rep("CP_OTHER", length(x))
  out[grepl("^GOBP_", x)] <- "GO:BP"
  out[grepl("^GOMF_", x)] <- "GO:MF"
  out[grepl("^GOCC_", x)] <- "GO:CC"
  out[grepl("^BIOCARTA_", x)] <- "BIOCARTA"
  out[grepl("^PID_", x)] <- "PID"
  out[grepl("^REACTOME_", x)] <- "REACTOME"
  out[grepl("^WP_", x)] <- "WIKIPATHWAYS"
  out[grepl("^KEGG_", x)] <- "KEGG_LEGACY"
  out[grepl("^KEGG_MEDICUS_|^N[0-9]", x)] <- "KEGG_MEDICUS"
  out
}

augment_category_map <- function(category_map, lisa_dict) {
  category_map <- category_map[!duplicated(category_map$category_id), , drop = FALSE]
  required <- unique(lisa_dict[, c("category_id", "category_display_name"), drop = FALSE])
  category_map <- category_map[
    as.character(category_map$category_id) %in%
      as.character(required$category_id),
    , drop = FALSE
  ]
  missing <- required[!required$category_id %in% category_map$category_id, , drop = FALSE]
  if (nrow(missing) == 0) return(category_map)
  missing$macrogroup_id <- ifelse(grepl("^MF_", missing$category_id), "MOLECULAR_FUNCTION",
                                  ifelse(grepl("^CC_", missing$category_id), "CELLULAR_COMPONENT", "UNCLASSIFIED"))
  missing$macrogroup_name <- ifelse(missing$macrogroup_id == "MOLECULAR_FUNCTION", "Molecular function",
                                    ifelse(missing$macrogroup_id == "CELLULAR_COMPONENT", "Cellular component", "Other or unclassified"))
  missing$macrogroup_order <- ifelse(missing$macrogroup_id == "MOLECULAR_FUNCTION", 12,
                                     ifelse(missing$macrogroup_id == "CELLULAR_COMPONENT", 13, 99))
  missing <- missing[order(missing$macrogroup_order, missing$category_id), , drop = FALSE]
  missing$category_order_within_macrogroup <- ave(
    seq_len(nrow(missing)),
    missing$macrogroup_id,
    FUN = seq_along
  )
  missing$notes <- "Auto-added from active LISA dictionary because this category is outside the GOBP-C2 macrogroup map."
  names(missing)[names(missing) == "category_display_name"] <- "display_name"
  common <- union(colnames(category_map), colnames(missing))
  for (cc in setdiff(common, colnames(category_map))) category_map[[cc]] <- NA
  for (cc in setdiff(common, colnames(missing))) missing[[cc]] <- NA
  rbind(category_map[, common, drop = FALSE], missing[, common, drop = FALSE])
}

prepare_registered_category_map <- function(category_map, lisa_dict) {
  active_categories <- unique(as.character(lisa_dict$category_id))
  category_map <- category_map[
    as.character(category_map$category_id) %in% active_categories,
    , drop = FALSE
  ]
  category_map <- category_map[order(
    as.numeric(category_map$macrogroup_order),
    as.numeric(category_map$category_order_within_macrogroup),
    as.character(category_map$category_id)
  ), , drop = FALSE]
  rownames(category_map) <- NULL
  category_map
}

build_lisa_palette <- function(category_map, palette) {
  if (palette != "lisa_default") warning("Only palette = 'lisa_default' is implemented in this first version.", call. = FALSE)
  if (!nrow(category_map)) return(character())
  macrogroup_id <- as.character(category_map$macrogroup_id)
  macrogroup_id[is.na(macrogroup_id) | !nzchar(macrogroup_id)] <-
    "UNCLASSIFIED"
  macrogroup_order <- if ("macrogroup_order" %in% names(category_map)) {
    suppressWarnings(as.numeric(category_map$macrogroup_order))
  } else {
    rep(Inf, nrow(category_map))
  }
  group_rows <- split(seq_len(nrow(category_map)), macrogroup_id)
  group_order <- vapply(group_rows, function(index) {
    value <- macrogroup_order[index]
    if (all(is.na(value))) Inf else min(value, na.rm = TRUE)
  }, numeric(1))
  group_ids <- names(group_rows)[order(group_order, names(group_rows))]
  group_colours <- stats::setNames(
    grDevices::hcl.colors(length(group_ids), palette = "Dark 3"),
    group_ids
  )
  cols <- character(nrow(category_map))
  for (mg in group_ids) {
    idx <- group_rows[[mg]]
    if (length(idx) == 0) next
    category_order <- if ("category_order_within_macrogroup" %in%
                          names(category_map)) {
      suppressWarnings(as.numeric(
        category_map$category_order_within_macrogroup[idx]
      ))
    } else {
      rep(Inf, length(idx))
    }
    category_id <- if ("category_id" %in% names(category_map)) {
      as.character(category_map$category_id[idx])
    } else {
      as.character(idx)
    }
    ordered_idx <- idx[order(category_order, category_id)]
    base <- unname(group_colours[[mg]])
    shades <- if (length(ordered_idx) == 1L) {
      base
    } else {
      grDevices::colorRampPalette(c(
        lighten_color(base, 0.22), base, darken_color(base, 0.18)
      ))(length(ordered_idx))
    }
    cols[ordered_idx] <- shades
  }
  cols
}

lighten_color <- function(color, amount = 0.25) {
  rgb <- grDevices::col2rgb(color) / 255
  rgb <- rgb + (1 - rgb) * amount
  grDevices::rgb(rgb[1], rgb[2], rgb[3])
}

darken_color <- function(color, amount = 0.25) {
  rgb <- grDevices::col2rgb(color) / 255
  rgb <- rgb * (1 - amount)
  grDevices::rgb(rgb[1], rgb[2], rgb[3])
}

prepare_category_plot_df <- function(summary_df, plot_order) {
  df <- summary_df
  if (plot_order %in% c("supracategory", "fixed")) {
    df <- order_category_summary(df)
  } else if (plot_order == "mean_NES" && "mean_NES" %in% colnames(df)) {
    df <- df[order(as.numeric(df$macrogroup_order), -abs(df$mean_NES), df$category_id), , drop = FALSE]
  } else if (plot_order == "n_genesets" && "n_genesets" %in% colnames(df)) {
    df <- df[order(as.numeric(df$macrogroup_order), -df$n_genesets, df$category_id), , drop = FALSE]
  } else if (plot_order == "consistency" && "consistency" %in% colnames(df)) {
    df <- df[order(as.numeric(df$macrogroup_order), -df$consistency, df$category_id), , drop = FALSE]
  }
  if ("macrogroup_name" %in% colnames(df)) {
    df$macrogroup_name <- factor(df$macrogroup_name, levels = unique(df$macrogroup_name))
  }
  df$category_plot_label <- factor(df$category_display_name, levels = rev(unique(df$category_display_name)))
  df
}

plot_lisa_gsea_lollipop <- function(summary_df, group_by_supracategory, plot_order) {
  lisa_assert_classified_plot_rows(summary_df, "GSEA category lollipop")
  df <- prepare_category_plot_df(summary_df, plot_order)
  df <- df[df$n_genesets > 0, , drop = FALSE]
  if (nrow(df) == 0) {
    return(
      ggplot2::ggplot(data.frame(x = 0, y = 1, label = "No significant LISA-mapped gene sets"),
                      ggplot2::aes(x = x, y = y, label = label)) +
        ggplot2::geom_text(size = 4, color = "grey35") +
        ggplot2::theme_void() +
        ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
    )
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_NES, y = category_plot_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = mean_NES, yend = category_plot_label), color = "grey70", linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(size = n_genesets, fill = color, alpha = n_genesets > 0), shape = 21, color = "grey25", stroke = 0.25) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.20), guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.5, 7), breaks = c(0, 1, 5, 10, 20)) +
    ggplot2::labs(x = "Mean NES across significant LISA-mapped gene sets", y = NULL, size = "n gene sets") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "right",
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  if (isTRUE(group_by_supracategory) && "macrogroup_name" %in% colnames(df)) {
    p <- p +
      ggplot2::facet_grid(macrogroup_name ~ ., scales = "free_y", space = "free_y", switch = "y") +
      ggplot2::theme(
        strip.placement = "outside",
        strip.background.y = ggplot2::element_rect(fill = "grey94", color = "grey80"),
        strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 8, face = "bold"),
        panel.spacing.y = grid::unit(0.25, "lines")
      )
  }
  p
}

lisa_apply_plot_context <- function(plot, display_title = NULL,
                                    comparison_subtitle = NULL,
                                    positive_direction = NULL) {
  nonempty <- function(x) {
    !is.null(x) && length(x) > 0L && !is.na(x[[1]]) &&
      nzchar(trimws(as.character(x[[1]])))
  }
  subtitle <- c(
    if (nonempty(comparison_subtitle)) as.character(comparison_subtitle[[1]]),
    if (nonempty(positive_direction)) as.character(positive_direction[[1]])
  )
  additions <- list()
  if (nonempty(display_title)) additions$title <- as.character(display_title[[1]])
  if (length(subtitle)) additions$subtitle <- paste(subtitle, collapse = "\n")
  if (!length(additions)) return(plot)
  plot +
    do.call(ggplot2::labs, additions) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "grey30")
    )
}

lisa_plot_device_text <- function(x) {
  # Base PDF and some bitmap devices cannot encode typographic dash glyphs
  # with their default font. Preserve source metadata and normalize only the
  # text object sent to the graphics device.
  gsub("[\u2013\u2014\u2212]", "-", as.character(x), perl = TRUE)
}

plot_lisa_gsea_direction_lollipop <- function(summary_df, group_by_supracategory, plot_order) {
  lisa_assert_classified_plot_rows(summary_df, "GSEA direction lollipop")
  df <- prepare_category_plot_df(summary_df, plot_order)
  df <- df[df$n_genesets > 0, , drop = FALSE]
  if (nrow(df) == 0) {
    return(
      ggplot2::ggplot(data.frame(x = 0, y = 1, label = "No significant LISA-mapped gene sets"),
                      ggplot2::aes(x = x, y = y, label = label)) +
        ggplot2::geom_text(size = 4, color = "grey35") +
        ggplot2::theme_void() +
        ggplot2::theme(plot.background = ggplot2::element_rect(fill = "white", color = NA))
    )
  }
  if (!"same_direction_pct" %in% colnames(df)) {
    df$same_direction_pct <- ifelse(df$n_genesets > 0, round(100 * df$consistency, 1), NA_real_)
  }
  df$direction_stats_label <- ifelse(
    df$n_genesets > 0,
    paste0(round(df$same_direction_pct, 0), "% (n=", df$n_genesets, ")"),
    ""
  )
  df$label_hjust <- ifelse(df$mean_NES >= 0, -0.15, 1.15)
  label_df <- df[df$n_genesets > 0, , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_NES, y = category_plot_label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = mean_NES, yend = category_plot_label), color = "grey70", linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(size = n_genesets, fill = color, alpha = n_genesets > 0), shape = 21, color = "grey25", stroke = 0.25) +
    ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(label = direction_stats_label, hjust = label_hjust),
      size = 2.7,
      fontface = "bold",
      color = "grey20"
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.20), guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.5, 7), breaks = c(0, 1, 5, 10, 20)) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.35, 0.35))) +
    ggplot2::labs(x = "Mean NES across significant LISA-mapped gene sets", y = NULL, size = "n gene sets") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "right",
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  if (isTRUE(group_by_supracategory) && "macrogroup_name" %in% colnames(df)) {
    p <- p +
      ggplot2::facet_grid(macrogroup_name ~ ., scales = "free_y", space = "free_y", switch = "y") +
      ggplot2::theme(
        strip.placement = "outside",
        strip.background.y = ggplot2::element_rect(fill = "grey94", color = "grey80"),
        strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 8, face = "bold"),
        panel.spacing.y = grid::unit(0.25, "lines")
      )
  }
  p
}

plot_lisa_gsea_dumbbell <- function(summary_df, group_by_supracategory, plot_order = "supracategory") {
  plot_lisa_gsea_lollipop(summary_df, group_by_supracategory, plot_order)
}

plot_lisa_ora_barplot <- function(summary_df, group_by_supracategory, plot_order) {
  lisa_assert_classified_plot_rows(summary_df, "ORA category barplot")
  df <- prepare_category_plot_df(summary_df, plot_order)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = total_overlap, y = category_plot_label, fill = color)) +
    ggplot2::geom_col(width = 0.65) +
    ggplot2::scale_fill_identity() +
    ggplot2::labs(x = "Total ORA overlapping genes across LISA-mapped gene sets", y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  if (isTRUE(group_by_supracategory) && "macrogroup_name" %in% colnames(df)) {
    p <- p +
      ggplot2::facet_grid(macrogroup_name ~ direction, scales = "free_y", space = "free_y", switch = "y") +
      ggplot2::theme(
        strip.placement = "outside",
        strip.background.y = ggplot2::element_rect(fill = "grey94", color = "grey80"),
        strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 8, face = "bold"),
        panel.spacing.y = grid::unit(0.25, "lines")
      )
  } else {
    p <- p + ggplot2::facet_wrap(~ direction, nrow = 1)
  }
  p
}

lisa_category_pathway_expected_count <- function(category_summary) {
  if (!nrow(category_summary) || !"n_genesets" %in% names(category_summary)) return(0L)
  as.integer(sum(!is.na(category_summary$n_genesets) & as.numeric(category_summary$n_genesets) > 0))
}

write_lisa_gsea_category_pathway_plots <- function(gsea_all, category_summary, gsea_padj_cutoff,
                                                   plot_order, output_dir, plot_formats,
                                                   bg = "white", label_chars = 62,
                                                   comparison_subtitle = NULL,
                                                   positive_direction = NULL,
                                                   stable_names = FALSE,
                                                   source_data = TRUE,
                                                   recipes = TRUE,
                                                   code_ledger = NULL) {
  lisa_assert_classified_plot_rows(gsea_all, "category member-set plots")
  lisa_assert_classified_plot_rows(category_summary, "category member-set manifest")
  lisa_guarded_dir_create(output_dir, output_dir)
  old_plot_files <- list.files(
    output_dir,
    pattern = "\\.(png|pdf|svg|tiff|tif|jpeg|jpg)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(old_plot_files) > 0) lapply(old_plot_files, lisa_guarded_delete, run_root = output_dir)
  old_manifest <- file.path(output_dir, "plot_manifest.tsv")
  if (file.exists(old_manifest)) lisa_guarded_delete(old_manifest, run_root = output_dir)
  old_source_files <- list.files(output_dir, pattern = "_(source|recipe)[.](tsv|R)$",
    full.names = TRUE, ignore.case = TRUE)
  if (length(old_source_files) > 0) lapply(old_source_files, lisa_guarded_delete, run_root = output_dir)

  categories <- order_category_summary(category_summary)
  categories <- categories[categories$n_genesets > 0, , drop = FALSE]
  if (nrow(categories) == 0) {
    manifest <- data.frame(
      category_order = integer(),
      category_id = character(),
      category_display_name = character(),
      macrogroup_name = character(),
      n_significant_genesets = integer(),
      fdr_cutoff = numeric(),
      color = character(),
      plot_stem = character(),
      source_tsv = character(),
      recipe_r = character(),
      stringsAsFactors = FALSE
    )
    write_tsv_local(manifest, file.path(output_dir, "plot_manifest.tsv"))
    return(invisible(manifest))
  }
  categories$category_display_name[is.na(categories$category_display_name) | categories$category_display_name == ""] <-
    categories$category_id[is.na(categories$category_display_name) | categories$category_display_name == ""]
  categories$macrogroup_name[is.na(categories$macrogroup_name) | categories$macrogroup_name == ""] <- "Other or unclassified"
  categories$color[is.na(categories$color) | categories$color == ""] <- "#737373"

  df <- gsea_all[!is.na(gsea_all$category_id) & !is.na(gsea_all$padj) & gsea_all$padj <= gsea_padj_cutoff, , drop = FALSE]
  if (nrow(df) > 0) {
    df <- prepare_pathway_plot_df(df, plot_order)
    finite_nes <- df$NES[is.finite(df$NES)]
    max_abs_nes <- if (length(finite_nes) > 0) max(abs(finite_nes), na.rm = TRUE) else 1
    if (!is.finite(max_abs_nes) || max_abs_nes == 0) max_abs_nes <- 1
    x_limits <- c(-max_abs_nes, max_abs_nes) * 1.08
  } else {
    x_limits <- c(-1, 1)
  }

  manifest <- vector("list", nrow(categories))
  for (ii in seq_len(nrow(categories))) {
    meta <- categories[ii, , drop = FALSE]
    cat_df <- df[df$category_id == meta$category_id[1], , drop = FALSE]
    p <- plot_lisa_gsea_one_category_pathways(
      cat_df = cat_df,
      category_meta = meta,
      gsea_padj_cutoff = gsea_padj_cutoff,
      label_chars = label_chars,
      x_limits = x_limits,
      comparison_subtitle = comparison_subtitle,
      positive_direction = positive_direction
    )
    file_label <- if (isTRUE(stable_names)) {
      safe_file_label(meta$category_id[1])
    } else {
      safe_file_label(sprintf("%02d_%s", ii, meta$category_id[1]))
    }
    stem <- file.path(output_dir, file_label)
    # This is deliberately a per-figure data contract, rather than an index or
    # shared enrichment table.  It retains candidates excluded by the FDR
    # selection and records the exact transforms used by the displayed plot.
    candidates <- gsea_all[!is.na(gsea_all$category_id) &
      gsea_all$category_id == meta$category_id[1], , drop = FALSE]
    candidates$figure_id <- paste0("lisa_category_gene_sets__", file_label)
    candidates$figure_type <- "lisa_category_gene_sets"
    candidates$source_row_order <- seq_len(nrow(candidates))
    candidates$category_display_name <- meta$category_display_name[1]
    candidates$macrogroup_name <- meta$macrogroup_name[1]
    candidates$category_color <- meta$color[1]
    candidates$plot_subtitle <- paste(c(
      meta$macrogroup_name[1],
      if (!is.null(comparison_subtitle)) as.character(comparison_subtitle[[1]]) else NULL,
      if (!is.null(positive_direction)) as.character(positive_direction[[1]]) else NULL
    ), collapse = " | ")
    candidates$figure_width <- 12.5
    candidates$figure_height <- category_pathway_plot_height(nrow(cat_df))
    candidates$figure_dpi <- 300
    candidates$figure_background <- bg
    candidates$figure_label_chars <- label_chars
    candidates$figure_comparison_subtitle <- if (is.null(comparison_subtitle)) "" else comparison_subtitle[[1L]]
    candidates$figure_positive_direction <- if (is.null(positive_direction)) "" else positive_direction[[1L]]
    candidates$selected_for_plot <- !is.na(candidates$padj) & candidates$padj <= gsea_padj_cutoff
    candidates$highlighted <- candidates$selected_for_plot
    candidates$labelled <- candidates$selected_for_plot
    candidates$selection_fdr_cutoff <- gsea_padj_cutoff
    candidates$plot_x_nes <- candidates$NES
    candidates$plot_point_size_neg_log10_fdr <- -log10(pmax(candidates$padj, 1e-300))
    candidates$plot_x_min <- x_limits[[1]]
    candidates$plot_x_max <- x_limits[[2]]
    candidates$plotted_order <- NA_integer_
    candidates$pathway_label <- ""
    if (nrow(cat_df)) {
      plotted <- cat_df[order(as.numeric(cat_df$source_family_order), cat_df$pathway_order_key,
        cat_df$pathway), , drop = FALSE]
      plotted$pathway_label <- trim_label(plotted$pathway_display_label %||% plotted$pathway, label_chars)
      hit <- match(candidates$pathway, plotted$pathway)
      candidates$plotted_order <- match(candidates$pathway, plotted$pathway)
      candidates$pathway_label <- ifelse(is.na(hit), "", plotted$pathway_label[hit])
    }
    source_tsv <- if (isTRUE(source_data)) paste0(stem, "_source.tsv") else ""
    recipe_r <- if (isTRUE(recipes)) paste0(stem, "_recipe.R") else ""
    if (isTRUE(source_data)) write_tsv_local(candidates, source_tsv)
    if (isTRUE(recipes)) {
      lisa_install_figure_recipe(recipe_r, code_ledger = code_ledger)
    }
    plot_height <- category_pathway_plot_height(nrow(cat_df))
    save_plot_multi(p, stem, plot_formats, width = 12.5, height = plot_height, bg = bg)
    manifest[[ii]] <- data.frame(
      category_order = ii,
      category_id = meta$category_id[1],
      category_display_name = meta$category_display_name[1],
      macrogroup_name = meta$macrogroup_name[1],
      n_significant_genesets = nrow(cat_df),
      fdr_cutoff = gsea_padj_cutoff,
      color = meta$color[1],
      plot_stem = stem,
      source_tsv = source_tsv,
      recipe_r = recipe_r,
      stringsAsFactors = FALSE
    )
  }

  manifest <- do.call(rbind, manifest)
  write_tsv_local(manifest, file.path(output_dir, "plot_manifest.tsv"))
  invisible(manifest)
}

plot_lisa_gsea_one_category_pathways <- function(cat_df, category_meta, gsea_padj_cutoff,
                                                 label_chars, x_limits,
                                                 comparison_subtitle = NULL,
                                                 positive_direction = NULL) {
  lisa_assert_classified_plot_rows(cat_df, "single category member-set plot")
  lisa_assert_classified_plot_rows(category_meta, "single category metadata")
  category_name <- category_meta$category_display_name[1]
  macrogroup_name <- category_meta$macrogroup_name[1]
  category_color <- category_meta$color[1] %||% "#737373"
  if (is.na(category_color) || category_color == "") category_color <- "#737373"
  context <- c(
    macrogroup_name,
    if (!is.null(comparison_subtitle) && nzchar(as.character(comparison_subtitle[[1]]))) {
      as.character(comparison_subtitle[[1]])
    },
    if (!is.null(positive_direction) && nzchar(as.character(positive_direction[[1]]))) {
      as.character(positive_direction[[1]])
    }
  )
  context <- paste(context, collapse = "\n")

  if (nrow(cat_df) == 0) {
    return(
      ggplot2::ggplot(data.frame(x = 0, y = 1, label = sprintf("No significant gene sets at FDR <= %s", gsea_padj_cutoff)),
                      ggplot2::aes(x = x, y = y, label = label)) +
        ggplot2::geom_text(size = 4, color = "grey35") +
        ggplot2::scale_x_continuous(limits = x_limits) +
        ggplot2::labs(title = category_name, subtitle = context, x = "NES", y = NULL) +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::theme(
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(face = "bold", size = 12),
          plot.background = ggplot2::element_rect(fill = "white", color = NA),
          panel.background = ggplot2::element_rect(fill = "white", color = NA)
        )
    )
  }

  cat_df <- cat_df[order(as.numeric(cat_df$source_family_order), cat_df$pathway_order_key, cat_df$pathway), , drop = FALSE]
  cat_df$pathway_label <- trim_label(cat_df$pathway_display_label %||% cat_df$pathway, label_chars)
  cat_df$pathway_plot_id <- make.unique(paste(cat_df$category_id, cat_df$pathway, sep = "__"))
  cat_df$pathway_plot_id <- factor(cat_df$pathway_plot_id, levels = rev(unique(cat_df$pathway_plot_id)))
  label_values <- stats::setNames(cat_df$pathway_label, as.character(cat_df$pathway_plot_id))
  label_size <- category_pathway_axis_text_size(nrow(cat_df))

  ggplot2::ggplot(cat_df, ggplot2::aes(x = NES, y = pathway_plot_id)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = pathway_plot_id),
                          color = lighten_color(category_color, 0.18), linewidth = 0.32) +
    ggplot2::geom_point(ggplot2::aes(size = -log10(pmax(padj, 1e-300))),
                        shape = 21, fill = category_color, color = "grey20", stroke = 0.18, alpha = 0.92) +
    ggplot2::scale_y_discrete(labels = label_values) +
    ggplot2::scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0.02, 0.02))) +
    ggplot2::scale_size_continuous(range = c(1.6, 5.2), name = "-log10(FDR)") +
    ggplot2::labs(
      title = sprintf("%s (n=%s)", category_name, nrow(cat_df)),
      subtitle = context,
      x = "NES",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = label_size, color = "grey15"),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey35"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

category_pathway_plot_height <- function(n_pathways) {
  max(4.5, 2.8 + 0.155 * max(1, n_pathways))
}

category_pathway_axis_text_size <- function(n_pathways) {
  if (n_pathways <= 50) return(7.2)
  if (n_pathways <= 120) return(6.2)
  if (n_pathways <= 220) return(5.4)
  4.8
}

plot_lisa_pathway_dotplot <- function(gsea_all, gsea_padj_cutoff, group_by_supracategory, plot_order) {
  lisa_assert_classified_plot_rows(gsea_all, "GSEA pathway dotplot")
  df <- gsea_all[!is.na(gsea_all$category_id) & !is.na(gsea_all$padj) & gsea_all$padj <= gsea_padj_cutoff, , drop = FALSE]
  if (nrow(df) == 0) df <- gsea_all[!is.na(gsea_all$category_id), , drop = FALSE]
  df <- add_pathway_order_metadata(df)
  df <- df[order(df$padj, -abs(df$NES), df$pathway), , drop = FALSE]
  df <- head(df, 60)
  df <- prepare_pathway_plot_df(df, plot_order)
  df$label <- trim_label(df$pathway_display_label %||% df$pathway, 55)
  df$pathway_plot_id <- make.unique(paste(df$category_id, df$pathway, sep = "__"))
  df$pathway_plot_id <- factor(df$pathway_plot_id, levels = rev(unique(df$pathway_plot_id)))
  label_values <- setNames(df$label, as.character(df$pathway_plot_id))
  fill_map <- df[!duplicated(df$category_display_name), c("category_display_name", "color"), drop = FALSE]
  fill_values <- setNames(fill_map$color, fill_map$category_display_name)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = NES, y = pathway_plot_id)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
    ggplot2::geom_point(ggplot2::aes(size = -log10(pmax(padj, 1e-300)), fill = category_display_name), shape = 21, color = "grey20", stroke = 0.2) +
    ggplot2::scale_fill_manual(values = fill_values, breaks = unique(df$category_display_name)) +
    ggplot2::scale_y_discrete(labels = label_values) +
    ggplot2::labs(x = "NES", y = NULL, size = "-log10(FDR)", fill = "LISA category") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "right",
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
  if (isTRUE(group_by_supracategory) && "macrogroup_name" %in% colnames(df)) {
    p <- p +
      ggplot2::facet_grid(macrogroup_name ~ ., scales = "free_y", space = "free_y", switch = "y") +
      ggplot2::theme(
        strip.placement = "outside",
        strip.background.y = ggplot2::element_rect(fill = "grey94", color = "grey80"),
        strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 8, face = "bold"),
        panel.spacing.y = grid::unit(0.25, "lines")
      )
  }
  p
}

prepare_pathway_plot_df <- function(df, plot_order) {
  if (!"macrogroup_order" %in% colnames(df)) df$macrogroup_order <- 99
  if (!"category_order_within_macrogroup" %in% colnames(df)) df$category_order_within_macrogroup <- 99
  if (!"macrogroup_name" %in% colnames(df)) df$macrogroup_name <- "Other or unclassified"
  if (!"category_display_name" %in% colnames(df)) df$category_display_name <- df$category_id
  if (!"color" %in% colnames(df)) df$color <- "#737373"
  df$macrogroup_order[is.na(df$macrogroup_order)] <- 99
  df$category_order_within_macrogroup[is.na(df$category_order_within_macrogroup)] <- 99
  df$macrogroup_name[is.na(df$macrogroup_name) | df$macrogroup_name == ""] <- "Other or unclassified"
  df$category_display_name[is.na(df$category_display_name) | df$category_display_name == ""] <- df$category_id[is.na(df$category_display_name) | df$category_display_name == ""]
  df$color[is.na(df$color) | df$color == ""] <- "#737373"
  df <- add_pathway_order_metadata(df)
  if (plot_order %in% c("supracategory", "fixed")) {
    df <- df[order(as.numeric(df$macrogroup_order),
                   as.numeric(df$category_order_within_macrogroup),
                   df$category_id,
                   as.numeric(df$source_family_order),
                   df$pathway_order_key,
                   df$pathway), , drop = FALSE]
  } else if (plot_order == "mean_NES") {
    df <- df[order(as.numeric(df$macrogroup_order), -abs(df$NES), df$padj, df$pathway), , drop = FALSE]
  } else if (plot_order == "n_genesets") {
    df <- df[order(as.numeric(df$macrogroup_order), df$category_id, df$padj, df$pathway), , drop = FALSE]
  } else if (plot_order == "consistency") {
    df <- df[order(as.numeric(df$macrogroup_order), df$category_id, -abs(df$NES), df$padj, df$pathway), , drop = FALSE]
  }
  df$macrogroup_name <- factor(df$macrogroup_name, levels = unique(df$macrogroup_name))
  df
}

add_pathway_order_metadata <- function(df) {
  if (!"source_family" %in% colnames(df)) df$source_family <- infer_source_family(df$pathway)
  source_order <- c(
    `GO:BP` = 1,
    `GO:MF` = 1,
    `GO:CC` = 1,
    REACTOME = 2,
    WIKIPATHWAYS = 3,
    KEGG_LEGACY = 4,
    KEGG_MEDICUS = 5,
    PID = 6,
    BIOCARTA = 7,
    CP_OTHER = 8
  )
  df$source_family_order <- unname(source_order[df$source_family])
  df$source_family_order[is.na(df$source_family_order)] <- 99
  df$pathway_order_key <- stable_pathway_order_key(df$pathway)
  df$pathway_display_label <- format_pathway_display_label(df$pathway)
  df
}

stable_pathway_order_key <- function(x) {
  key <- gsub("^REACTOME_|^GOBP_|^GOMF_|^GOCC_|^KEGG_MEDICUS_|^KEGG_|^WP_|^BIOCARTA_|^PID_", "", as.character(x))
  key <- gsub("[^A-Za-z0-9]+", " ", key)
  trimws(tolower(key))
}

save_plot_multi <- function(plot, stem, formats, width, height, bg = "white") {
  lisa_guarded_dir_create(dirname(stem))
  formats <- intersect(
    unique(tolower(formats)),
    c("pdf", "svg", "png", "tiff", "tif", "jpeg", "jpg")
  )
  for (format in formats) {
    destination <- paste0(stem, ".", format)
    lisa_guarded_write(destination, function(target) {
      if (identical(format, "pdf")) {
        ggplot2::ggsave(
          target, plot, width = width, height = height, device = "pdf",
          useDingbats = FALSE, limitsize = FALSE, bg = bg
        )
      } else if (identical(format, "svg")) {
        ggplot2::ggsave(
          target, plot, width = width, height = height,
          device = grDevices::svg, limitsize = FALSE, bg = bg
        )
      } else {
        device <- if (format %in% c("tif", "tiff")) {
          "tiff"
        } else if (format %in% c("jpg", "jpeg")) {
          "jpeg"
        } else {
          "png"
        }
        args <- list(
          filename = target, plot = plot, width = width, height = height,
          dpi = 300, device = device, limitsize = FALSE, bg = bg
        )
        if (identical(device, "tiff")) args$compression <- "lzw"
        do.call(ggplot2::ggsave, args)
      }
    })
  }
  invisible(stem)
}

add_wb_sheet <- function(wb, sheet, df) {
  sheet <- substr(sheet, 1, 31)
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet, df)
}

write_run_report <- function(path, comparison_name, input_type, dictionary, term2gene, term2gene_resolved,
                             ranking_qc, gsea_rows,
                             gsea_universe_gene_sets,
                             gsea_unclassified_gene_sets,
                             ora_rows, outputs) {
  lines <- c(
    sprintf("LISA run report: %s", comparison_name),
    sprintf("Date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    sprintf("Input type: %s", input_type),
    sprintf("Dictionary: %s", dictionary),
    sprintf("TERM2GENE: %s", term2gene),
    sprintf("MSigDB mode: %s", term2gene_resolved$msigdb_mode),
    sprintf("MSigDB db_species: %s", term2gene_resolved$db_species),
    sprintf("MSigDB target species: %s", term2gene_resolved$target_species),
    sprintf("TERM2GENE generated this run: %s", term2gene_resolved$generated),
    "",
    sprintf("Input rows: %s", ranking_qc$n_input_rows),
    sprintf("Ranked symbols: %s", ranking_qc$n_ranked_symbols),
    sprintf("LISA dictionary rows: %s", ranking_qc$n_lisa_dictionary_rows),
    sprintf("LISA gene sets: %s", ranking_qc$n_lisa_gene_sets),
    sprintf("Available TERM2GENE gene sets: %s", ranking_qc$n_available_term2gene_gene_sets),
    sprintf("Missing TERM2GENE gene sets: %s", ranking_qc$n_missing_term2gene_gene_sets),
    sprintf("Canonical KEGG gene sets available: %s", ranking_qc$has_canonical_kegg_gene_sets),
    sprintf("GSEA universe gene sets: %s", gsea_universe_gene_sets),
    sprintf("GSEA unclassified gene sets: %s", gsea_unclassified_gene_sets),
    sprintf("GSEA annotated rows: %s", gsea_rows),
    sprintf("ORA annotated rows: %s", ora_rows),
    sprintf("fgsea mode: %s", ranking_qc$fgsea_mode),
    sprintf("fgsea_nperm: %s", ifelse(is.na(ranking_qc$fgsea_nperm), "NULL (fgseaMultilevel adaptive mode)", ranking_qc$fgsea_nperm)),
    sprintf("GSEA adjusted-P cutoff: %s", ranking_qc$gsea_padj_cutoff),
    "",
    "Output directories:",
    sprintf("inputs: %s", outputs$inputs),
    sprintf("enrichment: %s", outputs$enrichment),
    sprintf("lisa_tables: %s", outputs$lisa_tables),
    sprintf("qc: %s", outputs$qc),
    sprintf("plots: %s", outputs$plots)
  )
  lisa_guarded_write(path, function(target) writeLines(lines, target))
}

infer_comparison_name <- function(input) {
  if (is.character(input) && length(input) == 1) return(tools::file_path_sans_ext(basename(input)))
  "LISA_DE_comparison"
}

safe_file_label <- function(x) {
  lisa_safe_id(as.character(x), "file label")
}

trim_label <- function(x, n = 60) {
  x <- format_pathway_display_label(x)
  ifelse(nchar(x) > n, paste0(substr(x, 1, n - 1), "..."), x)
}

format_pathway_display_label <- function(x) {
  x <- as.character(x)
  x <- gsub("_", " ", x)
  x <- gsub("^KEGG MEDICUS ", "KEGG_MEDICUS ", x)
  x
}
