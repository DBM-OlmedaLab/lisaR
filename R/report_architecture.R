# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Versioned report policy, output estimation and immutable report extensions.

lisa_output_policy_table <- function(plot_formats = "png", source_data = TRUE,
                                     recipes = FALSE) {
  formats <- unique(tolower(as.character(plot_formats)))
  formats <- formats[nzchar(formats)]
  if (any(!formats %in% c("png", "svg", "pdf"))) {
    stop("LISA-REPORT-POLICY-001 image formats must be png, svg or pdf.", call. = FALSE)
  }
  flags <- c(png = "png" %in% formats, svg = "svg" %in% formats,
             pdf = "pdf" %in% formats,
             source_data = lisa_report_bool(source_data, "report.source_data", TRUE),
             recipes = lisa_report_bool(recipes, "report.recipes", FALSE))
  data.frame(product = names(flags), requested = unname(flags),
    status = ifelse(flags, "requested", "not_requested"), stringsAsFactors = FALSE)
}

# One source of truth for the image formats of a run: the report output policy
# written by run_lisa_pipeline(). Every producer of figures for that run (the
# canonical pipeline, the FULL product core, the missing-product completion and
# the report generator) must derive its formats from here, never hard-code them.
lisa_requested_plot_formats <- function(project_dir, default = "png") {
  path <- file.path(project_dir, "report_output_policy.tsv")
  if (!file.exists(path)) return(default)
  policy <- read_lisa_tsv(path)
  if (!all(c("product", "requested") %in% names(policy))) return(default)
  requested <- tolower(as.character(policy$requested)) %in% c("true", "t", "1", "yes")
  formats <- intersect(c("png", "svg", "pdf"),
    as.character(policy$product)[requested])
  formats
}

# One representative per figure identity; PNG is only a preferred preview,
# never the identity of a product. Old callers may supply a PNG suffix pattern.
lisa_plot_files <- function(path, pattern = "[.]png$", formats = c("png", "svg", "pdf")) {
  if (!dir.exists(path) || !length(formats)) return(character())
  pattern <- sub("png\\$$", "(png|svg|pdf)$", pattern)
  paths <- sort(list.files(path, pattern = pattern, recursive = TRUE, full.names = TRUE))
  ext <- tolower(tools::file_ext(paths))
  paths <- paths[ext %in% formats]
  priority <- match(tolower(tools::file_ext(paths)), c("png", "svg", "pdf"))
  paths <- paths[order(priority, paths)]
  paths <- paths[!duplicated(tools::file_path_sans_ext(paths))]
  sort(paths)
}

# Collector file patterns derived from the requested formats: a generator that
# writes an SVG must never have it dropped by a hard-coded png|pdf pattern.
lisa_format_pattern <- function(formats = c("png", "pdf", "svg"), suffix = "") {
  paste0(suffix, "[.](", paste(formats, collapse = "|"), ")$")
}

lisa_report_contract_from_config <- function(config) {
  config_dir <- normalizePath(getwd(), mustWork = TRUE)
  cfg <- if (is.character(config) && length(config) == 1L) {
    config_path <- normalizePath(config, mustWork = TRUE)
    config_dir <- dirname(config_path)
    read_lisa_pipeline_config(config_path)
  } else config
  if (!is.list(cfg)) stop("LISA-REPORT-PLAN-001 config must be a YAML/JSON path or named list.", call. = FALSE)
  list(config = cfg, config_dir = config_dir,
    contract = lisa_validate_pipeline_config(cfg))
}

#' Estimate lisaR report outputs before execution
#'
#' @param config YAML/JSON path or equivalent named configuration list.
#' @param category_counts Optional named integer vector with the expected
#'   category count per collection. Unspecified collections use 100.
#' @param bytes_per_figure Conservative bytes per rendered figure format.
#' @return A list with summary, collection, product and output-policy tables.
#'   In `standard` mode the estimate includes the GSEA category overview
#'   families, one evidence matrix and one member-set chart per applicable LISA category, plus
#'   a linked SVG overview per individual analysis/semantic collection.
#'   It includes requested NES variants, the existing contrast overview families
#'   and paired evidence matrices for semantic contrasts. Direct HALLMARKS
#'   are not interpreted as LISA categories.
#'   With `report.category_evidence = FALSE`, legacy member-gene-set figures
#'   replace the evidence matrices.
#' @export
#'
#' @examples
#' config <- system.file("examples", "minimal-study.json", package = "lisaR")
#' if (nzchar(config)) {
#'   estimate <- plan_lisa_outputs(config, category_counts = c(`GOBP-C2` = 25))
#'   estimate$summary
#' }
plan_lisa_outputs <- function(config, category_counts = NULL,
                              bytes_per_figure = 1024^2) {
  resolved <- lisa_resolve_config(config)
  cfg <- resolved$cfg
  report <- resolved$contract$report
  if (identical(report$mode, "selected")) {
    stop("LISA-REPORT-PLAN-004 selected does not start from study inputs. Use plan_lisa_extension() with a completed standard run.", call. = FALSE)
  }
  pipeline <- lisa_config_get(cfg, "pipeline", list())
  analyses <- resolved$de_index
  contrasts <- resolved$contrast_index
  n_analyses <- nrow(analyses)
  n_contrasts <- nrow(contrasts)
  collections <- resolved$contract$collections
  registry <- lisa_collection_registry(collections)
  counts <- stats::setNames(rep.int(100L, length(collections)), collections)
  if (!is.null(category_counts)) {
    supplied <- as.integer(category_counts); names(supplied) <- names(category_counts)
    if (is.null(names(supplied)) || any(!names(supplied) %in% collections) ||
        anyNA(supplied) || any(supplied < 0L)) {
      stop("LISA-REPORT-PLAN-003 category_counts must be a non-negative named vector for configured collections.", call. = FALSE)
    }
    counts[names(supplied)] <- supplied
  }
  product_table <- if (identical(report$mode, "standard")) {
    data.frame(unit_type = "single_de", product = c("lisa_overview", "lisa_direction_overview", "gene_set_overview",
      if (!isTRUE(report$category_evidence)) "member_gene_sets"), stringsAsFactors = FALSE)
  } else {
    rbind(
      # The full extension includes KEGG gene-set figures, independently of
      # the unsupported painted-pathway-map option in the study config.
      data.frame(unit_type = "single_de", product = c("member_gene_sets", "gene_cards", "volcano", "heatmap",
        "kegg"), stringsAsFactors = FALSE),
      data.frame(unit_type = "contrast", product = c("contrast_profile", "contrast_heatmap"), stringsAsFactors = FALSE)
    )
  }
  selected_categories <- if (identical(report$mode, "standard")) 0L else sum(counts)
  standard_member_categories <- if (identical(report$mode, "standard") && !isTRUE(report$category_evidence)) {
    sum(counts[registry$analysis_collection[registry$run_category_pathway_plots]])
  } else 0L
  units <- if (identical(report$mode, "standard")) {
    n_analyses * (sum(ifelse(registry$post_lisa_profile == "hallmarks_minimal", 1L, 3L)) +
      sum(registry$run_ora & isTRUE(resolved$contract$ora$value)) + standard_member_categories)
  } else
    n_analyses * selected_categories * sum(product_table$unit_type == "single_de") +
      n_contrasts * selected_categories * sum(product_table$unit_type == "contrast")
  format_multiplier <- sum(report$formats)
  evidence_categories <- if (isTRUE(report$category_evidence)) {
    n_analyses * sum(counts[names(counts) != "HALLMARKS"])
  } else 0L
  if (evidence_categories > 0L) {
    # The member-set view stays in its category sheet, but is also rendered
    # as an ordinary downloadable figure; count it rather than hiding cost.
    units <- units + 2L * evidence_categories
    product_table <- rbind(product_table,
      data.frame(unit_type = "single_de", product = c("category_evidence", "category_member_evidence"), stringsAsFactors = FALSE))
  }
  semantic_scopes <- n_analyses * sum(collections != "HALLMARKS")
  all_contrast_scopes <- n_contrasts * sum(registry$run_contrasts)
  contrast_scopes <- n_contrasts * sum(registry$run_contrasts & registry$analysis_collection != "HALLMARKS")
  nes_variant_figures <- (semantic_scopes + 3L * contrast_scopes) * length(report$category_nes_variants)
  # Standard contrasts already render six profile views per supported scope.
  standard_contrast_figures <- if (identical(report$mode, "standard")) all_contrast_scopes * 6L else 0L
  contrast_evidence_categories <- if (isTRUE(report$category_evidence)) {
    n_contrasts * sum(counts[registry$analysis_collection[registry$run_contrasts & registry$analysis_collection != "HALLMARKS"]])
  } else 0L
  units <- units + nes_variant_figures + standard_contrast_figures + contrast_evidence_categories
  product_table <- rbind(product_table,
    data.frame(unit_type = "single_de", product = paste0("category_nes_", report$category_nes_variants)),
    if (contrast_scopes > 0L) data.frame(unit_type = "contrast", product = paste0("category_nes_", report$category_nes_variants)),
    if (standard_contrast_figures > 0L) data.frame(unit_type = "contrast", product = "contrast_profile"),
    if (contrast_evidence_categories > 0L) data.frame(unit_type = "contrast", product = "contrast_category_evidence"))
  interactive_overviews <- if (isTRUE(report$category_evidence)) semantic_scopes + contrast_scopes else 0L
  gene_explorer_pages <- as.integer(interactive_overviews > 0L)
  if (interactive_overviews > 0L) product_table <- rbind(product_table,
    data.frame(unit_type = "single_de", product = "category_navigation", stringsAsFactors = FALSE),
    data.frame(unit_type = "single_de", product = "gene_evidence", stringsAsFactors = FALSE))
  estimate <- (units * format_multiplier + interactive_overviews) * as.numeric(bytes_per_figure) +
    (if (report$source_data) units * 256 * 1024 else 0) +
    (if (report$recipes) units * 12 * 1024 else 0)
  list(
    summary = data.frame(mode = report$mode, analyses = n_analyses,
      contrasts = n_contrasts, collections = length(collections),
      selected_categories = selected_categories,
      member_gene_set_categories = standard_member_categories,
      category_evidence_categories = evidence_categories,
      category_member_evidence_categories = evidence_categories,
      contrast_evidence_categories = contrast_evidence_categories,
      category_nes_variant_figures = nes_variant_figures,
      standard_contrast_profile_figures = standard_contrast_figures,
      figure_families = nrow(product_table),
      format_multiplier = format_multiplier,
      expected_figure_files = units * format_multiplier + interactive_overviews,
      interactive_category_overviews = interactive_overviews,
      shared_gene_explorer_pages = gene_explorer_pages,
      conservative_bytes = estimate, conservative_gib = round(estimate / 1024^3, 3),
      warning = if (report$mode == "standard") {
        if (isTRUE(report$category_evidence)) {
          "Standard preserves category overview families and includes one evidence matrix per applicable LISA category and one downloadable member-set chart within each individual evidence sheet, plus linked SVG navigation and aligned contrast evidence. The unified category viewer offers selected NES variants; contrasts include all, same-direction and opposite-direction subsets. The shared gene explorer size depends on full membership and DE cardinality, not this figure estimate."
        } else "Standard includes one member-gene-set figure per applicable LISA category; review the estimated figure count and storage before execution."
      } else {
        "Selected/full expansion may create thousands of figures and use gigabytes or tens of gigabytes."
      },
      stringsAsFactors = FALSE),
    collections = data.frame(collection = collections,
      expected_categories = unname(counts), stringsAsFactors = FALSE),
    products = product_table,
    policy = lisa_output_policy_table(names(report$formats)[report$formats],
      report$source_data, report$recipes)
  )
}

lisa_extension_read_selection <- function(selection) {
  if (is.character(selection) && length(selection) == 1L && identical(tolower(selection), "full")) {
    return(list(mode = "full", filters = list(), report = list(mode = "full")))
  }
  object <- if (is.character(selection) && length(selection) == 1L) {
    read_lisa_pipeline_config(normalizePath(selection, mustWork = TRUE))
  } else selection
  if (!is.list(object)) stop("LISA-EXTENSION-003 selection must be 'full', a YAML/JSON path or a named list.", call. = FALSE)
  filters <- lisa_config_get(object, "selection", object)
  report <- lisa_config_get(object, "report", list())
  if (!is.list(filters)) stop("LISA-EXTENSION-004 selection must be a named object.", call. = FALSE)
  filters$report <- NULL; filters$from_run <- NULL
  mode <- tolower(as.character(lisa_config_get(filters, "mode", "selected")[[1]]))
  filters$mode <- NULL
  if (!mode %in% c("selected", "full")) stop("LISA-EXTENSION-005 extension mode must be selected or full.", call. = FALSE)
  allowed <- c("analyses", "contrasts", "collections", "categories", "supercategories",
    "products", "entities", "variants")
  unknown <- setdiff(names(filters), allowed)
  if (length(unknown)) stop("LISA-EXTENSION-006 unknown selection field(s): ", paste(unknown, collapse = ", "), call. = FALSE)
  filters <- lapply(filters, function(value) unique(as.character(unlist(value, use.names = FALSE))))
  if (is.null(report$mode)) report$mode <- mode
  list(mode = mode, filters = filters, report = report)
}

lisa_extension_nearest <- function(value, valid) {
  valid <- sort(unique(as.character(valid)))
  if (!length(valid)) return("<none>")
  d <- utils::adist(value, valid, ignore.case = TRUE)
  paste(utils::head(valid[order(d[1, ], valid)], 3L), collapse = ", ")
}

lisa_extension_find_one <- function(root, pattern) {
  hits <- sort(list.files(root, pattern = pattern, full.names = TRUE))
  if (length(hits) == 1L) hits[[1]] else ""
}

lisa_extension_gsea_padj_cutoff <- function(source_run) {
  contract_path <- file.path(source_run, "contract_manifest.tsv")
  if (file.exists(contract_path)) {
    contract <- read_lisa_tsv(contract_path)
    if (all(c("key", "value") %in% names(contract))) {
      hit <- contract$value[as.character(contract$key) == "gsea_padj_cutoff"]
      if (length(hit) == 1L) {
        return(lisa_validate_gsea_padj_cutoff(
          suppressWarnings(as.numeric(hit)),
          field = "source_run.contract_manifest.gsea_padj_cutoff"
        ))
      }
    }
  }
  tables <- c(
    list.files(file.path(source_run, "outputs", "single_de"),
      pattern = "_GSEA_category_summary[.]tsv$", recursive = TRUE, full.names = TRUE),
    list.files(file.path(source_run, "outputs", "category_contrasts"),
      pattern = "_GSEA_category_contrast[.]tsv$", recursive = TRUE, full.names = TRUE)
  )
  observed <- unique(unlist(lapply(tables, function(path) {
    x <- read_lisa_tsv(path)
    if (!"gsea_padj_cutoff" %in% names(x)) return(numeric())
    suppressWarnings(as.numeric(x$gsea_padj_cutoff))
  }), use.names = FALSE))
  observed <- observed[is.finite(observed)]
  if (length(observed) == 1L) {
    return(lisa_validate_gsea_padj_cutoff(
      observed,
      field = "source_run.gsea_padj_cutoff"
    ))
  }
  if (!length(observed)) return(0.25)
  stop("LISA-EXTENSION-027 source run contains inconsistent gsea_padj_cutoff values.", call. = FALSE)
}

lisa_extension_relative <- function(path, root) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root, "/")
  if (!startsWith(path, prefix)) stop("LISA-EXTENSION-007 source input escaped the canonical run.", call. = FALSE)
  substring(path, nchar(prefix) + 1L)
}

lisa_extension_matrix_columns <- function() {
  c(
    "expression_matrix", "expression_matrix_path", "counts_matrix",
    "counts_matrix_path", "vst_matrix", "vst_matrix_path",
    "normalized_matrix", "normalized_matrix_path", "tpm_matrix",
    "tpm_matrix_path", "sample_counts_path", "sample_counts", "matrix_path"
  )
}

# Resolve the matrix selected by the saved de_index through the run's own
# effective-input receipt. The receipt may retain an absolute build-host path;
# only its canonical `effective_inputs/<file>` suffix is used, and that suffix
# must identify exactly one manifest row with the same recorded and actual hash.
lisa_extension_effective_matrix <- function(source_run, analysis_row) {
  if (!nrow(analysis_row)) return(NULL)
  if (nrow(analysis_row) != 1L) {
    stop("LISA-EXTENSION-036 matrix selection requires exactly one de_index row.",
         call. = FALSE)
  }
  columns <- intersect(lisa_extension_matrix_columns(), names(analysis_row))
  selected <- columns[vapply(columns, function(column) {
    value <- as.character(analysis_row[[column]][[1L]])
    !is.na(value) && nzchar(trimws(value))
  }, logical(1L))]
  if (!length(selected)) return(NULL)
  source_column <- selected[[1L]]

  manifest_path <- file.path(source_run, "run_manifest.tsv")
  if (isTRUE(getOption("lisaR.full_products_in_run", FALSE)) &&
      !file.exists(manifest_path)) {
    # run_lisa(mode = "full") renders inside the unpromoted staging tree, which
    # has no manifest yet. The effective-input receipt written by this same
    # run binds the matrix bytes; the one manifest seals them at promotion.
    return(lisa_extension_effective_matrix_in_run(source_run, analysis_row,
      source_column))
  }
  manifest <- read_lisa_tsv(manifest_path)
  if (!all(c("path", "sha256") %in% names(manifest))) {
    stop("LISA-EXTENSION-036 the canonical manifest cannot bind an effective matrix.",
         call. = FALSE)
  }
  receipt_relative <- "effective_inputs.tsv"
  receipt_index <- which(as.character(manifest$path) == receipt_relative)
  if (length(receipt_index) != 1L) {
    stop("LISA-EXTENSION-036 effective_inputs.tsv is not uniquely declared in the canonical manifest.",
         call. = FALSE)
  }
  receipt_path <- file.path(source_run, receipt_relative)
  if (!file.exists(receipt_path) || dir.exists(receipt_path) ||
      !identical(lisa_extension_relative(receipt_path, source_run),
                 receipt_relative) ||
      !identical(lisa_sha256_file(receipt_path),
                 as.character(manifest$sha256[[receipt_index]]))) {
    stop("LISA-EXTENSION-036 the effective-input receipt is missing, outside, or changed.",
         call. = FALSE)
  }
  receipt <- read_lisa_tsv(receipt_path)
  required <- c("analysis_id", "input_class", "effective_path", "sha256")
  if (!all(required %in% names(receipt))) {
    stop("LISA-EXTENSION-036 effective_inputs.tsv lacks the matrix receipt schema.",
         call. = FALSE)
  }
  analysis <- as.character(analysis_row$analysis_id[[1L]])
  matches <- which(as.character(receipt$analysis_id) == analysis &
                   as.character(receipt$input_class) == source_column)
  if (length(matches) != 1L) {
    stop("LISA-EXTENSION-036 the selected matrix has no unique effective-input receipt.",
         call. = FALSE)
  }
  effective_path <- gsub("\\\\", "/",
                         as.character(receipt$effective_path[[matches]]))
  parts <- strsplit(effective_path, "/", fixed = TRUE)[[1L]]
  marker <- which(parts == "effective_inputs")
  if (length(marker) != 1L || marker[[1L]] == length(parts) ||
      length(parts) - marker[[1L]] != 1L) {
    stop("LISA-EXTENSION-036 the selected matrix receipt has no contained effective_inputs path.",
         call. = FALSE)
  }
  relative <- paste(parts[marker[[1L]]:length(parts)], collapse = "/")
  expected_name <- paste0(analysis, "_", source_column, ".tsv")
  if (!identical(basename(relative), expected_name)) {
    stop("LISA-EXTENSION-036 the selected matrix receipt does not name its original effective input.",
         call. = FALSE)
  }
  matrix_index <- which(as.character(manifest$path) == relative)
  if (length(matrix_index) != 1L) {
    stop("LISA-EXTENSION-036 the selected matrix is not uniquely declared in the canonical manifest.",
         call. = FALSE)
  }
  recorded <- as.character(receipt$sha256[[matches]])
  declared <- as.character(manifest$sha256[[matrix_index]])
  if (!lisa_sha256_all_valid(c(recorded, declared), 2L) ||
      !identical(recorded, declared)) {
    stop("LISA-EXTENSION-036 the selected matrix receipt and manifest disagree.",
         call. = FALSE)
  }
  matrix_path <- file.path(source_run, relative)
  if (!file.exists(matrix_path) || dir.exists(matrix_path) ||
      !identical(lisa_extension_relative(matrix_path, source_run), relative) ||
      !identical(lisa_sha256_file(matrix_path), declared)) {
    stop("LISA-EXTENSION-036 the selected effective matrix is missing, outside, or changed.",
         call. = FALSE)
  }
  list(path = relative, source_column = source_column)
}

lisa_extension_effective_matrix_in_run <- function(source_run, analysis_row,
                                                   source_column) {
  receipt <- read_lisa_tsv(file.path(source_run, "effective_inputs.tsv"))
  analysis <- as.character(analysis_row$analysis_id[[1L]])
  matches <- which(as.character(receipt$analysis_id) == analysis &
                   as.character(receipt$input_class) == source_column)
  if (length(matches) != 1L) {
    stop("LISA-EXTENSION-036 the selected matrix has no unique effective-input receipt.",
         call. = FALSE)
  }
  parts <- strsplit(gsub("\\\\", "/", as.character(receipt$effective_path[[matches]])),
    "/", fixed = TRUE)[[1L]]
  marker <- which(parts == "effective_inputs")
  if (length(marker) != 1L || length(parts) - marker[[1L]] != 1L) {
    stop("LISA-EXTENSION-036 the selected matrix receipt has no contained effective_inputs path.",
         call. = FALSE)
  }
  relative <- paste(parts[marker[[1L]]:length(parts)], collapse = "/")
  matrix_path <- file.path(source_run, relative)
  if (!file.exists(matrix_path) ||
      !identical(lisa_sha256_file(matrix_path), as.character(receipt$sha256[[matches]]))) {
    stop("LISA-EXTENSION-036 the selected effective matrix is missing or changed.",
         call. = FALSE)
  }
  list(path = relative, source_column = source_column)
}

lisa_extension_exact_catalog <- function(catalog, exact_products) {
  if (is.null(exact_products)) return(catalog)
  if (!is.list(exact_products)) {
    stop("LISA-EXTENSION-026 exact_products must be a named list.", call. = FALSE)
  }
  variants <- exact_products$heatmap_variants
  if (length(variants)) {
    variants <- unique(as.character(variants))
    is_heatmap <- catalog$product == "heatmap"
    if (any(is_heatmap)) {
      base_rows <- catalog[is_heatmap, , drop = FALSE]
      expanded <- do.call(rbind, lapply(variants, function(scale) {
        rows <- base_rows
        rows$variant <- scale
        rows
      }))
      catalog <- rbind(catalog[!is_heatmap, , drop = FALSE], expanded)
    }
  }
  maps <- exact_products$kegg_maps
  if (is.data.frame(maps) && nrow(maps)) {
    required <- c("analysis_id", "collection", "category_id", "kegg_id")
    missing <- setdiff(required, names(maps))
    if (length(missing)) {
      stop("LISA-EXTENSION-027 the prepared KEGG pathway index is missing ",
        "column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
    }
    # Every native map row inherits the scientific input paths of its own
    # analysis/collection, so `lisa_extension_source_contract()` still verifies
    # them against the canonical manifest. A pathway whose owning collection is
    # not in the catalogue is dropped rather than invented.
    owners <- unique(catalog[catalog$unit_type == "single_de",
      c("analysis_id", "collection", "category_id", "supercategory_id",
        "summary_path", "gsea_path", "de_path"), drop = FALSE])
    joined <- merge(maps[, unique(c(required, intersect(
      c("kegg_title", "rank", "kegg_cache_root", "kegg_snapshot_id",
        "kegg_species", "max_abs_log2fc", "color_power",
        "kegg_resource_digest"), names(maps)))),
      drop = FALSE], owners,
      by = c("analysis_id", "collection", "category_id"))
    if (nrow(joined)) {
      rows <- data.frame(
        unit_type = "single_de", analysis_id = joined$analysis_id,
        contrast_id = "", collection = joined$collection,
        category_id = joined$category_id,
        supercategory_id = joined$supercategory_id,
        product = "kegg_pathway_map",
        summary_path = joined$summary_path, gsea_path = joined$gsea_path,
        de_path = joined$de_path, contrast_path = "",
        entity = as.character(joined$kegg_id), variant = "",
        stringsAsFactors = FALSE)
      for (column in c("kegg_cache_root", "kegg_snapshot_id", "kegg_species",
                       "max_abs_log2fc", "color_power", "kegg_resource_digest",
                       "kegg_title", "rank")) {
        value <- if (column %in% names(joined)) as.character(joined[[column]]) else ""
        rows[[column]] <- value
        if (!column %in% names(catalog)) catalog[[column]] <- ""
      }
      for (column in setdiff(names(catalog), names(rows))) rows[[column]] <- ""
      rows <- rows[, names(catalog), drop = FALSE]
      catalog <- rbind(catalog, rows)
    }
  }
  catalog <- lisa_extension_native_catalog(catalog, exact_products)
  catalog
}

lisa_extension_native_catalog <- function(catalog, exact_products) {
  if (is.null(exact_products)) return(catalog)
  blank_row <- function(rows) {
    for (column in setdiff(names(catalog), names(rows))) rows[[column]] <- ""
    rows[, names(catalog), drop = FALSE]
  }
  ensure <- function(columns) {
    for (column in columns) if (!column %in% names(catalog)) catalog[[column]] <<- ""
  }
  ensure(c("preparation_digest", "prepared_dir", "de_a_path", "de_b_path"))

  if (isTRUE(exact_products$native_single_de)) {
    # Derived from the rows that already established gene-level applicability,
    # so a collection whose evidence cannot support gene-level products offers
    # no recurrent screen either.
    owners <- unique(catalog[catalog$unit_type == "single_de" &
      catalog$product == "gene_cards",
      c("analysis_id", "collection", "summary_path", "gsea_path", "de_path"),
      drop = FALSE])
    if (nrow(owners)) {
      rows <- data.frame(
        unit_type = "single_de", analysis_id = owners$analysis_id,
        contrast_id = "", collection = owners$collection,
        # No category, and no placeholder standing in for one.
        category_id = "", supercategory_id = "",
        product = "de_recurrent_genes",
        summary_path = owners$summary_path, gsea_path = owners$gsea_path,
        de_path = owners$de_path, contrast_path = "",
        entity = "", variant = "", stringsAsFactors = FALSE)
      catalog <- rbind(catalog, blank_row(rows))
    }
  }

  prepared <- exact_products$contrast_products
  if (!is.null(prepared) && is.data.frame(prepared$units) && nrow(prepared$units)) {
    # Every prepared bundle must name a contrast/collection the source run
    # really offers, so a bundle left over from another run cannot invent a row.
    owners <- unique(catalog[catalog$unit_type == "contrast",
      c("contrast_id", "collection", "contrast_path"), drop = FALSE])
    endpoints <- unique(catalog[catalog$unit_type == "single_de",
      c("analysis_id", "collection", "de_path"), drop = FALSE])
    units <- merge(prepared$units, owners, by = c("contrast_id", "collection"))
    contrast_row <- function(unit, product, category_id = "", entity = "",
                             variant = "", extra = NULL) {
      rows <- data.frame(
        unit_type = "contrast", analysis_id = "",
        contrast_id = unit$contrast_id, collection = unit$collection,
        category_id = category_id, supercategory_id = "", product = product,
        summary_path = "", gsea_path = "", de_path = "",
        contrast_path = unit$contrast_path, entity = entity, variant = variant,
        preparation_digest = unit$preparation_digest,
        prepared_dir = unit$prepared_dir,
        # `config/contrast_index.tsv` is deliberately NOT bound as a row input.
        # `render_lisa_categories()` already stages both config indexes into the
        # work root before any unit renders, and `lisa_extension_stage_files()`
        # copies with `overwrite = FALSE`; naming it here made every contrast
        # render fail with LISA-EXTENSION-016 on a file that was already there.
        # Preparation still verifies it against the run manifest.
        de_a_path = "", de_b_path = "", stringsAsFactors = FALSE)
      if (!is.null(extra)) for (column in names(extra)) rows[[column]] <- extra[[column]]
      rows
    }
    new_rows <- list()
    for (index in seq_len(nrow(units))) {
      unit <- units[index, , drop = FALSE]
      new_rows[[length(new_rows) + 1L]] <- contrast_row(unit,
        "contrast_paired_heatmap")
      new_rows[[length(new_rows) + 1L]] <- contrast_row(unit,
        "contrast_gene_category_network", variant = "all_categories_top7")
      if (is.data.frame(prepared$categories)) {
        own <- prepared$categories[
          prepared$categories$contrast_id == unit$contrast_id[[1L]] &
          prepared$categories$collection == unit$collection[[1L]], , drop = FALSE]
        if (nrow(own)) {
          new_rows[[length(new_rows) + 1L]] <- contrast_row(unit,
            "contrast_gene_card", category_id = as.character(own$category_id))
        }
      }
      if (is.data.frame(prepared$maps)) {
        own <- prepared$maps[
          prepared$maps$contrast_id == unit$contrast_id[[1L]] &
          prepared$maps$collection == unit$collection[[1L]], , drop = FALSE]
        # A contrast map needs BOTH endpoints' standardized DE tables to paint
        # its A/B node halves, so they are bound as inputs of the row itself.
        sides <- lisa_extension_contrast_endpoint_paths(unit, endpoints,
                                                        prepared$units)
        if (nrow(own) && !is.null(sides)) {
          settings <- intersect(c("kegg_cache_root", "kegg_snapshot_id",
            "kegg_species", "max_abs_log2fc", "color_power",
            "kegg_resource_digest", "kegg_title", "rank"), names(own))
          extra <- c(list(de_a_path = sides$a, de_b_path = sides$b),
                     lapply(stats::setNames(settings, settings),
                            function(column) as.character(own[[column]])))
          ensure(settings)
          new_rows[[length(new_rows) + 1L]] <- contrast_row(unit,
            "contrast_kegg_map", entity = as.character(own$kegg_id),
            extra = extra)
        }
      }
    }
    if (length(new_rows)) {
      catalog <- rbind(catalog,
                       do.call(rbind, lapply(new_rows, blank_row)))
    }
  }
  catalog
}

# The two endpoints' standardized DE paths for one prepared contrast bundle.
# Read from the bundle's own receipt rather than re-derived, so a row can never
# bind an endpoint the preparation did not actually pair.
lisa_extension_contrast_endpoint_paths <- function(unit, endpoints, units) {
  row <- units[units$contrast_id == unit$contrast_id[[1L]] &
               units$collection == unit$collection[[1L]], , drop = FALSE]
  if (!nrow(row) || !all(c("analysis_a", "analysis_b") %in% names(row))) return(NULL)
  pick <- function(analysis) {
    hit <- endpoints[endpoints$analysis_id == as.character(analysis) &
                     endpoints$collection == unit$collection[[1L]], , drop = FALSE]
    if (nrow(hit) != 1L) NULL else as.character(hit$de_path[[1L]])
  }
  a <- pick(row$analysis_a[[1L]]); b <- pick(row$analysis_b[[1L]])
  if (is.null(a) || is.null(b)) return(NULL)
  list(a = a, b = b)
}

lisa_extension_discover_catalog <- function(source_run, exact_products = NULL) {
  gsea_padj_cutoff <- lisa_extension_gsea_padj_cutoff(source_run)
  root <- file.path(source_run, "outputs", "single_de")
  summaries <- if (dir.exists(root)) list.files(root,
    pattern = "_GSEA_category_summary[.]tsv$", recursive = TRUE, full.names = TRUE) else character()
  de_index <- read_lisa_tsv(file.path(source_run, "config", "de_index.tsv"))
  rows <- list()
  for (summary_path in sort(summaries)) {
    rel <- lisa_extension_relative(summary_path, source_run)
    match <- regexec("^outputs/single_de/([^/]+)/collection_([^/]+)/lisa_tables/", rel)
    parts <- regmatches(rel, match)[[1]]
    if (length(parts) != 3L) next
    analysis <- parts[[2]]; collection <- parts[[3]]
    collection_root <- dirname(dirname(summary_path))
    gsea <- lisa_extension_find_one(file.path(collection_root, "enrichment"),
      paste0("^", analysis, "_GSEA_.*_annotated[.]tsv$"))
    de <- file.path(collection_root, "inputs", paste0(analysis, "_standardized_DE.tsv"))
    if (!nzchar(gsea) || !file.exists(de)) next
    summary <- read_lisa_tsv(summary_path)
    if (!all(c("category_id", "macrogroup_id") %in% names(summary))) next
    if ("n_genesets" %in% names(summary)) summary <- summary[as.numeric(summary$n_genesets) > 0, , drop = FALSE]
    if (!nrow(summary)) next
    analysis_row <- de_index[as.character(de_index$analysis_id) == analysis, , drop = FALSE]
    matrix <- lisa_extension_effective_matrix(source_run, analysis_row)
    has_matrix <- !is.null(matrix)
    gsea_columns <- names(read_lisa_tsv(gsea))
    de_columns <- names(read_lisa_tsv(de))
    gene_level_applicable <- "leadingEdge" %in% gsea_columns &&
      "symbol" %in% de_columns && "log2FoldChange" %in% de_columns &&
      any(c("padj", "pvalue") %in% de_columns)
    products <- c("member_gene_sets",
      if (gene_level_applicable) c("gene_cards", "volcano"),
      if (gene_level_applicable && has_matrix) "heatmap")
    kegg_categories <- character()
    gsea_table <- read_lisa_tsv(gsea)
    gsea_table <- lisa_classified_gsea_rows(gsea_table)
    if ("pathway" %in% names(gsea_table)) {
      is_kegg <- grepl("^KEGG_", as.character(gsea_table$pathway))
      if ("padj" %in% names(gsea_table)) is_kegg <- is_kegg &
        !is.na(as.numeric(gsea_table$padj)) & as.numeric(gsea_table$padj) <= gsea_padj_cutoff
      kegg_categories <- unique(as.character(gsea_table$category_id[is_kegg]))
    }
    for (product in products) {
      rows[[length(rows) + 1L]] <- data.frame(
        unit_type = "single_de", analysis_id = analysis, contrast_id = "",
        collection = collection,
        category_id = as.character(summary$category_id),
        supercategory_id = as.character(summary$macrogroup_id),
        product = product,
        summary_path = lisa_extension_relative(summary_path, source_run),
        gsea_path = lisa_extension_relative(gsea, source_run),
        de_path = lisa_extension_relative(de, source_run),
        contrast_path = "",
        matrix_path = if (identical(product, "heatmap")) matrix$path else "",
        matrix_source_column = if (identical(product, "heatmap"))
          matrix$source_column else "",
        stringsAsFactors = FALSE)
    }
    if (length(kegg_categories)) {
      kegg_summary <- summary[as.character(summary$category_id) %in% kegg_categories, , drop = FALSE]
      rows[[length(rows) + 1L]] <- data.frame(
        unit_type = "single_de", analysis_id = analysis, contrast_id = "",
        collection = collection, category_id = as.character(kegg_summary$category_id),
        supercategory_id = as.character(kegg_summary$macrogroup_id), product = "kegg",
        summary_path = lisa_extension_relative(summary_path, source_run),
        gsea_path = lisa_extension_relative(gsea, source_run),
        de_path = lisa_extension_relative(de, source_run), contrast_path = "",
        matrix_path = "", matrix_source_column = "",
        stringsAsFactors = FALSE)
    }
  }

  contrast_root <- file.path(source_run, "outputs", "category_contrasts")
  contrast_tables <- if (dir.exists(contrast_root)) list.files(contrast_root,
    pattern = "_GSEA_category_contrast[.]tsv$", recursive = TRUE, full.names = TRUE) else character()
  for (contrast_path in sort(contrast_tables)) {
    rel <- lisa_extension_relative(contrast_path, source_run)
    match <- regexec("^outputs/category_contrasts/([^/]+)/collection_([^/]+)/lisa_tables/", rel)
    parts <- regmatches(rel, match)[[1]]
    if (length(parts) != 3L) next
    contrast_id <- parts[[2]]; collection <- parts[[3]]
    contrast <- read_lisa_tsv(contrast_path)
    required <- c("category_id", "mean_NES_A", "mean_NES_B")
    if (!all(required %in% names(contrast)) || !nrow(contrast)) next
    # Keep the complete category map. A category with no significant member set
    # on either side is a deliberate blank contrast row, not an absent product.
    macro <- if ("macrogroup_id" %in% names(contrast)) as.character(contrast$macrogroup_id) else "CONTRAST"
    for (product in c("contrast_profile", "contrast_heatmap")) {
      rows[[length(rows) + 1L]] <- data.frame(
        unit_type = "contrast", analysis_id = "", contrast_id = contrast_id,
        collection = collection, category_id = as.character(contrast$category_id),
        supercategory_id = macro, product = product,
        summary_path = "", gsea_path = "", de_path = "",
        contrast_path = rel, matrix_path = "", matrix_source_column = "",
        stringsAsFactors = FALSE)
    }
  }
  if (!length(rows)) stop("LISA-EXTENSION-009 no renderable category units were discovered.", call. = FALSE)
  catalog <- do.call(rbind, rows)
  catalog$entity <- ""
  catalog$variant <- ""
  catalog <- lisa_extension_exact_catalog(catalog, exact_products)
  catalog <- unique(catalog)
  catalog[order(catalog$unit_type, catalog$analysis_id, catalog$contrast_id,
    catalog$collection, catalog$supercategory_id,
    catalog$category_id, catalog$product, catalog$entity, catalog$variant),
    , drop = FALSE]
}

lisa_extension_source_contract <- function(source_run, catalog) {
  source_run <- lisa_assert_run_tree_safe(source_run)
  check <- verify_run(source_run)
  if (!identical(check$gate, "PASS")) stop("LISA-EXTENSION-011 source run validation failed: ",
    paste(check$findings, collapse = "; "), call. = FALSE)
  manifest_path <- file.path(source_run, "run_manifest.tsv")
  manifest <- read_lisa_tsv(manifest_path)
  path_columns <- grep("_path$", names(catalog), value = TRUE)
  inputs <- sort(unique(unlist(catalog[path_columns], use.names = FALSE)))
  inputs <- inputs[nzchar(inputs)]
  matched <- match(inputs, manifest$path)
  if (anyNA(matched)) stop("LISA-EXTENSION-012 scientific input is absent from the canonical manifest: ",
    paste(inputs[is.na(matched)], collapse = ", "), call. = FALSE)
  actual <- vapply(file.path(source_run, inputs), lisa_sha256_file, character(1))
  if (!identical(unname(actual), unname(as.character(manifest$sha256[matched])))) {
    stop("LISA-EXTENSION-013 scientific input hash differs from the canonical manifest.", call. = FALSE)
  }
  list(root = source_run, manifest_hash = lisa_sha256_file(manifest_path), inputs = inputs)
}

lisa_extension_filter <- function(catalog, selection) {
  selected <- catalog
  if (selection$mode == "full") {
    # Prioritized-gene cards are part of FULL, alongside category evidence.
    return(selected)
  }
  mapping <- c(analyses = "analysis_id", contrasts = "contrast_id", collections = "collection",
    categories = "category_id", supercategories = "supercategory_id", products = "product",
    entities = "entity", variants = "variant")
  owner_filters <- intersect(c("analyses", "contrasts"), names(selection$filters))
  if (length(owner_filters) == 2L) {
    owner_selected <- rep(FALSE, nrow(selected))
    for (name in owner_filters) {
      wanted <- selection$filters[[name]]
      field <- mapping[[name]]
      valid <- unique(catalog[[field]])
      valid <- valid[nzchar(valid)]
      if (length(wanted) && !(length(wanted) == 1L && tolower(wanted) == "all")) {
        unknown <- setdiff(wanted, valid)
        if (length(unknown)) {
          suggestions <- vapply(unknown, lisa_extension_nearest, character(1), valid = valid)
          stop(sprintf("LISA-EXTENSION-014 unknown %s selection: %s. Nearest valid ID(s): %s.",
            name, paste(unknown, collapse = ", "), paste(suggestions, collapse = "; ")), call. = FALSE)
        }
      }
      owner_kind <- if (identical(name, "analyses")) "single_de" else "contrast"
      keep <- selected$unit_type == owner_kind
      if (length(wanted) && !(length(wanted) == 1L && tolower(wanted) == "all")) {
        keep <- keep & selected[[field]] %in% wanted
      }
      owner_selected <- owner_selected | keep
    }
    selected <- selected[owner_selected, , drop = FALSE]
  }
  remaining_filters <- setdiff(names(selection$filters), owner_filters)
  if (length(owner_filters) < 2L) remaining_filters <- names(selection$filters)
  for (name in remaining_filters) {
    wanted <- selection$filters[[name]]
    if (!length(wanted) || (length(wanted) == 1L && tolower(wanted) == "all")) next
    field <- mapping[[name]]
    unknown <- setdiff(wanted, unique(catalog[[field]]))
    if (length(unknown)) {
      suggestions <- vapply(unknown, lisa_extension_nearest, character(1), valid = catalog[[field]])
      stop(sprintf("LISA-EXTENSION-014 unknown %s selection: %s. Nearest valid ID(s): %s.",
        name, paste(unknown, collapse = ", "), paste(suggestions, collapse = "; ")), call. = FALSE)
    }
    selected <- selected[selected[[field]] %in% wanted, , drop = FALSE]
  }
  if (!nrow(selected)) stop("LISA-EXTENSION-015 selection resolved to zero renderable units.", call. = FALSE)
  selected
}

#' Plan additional selected or full report figures
#'
#' @param exact_products Optional named list of exact heatmap variants, KEGG
#'   maps or native products used by the exploration interface. Leave NULL
#'   for the standard selected or full report.
#' @param source_run A completed and checked lisaR run.
#' @param selection `"full"`, a YAML/JSON selection path or named list.
#' @return A plan listing the requested figures and expected files.
#' @export
#'
#' @examples
#' # Planning an extension requires a verified scientific run. The installed
#' # sample project is offline; the prerequisite run is kept interactive.
#' if (interactive()) {
#'   project <- lisa_init_project(tempfile("lisa-extension-project-"))
#'   result <- run_lisa(file.path(project, "study.yml"))
#'   selection <- list(
#'     selection = list(
#'       mode = "selected",
#'       categories = "SYN_SIGNAL",
#'       products = "member_gene_sets"
#'     ),
#'     report = list(
#'       mode = "selected",
#'       formats = list(png = TRUE, svg = FALSE, pdf = FALSE),
#'       source_data = TRUE,
#'       recipes = FALSE
#'     )
#'   )
#'   extension_plan <- plan_lisa_extension(result$output_dir, selection)
#'   extension_plan$expected_files
#' }
plan_lisa_extension <- function(source_run, selection, exact_products = NULL) {
  source_run <- lisa_assert_run_tree_safe(source_run)
  parsed <- lisa_extension_read_selection(selection)
  catalog <- lisa_extension_discover_catalog(source_run, exact_products)
  source <- lisa_extension_source_contract(source_run, catalog)
  gsea_padj_cutoff <- lisa_extension_gsea_padj_cutoff(source_run)
  selected <- lisa_extension_filter(catalog, parsed)
  report <- lisa_validate_report_config(list(report = parsed$report), allow_selected = TRUE)
  policy <- lisa_output_policy_table(names(report$formats)[report$formats], report$source_data, report$recipes)
  # The exact-selector dimensions join the fingerprint only when some selected
  # row actually carries one. Appending them unconditionally would change
  # `extension_id` for every existing FULL selection -- and `extension_id` is
  # written into the FULL receipt -- so the condition is what keeps the default
  # command-line behaviour byte-identical while still separating two variants of
  # the same figure.
  exact <- paste(selected$entity, selected$variant, sep = ":")
  exact_fingerprint <- if (any(nzchar(selected$entity) | nzchar(selected$variant)))
    exact else character()
  fingerprint <- paste(c(source$manifest_hash, format(gsea_padj_cutoff, digits = 17),
    apply(selected[c("unit_type", "analysis_id", "contrast_id", "collection", "category_id", "product")], 1L, paste, collapse = ":"),
    exact_fingerprint,
    policy$product, policy$status), collapse = "\n")
  extension_id <- paste0("lisa-extension-", substr(lisa_sha256_text(fingerprint), 1L, 16L))
  per_unit <- sum(report$formats) + as.integer(report$source_data) + as.integer(report$recipes)
  list(extension_id = extension_id, mode = parsed$mode, source_run = source$root,
    source_manifest_hash = source$manifest_hash,
    gsea_padj_cutoff = gsea_padj_cutoff,
    selected_categories = sort(unique(selected$category_id)), units = selected,
    policy = policy, report = report,
    expected_files = nrow(selected) * per_unit,
    warning = if (parsed$mode == "full")
      "Full expansion may create thousands of figures and use gigabytes or tens of gigabytes." else "")
}

lisa_extension_stage_files <- function(source_run, work_root, relative_paths) {
  managed_root <- getOption("lisaR.run_root", NULL)
  if (is.null(managed_root)) {
    stop("LISA-EXTENSION-016 staging requires an authorized run root.",
         call. = FALSE)
  }
  relative_paths <- sort(unique(as.character(relative_paths)))
  for (rel in relative_paths) {
    source <- file.path(source_run, rel)
    if (!file.exists(source) || dir.exists(source)) {
      stop("LISA-EXTENSION-016 canonical staging input is missing: ", rel, call. = FALSE)
    }
    # Re-resolve every path against the immutable source root. This deliberately
    # stages only manifest-bound scientific tables, never an entire collection.
    if (!identical(lisa_extension_relative(source, source_run), rel)) {
      stop("LISA-EXTENSION-016 canonical staging input escaped the source run: ", rel, call. = FALSE)
    }
    target <- file.path(work_root, rel)
    # Several collections can reuse one saved matrix. Reuse only the exact
    # already-staged bytes, never overwrite a conflicting file or follow links.
    if (lisa_path_entry_exists(target)) {
      lisa_assert_regular_managed_file(target, managed_root)
      if (!identical(lisa_sha256_file(source), lisa_sha256_file(target)))
        stop("LISA-EXTENSION-016 staged input differs from its canonical source: ", rel,
             call. = FALSE)
      next
    }
    copy_error <- ""
    copied <- tryCatch({
      lisa_guarded_copy(source, target, overwrite = FALSE,
                        run_root = managed_root)
      TRUE
    }, error = function(error) {
      if (inherits(error, "lisa_filesystem_recovery_error")) stop(error)
      copy_error <<- conditionMessage(error)
      FALSE
    })
    if (!copied) {
      stop("LISA-EXTENSION-016 failed to stage canonical input: ", rel, ": ", copy_error, call. = FALSE)
    }
  }
  invisible(relative_paths)
}

lisa_extension_file_label_prefix <- function(gsea_path, analysis_id) {
  name <- basename(gsea_path)
  sub(paste0("^", analysis_id, "_GSEA_(.*)_annotated[.]tsv$"), "\\1", name)
}

lisa_extension_collect <- function(source_dir, destination, category_ids) {
  if (!dir.exists(source_dir)) return(invisible(FALSE))
  managed_root <- getOption("lisaR.run_root", NULL)
  if (is.null(managed_root)) {
    stop("LISA-EXTENSION-017 collection requires an authorized run root.",
         call. = FALSE)
  }
  tree <- lisa_scan_run_tree(source_dir)
  files <- sort(tree$path[!tree$isdir])
  # A category attachment is declared by its own source/matrix sidecar, not by
  # a category-looking filename.  This keeps an ID collision from silently
  # moving a figure to the wrong card.  Rendering occurs before the output
  # policy removes optional source data, so the contract is always available
  # at collection time.
  sidecars <- files[grepl("_(source|matrix)[.]tsv$", files, ignore.case = TRUE)]
  families <- vapply(sidecars, lisa_extension_figure_family, character(1L))
  attached_families <- vapply(seq_along(sidecars), function(i) {
    tab <- tryCatch(read_lisa_tsv(sidecars[[i]]), error = function(error) NULL)
    if (is.null(tab) || !"category_id" %in% names(tab)) return("")
    ids <- unique(as.character(tab$category_id))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    if (length(ids) != 1L || !ids[[1L]] %in% as.character(category_ids)) return("")
    families[[i]]
  }, character(1L))
  attached_families <- unique(attached_families[nzchar(attached_families)])
  files <- files[lisa_extension_figure_family(files) %in% attached_families]
  for (path in files) {
    rel <- substring(path, nchar(source_dir) + 2L)
    target <- file.path(destination, rel)
    collection_error <- NULL
    copied <- tryCatch({
      lisa_guarded_copy(path, target, overwrite = FALSE,
                        run_root = managed_root)
      TRUE
    }, error = function(error) {
      if (inherits(error, "lisa_filesystem_recovery_error")) stop(error)
      collection_error <<- error
      FALSE
    })
    if (!copied) {
      stop(
        "LISA-EXTENSION-017 failed to collect rendered artifact: ",
        basename(path), "; ", conditionMessage(collection_error),
        call. = FALSE
      )
    }
  }
  invisible(TRUE)
}

lisa_extension_collect_inventory <- function(index_path, source_dir, destination,
                                             columns = c("output_png", "output_pdf"),
                                             stem_pattern = "_painted[.](png|pdf|svg)$",
                                             companions = c("_painted.png",
                                               "_painted.pdf", "_painted.svg",
                                               "_painted_source.tsv",
                                               "_painted_recipe.R",
                                               "_kegg_base.png")) {
  if (!file.exists(index_path)) {
    stop("LISA-EXTENSION-028 the generator wrote no product inventory at ",
      basename(index_path), "; refusing to guess which files were produced.",
      call. = FALSE)
  }
  managed_root <- getOption("lisaR.run_root", NULL)
  if (is.null(managed_root)) {
    stop("LISA-EXTENSION-028 collection requires an authorized run root.",
      call. = FALSE)
  }
  index <- read_lisa_tsv(index_path)
  named <- character()
  for (column in columns) {
    if (column %in% names(index)) named <- c(named, as.character(index[[column]]))
  }
  named <- named[nzchar(named) & !is.na(named)]
  if (!length(named)) {
    stop("LISA-EXTENSION-029 the product inventory names no generated figure; ",
      "the request produced nothing to return.", call. = FALSE)
  }
  # The inventory names the image files. Their source table, recipe and preserved
  # base diagram are the same product in the sense of the presentation policy, so
  # they are collected with the image and nothing else in the directory is.
  stems <- unique(sub(stem_pattern, "", basename(named)))
  wanted <- unlist(lapply(stems, function(stem) paste0(stem, companions)),
                   use.names = FALSE)
  wanted <- unique(c(wanted, basename(index_path)))
  collected <- character()
  for (name in wanted) {
    path <- file.path(source_dir, name)
    if (!file.exists(path)) next
    lisa_guarded_copy(path, file.path(destination, name), overwrite = FALSE,
      run_root = managed_root)
    collected <- c(collected, name)
  }
  missing_images <- setdiff(basename(named), collected)
  if (length(missing_images)) {
    stop("LISA-EXTENSION-030 the product inventory names file(s) the generator ",
      "did not write: ", paste(utils::head(missing_images, 3L), collapse = ", "),
      call. = FALSE)
  }
  invisible(collected)
}

lisa_extension_planned_scripts <- function(units) {
  products <- unique(as.character(units$product))
  scripts <- character()
  if ("gene_cards" %in% products) {
    scripts <- c(scripts, "build_single_de_category_gene_cards.R")
  }
  if ("volcano" %in% products) {
    scripts <- c(scripts, "build_single_de_category_volcano_overlays.R")
  }
  if ("heatmap" %in% products) {
    scripts <- c(scripts, "build_single_de_leading_edge_gene_heatmaps.R")
  }
  if ("kegg_pathway_map" %in% products) {
    scripts <- c(scripts, "build_single_de_kegg_pathway_painter.R")
  }
  if ("de_recurrent_genes" %in% products) {
    scripts <- c(scripts, "build_single_de_recurrent_gene_screen.R")
  }
  if (any(c("contrast_gene_card", "contrast_paired_heatmap") %in% products)) {
    scripts <- c(scripts, "build_contrast_gene_level_product.R")
  }
  if ("contrast_gene_category_network" %in% products) {
    scripts <- c(scripts, "build_contrast_gene_category_network.R")
  }
  if ("contrast_kegg_map" %in% products) {
    scripts <- c(scripts, "build_contrast_kegg_pathway_painter.R")
  }
  sort(unique(scripts))
}

lisa_extension_write_recipe <- function(stem, package_dir, code_ledger) {
  recipe_identity <- lisa_verify_code_recipe_ledger(
    code_ledger, package_dir, "reproduce_lisa_figure.R"
  )
  lisa_copy_verified_code_file(
    recipe_identity$path, paste0(stem, "_recipe.R"),
    recipe_identity$sha256
  )
}

lisa_extension_render_contrast <- function(unit, plan, work_root, report_root,
                                           package_dir, code_ledger) {
  side <- NULL
  contrast_id <- unit$contrast_id[[1]]; collection <- unit$collection[[1]]
  categories <- sort(unique(unit$category_id)); products <- sort(unique(unit$product))
  contrast <- read_lisa_tsv(file.path(work_root, unit$contrast_path[[1]]))
  contrast <- contrast[as.character(contrast$category_id) %in% categories, , drop = FALSE]
  formats <- names(plan$report$formats)[plan$report$formats]
  render_formats <- if (length(formats)) formats else "png"
  base <- file.path(report_root, "artifacts", "contrasts", contrast_id, collection)
  lisa_guarded_dir_create(base)
  for (i in seq_len(nrow(contrast))) {
    row <- contrast[i, , drop = FALSE]
    row$gsea_padj_cutoff <- plan$gsea_padj_cutoff
    category <- gsub("^_|_$", "", gsub("_+", "_", gsub("[^A-Za-z0-9._-]+", "_", row$category_id[[1]])))
    label <- if ("category_display_name" %in% names(row) && nzchar(as.character(row$category_display_name[[1]])))
      as.character(row$category_display_name[[1]]) else as.character(row$category_id[[1]])
    label_a <- if ("contrast_a_label" %in% names(row)) as.character(row$contrast_a_label[[1]]) else "A"
    label_b <- if ("contrast_b_label" %in% names(row)) as.character(row$contrast_b_label[[1]]) else "B"
    n_a <- if ("n_genesets_A" %in% names(row)) as.numeric(row$n_genesets_A[[1]]) else as.integer(is.finite(as.numeric(row$mean_NES_A[[1]])))
    n_b <- if ("n_genesets_B" %in% names(row)) as.numeric(row$n_genesets_B[[1]]) else as.integer(is.finite(as.numeric(row$mean_NES_B[[1]])))
    support_a <- if ("has_significant_support_A" %in% names(row)) isTRUE(as.logical(row$has_significant_support_A[[1]])) else n_a > 0
    support_b <- if ("has_significant_support_B" %in% names(row)) isTRUE(as.logical(row$has_significant_support_B[[1]])) else n_b > 0
    has_signal <- support_a || support_b
    display_a <- if ("display_mean_NES_A" %in% names(row)) as.numeric(row$display_mean_NES_A[[1]]) else if (has_signal) as.numeric(row$mean_NES_A[[1]]) else NA_real_
    display_b <- if ("display_mean_NES_B" %in% names(row)) as.numeric(row$display_mean_NES_B[[1]]) else if (has_signal) as.numeric(row$mean_NES_B[[1]]) else NA_real_
    values <- data.frame(
      side = factor(c(label_a, label_b), levels = c(label_a, label_b)),
      mean_NES = c(display_a, display_b),
      has_significant_support = c(support_a, support_b),
      stringsAsFactors = FALSE
    )
    if ("contrast_profile" %in% products) {
      product_dir <- file.path(base, "contrast_profile")
      lisa_guarded_dir_create(product_dir)
      stem <- file.path(product_dir, category)
      p <- plot_lisa_contrast_dumbbell(
        row, label_a, label_b, label_a, label_b,
        title = label,
        subtitle = paste(label_a, "versus", label_b),
        group_by_supracategory = FALSE,
        annotate_gene_sets = FALSE
      )
      save_plot_multi(p, stem, render_formats, width = 7, height = 4.2, bg = "white")
      lisa_write_contrast_figure_source(
        row, paste0(stem, "_source.tsv"), contrast_id, label,
        paste(label_a, "versus", label_b), "contrast_profile", "plain",
        label_a, label_b, plan$gsea_padj_cutoff,
        group_by_supracategory = FALSE,
        figure_width = 7,
        figure_height = 4.2,
        figure_dpi = 300
      )
      lisa_extension_write_recipe(stem, package_dir, code_ledger)
    }
    if ("contrast_heatmap" %in% products) {
      product_dir <- file.path(base, "contrast_heatmap")
      lisa_guarded_dir_create(product_dir)
      stem <- file.path(product_dir, category)
      if (has_signal && all(is.finite(values$mean_NES))) {
        p <- ggplot2::ggplot(values, ggplot2::aes(x = side, y = label, fill = mean_NES,
          alpha = has_significant_support)) +
          ggplot2::geom_tile(color = "white", linewidth = 0.5) +
          ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f%s", mean_NES,
            ifelse(has_significant_support, "", " (NS)"))), size = 4) +
          ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
          ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
          ggplot2::labs(title = label, x = NULL, y = NULL, fill = "Mean NES") +
          ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.grid = ggplot2::element_blank())
      } else {
        p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
          ggplot2::geom_blank() +
          ggplot2::annotate("text", x = 1, y = 1,
            label = sprintf("No member gene set met GSEA FDR <= %s on either side", plan$gsea_padj_cutoff),
            color = "grey35", size = 3.8) +
          ggplot2::labs(title = label, x = NULL, y = NULL) + ggplot2::theme_void()
      }
      save_plot_multi(p, stem, render_formats, width = 5.5, height = 3.8, bg = "white")
      source <- row; source$figure_id <- paste0("contrast_heatmap__", category)
      source$figure_type <- "lisa_contrast_heatmap"; source$source_row_order <- 1L
      source$figure_width <- 5.5; source$figure_height <- 3.8; source$figure_dpi <- 300
      source$figure_title <- label; source$figure_subtitle <- paste(label_a, "versus", label_b)
      source$selected_for_plot <- TRUE; source$highlighted <- TRUE; source$labelled <- TRUE
      write_lisa_tsv(source, paste0(stem, "_source.tsv")); lisa_extension_write_recipe(stem, package_dir, code_ledger)
    }
  }
  invisible(TRUE)
}

lisa_extension_render_unit <- function(unit, plan, work_root, report_root,
                                       package_dir, code_ledger) {
  path_columns <- grep("_path$", names(unit), value = TRUE)
  stage_paths <- unique(unlist(unit[path_columns], use.names = FALSE))
  stage_paths <- stage_paths[nzchar(stage_paths)]
  lisa_extension_stage_files(plan$source_run, work_root, stage_paths)
  if (identical(unit$unit_type[[1]], "contrast")) {
    if (any(unique(as.character(unit$product)) %in% lisa_explore_h3_products())) {
      return(lisa_extension_render_contrast_native(
        unit, plan, work_root, report_root, package_dir, code_ledger
      ))
    }
    return(lisa_extension_render_contrast(
      unit, plan, work_root, report_root, package_dir, code_ledger
    ))
  }
  analysis <- unit$analysis_id[[1]]; collection <- unit$collection[[1]]
  categories <- sort(unique(unit$category_id)); products <- sort(unique(unit$product))
  formats <- names(plan$report$formats)[plan$report$formats]
  render_formats <- if (length(formats)) formats else "png"
  category_arg <- paste(categories, collapse = ",")
  target_base <- file.path(report_root, "artifacts", analysis, collection)

  if ("member_gene_sets" %in% products) {
    member_categories <- sort(unique(unit$category_id[unit$product == "member_gene_sets"]))
    gsea <- read_lisa_tsv(file.path(work_root, unit$gsea_path[[1]]))
    gsea <- lisa_classified_gsea_rows(gsea)
    summary <- read_lisa_tsv(file.path(work_root, unit$summary_path[[1]]))
    summary <- summary[as.character(summary$category_id) %in% member_categories, , drop = FALSE]
    member_dir <- file.path(target_base, "member_gene_sets")
    write_lisa_gsea_category_pathway_plots(gsea, summary, gsea_padj_cutoff = plan$gsea_padj_cutoff,
      plot_order = "dictionary", output_dir = member_dir,
      plot_formats = render_formats, stable_names = TRUE,
      code_ledger = code_ledger)
    index <- file.path(member_dir, "plot_manifest.tsv")
    if (file.exists(index)) lisa_guarded_delete(index)
  }
  if ("kegg" %in% products) {
    kegg_categories <- sort(unique(unit$category_id[unit$product == "kegg"]))
    gsea <- read_lisa_tsv(file.path(work_root, unit$gsea_path[[1]]))
    gsea <- lisa_classified_gsea_rows(gsea)
    gsea <- gsea[grepl("^KEGG_", as.character(gsea$pathway)), , drop = FALSE]
    summary <- read_lisa_tsv(file.path(work_root, unit$summary_path[[1]]))
    summary <- summary[as.character(summary$category_id) %in% kegg_categories, , drop = FALSE]
    kegg_dir <- file.path(target_base, "kegg")
    write_lisa_gsea_category_pathway_plots(gsea, summary, gsea_padj_cutoff = plan$gsea_padj_cutoff,
      plot_order = "dictionary", output_dir = kegg_dir,
      plot_formats = render_formats, stable_names = TRUE,
      code_ledger = code_ledger)
    index <- file.path(kegg_dir, "plot_manifest.tsv")
    if (file.exists(index)) lisa_guarded_delete(index)
  }

  gene_products <- intersect(products, c("gene_cards", "volcano", "heatmap",
    "kegg_pathway_map", "de_recurrent_genes"))
  if (length(gene_products)) {
    prefix <- lisa_extension_file_label_prefix(unit$gsea_path[[1]], analysis)
    built <- lisa_build_single_gene_level_tables(analysis, collection,
      file.path(work_root, "outputs"), prefix)
    if (!identical(built$status, "completed")) stop("LISA-EXTENSION-018 gene-level derivative build failed: ", built$message, call. = FALSE)
  }
  run_builder <- function(script, args) {
    result <- lisa_run_post_script(package_dir, script,
      c("--project-dir", work_root, "--analysis-id", analysis,
        "--universe", collection, args), trusted_run_root = work_root,
      code_ledger = code_ledger)
    if (!identical(result$status, "completed")) stop("LISA-EXTENSION-019 builder failed: ", script, ": ", result$message, call. = FALSE)
  }
  gene_root <- file.path(work_root, "outputs", "gene_level", "single_de", analysis, paste0("collection_", collection))
  if ("gene_cards" %in% products) {
    run_builder("build_single_de_category_gene_cards.R", c("--categories", category_arg,
      "--formats", paste(render_formats, collapse = ",")))
    lisa_extension_collect(file.path(gene_root, "category_gene_cards"),
      file.path(target_base, "gene_cards"), categories)
  }
  if ("volcano" %in% products) {
    run_builder("build_single_de_category_volcano_overlays.R", c("--categories", category_arg,
      "--formats", paste(render_formats, collapse = ",")))
    lisa_extension_collect(file.path(gene_root, "category_volcano_overlays"),
      file.path(target_base, "volcano"), categories)
  }
  if ("heatmap" %in% products) {
    scale <- lisa_extension_unit_scalar(unit, "variant", "heatmap", "zscore")
    matrix_path <- lisa_extension_unit_scalar(unit, "matrix_path", "heatmap")
    matrix_source_column <- lisa_extension_unit_scalar(
      unit, "matrix_source_column", "heatmap"
    )
    staged_matrix <- file.path(work_root, matrix_path)
    run_builder("build_single_de_leading_edge_gene_heatmaps.R", c(
      "--categories", category_arg,
      "--formats", paste(render_formats, collapse = ","), "--scale", scale,
      "--expression-matrix", staged_matrix,
      "--expression-matrix-column", matrix_source_column
    ))
    lisa_extension_collect(file.path(gene_root, "leading_edge_gene_heatmaps"),
      file.path(target_base, "heatmap"), categories)
  }
  if ("kegg_pathway_map" %in% products) {
    lisa_extension_render_kegg_map(unit, gene_root, target_base, run_builder)
  }
  if ("de_recurrent_genes" %in% products) {
    run_builder("build_single_de_recurrent_gene_screen.R", character())
    recurrent_dir <- file.path(gene_root, "recurrent_gene_screen")
    index_path <- lisa_extension_find_one(recurrent_dir,
      "_recurrent_gene_screen_index[.]tsv$")
    if (!nzchar(index_path)) {
      stop("LISA-EXTENSION-035 the recurrent gene screen wrote no product ",
        "inventory; refusing to guess which files were produced.", call. = FALSE)
    }
    lisa_extension_collect_inventory(index_path, recurrent_dir,
      file.path(target_base, "de_recurrent_genes"),
      columns = c("plot_png", "plot_pdf", "plot_svg"),
      stem_pattern = lisa_format_pattern(),
      companions = c(".png", ".pdf", ".svg", "_source.tsv", "_recipe.R"))
  }
  invisible(TRUE)
}

lisa_extension_render_contrast_native <- function(unit, plan, work_root,
                                                  report_root, package_dir,
                                                  code_ledger) {
  products <- intersect(unique(as.character(unit$product)),
                        lisa_explore_h3_products())
  if (length(products) != 1L) {
    stop("LISA-EXTENSION-036 a native contrast render unit must name exactly ",
      "one product; it named ", length(products), ".", call. = FALSE)
  }
  product <- products[[1L]]
  contrast_id <- unit$contrast_id[[1]]
  collection <- unit$collection[[1]]
  prepared <- lisa_extension_stage_prepared_contrast(unit, plan, work_root,
                                                     product)
  base <- file.path(report_root, "artifacts", "contrasts", contrast_id,
                    collection)
  lisa_guarded_dir_create(base)
  out_root <- file.path(work_root, "outputs", "gene_level", "category_contrasts",
                        contrast_id, paste0("collection_", collection))
  prefix <- paste(prepared$output_id, collection, "contrast_gene_level",
                  sep = "_")
  run_builder <- function(script, args) {
    result <- lisa_run_post_script(package_dir, script,
      c("--project-dir", work_root, "--contrast-id", contrast_id,
        "--universe", collection,
        "--paired-input", prepared$paired, "--summary-input", prepared$summary,
        args), trusted_run_root = work_root, code_ledger = code_ledger)
    if (!identical(result$status, "completed")) {
      stop("LISA-EXTENSION-019 builder failed: ", script, ": ", result$message,
        call. = FALSE)
    }
  }
  collect <- function(index_path, source_dir, stem_pattern, companions) {
    if (!nzchar(index_path)) {
      stop("LISA-EXTENSION-037 the ", product, " generator wrote no product ",
        "inventory; refusing to guess which files were produced.", call. = FALSE)
    }
    lisa_extension_collect_inventory(index_path, source_dir,
      file.path(base, product), stem_pattern = stem_pattern,
      companions = companions)
  }
  plain <- c(".png", ".pdf", ".svg", "_source.tsv", "_recipe.R")

  if (identical(product, "contrast_gene_card")) {
    category <- lisa_extension_unit_scalar(unit, "category_id", product)
    if (!nzchar(category)) {
      stop("LISA-EXTENSION-038 a contrast gene card request must name one ",
        "category.", call. = FALSE)
    }
    run_builder("build_contrast_gene_level_product.R",
      c("--exact-product", "contrast_gene_card", "--category-id", category))
    cards_dir <- file.path(out_root, "contrast_category_cards")
    collect(lisa_extension_find_one(cards_dir,
      "_contrast_gene_card_index[.]tsv$"), cards_dir, "[.](png|pdf|svg)$", plain)
  } else if (identical(product, "contrast_paired_heatmap")) {
    run_builder("build_contrast_gene_level_product.R",
      c("--exact-product", "contrast_paired_heatmap"))
    collect(lisa_extension_find_one(out_root,
      "_paired_gene_heatmap_index[.]tsv$"), out_root, "[.](png|pdf|svg)$", plain)
  } else if (identical(product, "contrast_gene_category_network")) {
    variant <- lisa_extension_unit_scalar(unit, "variant", product,
                                          "all_categories_top7")
    run_builder("build_contrast_gene_category_network.R",
      c("--variant", variant))
    network_dir <- file.path(out_root, "contrast_gene_category_network")
    collect(lisa_extension_find_one(network_dir,
      "_contrast_gene_category_network_index[.]tsv$"), network_dir,
      "[.](png|pdf|svg)$", plain)
  } else {
    lisa_extension_render_contrast_kegg_map(unit, out_root, base, run_builder,
                                            collect)
  }
  invisible(TRUE)
}

# The native contrast pathway map: exactly one map, selected by pathway id,
# painted with the complete contrast-wide ranking and paint context of its
# accepted generator. `--kegg-id` is applied inside the painter AFTER the whole
# contrast has been ranked, so the map's `contrast_kegg_rank`, title, A/B node
# halves, thresholds and colour scale are the ones the default generator would
# have produced for that pathway. Resources are resolved `cache_only`: a missing
# snapshot is an error naming what is absent, never a download or an install.
lisa_extension_render_contrast_kegg_map <- function(unit, out_root, base,
                                                    run_builder, collect) {
  kegg_id <- lisa_extension_unit_scalar(unit, "entity", "contrast_kegg_map")
  if (!nzchar(kegg_id)) {
    stop("LISA-EXTENSION-039 a native contrast KEGG map request must name one ",
      "pathway.", call. = FALSE)
  }
  resource <- function(column, label) {
    value <- lisa_extension_unit_scalar(unit, column, "contrast_kegg_map")
    if (!nzchar(value)) {
      stop("LISA-EXTENSION-033 native KEGG maps need ", label,
        "; none is declared for this run. Declare the cached snapshot ",
        "explicitly -- lisaR never downloads or installs one.", call. = FALSE)
    }
    value
  }
  species <- resource("kegg_species", "a species")
  cache_root <- resource("kegg_cache_root", "a KEGG snapshot cache root")
  snapshot_id <- resource("kegg_snapshot_id", "a KEGG snapshot id")
  recorded_digest <- resource("kegg_resource_digest",
                              "a validated KEGG resource digest")
  organism <- lisa_species_contract(species)$kegg_code
  observed_digest <- lisa_explore_kegg_pathway_resource_digest(
    cache_root, snapshot_id, organism, kegg_id)
  if (!identical(recorded_digest, observed_digest)) {
    stop("LISA-EXTENSION-034 the validated KEGG resources for ", kegg_id,
      " no longer match the prepared index. Prepare the contrast products ",
      "again against a new immutable snapshot before painting.", call. = FALSE)
  }
  args <- c("--species", species, "--kegg-cache-root", cache_root,
    "--kegg-snapshot-id", snapshot_id, "--kegg-access-mode", "cache_only",
    "--kegg-id", kegg_id)
  for (setting in c("max_abs_log2fc", "color_power")) {
    value <- lisa_extension_unit_scalar(unit, setting, "contrast_kegg_map")
    if (nzchar(value)) args <- c(args, paste0("--", gsub("_", "-", setting)), value)
  }
  run_builder("build_contrast_kegg_pathway_painter.R", args)
  painter_dir <- file.path(out_root, "contrast_kegg_pathway_painter")
  # The painter names its exact-mode inventory after the pathway, so it can
  # never overwrite or be confused with a default full-run index.
  collect(file.path(painter_dir,
    paste0(kegg_id, "_contrast_kegg_pathway_painter_index.tsv")), painter_dir,
    "_contrast_painted[.](png|pdf|svg)$",
    c("_contrast_painted.png", "_contrast_painted.pdf", "_contrast_painted.svg",
      "_contrast_painted_source.tsv", "_contrast_painted_recipe.R",
      "_kegg_base.png"))
}

lisa_extension_stage_prepared_contrast <- function(unit, plan, work_root,
                                                   product) {
  prepared_dir <- lisa_extension_unit_scalar(unit, "prepared_dir", product)
  if (!nzchar(prepared_dir) || !dir.exists(prepared_dir)) {
    stop("LISA-EXTENSION-040 the native contrast product ", product,
      " needs a prepared bundle; none is recorded for this request. Run ",
      "lisa_explore_prepare_contrast_products() first.", call. = FALSE)
  }
  if (lisa_path_within(prepared_dir, plan$source_run)) {
    stop("LISA-EXTENSION-041 a prepared contrast bundle may never live inside ",
      "the immutable source run.", call. = FALSE)
  }
  receipt <- lisa_explore_contrast_read_receipt(prepared_dir)
  if (is.null(receipt)) {
    stop("LISA-EXTENSION-042 the prepared contrast bundle at ", prepared_dir,
      " has no readable preparation receipt.", call. = FALSE)
  }
  expected <- lisa_extension_unit_scalar(unit, "preparation_digest", product)
  observed <- lisa_explore_contrast_preparation_digest(receipt)
  if (nzchar(expected) && !identical(expected, observed)) {
    stop("LISA-EXTENSION-043 the prepared contrast bundle changed after this ",
      "request was planned; prepare again before rendering.", call. = FALSE)
  }
  staged <- file.path(work_root, "prepared_contrast")
  lisa_guarded_dir_create(staged)
  resolved <- list(output_id = as.character(receipt$output_id))
  for (file in receipt$files) {
    relative <- lisa_explore_presentation_path(as.character(file$path),
                                               "prepared file")
    source <- file.path(prepared_dir, relative)
    if (!file.exists(source) ||
        !identical(lisa_sha256_file(source), as.character(file$sha256))) {
      stop("LISA-EXTENSION-044 a prepared contrast file is missing or changed ",
        "on disk: ", relative, ". Prepare the contrast products again.",
        call. = FALSE)
    }
    target <- file.path(staged, relative)
    lisa_guarded_copy(source, target, overwrite = FALSE, run_root = work_root)
  }
  bundle_names <- lisa_explore_contrast_bundle_files()
  resolved$paired <- file.path(staged, bundle_names[["paired"]])
  resolved$summary <- file.path(staged, bundle_names[["summary"]])
  if (!file.exists(resolved$paired) || !file.exists(resolved$summary)) {
    stop("LISA-EXTENSION-045 the prepared contrast bundle does not carry both ",
      "the paired evidence and its ranked summary.", call. = FALSE)
  }
  resolved
}

# One exact value for one product inside a render unit. A unit groups rows by
# owner and collection, so it can legitimately carry several products at once;
# this refuses rather than guesses when one product's rows disagree, which is how
# a wrong scale or a wrong map is kept from being rendered silently.
lisa_extension_unit_scalar <- function(unit, column, product, default = "") {
  if (!column %in% names(unit)) return(default)
  values <- unique(as.character(unit[[column]][unit$product == product]))
  values <- values[!is.na(values) & nzchar(values)]
  if (!length(values)) return(default)
  if (length(values) != 1L) {
    stop("LISA-EXTENSION-031 the render unit names ", length(values), " ", column,
      " values for product ", product, " (",
      paste(utils::head(values, 5L), collapse = ", "),
      "); exactly one is required.", call. = FALSE)
  }
  values[[1L]]
}

# The native KEGG pathway map: exactly one map, selected by pathway id, painted
# with the collection-wide ranking and paint context of its accepted generator.
#
# `--kegg-id` is applied inside the painter AFTER the whole collection has been
# ranked, so the map's rank, title, node values, colour scale and audit table are
# the ones the default generator would have produced for that pathway. Resources
# are resolved `cache_only`: a missing snapshot is an error naming what is
# absent, never a download or an install.
lisa_extension_render_kegg_map <- function(unit, gene_root, target_base, run_builder) {
  kegg_id <- lisa_extension_unit_scalar(unit, "entity", "kegg_pathway_map")
  if (!nzchar(kegg_id)) {
    stop("LISA-EXTENSION-032 a native KEGG map request must name one pathway.",
      call. = FALSE)
  }
  resource <- function(column, label) {
    value <- lisa_extension_unit_scalar(unit, column, "kegg_pathway_map")
    if (!nzchar(value)) {
      stop("LISA-EXTENSION-033 native KEGG maps need ", label,
        "; none is declared for this run. Declare the cached snapshot ",
        "explicitly -- lisaR never downloads or installs one.", call. = FALSE)
    }
    value
  }
  species <- resource("kegg_species", "a species")
  cache_root <- resource("kegg_cache_root", "a KEGG snapshot cache root")
  snapshot_id <- resource("kegg_snapshot_id", "a KEGG snapshot id")
  recorded_digest <- resource("kegg_resource_digest",
                              "a validated KEGG resource digest")
  organism <- lisa_species_contract(species)$kegg_code
  observed_digest <- lisa_explore_kegg_pathway_resource_digest(
    cache_root, snapshot_id, organism, kegg_id)
  if (!identical(recorded_digest, observed_digest)) {
    stop("LISA-EXTENSION-034 the validated KEGG resources for ", kegg_id,
      " no longer match the prepared index. Prepare the pathway index again ",
      "against a new immutable snapshot before painting.", call. = FALSE)
  }
  args <- c("--species", species, "--kegg-cache-root", cache_root,
    "--kegg-snapshot-id", snapshot_id, "--kegg-access-mode", "cache_only",
    "--kegg-id", kegg_id)
  for (setting in c("max_abs_log2fc", "color_power")) {
    value <- lisa_extension_unit_scalar(unit, setting, "kegg_pathway_map")
    if (nzchar(value)) {
      args <- c(args, paste0("--", gsub("_", "-", setting)), value)
    }
  }
  run_builder("build_single_de_kegg_pathway_painter.R", args)
  painter_dir <- file.path(gene_root, "kegg_painter")
  index_path <- file.path(painter_dir,
    paste0(basename(gene_root), "_kegg_pathway_painter_index.tsv"))
  # The painter names its exact-mode inventory after the pathway so it can never
  # overwrite or be confused with a default full-run index.
  exact_index <- file.path(painter_dir,
    paste0(kegg_id, "_kegg_pathway_painter_index.tsv"))
  if (file.exists(exact_index)) index_path <- exact_index
  lisa_extension_collect_inventory(index_path, painter_dir,
    file.path(target_base, "kegg_pathway_map"))
}

lisa_extension_apply_file_policy <- function(root, report) {
  root <- lisa_existing_run_root(root)
  tree <- lisa_scan_run_tree(root)
  files <- tree$path[!tree$isdir]
  remove <- rep(FALSE, length(files))
  if (!report$formats[["png"]]) {
    remove <- remove | grepl("[.]png$", files, ignore.case = TRUE)
  }
  if (!report$formats[["svg"]]) {
    remove <- remove | grepl("[.]svg$", files, ignore.case = TRUE)
  }
  if (!report$formats[["pdf"]]) {
    remove <- remove | grepl("[.]pdf$", files, ignore.case = TRUE)
  }
  if (!report$source_data) {
    remove <- remove |
      grepl("_(source|matrix)[.]tsv$", files, ignore.case = TRUE)
  }
  if (!report$recipes) {
    remove <- remove | grepl("_recipe[.]R$", files, ignore.case = TRUE)
  }
  for (path in sort(files[remove])) {
    lisa_guarded_delete(path, run_root = root)
  }
  invisible(TRUE)
}

lisa_extension_validate_manifest <- function(report_root,
                                             manifest_path = file.path(
                                               report_root,
                                               "extension_manifest.tsv"
                                             )) {
  fail <- function(detail) {
    stop("LISA-EXTENSION-023 staging validation failed: ", detail,
         call. = FALSE)
  }
  report_root <- lisa_existing_run_root(report_root)
  manifest_path <- lisa_guarded_path(manifest_path, report_root)
  tryCatch(
    lisa_assert_regular_managed_file(
      manifest_path, report_root, "Extension manifest"
    ),
    error = function(error) fail(conditionMessage(error))
  )
  manifest <- tryCatch(read_lisa_tsv(manifest_path), error = function(error) {
    fail(paste0("the manifest could not be read: ", conditionMessage(error)))
  })
  required <- c("path", "bytes", "sha256")
  if (!identical(names(manifest), required)) {
    fail("the manifest schema is not exactly path/bytes/sha256")
  }
  paths <- as.character(manifest$path)
  if (anyNA(paths) || any(!nzchar(paths)) || anyDuplicated(paths) ||
      any(grepl("[[:cntrl:]\\\\]", paths)) ||
      any(grepl("^/|^[A-Za-z]:/", paths))) {
    fail("manifest paths must be unique, relative, and portable")
  }
  unsafe_components <- vapply(strsplit(paths, "/", fixed = TRUE), function(x) {
    any(!nzchar(x) | x %in% c(".", ".."))
  }, logical(1))
  if (any(unsafe_components)) fail("manifest paths contain traversal")
  if (anyNA(manifest$bytes) ||
      any(!is.finite(as.numeric(manifest$bytes))) ||
      any(as.numeric(manifest$bytes) < 0) ||
      any(as.numeric(manifest$bytes) != floor(as.numeric(manifest$bytes)))) {
    fail("manifest byte counts are invalid")
  }
  hashes <- as.character(manifest$sha256)
  if (anyNA(hashes) || any(!grepl("^[0-9a-f]{64}$", hashes))) {
    fail("manifest SHA-256 values are invalid")
  }

  tree <- lisa_scan_run_tree(report_root)
  actual <- substring(tree$path[!tree$isdir], nchar(report_root) + 2L)
  manifest_relative <- substring(manifest_path, nchar(report_root) + 2L)
  actual <- sort(setdiff(actual, manifest_relative))
  if (!identical(paths, sort(paths)) || !identical(paths, actual)) {
    fail("declared and observed files differ")
  }
  declared_files <- file.path(report_root, paths)
  if (length(declared_files)) {
    declared_files <- vapply(
      declared_files, lisa_guarded_path, character(1), run_root = report_root
    )
  }
  observed_bytes <- if (length(declared_files)) {
    as.numeric(file.info(declared_files)$size)
  } else {
    numeric()
  }
  observed_hashes <- vapply(
    declared_files, lisa_sha256_file, character(1), USE.NAMES = FALSE
  )
  if (!identical(observed_bytes, as.numeric(manifest$bytes)) ||
      !identical(observed_hashes, hashes)) {
    fail("a declared file changed after the manifest was written")
  }
  invisible(manifest)
}

lisa_extension_html_escape <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

lisa_extension_url_path <- function(x) {
  normalized <- gsub("\\\\", "/", as.character(x))
  vapply(strsplit(normalized, "/", fixed = TRUE), function(parts) {
    paste(vapply(parts, utils::URLencode, character(1), reserved = TRUE,
      repeated = TRUE), collapse = "/")
  }, character(1), USE.NAMES = FALSE)
}

# Format variants share one figure card. Heatmap source/recipe stems omit the
# display-scale suffix used by the image; all links retain their actual names.
#
# `raw` is in this list because the builder has always accepted `--scale raw`
# while this pattern did not know about it: a raw-scale image's family would have
# been `..._leading_edge_gene_heatmap_raw`, which matches no `_matrix.tsv`
# sidecar family, so `lisa_extension_collect()` would have silently dropped the
# very figure that was requested. FULL only ever renders `zscore`, so no existing
# file name contains `_raw` and adding it changes nothing for accepted runs.
lisa_extension_figure_family <- function(path) {
  stem <- sub("_(source|matrix)[.]tsv$|_settings[.]json$|_recipe[.]R$|[.](png|svg|pdf)$", "",
    path, ignore.case = TRUE)
  sub("(_leading_edge_gene_heatmap)_(zscore|raw|log2|none)$", "\\1", stem)
}

lisa_extension_index_groups <- function(plan, paths) {
  groups <- list()
  units <- plan$units
  required <- c("unit_type", "analysis_id", "contrast_id", "collection", "product")
  if (is.data.frame(units) && nrow(units) && all(required %in% names(units))) {
    keys <- unique(units[required])
    for (i in seq_len(nrow(keys))) {
      key <- keys[i, , drop = FALSE]
      contrast <- identical(as.character(key$unit_type[[1L]]), "contrast")
      owner <- as.character(if (contrast) key$contrast_id[[1L]] else key$analysis_id[[1L]])
      collection <- as.character(key$collection[[1L]])
      product <- as.character(key$product[[1L]])
      prefix <- paste(c("artifacts", if (contrast) "contrasts", owner, collection, product), collapse = "/")
      selected <- rep(TRUE, nrow(units))
      for (name in required) selected <- selected & as.character(units[[name]]) == as.character(key[[name]][[1L]])
      groups[[length(groups) + 1L]] <- list(kind = if (contrast) "contrast" else "analysis",
        owner = owner, collection = collection, product = product,
        units = units[selected, , drop = FALSE], paths = paths[startsWith(paths, paste0(prefix, "/"))])
    }
  }
  assigned <- unlist(lapply(groups, `[[`, "paths"), use.names = FALSE)
  extra <- setdiff(paths, assigned)
  # Retain unclassified artifacts and the legacy index-only helper contract.
  if (length(extra) || !length(groups)) groups[[length(groups) + 1L]] <- list(
    kind = "analysis", owner = "Saved results", collection = "Other saved files",
    product = "saved_files", units = data.frame(), paths = extra)
  groups
}

lisa_extension_category_attachments <- function(report_root, group) {
  units <- group$units
  required <- c("category_id", "supercategory_id")
  if (!nrow(units) || !all(required %in% names(units))) {
    return(list(category_products = list(), unclassified = sort(unique(group$paths))))
  }
  paths <- sort(unique(group$paths))
  sidecars <- paths[grepl("_(source|matrix)[.]tsv$", paths, ignore.case = TRUE)]
  by_family <- split(paths, lisa_extension_figure_family(paths))
  attached <- list(); classified <- character()
  for (sidecar in sidecars) {
    absolute <- file.path(report_root, sidecar)
    source <- tryCatch(read_lisa_tsv(absolute), error = function(error) NULL)
    if (is.null(source) || !"category_id" %in% names(source)) next
    ids <- unique(as.character(source$category_id))
    ids <- ids[!is.na(ids) & nzchar(ids)]
    # The source must state exactly one ID and that ID must occur exactly once
    # in this already exact context/collection/product unit group.
    if (length(ids) != 1L) next
    row <- units[as.character(units$category_id) == ids[[1L]], , drop = FALSE]
    if (nrow(row) != 1L) next
    family <- lisa_extension_figure_family(sidecar)
    assets <- by_family[[family]]
    if (!length(assets)) next
    key <- paste(ids[[1L]], group$product, sep = "\r")
    if (is.null(attached[[key]])) {
      attached[[key]] <- list(
        category_id = ids[[1L]],
        supercategory_id = if (is.na(row$supercategory_id[[1L]])) "" else as.character(row$supercategory_id[[1L]]),
        product = group$product,
        label = basename(family), assets = assets
      )
    } else attached[[key]]$assets <- sort(unique(c(attached[[key]]$assets, assets)))
    classified <- c(classified, assets)
  }
  list(category_products = unname(attached),
    unclassified = sort(setdiff(paths, unique(classified))))
}

# C4: a PNG product whose SVG twin is absent is flagged whenever SVG was
# requested. Native KEGG diagrams paint over the KEGG raster; they are the
# documented exception and are labelled as such rather than as "missing".
lisa_svg_status_pill <- function(formats, product, files, svg_requested) {
  formats <- tolower(as.character(formats))
  if (!isTRUE(svg_requested) || !"png" %in% formats || "svg" %in% formats) return("")
  label <- if (lisa_svg_exception_product(product, files)) {
    "SVG not available for this product type"
  } else "SVG missing"
  paste0(' <span class="pill status-warn lisa-svg-status">', label, '</span>')
}

lisa_svg_exception_product <- function(product, files = character()) {
  any(grepl("^(kegg_pathway_map|contrast_kegg_map)$", as.character(product))) ||
    any(grepl("_painted[.]png$|_kegg_base[.]png$", as.character(files)))
}

lisa_extension_policy_formats <- function(report_root) {
  path <- file.path(report_root, "output_policy.tsv")
  if (!file.exists(path)) return(character())
  policy <- read_lisa_tsv(path)
  if (!all(c("product", "requested") %in% names(policy))) return(character())
  requested <- tolower(as.character(policy$requested)) %in% c("true", "t", "1", "yes")
  intersect(c("png", "svg", "pdf"), as.character(policy$product)[requested])
}

lisa_extension_gallery <- function(report_root, group, section_id,
                                   svg_requested = NA) {
  esc <- lisa_extension_html_escape
  if (is.na(svg_requested)) {
    svg_requested <- "svg" %in% lisa_extension_policy_formats(report_root)
  }
  link <- function(path, label) paste0('<a download href="',
    esc(lisa_extension_url_path(path)), '">', esc(label), '</a>')
  groups <- if (!is.null(group$paths)) list(group) else group
  attachments <- lapply(groups, function(one) lisa_extension_category_attachments(report_root, one))
  products <- unlist(lapply(attachments, `[[`, "category_products"), recursive = FALSE)
  unclassified <- sort(unique(unlist(lapply(attachments, `[[`, "unclassified"), use.names = FALSE)))
  by_category <- split(products, vapply(products, `[[`, character(1L), "category_id"))
  cards <- vapply(by_category, function(category_products) {
    category_id <- category_products[[1L]]$category_id
    supercategory <- category_products[[1L]]$supercategory_id
    item_html <- vapply(category_products, function(product) {
      files <- product$assets; formats <- tolower(tools::file_ext(files))
      images <- files[formats %in% c("png", "svg")]
      if (length(images)) images <- images[order(match(tolower(tools::file_ext(images)), c("png", "svg")))]
      labels <- ifelse(grepl("_source[.]tsv$", files), "Source TSV",
        ifelse(grepl("_matrix[.]tsv$", files), "Matrix TSV", ifelse(grepl("_recipe[.]R$", files), "Recipe R",
          ifelse(grepl("_settings[.]json$", files), "Settings JSON", toupper(formats)))))
      svg_note <- lisa_svg_status_pill(formats, product$product, files, svg_requested)
      paste0('<section class="lisa-extension-product"><h4>', esc(product$product), '</h4><div class="lisa-extension-downloads" aria-label="Download this figure">',
        paste(mapply(link, files, labels, USE.NAMES = FALSE), collapse = " &middot; "), svg_note, '</div>',
        if (length(images)) paste0('<a href="', esc(lisa_extension_url_path(images[[1L]])), '"><img src="',
          esc(lisa_extension_url_path(images[[1L]])), '" alt="', esc(category_id), '" loading="lazy" decoding="async"></a>') else
          '<p class="lisa-extension-empty">No inline image was requested. Use the available downloads above.</p>', '</section>')
    }, character(1L), USE.NAMES = FALSE)
    paste0('<article class="lisa-extension-card" data-extension-card data-category="', esc(category_id),
      '" data-supercategory="', esc(supercategory), '"><h3>', esc(category_id), '</h3>',
      if (nzchar(supercategory)) paste0('<p class="lisa-extension-super">', esc(supercategory), '</p>') else "",
      paste(item_html, collapse = "\n"), '</article>')
  }, character(1L), USE.NAMES = FALSE)
  supers <- sort(unique(vapply(products, `[[`, character(1L), "supercategory_id")))
  supers <- supers[nzchar(supers)]
  paste0('<div data-extension-gallery><div class="lisa-extension-controls">',
    '<label for="', section_id, '-search">Find a category <input type="search" id="', section_id,
    '-search" data-extension-search placeholder="Category name or identifier"></label>',
    if (length(supers)) paste0('<label for="', section_id, '-super">Supercategory <select id="',
      section_id, '-super" data-extension-super><option value="">All supercategories</option>',
      paste(vapply(supers, function(value) paste0('<option value="', esc(value), '">', esc(value), '</option>'),
        character(1L)), collapse = ""), '</select></label>') else "",
    '<button type="button" data-extension-prev>Previous</button>',
    '<button type="button" data-extension-next>Next</button>',
    '<span data-extension-count role="status" aria-live="polite"></span></div>',
    '<div class="lisa-extension-grid">', paste(cards, collapse = "\n"), '</div>',
    if (!length(cards)) '<p>No category product had an exact source metadata attachment. Saved files remain listed below.</p>' else "",
    if (length(unclassified)) paste0('<details class="lisa-extension-files"><summary>Unclassified saved files (',
      length(unclassified), ')</summary><p>These files lack one exact category attachment in their source metadata and are not assigned to a category card.</p><ul>', paste(vapply(unclassified, function(path) paste0('<li>',
        link(path, basename(path)), '</li>'), character(1L)), collapse = ""), '</ul></details>') else "",
    '</div>')
}

# An extension may reuse evidence already rendered with the immutable source
# run.  Do not regenerate its statistical inputs: copy only its saved tables
# and figures, replace the viewer chrome from the installed package, and add
# exact-ID references to the canonical extension artifacts. Legacy source runs
# legitimately lack these pages; callers then retain the older gallery surface.
lisa_extension_read_payload <- function(path, data_id) {
  if (!file.exists(path)) return(NULL)
  html <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  marker <- paste0('<script type="application/json" id="', data_id, '">')
  start <- regexpr(marker, html, fixed = TRUE)[[1L]]
  if (start < 1L) return(NULL)
  from <- start + nchar(marker); tail <- substr(html, from, nchar(html))
  end <- regexpr("</script>", tail, fixed = TRUE)[[1L]]
  if (end < 1L) return(NULL)
  tryCatch(jsonlite::fromJSON(substr(tail, 1L, end - 1L), simplifyVector = FALSE), error = function(e) NULL)
}

lisa_extension_payload_scope <- function(payload, owner, collection, contrast, selected_ids) {
  if (is.null(payload) || !is.list(payload) || is.null(payload$metadata) || is.null(payload$categories)) return(NULL)
  key <- if (contrast) "contrast_id" else "analysis_id"
  if (!identical(as.character(payload$metadata[[key]]), as.character(owner)) ||
      !identical(as.character(payload$metadata$collection), as.character(collection))) return(NULL)
  ids <- vapply(payload$categories, function(x) {
    if (!is.list(x) || length(x$category_id) != 1L) return(NA_character_)
    as.character(x$category_id)
  }, character(1L))
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) || any(!selected_ids %in% ids)) return(NULL)
  # A selected extension is a true subset, not a full source page with a
  # hidden selector. Filter only records that explicitly carry category_id;
  # context-wide tables (for example genes) remain intact by contract.
  keep <- function(rows) {
    if (!is.list(rows) || !length(rows) || !all(vapply(rows, is.list, logical(1L)))) return(rows)
    has_id <- vapply(rows, function(row) !is.null(row$category_id), logical(1L))
    if (!any(has_id)) return(rows)
    rows[vapply(rows, function(row) is.null(row$category_id) || as.character(row$category_id) %in% selected_ids, logical(1L))]
  }
  for (name in names(payload)) payload[[name]] <- keep(payload[[name]])
  payload
}

lisa_extension_copy_file <- function(from, to) {
  lisa_guarded_dir_create(dirname(to))
  if (!file.copy(from, to, overwrite = TRUE)) stop("Cannot stage saved extension evidence file.", call. = FALSE)
  if (!identical(lisa_sha256_file(from), lisa_sha256_file(to))) stop("Staged extension evidence hash mismatch.", call. = FALSE)
  invisible(to)
}

lisa_category_product_reference_asset <- function(
    report_root, evidence_dir, path, inventory = NULL,
    allowed_roots = c("artifacts", "outputs"),
    label = "Category product") {
  report_root <- normalizePath(report_root, winslash = "/", mustWork = TRUE)
  evidence_dir <- lisa_guarded_path(evidence_dir, report_root)
  evidence_dir <- normalizePath(
    evidence_dir, winslash = "/", mustWork = FALSE
  )
  path <- lisa_assert_regular_managed_file(path, report_root, label)
  source_path <- gsub(
    "\\\\", "/",
    lisa_presentation_relative(path, report_root)
  )
  parts <- strsplit(source_path, "/", fixed = TRUE)[[1L]]
  allowed_roots <- unique(as.character(allowed_roots))
  allowed_roots <- allowed_roots[!is.na(allowed_roots) & nzchar(allowed_roots)]
  if (!length(parts) || !length(allowed_roots) ||
      !parts[[1L]] %in% allowed_roots ||
      any(!nzchar(parts) | parts %in% c(".", "..")) ||
      grepl("^[A-Za-z]:|^/|[[:cntrl:]\\\\]", source_path)) {
    stop(label, " is outside the allowed canonical artifact roots.",
         call. = FALSE)
  }
  sha256 <- lisa_sha256_file(path)
  if (!is.null(inventory)) {
    required <- c("path", "sha256")
    if (!is.data.frame(inventory) || !all(required %in% names(inventory))) {
      stop(label, " reference requires an artifact path/SHA-256 inventory.",
           call. = FALSE)
    }
    inventory_path <- gsub("\\\\", "/", as.character(inventory$path))
    matched <- which(inventory_path == source_path)
    if (length(matched) != 1L ||
        !identical(as.character(inventory$sha256[[matched]]), sha256)) {
      stop(label, " is not uniquely bound to its artifact inventory: ",
           source_path, call. = FALSE)
    }
    sha256 <- as.character(inventory$sha256[[matched]])
  }
  relative <- lisa_presentation_relative(path, evidence_dir)
  list(
    format = tolower(tools::file_ext(source_path)),
    href = lisa_extension_url_path(relative),
    source_path = source_path,
    source_name = basename(source_path),
    sha256 = sha256
  )
}

lisa_extension_reference_products <- function(report_root, evidence_dir, products, inventory) {
  required <- c("path", "sha256")
  if (!is.data.frame(inventory) || !all(required %in% names(inventory))) {
    stop("Extension category-product references require the artifact path/SHA-256 inventory.", call. = FALSE)
  }
  inventory$path <- gsub("\\\\", "/", as.character(inventory$path))
  out <- lapply(products, function(product) {
    assets <- lapply(product$assets, function(path) {
      path <- gsub("\\\\", "/", as.character(path))
      parts <- strsplit(path, "/", fixed = TRUE)[[1L]]
      if (length(path) != 1L || is.na(path) || !nzchar(path) || length(parts) < 2L ||
          !identical(parts[[1L]], "artifacts") || any(!nzchar(parts) | parts %in% c(".", "..")) ||
          grepl("^[A-Za-z]:|^/|[[:cntrl:]\\\\]", path)) {
        stop("Extension category product is not one canonical contained artifact path.", call. = FALSE)
      }
      artifact <- file.path(report_root, path)
      lisa_category_product_reference_asset(
        report_root, evidence_dir, artifact, inventory,
        allowed_roots = "artifacts",
        label = "Extension category artifact"
      )
    })
    list(category_id = product$category_id, product = product$product, label = product$label, assets = assets)
  })
  unname(out)
}

lisa_extension_validate_product_references <- function(report_root, expected_by_page, inventory) {
  if (!length(expected_by_page)) return(invisible(TRUE))
  required <- c("path", "sha256")
  if (!is.data.frame(inventory) || !all(required %in% names(inventory))) {
    stop("Extension product-reference validation requires the artifact path/SHA-256 inventory.", call. = FALSE)
  }
  inventory$path <- gsub("\\\\", "/", as.character(inventory$path))
  for (page in names(expected_by_page)) {
    page_file <- file.path(report_root, page)
    lisa_assert_regular_managed_file(page_file, report_root, "Extension evidence page")
    data_id <- if (startsWith(page, "contrast_evidence/")) "contrast-evidence-data" else "evidence-data"
    payload <- lisa_extension_read_payload(page_file, data_id)
    if (is.null(payload) || !is.list(payload$category_products)) {
      stop("Extension evidence page has no readable category-product payload: ", page, call. = FALSE)
    }
    products <- payload$category_products
    assets <- unlist(lapply(products, function(product) {
      if (!is.list(product) || !is.list(product$assets)) return(list())
      product$assets
    }), recursive = FALSE)
    fields <- c("format", "href", "source_path", "source_name", "sha256")
    valid <- vapply(assets, function(asset) is.list(asset) && all(fields %in% names(asset)), logical(1L))
    if (!all(valid)) stop("Extension category-product payload has an incomplete asset record: ", page, call. = FALSE)
    source_paths <- vapply(assets, function(asset) gsub("\\\\", "/", as.character(asset$source_path)), character(1L))
    expected <- sort(unique(gsub("\\\\", "/", as.character(expected_by_page[[page]]))))
    if (anyDuplicated(source_paths) || !identical(sort(source_paths), expected)) {
      stop("Extension category-product payload does not match its classified artifact paths: ", page, call. = FALSE)
    }
    page_dir <- dirname(page_file)
    for (i in seq_along(assets)) {
      asset <- assets[[i]]; source_path <- source_paths[[i]]
      matched <- which(inventory$path == source_path)
      if (length(matched) != 1L || !startsWith(source_path, "artifacts/")) {
        stop("Extension category-product payload is not inventory-bound: ", source_path, call. = FALSE)
      }
      artifact <- file.path(report_root, source_path)
      lisa_assert_regular_managed_file(artifact, report_root, "Extension category artifact")
      relative <- lisa_presentation_relative(artifact, page_dir)
      expected_href <- lisa_extension_url_path(relative)
      if (!startsWith(expected_href, "../../../artifacts/") ||
          !identical(as.character(asset$href), expected_href) ||
          !identical(as.character(asset$source_name), basename(source_path)) ||
          !identical(tolower(as.character(asset$format)), tolower(tools::file_ext(source_path))) ||
          !identical(as.character(asset$sha256), as.character(inventory$sha256[[matched]]))) {
        stop("Extension category-product payload reference failed reconciliation: ", source_path, call. = FALSE)
      }
    }
  }
  invisible(TRUE)
}

lisa_extension_write_saved_evidence <- function(report_root, plan, group, products, inventory) {
  contrast <- identical(group$kind, "contrast")
  owner <- group$owner; collection <- group$collection
  exact_owner <- owner
  if (contrast) {
    index_path <- file.path(plan$source_run, "config", "contrast_index.tsv")
    index <- if (file.exists(index_path)) tryCatch(read_lisa_tsv(index_path), error = function(e) NULL) else NULL
    if (!is.null(index) && all(c("contrast_id", "output_id") %in% names(index))) {
      matched <- which(paste(index$contrast_id, index$output_id, sep = "_") == owner)
      if (length(matched) == 1L) exact_owner <- as.character(index$contrast_id[[matched]])
    }
  }
  source_dir <- file.path(plan$source_run, "report_pages", if (contrast) "contrast_evidence" else "evidence", owner, collection)
  source_index <- file.path(source_dir, "index.html")
  data_id <- if (contrast) "contrast-evidence-data" else "evidence-data"
  payload <- lisa_extension_payload_scope(lisa_extension_read_payload(source_index, data_id), exact_owner, collection,
    contrast, sort(unique(as.character(group$units$category_id))))
  if (is.null(payload)) return(NULL)
  evidence_dir <- file.path(report_root, if (contrast) "contrast_evidence" else "evidence", owner, collection)
  # Saved tables/figures are immutable presentation inputs.  Do not copy the
  # old HTML/assets: installed templates and JS/CSS are deliberately refreshed.
  selected_ids <- sort(unique(as.character(group$units$category_id)))
  for (part in c("tables", "figures")) if (dir.exists(file.path(source_dir, part))) {
    files <- list.files(file.path(source_dir, part), recursive = TRUE, full.names = TRUE)
    for (from in files[!dir.exists(files)]) {
      to <- file.path(evidence_dir, part, substring(from, nchar(file.path(source_dir, part)) + 2L))
      if (part == "tables" && grepl("[.]tsv$", from, ignore.case = TRUE)) {
        tab <- tryCatch(read_lisa_tsv(from), error = function(e) NULL)
        if (!is.null(tab) && "category_id" %in% names(tab)) {
          tab <- tab[as.character(tab$category_id) %in% selected_ids, , drop = FALSE]
          lisa_guarded_dir_create(dirname(to)); write_lisa_tsv(tab, to)
        } else lisa_extension_copy_file(from, to)
      } else lisa_extension_copy_file(from, to)
    }
  }
  recipe <- if (contrast) "reproduce_contrast_evidence.R" else "reproduce_category_evidence.R"
  if (file.exists(file.path(source_dir, recipe)))
    lisa_extension_copy_file(file.path(source_dir, recipe), file.path(evidence_dir, recipe))
  payload$category_products <- lisa_extension_reference_products(report_root, evidence_dir, products, inventory)
  template_root <- system.file(if (contrast) "contrast-evidence" else "category-evidence", package = "lisaR")
  if (!nzchar(template_root)) stop("Installed evidence template is missing.", call. = FALSE)
  for (extension in c("css", "js")) lisa_extension_copy_file(file.path(template_root, paste0("viewer.", extension)),
    file.path(evidence_dir, "assets", paste0(if (contrast) "contrast-evidence" else "category-evidence", ".", extension)))
  html <- paste(readLines(file.path(template_root, "viewer.html"), warn = FALSE), collapse = "\n")
  json <- gsub("<", "\\u003c", jsonlite::toJSON(payload, auto_unbox = TRUE, dataframe = "rows", na = "null", null = "null", digits = 17), fixed = TRUE)
  if (contrast) {
    html <- sub("<!-- CONTRAST_DATA -->", paste0('<script type="application/json" id="contrast-evidence-data">', json, "</script>"), html, fixed = TRUE)
    if (!file.exists(file.path(evidence_dir, recipe))) html <- sub(
      '<a href="reproduce_contrast_evidence.R" download>R script</a>', "", html, fixed = TRUE)
  } else {
    html <- sub("<!-- EVIDENCE_STYLE -->", '<link rel="stylesheet" href="assets/category-evidence.css">', html, fixed = TRUE)
    html <- sub("<!-- EVIDENCE_DATA -->", paste0('<script type="application/json" id="evidence-data">', json, "</script>"), html, fixed = TRUE)
    html <- sub("<!-- EVIDENCE_SCRIPT -->", '<script src="assets/category-evidence.js"></script>', html, fixed = TRUE)
  }
  # This gives the evidence page current shell assets without attempting to
  # derive or modify scientific evidence. Extension navigation is injected once
  # the complete inventory is known below.
  html <- lisa_present_evidence_html(html, evidence_dir, payload$metadata, if (contrast) "contrasts" else "analyses")
  lisa_guarded_write(file.path(evidence_dir, "index.html"), function(target) writeLines(html, target, useBytes = TRUE))
  file.path(if (contrast) "contrast_evidence" else "evidence", owner, collection, "index.html")
}

lisa_extension_write_index <- function(report_root, plan, inventory) {
  esc <- lisa_extension_html_escape
  paths <- gsub("\\\\", "/", as.character(inventory$path))
  invalid <- anyNA(paths) || any(!nzchar(paths)) || anyDuplicated(paths) ||
    any(grepl("^[A-Za-z]:|^/|[[:cntrl:]]", paths)) ||
    any(vapply(strsplit(paths, "/", fixed = TRUE), function(x) any(x %in% c("", ".", "..")), logical(1L)))
  if (invalid) stop("Extension index requires unique contained relative inventory paths.", call. = FALSE)
  groups <- lisa_extension_index_groups(plan, paths)
  # Build one selected evidence surface per exact context/collection before the
  # overview index.  Every product is attached to its selected category card;
  # the legacy gallery remains only when immutable source evidence is absent.
  scope_key <- function(group) paste(group$kind, group$owner, group$collection, sep = "\r")
  evidence_pages <- list()
  evidence_assets <- list()
  for (key in unique(vapply(groups, scope_key, character(1L)))) {
    scoped <- groups[vapply(groups, function(group) identical(scope_key(group), key), logical(1L))]
    if (!length(scoped) || !nrow(scoped[[1L]]$units)) next
    products <- unlist(lapply(scoped, function(group) lisa_extension_category_attachments(report_root, group)$category_products), recursive = FALSE)
    if (length(products)) {
      evidence_group <- scoped[[1L]]
      evidence_group$units <- unique(do.call(rbind, lapply(scoped, `[[`, "units")))
      page <- lisa_extension_write_saved_evidence(
        report_root, plan, evidence_group, products, inventory
      )
      if (!is.null(page)) {
        evidence_pages[[key]] <- page
        evidence_assets[[page]] <- sort(unique(unlist(lapply(products, `[[`, "assets"), use.names = FALSE)))
      }
    }
  }
  product_labels <- c(member_gene_sets = "Member gene sets", gene_cards = "GeneCards",
    volcano = "Volcano overlays", heatmap = "Gene heatmaps", kegg = "KEGG gene sets",
    contrast_profile = "Category profiles", contrast_heatmap = "Category heatmaps",
    saved_files = "Saved files")
  context_keys <- unique(vapply(groups, function(g) paste(g$kind, g$owner, sep = "\r"), character(1L)))
  contexts <- list(); body <- character()
  for (i in seq_along(context_keys)) {
    selected <- groups[vapply(groups, function(g) identical(paste(g$kind, g$owner, sep = "\r"), context_keys[[i]]), logical(1L))]
    kind <- selected[[1L]]$kind; owner <- selected[[1L]]$owner
    context_id <- paste0("extension-", kind, "-", substr(lisa_sha256_text(owner), 1L, 16L))
    label_key <- paste(kind, owner, sep = "\r")
    label <- if (!is.null(plan$context_labels[[label_key]])) {
      as.character(plan$context_labels[[label_key]])
    } else {
      lisa_presentation_label(
        plan$source_run, owner, contrast = identical(kind, "contrast")
      )
    }
    context <- list(id = context_id, scientific_id = owner, label = label, kind = kind, collections = list())
    context_html <- character()
    collections <- unique(vapply(selected, `[[`, character(1L), "collection"))
    for (j in seq_along(collections)) {
      collection <- collections[[j]]
      sections <- list(); collection_html <- character()
      products <- selected[vapply(selected, function(g) identical(g$collection, collection), logical(1L))]
      page_key <- scope_key(products[[1L]])
      evidence_page <- evidence_pages[[page_key]]
      id <- paste0("extension-section-", substr(lisa_sha256_text(paste(context_id, collection, "category-products", sep = "\r")), 1L, 16L))
      attachments <- lapply(products, function(group) lisa_extension_category_attachments(report_root, group))
      category_products <- unlist(lapply(attachments, `[[`, "category_products"), recursive = FALSE)
      category_products <- lapply(category_products, function(product) list(
          category_id = product$category_id, product = product$product, label = product$label,
          assets = lapply(product$assets, function(path) list(
            format = tolower(tools::file_ext(path)), href = lisa_extension_url_path(path)
          ))
      ))
      if (!is.null(evidence_page)) {
        # Evidence is the single category-specific surface.  The viewer payload
        # contains exact-ID canonical product references; do not duplicate them in a
        # second gallery.  Collection-level products can be added after this
        # section only when explicitly declared in a future source contract.
        sections[[1L]] <- list(id = id, label = "Evidence", kind = if (kind == "contrast") "contrast-evidence" else "category-evidence",
          href = evidence_page, category_products = category_products, global_products = list())
        collection_html <- c(collection_html, paste0('<section class="lisa-extension-section" id="', id,
          '"><h2>Evidence</h2><p>Inspect selected categories and their exact saved products in one evidence sheet.</p><a class="file-link" href="',
          esc(lisa_extension_url_path(evidence_page)), '">Open evidence</a></section>'))
      } else {
        # Compatible fallback for older immutable runs without saved evidence.
        sections[[1L]] <- list(id = id, label = "Saved category products", kind = "category-products",
          href = paste0("#", id), category_products = category_products, global_products = list())
        collection_html <- c(collection_html, paste0('<details class="lisa-extension-section" id="', id,
          '" open><summary>Saved category products</summary>',
          lisa_extension_gallery(report_root, products, id), '</details>'))
      }
      context$collections[[length(context$collections) + 1L]] <- list(id = collection, label = collection, sections = sections)
      scientific_context <- jsonlite::toJSON(list(analysis = label, collection = collection,
        exact_ids = stats::setNames(list(owner, collection), c(if (kind == "contrast") "contrast_id" else "analysis_id", "collection"))),
        auto_unbox = TRUE, null = "null")
      context_html <- c(context_html, paste0('<details class="lisa-extension-collection" data-lisa-nav-collection="', esc(collection),
        '" data-lisa-context-json="', esc(scientific_context), '"', if (i == 1L && j == 1L) ' open' else "", '><summary>',
        esc(collection), '</summary>', paste(collection_html, collapse = "\n"), '</details>'))
    }
    contexts[[length(contexts) + 1L]] <- context
    body <- c(body, paste0('<details class="lisa-extension-context" data-lisa-nav-context="', context_id, '"',
      if (i == 1L) ' open' else "", '><summary>', esc(if (kind == "contrast") "Contrast: " else "Analysis: "),
      esc(label), '</summary>', paste(context_html, collapse = "\n"), '</details>'))
  }
  navigation <- list(version = 1L, page = "extension", contexts = contexts)
  navigation_json <- jsonlite::toJSON(navigation, auto_unbox = TRUE, null = "null", digits = NA)
  navigation_json <- gsub("<", "\\u003c", navigation_json, fixed = TRUE)
  # Evidence pages share the same complete context/collection navigation.  A
  # nested page needs relative routes, whereas the root extension index does
  # not; calculate those paths rather than retaining stale source-run links.
  for (page in unname(unlist(evidence_pages, use.names = FALSE))) {
    page_file <- file.path(report_root, page); page_dir <- dirname(page_file)
    local_nav <- navigation
    for (i in seq_along(local_nav$contexts)) for (j in seq_along(local_nav$contexts[[i]]$collections)) {
      sections <- local_nav$contexts[[i]]$collections[[j]]$sections
      for (k in seq_along(sections)) {
        href <- sections[[k]]$href
        sections[[k]]$href <- if (startsWith(href, "#")) paste0(
          lisa_extension_url_path(lisa_presentation_relative(file.path(report_root, "index.html"), page_dir)), href) else
          lisa_extension_url_path(lisa_presentation_relative(file.path(report_root, href), page_dir))
        for (p in seq_along(sections[[k]]$category_products)) {
          assets <- sections[[k]]$category_products[[p]]$assets
          for (a in seq_along(assets)) assets[[a]]$href <- lisa_extension_url_path(lisa_presentation_relative(
            file.path(report_root, utils::URLdecode(assets[[a]]$href)), page_dir))
          sections[[k]]$category_products[[p]]$assets <- assets
        }
        if (identical(href, page)) local_nav$current <- list(context = local_nav$contexts[[i]]$id,
          collection = local_nav$contexts[[i]]$collections[[j]]$id, section = sections[[k]]$id)
      }
      local_nav$contexts[[i]]$collections[[j]]$sections <- sections
    }
    local_nav$page <- if (startsWith(page, "contrast_evidence/")) "contrasts" else "analyses"
    local_nav$mode <- "external"
    encoded <- gsub("<", "\\u003c", jsonlite::toJSON(local_nav, auto_unbox = TRUE, null = "null", digits = NA), fixed = TRUE)
    html <- paste(readLines(page_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
    html <- gsub('(?s)<script[^>]*id="lisa-report-navigation"[^>]*>.*?</script>', "", html, perl = TRUE)
    marker <- regexpr("(?i)</body\\s*>", html, perl = TRUE)
    if (marker[[1L]] < 1L) stop("Saved extension evidence page has no body.", call. = FALSE)
    html <- paste0(substr(html, 1L, marker[[1L]] - 1L), '<script id="lisa-report-navigation" type="application/json">',
      encoded, "</script>", substr(html, marker[[1L]], nchar(html)))
    lisa_guarded_write(page_file, function(target) writeLines(html, target, useBytes = TRUE))
  }
  lisa_extension_validate_product_references(report_root, evidence_assets, inventory)
  routes <- list(overview = "index.html", methods = "#extension-methods")
  for (kind in c("analysis", "contrast")) {
    found <- which(vapply(contexts, function(x) identical(x$kind, kind), logical(1L)))
    if (length(found)) routes[[if (kind == "analysis") "analyses" else "contrasts"]] <-
      contexts[[found[[1L]]]]$collections[[1L]]$sections[[1L]]$href
  }
  shell <- .lisa_report_shell(routes = routes, active = "overview",
    context = list(study = basename(plan$source_run), selection = paste(plan$mode, "figures")),
    assets = .lisa_copy_report_shell_assets(report_root, asset_subdir = "report_assets/lisa-shell"))
  gallery_script <- paste0('<script>(function(){"use strict";',
    'document.querySelectorAll("[data-extension-gallery]").forEach(function(root){',
    'var cards=Array.from(root.querySelectorAll(".lisa-extension-card")),search=root.querySelector("[data-extension-search]"),',
    'superSelect=root.querySelector("[data-extension-super]"),prev=root.querySelector("[data-extension-prev]"),',
    'next=root.querySelector("[data-extension-next]"),count=root.querySelector("[data-extension-count]"),page=0,size=12;',
    'function draw(reset){if(reset)page=0;var term=search.value.trim().toLowerCase(),sup=superSelect?superSelect.value:"";',
    'var found=cards.filter(function(card){return(!term||card.dataset.category.toLowerCase().indexOf(term)>=0)',
    '&&(!sup||card.dataset.supercategory===sup);});page=Math.max(0,Math.min(page,Math.ceil(found.length/size)-1));',
    'var start=page*size,visible=new Set(found.slice(start,start+size));cards.forEach(function(card){card.hidden=!visible.has(card);});',
    'prev.disabled=page===0;next.disabled=start+size>=found.length;',
    'count.textContent=found.length?"Showing "+(start+1)+"-"+Math.min(start+size,found.length)+" of "+found.length+" figures":"No matching figures";}',
    'search.addEventListener("input",function(){draw(true);});if(superSelect)superSelect.addEventListener("change",function(){draw(true);});',
    'prev.addEventListener("click",function(){page--;draw(false);});next.addEventListener("click",function(){page++;draw(false);});draw(true);',
    '});})();</script>')
  html <- paste0('<!doctype html><html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width,initial-scale=1"><title>LISA ', esc(plan$mode), ' report</title>',
    shell$head,
    '<style>.lisa-extension-context,.lisa-extension-collection,.lisa-extension-section{margin:14px 0;border:1px solid #d9dcdf;border-radius:8px;background:#fff;padding:12px}',
    '.lisa-shell-main{max-width:1450px;margin:auto;padding:24px}.lisa-extension-context>summary{font-size:1.22rem;font-weight:650}',
    '.lisa-extension-collection>summary,.lisa-extension-section>summary{font-size:1.05rem;font-weight:600}summary{cursor:pointer}',
    '.lisa-extension-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(min(100%,360px),1fr));gap:16px}',
    '.lisa-extension-card{min-width:0;border:1px solid #e3e5e7;border-radius:6px;padding:14px}.lisa-extension-card[hidden]{display:none}',
    '.lisa-extension-card h3{font-size:1rem;overflow-wrap:anywhere;margin:0 0 8px}.lisa-extension-super{font-size:.85rem;color:#50565e}',
    '.lisa-extension-card img{width:100%;height:310px;object-fit:contain;background:#fff}.lisa-extension-downloads{line-height:1.8;font-size:.9rem;margin-bottom:10px}',
    '.lisa-extension-controls{display:flex;align-items:end;flex-wrap:wrap;gap:12px;margin:16px 0}.lisa-extension-controls label{display:grid;gap:5px}',
    '.lisa-extension-controls input,.lisa-extension-controls select,.lisa-extension-controls button{font:inherit;min-height:38px;max-width:100%}',
    '.lisa-extension-controls [role=status]{font-size:.9rem;align-self:center}.lisa-extension-files{overflow-wrap:anywhere}',
    '.lisa-extension-section,#extension-methods{scroll-margin-top:var(--lisa-shell-height,240px)}.lisa-extension-empty{color:#50565e}',
    '@media(max-width:600px){.lisa-shell-main{padding:12px}.lisa-extension-context,.lisa-extension-collection,.lisa-extension-section{padding:8px}}',
    '</style></head><body class="', shell$body_class, '">', shell$header, shell$main_open,
    '<h1>LISA ', esc(plan$mode), ' report</h1><p>Choose an analysis or contrast, collection and figure type using the fixed navigation. ',
    'Search and supercategory filters below apply to the saved category figures; downloads always match their figure.</p>',
    paste(body, collapse = "\n"), '<section id="extension-methods"><h2>Data and methods</h2>',
    '<p>These figures were created from the verified saved results of <strong>', esc(basename(plan$source_run)),
    '</strong>. Enrichment was not repeated and the source run was not changed.</p>',
    '<p><a download href="extension_inventory.tsv">File inventory</a> &middot; ',
    '<a download href="extension_receipt.tsv">Report receipt</a></p></section>', shell$main_close,
    '<script id="lisa-report-navigation" type="application/json">', navigation_json, '</script>', gallery_script,
    '</body></html>')
  lisa_guarded_write(file.path(report_root, "navigation_inventory.json"),
    function(path) writeLines(navigation_json, path, useBytes = TRUE))
  lisa_guarded_write(file.path(report_root, "index.html"),
    function(path) writeLines(html, path, useBytes = TRUE))
  invisible(file.path(report_root, "index.html"))
}

#' Create additional selected or full report figures
#'
#' Creates the requested figures from the saved scientific tables. It does not
#' repeat enrichment or change the original run.
#'
#' @param exact_products Optional named list of exact heatmap variants, KEGG
#'   maps or native products used by the exploration interface. Leave NULL
#'   for the standard selected or full report.
#' @param source_run A completed and checked lisaR run.
#' @param selection `"full"`, a selection file or named list.
#' @param output_dir A new folder for the additional report.
#' @return Invisibly, a summary of the created and checked report.
#' @export
#'
#' @examples
#' # Rendering uses only a verified saved run and does not repeat enrichment.
#' # The offline demonstration run is kept interactive because it creates figures.
#' if (interactive()) {
#'   project <- lisa_init_project(tempfile("lisa-render-project-"))
#'   result <- run_lisa(file.path(project, "study.yml"))
#'   selection <- list(
#'     selection = list(
#'       mode = "selected",
#'       categories = "SYN_SIGNAL",
#'       products = "member_gene_sets"
#'     ),
#'     report = list(
#'       mode = "selected",
#'       formats = list(png = TRUE, svg = FALSE, pdf = FALSE),
#'       source_data = TRUE,
#'       recipes = FALSE
#'     )
#'   )
#'   rendered <- render_lisa_categories(
#'     result$output_dir,
#'     selection,
#'     tempfile("lisa-extension-example-")
#'   )
#'   rendered$validation_status
#' }
render_lisa_categories <- function(source_run, selection, output_dir,
                                   exact_products = NULL) {
  plan <- plan_lisa_extension(source_run, selection, exact_products)
  destination <- lisa_managed_destination(output_dir, create_parent = TRUE)
  output_dir <- destination$path
  if (lisa_path_within(output_dir, plan$source_run)) {
    stop("LISA-EXTENSION-020 output_dir must be outside the immutable source run.", call. = FALSE)
  }
  parent <- destination$parent
  if (lisa_path_entry_exists(output_dir)) stop("LISA-EXTENSION-021 output_dir already exists and extensions are immutable.", call. = FALSE)
  staging <- lisa_short_staging_path(output_dir, plan$extension_id)
  if (lisa_path_entry_exists(staging)) stop("LISA-EXTENSION-022 deterministic staging directory already exists.", call. = FALSE)
  staging <- lisa_run_root(staging)
  old_run_root <- getOption("lisaR.run_root", NULL)
  options(lisaR.run_root = staging)
  on.exit(options(lisaR.run_root = old_run_root), add = TRUE)
  promoted <- FALSE
  recovery_required <- FALSE
  on.exit(if (!promoted && !recovery_required &&
              lisa_path_entry_exists(staging)) {
    lisa_guarded_delete(staging, recursive = TRUE, run_root = NULL)
  }, add = TRUE)
  withCallingHandlers({
  work_root <- file.path(staging, "work"); report_root <- file.path(staging, "report")
  lisa_guarded_dir_create(work_root, staging)
  lisa_guarded_dir_create(report_root, staging)
  config_files <- c("config/de_index.tsv", "config/contrast_index.tsv")
  config_files <- config_files[file.exists(file.path(plan$source_run, config_files))]
  lisa_extension_stage_files(plan$source_run, work_root, config_files)
  package_dir <- lisa_resolve_package_dir()
  code_ledger <- lisa_write_code_identity_ledger(
    package_dir, lisa_extension_planned_scripts(plan$units),
    file.path(report_root, "code_identity.tsv"),
    executable_recipes = "reproduce_lisa_figure.R"
  )
  lisa_extension_render_units(plan, work_root, report_root, package_dir,
    code_ledger)
  lisa_assert_run_tree_safe(staging)
  lisa_extension_apply_file_policy(report_root, plan$report)
  lisa_guarded_delete(work_root, recursive = TRUE, run_root = staging)
  write_lisa_tsv(plan$policy, file.path(report_root, "output_policy.tsv"))
  artifact_root <- file.path(report_root, "artifacts")
  artifact_tree <- lisa_scan_run_tree(artifact_root)
  artifact_files <- sort(substring(
    artifact_tree$path[!artifact_tree$isdir], nchar(artifact_root) + 2L
  ))
  inventory <- data.frame(path = file.path("artifacts", artifact_files),
    bytes = as.numeric(file.info(file.path(report_root, "artifacts", artifact_files))$size),
    sha256 = vapply(file.path(report_root, "artifacts", artifact_files), lisa_sha256_file, character(1)),
    stringsAsFactors = FALSE)
  write_lisa_tsv(inventory, file.path(report_root, "extension_inventory.tsv"))
  lisa_extension_write_index(report_root, plan, inventory)
  receipt <- data.frame(extension_id = plan$extension_id, mode = plan$mode,
    source_run = basename(plan$source_run), source_manifest_hash = plan$source_manifest_hash,
    gsea_padj_cutoff = plan$gsea_padj_cutoff,
    selected_analyses = paste(sort(unique(plan$units$analysis_id[nzchar(plan$units$analysis_id)])), collapse = ";"),
    selected_contrasts = paste(sort(unique(plan$units$contrast_id[nzchar(plan$units$contrast_id)])), collapse = ";"),
    selected_categories = paste(plan$selected_categories, collapse = ";"),
    products = paste(sort(unique(plan$units$product)), collapse = ";"),
    created_at = as.character(file.info(file.path(plan$source_run, "run_manifest.tsv"))$mtime),
    warnings = plan$warning, validation_status = "PASS", stringsAsFactors = FALSE)
  write_lisa_tsv(receipt, file.path(report_root, "extension_receipt.tsv"))
  report_tree <- lisa_scan_run_tree(report_root)
  files <- substring(report_tree$path[!report_tree$isdir], nchar(report_root) + 2L)
  files <- sort(files)
  manifest <- data.frame(path = files,
    bytes = as.numeric(file.info(file.path(report_root, files))$size),
    sha256 = vapply(file.path(report_root, files), lisa_sha256_file, character(1)),
    stringsAsFactors = FALSE)
  write_lisa_tsv(manifest, file.path(report_root, "extension_manifest.tsv"))
  lisa_extension_validate_manifest(report_root)
  source_after <- verify_run(plan$source_run)
  if (!identical(source_after$gate, "PASS") ||
      !identical(lisa_sha256_file(file.path(plan$source_run, "run_manifest.tsv")), plan$source_manifest_hash)) {
    stop("LISA-EXTENSION-024 canonical source changed during extension rendering.", call. = FALSE)
  }
  # Re-read, reconcile, and re-hash immediately before the one-rename
  # promotion. This is intentionally separate from the first validation.
  lisa_extension_validate_manifest(report_root)
  promotion_error <- NULL
  promoted_ok <- tryCatch({
    lisa_promote_managed_directory(report_root, output_dir)
    TRUE
  }, error = function(error) {
    promotion_error <<- error
    FALSE
  })
  if (!promoted_ok) {
    if (inherits(promotion_error, "lisa_filesystem_recovery_error")) {
      stop(promotion_error)
    }
    stop("LISA-EXTENSION-025 atomic promotion failed: ",
         conditionMessage(promotion_error), call. = FALSE)
  }
  lisa_guarded_delete(staging, recursive = TRUE, run_root = NULL)
  promoted <- TRUE
  invisible(c(as.list(receipt[1, , drop = FALSE]), list(output_dir = output_dir,
    manifest = file.path(output_dir, "extension_manifest.tsv"))))
  }, lisa_filesystem_recovery_error = function(error) {
    recovery_required <<- TRUE
  })
}


# The FULL product core. It renders every unit of a plan into
# `<report_root>/artifacts/<owner>/<collection>/<product>/`. It is shared by the
# on-demand extension (render_lisa_categories(), which keeps its own gallery)
# and by run_lisa(mode = "full"), which calls it inside the run's own staging
# tree so the FULL products belong to the same run, manifest and report.
lisa_extension_render_units <- function(plan, work_root, report_root,
                                        package_dir, code_ledger,
                                        workers = 1L) {
  group_key <- paste(plan$units$unit_type, plan$units$analysis_id,
    plan$units$contrast_id, plan$units$collection, sep = "\r")
  keys <- sort(unique(group_key))
  if (workers <= 1L || length(keys) <= 1L) {
    for (key in keys) {
      unit <- plan$units[group_key == key, , drop = FALSE]
      lisa_extension_render_unit(
        unit, plan, work_root, report_root, package_dir, code_ledger
      )
    }
    return(invisible(TRUE))
  }
  # Parallel units write disjoint owner/collection trees. Shared inputs (for
  # example one expression matrix used by several collections) are staged once,
  # serially, before any worker starts, so no two workers stage the same file.
  path_columns <- grep("_path$", names(plan$units), value = TRUE)
  stage_paths <- unique(unlist(plan$units[path_columns], use.names = FALSE))
  lisa_extension_stage_files(plan$source_run, work_root,
    stage_paths[!is.na(stage_paths) & nzchar(stage_paths)])
  task_ids <- paste0("full_unit:", seq_along(keys))
  values <- lisa_pipeline_map(task_ids, function(task_id) {
    key <- keys[[match(task_id, task_ids)]]
    unit <- plan$units[group_key == key, , drop = FALSE]
    lisa_extension_render_unit(
      unit, plan, work_root, report_root, package_dir, code_ledger
    )
    TRUE
  }, workers = workers)
  invisible(values)
}

# run_lisa(mode = "full"): render the FULL products of a run into its OWN
# staging tree before the one report generator and the one promotion. Nothing
# is written outside `run_dir`; the scientific tables are only read.
lisa_render_run_full_products <- function(run_dir, report, workers = 1L) {
  run_dir <- normalizePath(run_dir, winslash = "/", mustWork = TRUE)
  managed_root <- getOption("lisaR.run_root", NULL)
  if (is.null(managed_root)) {
    stop("LISA-FULL-001 FULL products require the authorized run staging root.",
      call. = FALSE)
  }
  artifacts <- file.path(run_dir, "artifacts")
  work_root <- file.path(run_dir, "full_work")
  if (lisa_path_entry_exists(artifacts) || lisa_path_entry_exists(work_root)) {
    stop("LISA-FULL-002 FULL products were already rendered into this run.",
      call. = FALSE)
  }
  report <- lisa_validate_report_config(list(report = report),
    allow_selected = TRUE)
  old_in_run <- options(lisaR.full_products_in_run = TRUE)
  on.exit(options(old_in_run), add = TRUE)
  catalog <- lisa_extension_discover_catalog(run_dir)
  units <- lisa_extension_filter(catalog, list(mode = "full", filters = list()))
  plan <- list(extension_id = "run-full-products", mode = "full",
    source_run = run_dir, source_manifest_hash = "",
    gsea_padj_cutoff = lisa_extension_gsea_padj_cutoff(run_dir),
    selected_categories = sort(unique(units$category_id)), units = units,
    policy = lisa_output_policy_table(names(report$formats)[report$formats],
      report$source_data, report$recipes),
    report = report, warning = "")
  lisa_guarded_dir_create(work_root, managed_root)
  lisa_guarded_dir_create(artifacts, managed_root)
  config_files <- c("config/de_index.tsv", "config/contrast_index.tsv")
  config_files <- config_files[file.exists(file.path(run_dir, config_files))]
  lisa_extension_stage_files(run_dir, work_root, config_files)
  package_dir <- lisa_resolve_package_dir()
  code_ledger <- lisa_write_code_identity_ledger(
    package_dir, lisa_extension_planned_scripts(plan$units),
    file.path(artifacts, "code_identity.tsv"),
    executable_recipes = "reproduce_lisa_figure.R"
  )
  lisa_extension_render_units(plan, work_root, run_dir, package_dir,
    code_ledger, workers = workers)
  lisa_guarded_delete(work_root, recursive = TRUE, run_root = managed_root)
  # Only the image formats are enforced here: the integrated report attaches
  # every FULL figure to its category through its exact source sidecar, so the
  # sidecars stay with the figure (as in the standard figure cards).
  lisa_apply_format_policy(artifacts, report$formats)
  lisa_write_full_products_inventory(run_dir, units)
  invisible(list(artifacts = artifacts, units = units))
}

lisa_apply_format_policy <- function(root, formats) {
  if (!dir.exists(root)) return(invisible(FALSE))
  managed_root <- getOption("lisaR.run_root", NULL)
  tree <- lisa_scan_run_tree(root)
  files <- tree$path[!tree$isdir]
  remove <- rep(FALSE, length(files))
  for (format in c("png", "svg", "pdf")) {
    if (!isTRUE(formats[[format]])) {
      # The KEGG painter's preserved base diagram is an input witness, not a
      # figure format, and is kept whatever the requested formats are.
      remove <- remove | (grepl(paste0("[.]", format, "$"), files, ignore.case = TRUE) &
        !grepl("_kegg_base[.]png$", files))
    }
  }
  for (path in sort(files[remove])) {
    lisa_guarded_delete(path, run_root = managed_root)
  }
  invisible(TRUE)
}

lisa_write_full_products_inventory <- function(run_dir, units = NULL) {
  roots <- c("artifacts", file.path("outputs", "gene_level"))
  rows <- list()
  for (root in roots) {
    absolute <- file.path(run_dir, root)
    if (!dir.exists(absolute)) next
    tree <- lisa_scan_run_tree(absolute)
    files <- sort(tree$path[!tree$isdir])
    files <- files[basename(files) != "full_products_inventory.tsv"]
    if (!length(files)) next
    rows[[length(rows) + 1L]] <- data.frame(
      path = substring(files, nchar(run_dir) + 2L),
      bytes = as.numeric(file.info(files)$size),
      format = tolower(tools::file_ext(files)), stringsAsFactors = FALSE)
  }
  inventory <- if (length(rows)) do.call(rbind, rows) else
    data.frame(path = character(), bytes = numeric(), format = character())
  path <- file.path(run_dir, "artifacts", "full_products_inventory.tsv")
  if (!dir.exists(dirname(path))) lisa_guarded_dir_create(dirname(path))
  write_lisa_tsv(inventory, path)
  invisible(path)
}
