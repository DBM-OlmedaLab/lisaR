# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé. GPL-3 or later.

lisa_support_save_lollipop <- function(source, path) {
  rendered <- lisa_lollipop_from_source(source)
  ext <- tolower(tools::file_ext(path))
  if (ext == "png") grDevices::png(path, width = rendered$width, height = rendered$height,
    units = "in", res = rendered$dpi, type = "cairo", bg = rendered$background) else
  if (ext == "svg") grDevices::svg(path, width = rendered$width, height = rendered$height, bg = rendered$background) else
  if (ext == "pdf") grDevices::cairo_pdf(path, width = rendered$width, height = rendered$height, bg = rendered$background) else
    stop("Support lollipops support png, svg and pdf.", call. = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(rendered$plot)
  invisible(path)
}

lisa_write_support_lollipops <- function(directory, categories, analysis_id, collection, plot_formats) {
  sources <- list.files(file.path(directory, "plots"), pattern = "_GSEA_lollipop(_direction_stats)?_source[.]tsv$", full.names = TRUE)
  written <- character()
  for (path in sources) {
    source <- lisa_read_figure_source_tsv(path)
    # Refuse incomplete legacy geometry/metadata rather than changing the
    # historical meanings or silently guessing plotting options.
    lisa_lollipop_from_source(source)
    decorated <- lisa_support_lollipop_source(source, categories, analysis_id, collection)
    stem <- sub("_source[.]tsv$", "_support", path)
    # Existing support-v1 files belong to the previous presentation. Preserve
    # their bytes when refreshing a run that already contains those variants.
    if (!file.exists(paste0(stem, "_source.tsv"))) lisa_write_inference_tsv(decorated, paste0(stem, "_source.tsv"))
    for (format in plot_formats) if (!file.exists(paste0(stem, ".", format)))
      lisa_support_save_lollipop(lisa_read_figure_source_tsv(paste0(stem, "_source.tsv")), paste0(stem, ".", format))
    written <- c(written, paste0(stem, "_source.tsv"))
    if (identical(unique(as.character(source$annotation_variant)), "plain")) {
      hommel <- lisa_hommel_lollipop_source(source, categories, analysis_id, collection)
      new_stem <- sub("_source[.]tsv$", "_hommel_support", path)
      lisa_write_inference_tsv(hommel, paste0(new_stem, "_source.tsv"))
      for (format in plot_formats) lisa_support_save_lollipop(hommel, paste0(new_stem, ".", format))
      written <- c(written, paste0(new_stem, "_source.tsv"))
    }
  }
  invisible(written)
}

lisa_support_presentation_metadata <- function() {
  list(schema_version = "1.0", ui_schema_version = "hommel-support-ui-v1", presentation_only = TRUE,
    bands = "S1:(0,10); S2:[10,25); S3:[25,50); S4:[50,75); S5:[75,100] percent",
    assignment = "Exact integer cross-products d*100 versus N*cutoff; never rounded percentages",
    percent_format = "Conservative decimal long division, two places extended until a positive decimal is visible; trailing zeros removed",
    zeros = "No grade/graph annotation for d=0 or N=0; separate table states",
    direction = "Observed NES signs among evaluable members; not a directional guarantee for the supported minimum",
    lollipop = "Original mean_NES/n_genesets/dictionary colours and significant-member same_direction_pct unchanged; neutral support column from the complete inference family",
    confidence = .95, helper_sha256 = lisa_cache_digest(body(lisa_support_grade)),
    renderer_sha256 = lisa_cache_digest(list(body(lisa_plot_support_lollipop), body(lisa_plot_category_inference), body(lisa_plot_hommel_support))),
    inference_recomputed = FALSE, difference_grades = FALSE, HALLMARKS_excluded = TRUE)
}

#' Refresh minimum-support presentation in a new copy of an inferred LISA run
#'
#' Copies a manifest-verified run and decorates its already calculated category
#' inference. No DE, fgsea, closed tests, P values, lower bounds or families are
#' recalculated. Historical lollipop sources and figures remain unchanged.
#' Adds a distinct Hommel support view and refreshes exact-key navigator and
#' category-evidence metadata. Repeated refreshes preserve previous variants.
#' Rebuild the new copy with the installed `scripts/build_LISA_report.R` to
#' create the STANDARD or FULL report and seal its new manifest.
#' @param source_run Complete run with saved category inference and manifest.
#' @param output_dir New destination outside `source_run`; must not exist.
#' @param plot_formats Any of `png`, `svg`, `pdf`.
#' @return Invisibly, the new run directory.
#' @export
refresh_lisa_support_grades <- function(source_run, output_dir,
    plot_formats = c("png", "svg", "pdf")) {
  source_run <- lisa_assert_run_tree_safe(source_run)
  proposed <- lisa_path_canonical(output_dir)
  if (identical(proposed, source_run) || startsWith(proposed, paste0(source_run, "/")))
    stop("Support destination must be outside the source run.", call. = FALSE)
  if (file.exists(output_dir) || dir.exists(output_dir)) stop("Support destination must be new.", call. = FALSE)
  if (any(!plot_formats %in% c("png", "svg", "pdf"))) stop("Invalid support plot formats.", call. = FALSE)
  manifest_path <- file.path(source_run, "run_manifest.tsv")
  if (!file.exists(manifest_path)) stop("Support refresh requires a complete source manifest.", call. = FALSE)
  source_check <- lisa_validate_built_run_manifest(source_run, read_lisa_tsv(manifest_path))
  if (!identical(source_check$gate, "PASS"))
    stop("Source run failed integrity verification: ", paste(source_check$findings, collapse = "; "), call. = FALSE)
  family_root <- file.path(source_run, "outputs", "category_inference")
  tables <- list.files(family_root, pattern = "^category_results[.]tsv$", recursive = TRUE, full.names = TRUE)
  if (!length(tables)) stop("Support refresh requires existing category inference; it never recomputes inference.", call. = FALSE)
  # Validate the scientific tables before allocating/copying a report tree.
  for (path in tables) lisa_category_support_presentation(read_lisa_tsv(path))
  lisa_managed_destination(output_dir, create_parent = TRUE)
  fs::dir_copy(source_run, output_dir)
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  archive_base <- file.path(output_dir, "source_provenance", "hommel_support_ui_v1")
  archive <- archive_base
  generation <- 1L
  while (dir.exists(archive)) { generation <- generation + 1L; archive <- paste0(archive_base, "_", generation) }
  lisa_guarded_dir_create(archive, output_dir)
  if (!file.rename(file.path(output_dir, "run_manifest.tsv"), file.path(archive, "parent_run_manifest.tsv")))
    stop("Cannot archive parent manifest before refreshing support presentation.", call. = FALSE)
  prior_receipt <- file.path(output_dir, "support_grades_derivation.tsv")
  if (file.exists(prior_receipt) && !file.copy(prior_receipt, file.path(archive, "parent_support_grades_derivation.tsv")))
    stop("Cannot archive parent support receipt.", call. = FALSE)
  lisa_rebase_run_text_paths(output_dir, source_run, output_dir)
  # Builder caches figure-source exports by deterministic path. Discard only
  # copied derived report surfaces so a fresh build cannot reuse old recipes.
  for (name in c("report_media", "report_files", "report_figure_data")) {
    path <- file.path(output_dir, name)
    if (dir.exists(path)) fs::dir_delete(path)
  }
  for (old_table in tables) {
    relative <- substring(old_table, nchar(source_run) + 2L)
    path <- file.path(output_dir, relative)
    original <- read_lisa_tsv(path)
    table <- lisa_category_support_presentation(original)
    analysis <- unique(as.character(table$analysis_id))
    if (length(analysis) != 1L) stop("A category family must have one analysis identity.", call. = FALSE)
    file.copy(path, file.path(archive, paste0(analysis, "_category_results.tsv")))
    prior_metadata <- file.path(dirname(path), "support_presentation.json")
    if (file.exists(prior_metadata) && !file.copy(prior_metadata, file.path(archive, paste0(analysis, "_support_presentation.json"))))
      stop("Cannot archive parent support presentation metadata.", call. = FALSE)
    lisa_write_inference_tsv(table, path)
    for (collection in unique(table$collection)) {
      directory <- file.path(output_dir, "outputs", "single_de", analysis, paste0("collection_", collection))
      selected <- table[table$collection == collection, , drop = FALSE]
      target <- lisa_inference_one_file(file.path(directory, "lisa_tables"), "_LISA_category_inference[.]tsv$")
      lisa_write_inference_tsv(selected, target)
      source_path <- lisa_inference_one_file(file.path(directory, "plots"), "_LISA_category_inference_source[.]tsv$")
      source <- lisa_read_figure_source_tsv(source_path)
      member_key <- lisa_support_key(source)
      selected_key <- lisa_support_key(selected, unique = TRUE)
      if (any(!member_key %in% selected_key)) stop("Inference figure contains an unknown category key.", call. = FALSE)
      # Preserve the exact prior source/figures under a versioned stem. A
      # second refresh starts from the S-grade presentation, not the original
      # inference figure, so it needs its own historical identity.
      stem <- sub("_source[.]tsv$", "", source_path)
      for (old in c(source_path, paste0(stem, ".", c("png", "svg", "pdf")))) {
        if (file.exists(old)) {
          suffix <- if ("support_schema_version" %in% names(source))
            "_LISA_category_inference_historical_support_v1" else "_LISA_category_inference_historical_v1"
          historical <- sub("_LISA_category_inference", suffix, old, fixed = TRUE)
          if (!file.exists(historical)) file.copy(old, historical)
        }
      }
      members <- source[!is.na(source$pathway) & nzchar(source$pathway),
        c("analysis_id", "collection", "category_id", "pathway", "NES", "pval", "hypothesis_id"), drop = FALSE]
      members$included_in_family <- TRUE
      members$raw_p_available <- TRUE
      refreshed <- lisa_category_inference_source(selected, members, source$figure_title[[1L]])
      # Refresh presentation only, including existing saved figure formats;
      # source inferential values are copied unchanged from selected.
      lisa_write_inference_tsv(refreshed, source_path)
      for (format in plot_formats)
        lisa_category_inference_save_plot(refreshed, paste0(stem, ".", format))
      lisa_write_support_lollipops(directory, selected, analysis, collection, plot_formats)
      lisa_refresh_hommel_views(output_dir, analysis, collection, selected)
    }
    # Inferential method.json, raw_results.tsv and hypothesis_family.tsv stay
    # byte-for-byte unchanged; presentation provenance is deliberately separate.
    jsonlite::write_json(lisa_support_presentation_metadata(), file.path(dirname(path), "support_presentation.json"),
      auto_unbox = TRUE, pretty = TRUE, digits = NA)
  }
  receipt <- data.frame(operation = "refresh_support_grades_without_recomputing_inference",
    status_at_creation = "support_added_report_rebuild_pending", source_integrity = source_check$gate,
    source_run = basename(source_run), source_manifest_sha256 = lisa_sha256_file(manifest_path),
    original_DE_GSEA_recomputed = FALSE, inference_recomputed = FALSE,
    support_schema_version = "1.0", ui_schema_version = "hommel-support-ui-v1", stringsAsFactors = FALSE)
  write_lisa_tsv(receipt, file.path(output_dir, "support_grades_derivation.tsv"))
  invisible(output_dir)
}

lisa_finalize_support_grades_derivation <- function(run_dir) {
  path <- file.path(run_dir, "support_grades_derivation.tsv")
  if (!file.exists(path)) return(invisible(FALSE))
  receipt <- read_lisa_tsv(path)
  if (nrow(receipt) != 1L || receipt$operation != "refresh_support_grades_without_recomputing_inference")
    stop("Invalid support presentation derivation receipt.", call. = FALSE)
  finished <- data.frame(operation = receipt$operation, status = "report_rebuilt_and_resealed",
    original_DE_GSEA_recomputed = FALSE, inference_recomputed = FALSE,
    support_schema_version = "1.0", lisaR_version = as.character(utils::packageVersion("lisaR")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE), stringsAsFactors = FALSE)
  write_lisa_tsv(finished, file.path(run_dir, "audit", "support_grades_report_receipt.tsv"))
  manifest <- lisa_run_manifest(run_dir, workers = 1L)
  write_lisa_tsv(manifest, file.path(run_dir, "run_manifest.tsv"))
  check <- lisa_validate_built_run_manifest(run_dir, manifest)
  if (!identical(check$gate, "PASS")) stop("Support report failed final manifest verification.", call. = FALSE)
  invisible(TRUE)
}
