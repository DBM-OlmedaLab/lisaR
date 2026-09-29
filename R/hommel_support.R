# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé. GPL-3 or later.

# Presentation only. N is the saved raw-P-evaluable count, never the viewer's
# NES/FDR-evaluable count or its individually significant-set count.
lisa_hommel_attach <- function(object, categories) {
  m <- object$metadata
  if (!inherits(object, "lisa_category_navigation") && !inherits(object, "lisa_category_evidence"))
    stop("Hommel support attaches only to single-analysis category views.", call. = FALSE)
  if (!is.null(m$contrast_id) || !is.null(m$analysis_a) || !is.null(m$analysis_b) ||
      any(c("mean_NES_A", "mean_NES_B", "delta_mean_NES") %in% names(object$categories)))
    stop("Hommel support is not inference on an analysis contrast.", call. = FALSE)
  categories <- lisa_category_support_presentation(categories)
  for (field in c("analysis_id", "collection")) {
    if (length(m[[field]]) != 1L || is.na(m[[field]]) || !nzchar(m[[field]]) ||
        any(categories[[field]] != m[[field]]))
      stop("Hommel support ", field, " scope mismatch.", call. = FALSE)
    if (field %in% names(object$categories) && any(is.na(object$categories[[field]]) |
        object$categories[[field]] != m[[field]]))
      stop("Category view ", field, " scope mismatch.", call. = FALSE)
    object$categories[[field]] <- rep(m[[field]], nrow(object$categories))
  }
  key <- lisa_support_key(object$categories, unique = TRUE)
  ix <- match(key, lisa_support_key(categories, unique = TRUE))
  if (anyNA(ix)) stop("Category view has no exact Hommel analysis/collection/category result.", call. = FALSE)
  fields <- c("category_p_adjusted", "minimum_enriched_sets", "minimum_enriched_pct",
    "n_sets_evaluable", "n_sets_total", "n_sets_missing_p", "n_sets_not_eligible_or_absent",
    "confidence_level", "significant", "status", names(lisa_support_grade(1, 1)),
    "observed_NES_pattern", "observed_NES_counts")
  lisa_require_columns(categories, fields, "saved Hommel category results")
  for (field in fields) object$categories[[field]] <- categories[[field]][ix]
  object$metadata$hommel_support_schema <- "hommel-support-ui-v1"
  object$metadata$hommel_support_scope <- "within_analysis_only; not an ON-PRE difference test"
  object$metadata$hommel_support_method <- "robust Hommel closed testing; hommel::hommel(p, simes=FALSE)"
  object
}

lisa_hommel_read_results <- function(collection_dir, analysis_id, collection) {
  if (identical(collection, "HALLMARKS")) return(NULL)
  files <- list.files(file.path(collection_dir, "lisa_tables"),
    pattern = "_LISA_category_inference[.]tsv$", full.names = TRUE)
  if (!length(files)) return(NULL)
  if (length(files) != 1L) stop("Ambiguous saved Hommel results.", call. = FALSE)
  x <- lisa_category_support_presentation(read_lisa_tsv(files[[1L]]))
  if (any(x$analysis_id != analysis_id | x$collection != collection))
    stop("Saved Hommel results have a different analysis or collection.", call. = FALSE)
  x
}

lisa_hommel_help_button <- function(label = "How to read category significance") {
  paste0('<button type="button" class="hommel-help-button" data-hommel-help aria-haspopup="dialog" ',
    'aria-controls="hommel-support-help">', lisa_navigation_escape(label), '</button>')
}

lisa_hommel_reference_html <- function() {
  paste0('<p class="hommel-method-reference">Method: <strong>robust Hommel closed testing</strong>. ',
    '<a href="https://doi.org/10.1093/biomet/asz041">Goeman et al. (2019), doi:10.1093/biomet/asz041</a>.</p>')
}

lisa_hommel_help_html <- function() {
  paste0('<!-- HOMMEL_HELP_START --><dialog id="hommel-support-help" class="hommel-support-help" ',
    'aria-labelledby="hommel-help-title"><div class="hommel-help-heading"><h2 id="hommel-help-title">Category significance and support</h2>',
    '<button type="button" data-hommel-close autofocus aria-label="Close category significance help">Close</button></div>',
    '<p><strong>Hommel is a multiple-testing correction</strong>: robust Hommel closed testing accounts for the many gene sets assessed together. ',
    'It limits false positives by controlling the chance of even one false positive within one analysis ',
    '(familywise error), assuming valid input P values.</p>',
    '<h3>What does the adjusted category P mean?</h3>',
    '<p>The <strong>adjusted category P produced by Hommel</strong> tests whether ',
    '<strong>none of the evaluable gene sets in this category is enriched</strong>. ',
    'A small adjusted P gives evidence that <strong>at least one</strong> is enriched. ',
    'The stars on the figure use this adjusted P: * for P &le; 0.05, ** for P &le; 0.01, ',
    '*** for P &le; 0.001, and **** for P &le; 0.0001. Only the strongest level reached is shown. ',
    'No stars means P &gt; 0.05 or not evaluable; the category card distinguishes these cases.</p>',
    '<h3>What does minimum support mean?</h3>',
    '<p>The category card also shows <strong>d/N</strong>: among N eligible gene sets with a usable raw P value, ',
    'at least d can be supported as enriched with <strong>95% simultaneous confidence</strong>. ',
    'Here is an invented example, not a study result: <strong>5/20</strong> means ',
    '<strong>at least 5 of the 20 sets</strong>, not exactly 5; more may be enriched. ',
    'The percentage is the same minimum expressed as a fraction, not an additional test. ',
    'If d=0, no positive minimum is supported; this does not mean no enrichment. ',
    'If N=0, the category cannot be evaluated. Check evaluated / total coverage in the card.</p>',
    '<p>The adjusted P and d/N answer different questions: evidence for at least one enriched set versus ',
    'a simultaneous lower bound on how many are enriched. Neither measures effect size, proves the whole biological ',
    'process is activated, establishes direction, or tests a difference between conditions. ',
    'The separate NES panel describes observed directions among individually significant sets. ',
    'Overlapping gene sets are not independent mechanisms.</p>',
    lisa_hommel_reference_html(),
    '</dialog><!-- HOMMEL_HELP_END -->')
}

# Presentation only: the saved adjusted P is never recomputed or rounded for
# threshold decisions. Missing P receives no mark, as does P > 0.05.
lisa_category_significance_stars <- function(p) {
  if (!is.numeric(p) || any(!is.na(p) & (!is.finite(p) | p < 0 | p > 1)))
    stop("Adjusted category P must be numeric in [0,1] or missing.", call. = FALSE)
  stars <- rep("", length(p))
  stars[!is.na(p) & p <= .05] <- "*"
  stars[!is.na(p) & p <= .01] <- "**"
  stars[!is.na(p) & p <= .001] <- "***"
  stars[!is.na(p) & p <= .0001] <- "****"
  stars
}

lisa_hommel_install_help <- function(html) {
  html <- paste(html, collapse = "\n")
  html <- gsub('(?s)<!-- HOMMEL_HELP_START -->.*?<!-- HOMMEL_HELP_END -->', "", html, perl = TRUE)
  at <- regexpr("(?i)</body\\s*>", html, perl = TRUE)[[1L]]
  if (at < 1L) stop("Hommel help requires a complete HTML page.", call. = FALSE)
  paste0(substr(html, 1L, at - 1L), lisa_hommel_help_html(), substr(html, at, nchar(html)))
}

lisa_hommel_lollipop_source <- function(source, categories, analysis_id, collection) {
  if (!identical(unique(as.character(source$annotation_variant)), "plain"))
    stop("Hommel support view requires the original plain lollipop source.", call. = FALSE)
  x <- lisa_support_lollipop_source(source, categories, analysis_id, collection)
  x$annotation_variant <- "hommel_support_v1"
  x$hommel_support_schema <- "hommel-support-ui-v1"
  x$figure_width <- 13
  x$figure_height <- max(5.6, 2.8 + .32 * nrow(x))
  x
}

lisa_plot_hommel_support <- function(summary_df, group_by_supracategory, plot_order) {
  lisa_assert_classified_plot_rows(summary_df, "Hommel support lollipop")
  x <- prepare_category_plot_df(lisa_category_support_presentation(summary_df), plot_order)
  wrap_label <- function(value, width) vapply(as.character(value), function(label)
    paste(strwrap(label, width = width), collapse = "\n"), character(1L))
  x$key <- lisa_support_key(x, unique = TRUE)
  x$category_plot_label <- factor(x$key, levels = rev(x$key))
  labels <- stats::setNames(wrap_label(x$category_display_name, 42L), x$key)
  panels <- c("Category significance", "Mean NES\n(individually significant sets)")
  support <- x; support$panel <- factor(panels[[1L]], levels = panels)
  support$hommel_text <- lisa_category_significance_stars(support$category_p_adjusted)
  geometry <- x[is.finite(x$mean_NES) & is.finite(x$n_genesets) & x$n_genesets > 0, , drop = FALSE]
  geometry$panel <- factor(rep(panels[[2L]], nrow(geometry)), levels = panels)
  rows <- x; rows$panel <- factor(panels[[2L]], levels = panels)
  p <- ggplot2::ggplot() +
    ggplot2::geom_blank(data = support, ggplot2::aes(x = 0, y = category_plot_label)) +
    ggplot2::geom_blank(data = support, ggplot2::aes(x = 30, y = category_plot_label)) +
    ggplot2::geom_blank(data = rows, ggplot2::aes(x = 0, y = category_plot_label)) +
    ggplot2::geom_text(data = support, ggplot2::aes(x = 4, y = category_plot_label, label = hommel_text),
      hjust = 0, color = "#172b3a", size = 6) +
    ggplot2::geom_vline(data = rows, ggplot2::aes(xintercept = 0), linetype = "dashed", color = "grey55", linewidth = .35) +
    ggplot2::geom_segment(data = geometry, ggplot2::aes(x = 0, xend = mean_NES, y = category_plot_label, yend = category_plot_label),
      color = "grey70", linewidth = .35) +
    ggplot2::geom_point(data = geometry, ggplot2::aes(x = mean_NES, y = category_plot_label, size = n_genesets, fill = color),
      alpha = .95, shape = 21, color = "grey25", stroke = .25) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_size_continuous(range = c(1.5, 7), breaks = c(0, 1, 5, 10, 20)) +
    ggplot2::scale_y_discrete(labels = labels, drop = TRUE) +
    ggplot2::scale_x_continuous(breaks = function(limits) pretty(limits, n = 5),
      expand = ggplot2::expansion(mult = c(.08, .08))) +
    ggplot2::labs(x = "Mean NES (right panel)", y = NULL, size = "n gene sets",
      subtitle = "Adjusted category P \u00b7 robust Hommel multiple-testing correction",
      caption = paste0("Stars: * P \u2264 0.05; ** P \u2264 0.01; *** P \u2264 0.001; **** P \u2264 0.0001 (adjusted category P). Blank: P > 0.05 or not evaluable; see category card.\n",
        "A significant category has evidence for at least one enriched set; stars do not measure effect size or extent.\n",
        "Minimum supported enrichment d/N, with 95% simultaneous confidence, is shown in Category evidence.\n",
        "Method: robust Hommel closed testing. Goeman et al. (2019), doi:10.1093/biomet/asz041.")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "right", panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold", size = 12),
      axis.text.y = ggplot2::element_text(size = 11, lineheight = .95),
      plot.caption = ggplot2::element_text(hjust = 0, size = 10),
      plot.title.position = "plot", plot.caption.position = "plot",
      plot.background = ggplot2::element_rect(fill = "white", color = NA))
  if (isTRUE(group_by_supracategory) && "macrogroup_name" %in% names(x)) {
    p <- p + ggplot2::facet_grid(macrogroup_name ~ panel, scales = "free", space = "free_y", switch = "y",
      labeller = ggplot2::labeller(macrogroup_name = function(value) wrap_label(value, 26L))) +
      ggplot2::theme(strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 10),
        panel.spacing.y = grid::unit(.25, "lines"))
  } else p <- p + ggplot2::facet_grid(. ~ panel, scales = "free_x")
  # The stars panel has layout coordinates, not a numerical measurement.
  # Select its scale by facet identity, never by a numeric range that NES
  # values might also reach. Retain the ordinary NES breaks and gridlines.
  original_facet <- p$facet
  p$facet <- ggplot2::ggproto(NULL, original_facet,
    init_scales = function(layout, x_scale = NULL, y_scale = NULL, params) {
      scales <- ggplot2::FacetGrid$init_scales(layout, x_scale, y_scale, params)
      for (i in unique(layout$SCALE_X[as.character(layout$panel) == "Category significance"])) {
        if (!is.null(scales$x)) {
          scales$x[[i]]$breaks <- NULL
          scales$x[[i]]$minor_breaks <- NULL
        }
      }
      scales
    })
  p
}

# Refresh only presentation metadata and HTML assets. Existing set/gene data,
# matrix figures, exports and enrichment results are never recalculated.
lisa_refresh_hommel_evidence_page <- function(page, categories) {
  if (!file.exists(page)) return(invisible(FALSE))
  html <- paste(readLines(page, warn = FALSE), collapse = "\n")
  marker <- '<script type="application/json" id="evidence-data">'
  start <- regexpr(marker, html, fixed = TRUE)[[1L]]
  if (start < 1L) stop("Category evidence has no exact data payload.", call. = FALSE)
  start <- start + nchar(marker)
  tail <- substr(html, start, nchar(html))
  end <- regexpr("</script>", tail, fixed = TRUE)[[1L]]
  if (end < 1L) stop("Category evidence payload is unterminated.", call. = FALSE)
  raw <- substr(tail, 1L, end - 1L)
  payload <- jsonlite::fromJSON(raw, simplifyVector = FALSE)
  flat <- jsonlite::fromJSON(raw, simplifyVector = TRUE)
  view <- structure(list(metadata = payload$metadata, categories = flat$categories), class = "lisa_category_evidence")
  view <- lisa_hommel_attach(view, categories)
  payload$metadata <- view$metadata
  payload$categories <- view$categories
  json <- jsonlite::toJSON(payload, auto_unbox = TRUE, dataframe = "rows", na = "null", null = "null", digits = 17)
  json <- gsub("<", "\\u003c", json, fixed = TRUE)
  html <- paste0(substr(html, 1L, start - 1L), json, substr(tail, end, nchar(tail)))
  html <- lisa_hommel_install_help(html)
  directory <- dirname(page)
  lisa_write_inference_tsv(view$categories, file.path(directory, "tables", "categories.tsv"))
  write_lisa_tsv(data.frame(key = names(view$metadata), value = vapply(view$metadata,
    lisa_evidence_text, character(1L)), stringsAsFactors = FALSE), file.path(directory, "tables", "metadata.tsv"))
  asset_dir <- file.path(directory, "assets"); lisa_guarded_dir_create(asset_dir)
  for (extension in c("css", "js")) {
    from <- system.file("category-evidence", paste0("viewer.", extension), package = "lisaR")
    if (!nzchar(from) || !file.copy(from, file.path(asset_dir, paste0("category-evidence.", extension)), overwrite = TRUE))
      stop("Cannot refresh verified category-evidence assets.", call. = FALSE)
  }
  .lisa_copy_report_shell_assets(directory)
  lisa_guarded_write(page, function(target) writeLines(html, target, useBytes = TRUE))
  invisible(TRUE)
}

lisa_refresh_hommel_views <- function(run_dir, analysis_id, collection, categories) {
  navigation_dir <- file.path(run_dir, "report_pages", "category_navigation", analysis_id, collection)
  metadata_path <- file.path(navigation_dir, "metadata.json")
  if (file.exists(metadata_path)) {
    navigation <- structure(list(metadata = jsonlite::read_json(metadata_path, simplifyVector = TRUE),
      categories = read_lisa_tsv(file.path(navigation_dir, "tables", "category_navigation_source.tsv"))),
      class = "lisa_category_navigation")
    navigation <- lisa_hommel_attach(navigation, categories)
    render_lisa_category_navigation(navigation, navigation_dir)
  }
  pages <- c(file.path(run_dir, "report_pages", "evidence", analysis_id, collection, "index.html"),
    file.path(run_dir, "outputs", "gene_level", "single_de", analysis_id, paste0("collection_", collection), "category_evidence", "index.html"))
  for (page in pages) lisa_refresh_hommel_evidence_page(page, categories)
  invisible(TRUE)
}
