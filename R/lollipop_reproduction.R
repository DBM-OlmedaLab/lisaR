# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Lollipop sources retain the complete plotting contract. Reproduction calls
# the same renderers as run_LISA_DE instead of approximating them with bars.
lisa_write_lollipop_figure_source <- function(summary_df, path,
    group_by_supracategory, plot_order, display_title = NULL,
    comparison_subtitle = NULL, positive_direction = NULL,
    annotation_variant = "plain", figure_width = 11, figure_height = 12,
    figure_dpi = 300, plot_bg = "white") {
  lisa_assert_classified_plot_rows(summary_df, "lollipop source")
  if (!nrow(summary_df)) stop("Lollipop source requires category rows.", call. = FALSE)
  x <- summary_df
  text <- function(value) {
    if (is.null(value) || !length(value) || is.na(value[[1L]]) || !nzchar(value[[1L]])) "%00" else
      lisa_encode_figure_source_text(as.character(value[[1L]]))
  }
  x$figure_type <- "lisa_lollipop"
  x$source_row_order <- seq_len(nrow(x))
  x$selected_for_plot <- is.finite(x$n_genesets) & x$n_genesets > 0
  x$group_by_supracategory <- group_by_supracategory
  x$plot_order <- plot_order
  x$annotation_variant <- annotation_variant
  x$figure_title_encoded <- text(display_title)
  x$figure_subtitle_encoded <- text(comparison_subtitle)
  x$positive_direction_encoded <- text(positive_direction)
  x$figure_width <- figure_width
  x$figure_height <- figure_height
  x$figure_dpi <- figure_dpi
  x$figure_background <- plot_bg
  write_tsv_local(x, path)
  invisible(x)
}

lisa_lollipop_from_source <- function(x) {
  required <- c("category_id", "category_display_name", "mean_NES", "n_genesets", "color")
  if (!nrow(x) || !all(required %in% names(x))) {
    stop("Lollipop source lacks category geometry.", call. = FALSE)
  }
  constant <- function(field, empty = FALSE) {
    values <- as.character(x[[field]])
    if (!field %in% names(x) || length(values) != nrow(x) || anyNA(values) ||
        length(unique(values)) != 1L || (!empty && any(!nzchar(values)))) {
      stop("Lollipop source lacks consistent metadata: ", field,
        ". Regenerate source data with the corrected package.", call. = FALSE)
    }
    values[[1L]]
  }
  flag <- tolower(constant("group_by_supracategory"))
  if (!flag %in% c("true", "false")) stop("Invalid lollipop grouping flag.", call. = FALSE)
  order <- constant("plot_order")
  if (!order %in% c("supracategory", "fixed", "mean_NES", "n_genesets", "consistency")) {
    stop("Invalid lollipop plot order.", call. = FALSE)
  }
  variant <- constant("annotation_variant")
  if (!variant %in% c("plain", "direction", "support_v1", "direction_support_v1", "hommel_support_v1")) stop("Invalid lollipop variant.", call. = FALSE)
  decode <- function(field) {
    value <- constant(field, empty = TRUE)
    if (identical(value, "%00")) return("")
    value <- gsub("%09", "\t", value, fixed = TRUE)
    value <- gsub("%0D", "\r", value, fixed = TRUE)
    value <- gsub("%0A", "\n", value, fixed = TRUE)
    gsub("%25", "%", value, fixed = TRUE)
  }
  dims <- vapply(c("figure_width", "figure_height", "figure_dpi"), function(field) {
    value <- suppressWarnings(as.numeric(constant(field)))
    if (!is.finite(value) || value <= 0) stop("Invalid lollipop dimensions.", call. = FALSE)
    value
  }, numeric(1))
  background <- constant("figure_background")
  renderer <- if (variant == "direction") plot_lisa_gsea_direction_lollipop else plot_lisa_gsea_lollipop
  plot <- if (variant == "hommel_support_v1") {
    lisa_plot_hommel_support(x, flag == "true", order)
  } else if (variant %in% c("support_v1", "direction_support_v1")) {
    if (!"support_schema_version" %in% names(x) || !identical(unique(as.character(x$support_schema_version)), "support-grades-v1"))
      stop("Unsupported support presentation schema.", call. = FALSE)
    lisa_plot_support_lollipop(x, flag == "true", order, direction = variant == "direction_support_v1")
  } else renderer(x, flag == "true", order)
  plot <- lisa_apply_plot_context(plot, decode("figure_title_encoded"),
    if (variant == "hommel_support_v1") "Adjusted category P \u00b7 robust Hommel multiple-testing correction" else
      decode("figure_subtitle_encoded"), decode("positive_direction_encoded"))
  list(plot = plot, width = dims[[1L]], height = dims[[2L]], dpi = dims[[3L]],
    background = background)
}
