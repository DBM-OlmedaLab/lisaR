#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop("Usage: Rscript reproduce_lisa_figure.R FIGURE_SOURCE.tsv[.gz] [OUTPUT.png]", call. = FALSE)
}

source_tsv <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
if (!grepl("[.]tsv(?:[.]gz)?$", source_tsv, ignore.case = TRUE, perl = TRUE)) {
  stop("Figure source must be a .tsv or .tsv.gz file.", call. = FALSE)
}
output_arg <- if (length(args) == 2L) {
  path.expand(args[[2]])
} else {
  sub("([.]tsv)([.]gz)?$", "_regenerated.png", source_tsv)
}
if (!grepl("[.]png$", output_arg, ignore.case = TRUE)) {
  stop("Figure output must be a .png file.", call. = FALSE)
}
output_parent <- normalizePath(dirname(output_arg), winslash = "/", mustWork = TRUE)
if (!dir.exists(output_parent)) {
  stop("Figure output parent must be an existing directory.", call. = FALSE)
}
output_png <- file.path(output_parent, basename(output_arg))
output_link <- Sys.readlink(output_png)
if (!is.na(output_link) && nzchar(output_link)) {
  stop("Figure output must not be a symbolic link.", call. = FALSE)
}
if (file.exists(output_png)) {
  # Refuse all existing destinations. Besides preserving prior output, this
  # closes hard-link and path-alias cases that cannot be distinguished
  # portably after a destructive graphics-device open.
  stop("Figure output already exists; choose a new .png path.", call. = FALSE)
}
if (identical(output_png, source_tsv)) {
  stop("Figure output must differ from the source TSV.", call. = FALSE)
}

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
}

one <- function(x, default = "") {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x)) x[[1]] else default
}

constant_metadata <- function(x, field, default = NULL, required = FALSE) {
  if (!field %in% names(x)) {
    if (isTRUE(required)) {
      stop(sprintf("Figure source TSV lacks required metadata '%s'.", field), call. = FALSE)
    }
    return(default)
  }
  values <- as.character(x[[field]])
  if (length(values) != nrow(x) || anyNA(values) || any(!nzchar(trimws(values))) ||
      length(unique(values)) != 1L) {
    stop(sprintf("Figure source TSV has inconsistent metadata '%s'.", field), call. = FALSE)
  }
  values[[1L]]
}

constant_flag_metadata <- function(x, field, default = NULL, required = FALSE) {
  value <- constant_metadata(x, field, default = default, required = required)
  if (is.null(value)) return(default)
  normalized <- tolower(trimws(as.character(value)))
  if (!normalized %in% c("true", "t", "1", "yes", "false", "f", "0", "no")) {
    stop(sprintf("Figure source TSV has invalid logical metadata '%s'.", field), call. = FALSE)
  }
  normalized %in% c("true", "t", "1", "yes")
}

positive_metadata_number <- function(x, field, default = NULL, required = FALSE) {
  value <- constant_metadata(x, field, default = default, required = required)
  if (is.null(value)) return(default)
  value <- suppressWarnings(as.numeric(value))
  if (length(value) != 1L || !is.finite(value) || value <= 0) {
    stop(sprintf("Figure source TSV has invalid positive metadata '%s'.", field), call. = FALSE)
  }
  value
}

decode_figure_source_text <- function(x) {
  x <- as.character(x)
  x <- gsub("%09", "\t", x, fixed = TRUE)
  x <- gsub("%0D", "\r", x, fixed = TRUE)
  x <- gsub("%0A", "\n", x, fixed = TRUE)
  gsub("%25", "%", x, fixed = TRUE)
}

num <- function(x) suppressWarnings(as.numeric(x))
truth <- function(x) tolower(as.character(x)) %in% c("true", "t", "1", "yes")

save_plot <- function(plot, width = 9, height = 7, dpi = 220) {
  ggplot2::ggsave(output_png, plot, width = width, height = height, dpi = dpi,
    bg = "white", limitsize = FALSE)
}

render_category_gene_sets <- function(x) {
  keep <- truth(x$selected_for_plot)
  d <- x[keep, , drop = FALSE]
  if (!nrow(d)) stop("No selected category gene sets in source data.", call. = FALSE)
  d <- d[order(num(d$plotted_order)), , drop = FALSE]
  # Replay the original plot implementation, including its axis limits,
  # typography and dimensions. Saved order/labels are authoritative.
  d$NES <- num(d$plot_x_nes)
  d$source_family_order <- 1L
  d$pathway_order_key <- seq_len(nrow(d))
  d$pathway_display_label <- d$pathway_label
  meta <- d[1L, , drop = FALSE]
  meta$color <- one(d$category_color, "#737373")
  positive <- one(x$figure_positive_direction)
  comparison <- one(x$figure_comparison_subtitle)
  p <- lisaR:::plot_lisa_gsea_one_category_pathways(d, meta,
    gsea_padj_cutoff = num(one(x$selection_fdr_cutoff, "0.05")),
    label_chars = num(one(x$figure_label_chars, "62")),
    x_limits = c(num(one(x$plot_x_min)), num(one(x$plot_x_max))),
    comparison_subtitle = if (nzchar(comparison)) comparison else NULL,
    positive_direction = if (nzchar(positive)) positive else NULL)
  ggplot2::ggsave(output_png, p,
    width = positive_metadata_number(x, "figure_width", 12.5),
    height = positive_metadata_number(x, "figure_height", lisaR:::category_pathway_plot_height(nrow(d))),
    dpi = positive_metadata_number(x, "figure_dpi", 300),
    bg = one(x$figure_background, "white"), limitsize = FALSE)
}

render_volcano <- function(x) {
  d <- x[truth(x$selected_for_plot), , drop = FALSE]
  native_style <- "native_volcano_style" %in% names(d) && all(truth(d$native_volcano_style))
  d$highlighted <- truth(d$highlighted)
  d$labelled <- truth(d$labelled)
  d$de_class <- as.character(d$de_class)
  if (!"plot_x_log2FC" %in% names(d)) d$plot_x_log2FC <- num(d$log2FC)
  if (!"plot_y_neg_log10_fdr" %in% names(d)) d$plot_y_neg_log10_fdr <- -log10(pmax(num(d$padj), 1e-300))
  color <- one(d$category_color, "#009E73")
  lfc <- num(one(d$threshold_abs_log2FC, 0.58))
  fdr <- num(one(d$threshold_de_fdr, 0.05))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = plot_x_log2FC, y = plot_y_neg_log10_fdr)) +
    ggplot2::geom_point(ggplot2::aes(color = de_class), alpha = 0.15, size = 0.55, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", color = "grey70", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(fdr), linetype = "dashed", color = "grey70", linewidth = 0.3) +
    ggplot2::geom_point(data = d[d$highlighted, , drop = FALSE], color = "grey10", fill = color,
      shape = 21, size = 1.55, alpha = 0.9, stroke = 0.22, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(`DE down` = "#4A72B8", `DE up` = "#C44E52", `not DE` = "#B8B8B8")) +
    ggplot2::labs(title = one(d$category_display_name, one(d$category_id)), subtitle = one(d$plot_subtitle),
      x = "log2FC", y = "-log10(FDR or p-value)", color = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
  if (native_style) p <- p + ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    plot.subtitle = ggplot2::element_text(size = 9, color = "grey35"),
    plot.background = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA))
  lab <- d[d$labelled, , drop = FALSE]
  if (nrow(lab)) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      repel <- list(data = lab, mapping = ggplot2::aes(label = symbol), size = 3.25,
        max.overlaps = Inf, box.padding = 0.42, point.padding = 0.16, seed = 1)
      if (native_style) repel <- c(repel, list(force = 2, force_pull = 0.15,
        min.segment.length = 0, segment.size = 0.22, segment.color = "grey55", na.rm = TRUE))
      p <- p + do.call(ggrepel::geom_text_repel, repel)
    } else p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(label = symbol), size = 2.5, check_overlap = TRUE)
  }
  save_plot(p, 8.5, 7.2, 240)
}

render_gene_card <- function(x) {
  d <- x[truth(x$selected_for_plot), , drop = FALSE]
  d <- d[order(num(d$plotted_order)), , drop = FALSE]
  if (!nrow(d)) stop("No selected genes in GeneCard source data.", call. = FALSE)
  d$gene <- factor(d$symbol, levels = rev(d$symbol))
  d$score <- num(d$gene_contribution_score)
  color <- one(d$category_color, one(d$color, "#525252"))
  p1 <- ggplot2::ggplot(d, ggplot2::aes(x = gene, y = 1, fill = num(log2FC))) +
    ggplot2::geom_tile(color = "white") + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", num(log2FC))), size = 2.5) +
    ggplot2::coord_flip() + ggplot2::scale_fill_gradient2(low = "#3B66B0", mid = "white", high = "#B73A3A", midpoint = 0, name = "log2FC") +
    ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_minimal(base_size = 9) + ggplot2::theme(axis.text.x = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
  p2 <- ggplot2::ggplot(d, ggplot2::aes(x = gene, y = score)) + ggplot2::geom_col(fill = color) +
    ggplot2::coord_flip() + ggplot2::labs(x = NULL, y = "DE-weighted contribution") + ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank())
  marks <- rbind(
    data.frame(gene = d$gene, mark = "FDR<0.05", value = truth(d$de_is_significant)),
    data.frame(gene = d$gene, mark = "LE", value = truth(d$in_gsea_leading_edge)),
    data.frame(gene = d$gene, mark = "ORA", value = truth(d$in_ora_overlap)))
  p3 <- ggplot2::ggplot(marks, ggplot2::aes(x = mark, y = gene)) +
    ggplot2::geom_point(ggplot2::aes(fill = value), shape = 21, size = 2.8) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#0AA36E", `FALSE` = "white")) +
    ggplot2::labs(x = NULL, y = NULL, title = "Evidence") + ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid = ggplot2::element_blank(), legend.position = "none")
  p4 <- ggplot2::ggplot(d, ggplot2::aes(x = num(category_recurrence), y = gene)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = num(category_recurrence), yend = gene), color = "#D2D7DD") +
    ggplot2::geom_point(color = "#4C6A92", fill = "#DDE8F7", shape = 21) +
    ggplot2::labs(x = "LISA categories", y = NULL) + ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), panel.grid.major.y = ggplot2::element_blank())
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Package 'patchwork' is required for GeneCard reproduction.", call. = FALSE)
  title <- sprintf("%s\n%s | %d supporting genes shown", one(d$category_display_name, one(d$category_id)), one(d$macrogroup_name), nrow(d))
  p <- (p1 | p2 | p3 | p4) + patchwork::plot_annotation(title = title)
  save_plot(p, 8.9, max(8.2, 5.3 + 0.34 * nrow(d)))
}

render_saved_native_gene <- function(x, type) {
  script <- if(type == "gene_card") "build_single_de_category_gene_cards.R" else "build_contrast_gene_level_product.R"
  builder <- system.file("scripts",script,package="lisaR")
  if(!nzchar(builder)) stop("Native gene figure builder is not installed.")
  env <- new.env(parent=globalenv());sys.source(builder,env)
  if(type == "gene_card") {
    n <- as.integer(one(x$selection_top_genes,15L))
    p <- env$build_card_plot(x,one(x$category_id),n)
    save_plot(p,8.9,max(8.2,5.3+0.34*n),220)
  } else {
    pdf <- tempfile(fileext=".pdf");on.exit(unlink(pdf),add=TRUE)
    cfg <- list(lfc_cap=num(one(x$lfc_cap,1.5)),
      top_genes_per_category=if(type=="contrast_gene_card") sum(truth(x$selected_for_plot)) else 2L*max(table(x$category_id)))
    if(type=="contrast_gene_card") env$plot_category_card(x,output_png,pdf,cfg)
    else env$plot_combined_heatmap(x,unique(as.character(x$category_id)),output_png,pdf,cfg)
  }
}

render_contrast_gene_card <- function(x) {
  d <- x[truth(x$selected_for_plot), , drop = FALSE]
  d <- d[order(num(d$plotted_order)), , drop = FALSE]
  d$symbol <- factor(d$symbol, levels = rev(unique(d$symbol)))
  long <- rbind(
    data.frame(symbol = d$symbol, analysis = "A", log2FC = num(d$log2FC_A)),
    data.frame(symbol = d$symbol, analysis = "B", log2FC = num(d$log2FC_B)))
  cap <- num(one(d$lfc_cap, 1.5))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = analysis, y = symbol,
      fill = pmax(pmin(log2FC, cap), -cap))) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(log2FC), "", sprintf("%.2f", log2FC))), size = 2.5) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-cap, cap), na.value = "grey92", name = "log2FC") +
    ggplot2::labs(title = paste(one(d$category_display_name), one(d$direction_class), sep = " | "),
      subtitle = sprintf("A = %s; B = %s", one(d$label_A), one(d$label_B)), x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.grid = ggplot2::element_blank())
  save_plot(p, 7.1, max(4.5, 1 + 0.22 * nrow(d)), 180)
}

render_contrast_heatmap <- function(x) {
  x$row_label <- paste(x$category_display_name, x$symbol, sep = " | ")
  x$row_label <- factor(x$row_label, levels = rev(unique(x$row_label)))
  long <- rbind(
    data.frame(row_label = x$row_label, analysis = "A", log2FC = num(x$log2FC_A)),
    data.frame(row_label = x$row_label, analysis = "B", log2FC = num(x$log2FC_B)))
  cap <- num(one(x$lfc_cap, 1.5))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = analysis, y = row_label,
      fill = pmax(pmin(log2FC, cap), -cap))) +
    ggplot2::geom_tile(color = "white", linewidth = 0.18) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-cap, cap), na.value = "grey92", name = "log2FC") +
    ggplot2::labs(title = "Paired heatmap", subtitle = sprintf("A = %s; B = %s", one(x$label_A), one(x$label_B)), x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 8.5) + ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none")
  save_plot(p, 3.9, max(5.8, 0.125 * nrow(x)), 180)
}

value_rgb <- function(value, max_abs, power = 1) {
  if (!is.finite(value)) return(c(0.82, 0.82, 0.82))
  z <- sign(value) * (min(abs(value), max_abs) / max_abs)^power
  white <- c(1, 1, 1)
  target <- if (z >= 0) c(0.76, 0.14, 0.16) else c(0.20, 0.38, 0.72)
  white * (1 - abs(z)) + target * abs(z)
}

paint_region <- function(img, x, y, width, height, rgb, side = "full") {
  h <- dim(img)[1]; w <- dim(img)[2]
  xmin <- max(1L, floor(x - width / 2)); xmax <- min(w, ceiling(x + width / 2))
  ymin <- max(1L, floor(y - height / 2)); ymax <- min(h, ceiling(y + height / 2))
  if (side == "left") xmax <- floor((xmin + xmax) / 2)
  if (side == "right") xmin <- ceiling((xmin + xmax) / 2)
  if (xmin > xmax || ymin > ymax) return(img)
  region <- img[ymin:ymax, xmin:xmax, , drop = FALSE]
  if (dim(region)[3] < 3L) return(img)
  alpha <- 0.72
  for (ch in 1:3) region[, , ch] <- region[, , ch] * (1 - alpha) + rgb[[ch]] * alpha
  img[ymin:ymax, xmin:xmax, 1:3] <- region[, , 1:3]
  img
}

render_kegg <- function(x, contrast = FALSE) {
  if (!requireNamespace("png", quietly = TRUE)) stop("Package 'png' is required for KEGG reproduction.", call. = FALSE)
  base_file <- file.path(dirname(source_tsv), one(x$base_image_file))
  if (!file.exists(base_file)) stop("Companion KEGG base image is missing: ", base_file, call. = FALSE)
  builder <- system.file("scripts", if (contrast) "build_contrast_kegg_pathway_painter.R"
    else "build_single_de_kegg_pathway_painter.R", package="lisaR")
  if (!nzchar(builder)) stop("Install the recorded lisaR native KEGG painter.")
  env <- new.env(parent=globalenv());sys.source(builder,env)
  pdf <- tempfile(fileext=".pdf");on.exit(unlink(pdf),add=TRUE)
  env$draw_painted_map(png::readPNG(base_file), x, one(x$figure_title),
    one(x$figure_subtitle), output_png, pdf, num(one(x$max_abs_log2fc,0.5)),
    num(one(x$color_power,1)))
}

render_heatmap <- function(x) {
  d <- x[truth(x$selected_for_plot), , drop = FALSE]
  row_cols <- intersect(c("symbol", "gene", "row", "category_display_name", "category_id"), names(d))
  if (!length(row_cols)) stop("Heatmap source data lacks a row identifier.", call. = FALSE)
  row_col <- row_cols[[1]]

  value_cols <- intersect(c("plot_value", "zscore", "scaled_value", "value", "expression"), names(d))
  col_cols <- intersect(c("sample", "column", "condition", "side"), names(d))
  if (length(value_cols) && length(col_cols)) {
    value_col <- value_cols[[1]]
    col_col <- col_cols[[1]]
    long <- d
  } else if (all(c("log2FC_A", "log2FC_B") %in% names(d))) {
    label_a <- one(d$label_A, "A")
    label_b <- one(d$label_B, "B")
    cap <- num(one(d$lfc_cap, NA_real_))
    a <- num(d$log2FC_A)
    b <- num(d$log2FC_B)
    if (is.finite(cap) && cap > 0) {
      a <- pmax(-cap, pmin(cap, a))
      b <- pmax(-cap, pmin(cap, b))
    }
    long <- rbind(
      data.frame(row_id = d[[row_col]], column_id = label_a, plot_value = a),
      data.frame(row_id = d[[row_col]], column_id = label_b, plot_value = b)
    )
    row_col <- "row_id"
    col_col <- "column_id"
    value_col <- "plot_value"
  } else {
    stop("Heatmap source data lacks supported long or paired A/B geometry.", call. = FALSE)
  }

  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data[[col_col]], y = .data[[row_col]], fill = num(.data[[value_col]]))) +
    ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#3B66B0", mid = "white", high = "#B73A3A", midpoint = 0) +
    ggplot2::labs(x = NULL, y = NULL, fill = value_col) + ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  save_plot(p, max(7, 0.35 * length(unique(long[[col_col]])) + 4), max(5, 0.22 * length(unique(long[[row_col]])) + 2))
}

lighten_color <- function(color, amount = 0.25) {
  rgb <- grDevices::col2rgb(color) / 255
  rgb <- rgb + (1 - rgb) * amount
  grDevices::rgb(rgb[1], rgb[2], rgb[3])
}

render_dumbbell <- function(x) {
  metadata <- x
  plot_set <- constant_metadata(metadata, "plot_set", required = TRUE)
  annotation_variant <- tolower(constant_metadata(
    metadata, "annotation_variant", required = TRUE
  ))
  if (!annotation_variant %in% c("plain", "annotated")) {
    stop("Dumbbell source TSV has unsupported annotation_variant.", call. = FALSE)
  }
  annotate_gene_sets <- identical(annotation_variant, "annotated")
  group_by_supracategory <- constant_flag_metadata(
    metadata, "group_by_supracategory", required = TRUE
  )
  figure_width <- positive_metadata_number(
    metadata, "figure_width", required = TRUE
  )
  figure_height <- positive_metadata_number(
    metadata, "figure_height", required = TRUE
  )
  figure_dpi <- positive_metadata_number(
    metadata, "figure_dpi", required = TRUE
  )
  figure_title <- decode_figure_source_text(constant_metadata(
    metadata, "figure_title", required = TRUE
  ))
  figure_subtitle <- decode_figure_source_text(constant_metadata(
    metadata, "figure_subtitle", required = TRUE
  ))
  contrast_a_label <- decode_figure_source_text(constant_metadata(
    metadata, "contrast_a_label", required = TRUE
  ))
  contrast_b_label <- decode_figure_source_text(constant_metadata(
    metadata, "contrast_b_label", required = TRUE
  ))
  if ("selected_for_plot" %in% names(x)) {
    x <- x[truth(x$selected_for_plot), , drop = FALSE]
  }
  if (!nrow(x)) {
    p <- ggplot2::ggplot(
      data.frame(x = 0, y = 1, label = "No categories match this filter"),
      ggplot2::aes(x = x, y = y, label = label)
    ) +
      ggplot2::geom_text(size = 4, color = "grey35") +
      ggplot2::labs(
        title = figure_title,
        subtitle = figure_subtitle,
        x = "Mean NES", y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        panel.background = ggplot2::element_rect(fill = "white", color = NA)
      )
    cat(paste(c(
      "LISA_DUMBBELL_REPRODUCTION",
      paste0("plot_set=", gsub("[[:cntrl:]]+", " ", plot_set)),
      paste0("variant=", annotation_variant),
      paste0("grouped=", tolower(as.character(group_by_supracategory))),
      "facets=0", "panel_rows=", "categories=0", "blank=0",
      "segments=0", "supported=0", "contextual=0", "labels=0",
      paste0("title=", gsub("[[:cntrl:]]+", " ", figure_title)),
      paste0("width=", format(figure_width, trim = TRUE)),
      paste0("height=", format(figure_height, trim = TRUE)),
      paste0("dpi=", format(figure_dpi, trim = TRUE))
    ), collapse = "\t"), "\n", sep = "")
    save_plot(p, figure_width, figure_height, figure_dpi)
    return(invisible(NULL))
  }

  d <- x
  if (!"n_genesets_A" %in% names(d)) {
    d$n_genesets_A <- as.integer(is.finite(num(d$mean_NES_A)))
  }
  if (!"n_genesets_B" %in% names(d)) {
    d$n_genesets_B <- as.integer(is.finite(num(d$mean_NES_B)))
  }
  if (!"same_direction_pct_A" %in% names(d)) d$same_direction_pct_A <- NA_real_
  if (!"same_direction_pct_B" %in% names(d)) d$same_direction_pct_B <- NA_real_
  if (!"has_significant_support_A" %in% names(d)) {
    d$has_significant_support_A <- num(d$n_genesets_A) > 0
  } else {
    d$has_significant_support_A <- truth(d$has_significant_support_A)
  }
  if (!"has_significant_support_B" %in% names(d)) {
    d$has_significant_support_B <- num(d$n_genesets_B) > 0
  } else {
    d$has_significant_support_B <- truth(d$has_significant_support_B)
  }
  if (!"plot_has_any_significant_support" %in% names(d)) {
    d$plot_has_any_significant_support <-
      d$has_significant_support_A | d$has_significant_support_B
  } else {
    d$plot_has_any_significant_support <- truth(d$plot_has_any_significant_support)
  }
  if (!"display_mean_NES_A" %in% names(d)) {
    contextual <- if ("mean_NES_contextual_A" %in% names(d)) d$mean_NES_contextual_A else d$mean_NES_A
    d$display_mean_NES_A <- ifelse(
      d$plot_has_any_significant_support,
      ifelse(d$has_significant_support_A, num(d$mean_NES_A), num(contextual)),
      NA_real_
    )
  } else {
    d$display_mean_NES_A <- num(d$display_mean_NES_A)
  }
  if (!"display_mean_NES_B" %in% names(d)) {
    contextual <- if ("mean_NES_contextual_B" %in% names(d)) d$mean_NES_contextual_B else d$mean_NES_B
    d$display_mean_NES_B <- ifelse(
      d$plot_has_any_significant_support,
      ifelse(d$has_significant_support_B, num(d$mean_NES_B), num(contextual)),
      NA_real_
    )
  } else {
    d$display_mean_NES_B <- num(d$display_mean_NES_B)
  }
  if (!"n_genesets_evaluable_A" %in% names(d)) d$n_genesets_evaluable_A <- d$n_genesets_A
  if (!"n_genesets_evaluable_B" %in% names(d)) d$n_genesets_evaluable_B <- d$n_genesets_B
  if (!"gsea_padj_cutoff" %in% names(d)) d$gsea_padj_cutoff <- NA_real_
  if (!"category_id" %in% names(d)) stop("Dumbbell source TSV lacks category_id.", call. = FALSE)
  if (!"category_display_name" %in% names(d)) d$category_display_name <- d$category_id
  if (!"macrogroup_name" %in% names(d)) d$macrogroup_name <- "Other or unclassified"
  if (!"color" %in% names(d)) d$color <- "#737373"
  missing_label <- is.na(d$category_display_name) | d$category_display_name == ""
  d$category_display_name[missing_label] <- d$category_id[missing_label]
  d$macrogroup_name[is.na(d$macrogroup_name) | d$macrogroup_name == ""] <- "Other or unclassified"
  d$color[is.na(d$color) | d$color == ""] <- "#737373"

  d$category_plot_id <- make.unique(paste(d$category_id, d$category_display_name, sep = "__"))
  d$category_plot_id <- factor(d$category_plot_id, levels = rev(unique(d$category_plot_id)))
  d$macrogroup_name <- factor(d$macrogroup_name, levels = unique(d$macrogroup_name))
  label_values <- stats::setNames(d$category_display_name, as.character(d$category_plot_id))
  finite_nes <- c(d$display_mean_NES_A, d$display_mean_NES_B)
  finite_nes <- finite_nes[is.finite(finite_nes)]
  max_abs_nes <- if (length(finite_nes)) max(abs(finite_nes), na.rm = TRUE) else 1
  if (!is.finite(max_abs_nes) || max_abs_nes == 0) max_abs_nes <- 1
  x_limit_multiplier <- if (annotate_gene_sets) 1.52 else 1.16
  x_limits <- c(-max_abs_nes, max_abs_nes) * x_limit_multiplier
  label_pad <- max_abs_nes * 0.055
  a_is_left_endpoint <- d$display_mean_NES_A <= d$display_mean_NES_B
  a_is_left_endpoint[is.na(a_is_left_endpoint)] <- TRUE

  point_df <- rbind(
    data.frame(
      category_plot_id = d$category_plot_id,
      macrogroup_name = d$macrogroup_name,
      mean_NES = d$display_mean_NES_A,
      n_genesets = d$n_genesets_A,
      n_genesets_evaluable = d$n_genesets_evaluable_A,
      same_direction_pct = d$same_direction_pct_A,
      has_significant_support = d$has_significant_support_A,
      category_has_signal = d$plot_has_any_significant_support,
      contrast = contrast_a_label,
      label_side = ifelse(a_is_left_endpoint, -1, 1),
      point_fill = vapply(d$color, lighten_color, character(1), amount = 0.55),
      stringsAsFactors = FALSE
    ),
    data.frame(
      category_plot_id = d$category_plot_id,
      macrogroup_name = d$macrogroup_name,
      mean_NES = d$display_mean_NES_B,
      n_genesets = d$n_genesets_B,
      n_genesets_evaluable = d$n_genesets_evaluable_B,
      same_direction_pct = d$same_direction_pct_B,
      has_significant_support = d$has_significant_support_B,
      category_has_signal = d$plot_has_any_significant_support,
      contrast = contrast_b_label,
      label_side = ifelse(a_is_left_endpoint, 1, -1),
      point_fill = d$color,
      stringsAsFactors = FALSE
    )
  )
  point_df$contrast <- factor(point_df$contrast, levels = c(contrast_a_label, contrast_b_label))
  point_df$n_genesets <- num(point_df$n_genesets)
  point_df$n_genesets_evaluable <- num(point_df$n_genesets_evaluable)
  point_df$same_direction_pct <- num(point_df$same_direction_pct)
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

  segment_df <- d[
    d$plot_has_any_significant_support &
      is.finite(d$display_mean_NES_A) & is.finite(d$display_mean_NES_B),
    , drop = FALSE
  ]
  visible_points <- point_df[point_df$endpoint_visible, , drop = FALSE]
  supported_points <- visible_points[visible_points$has_significant_support, , drop = FALSE]
  contextual_points <- visible_points[!visible_points$has_significant_support, , drop = FALSE]

  p <- ggplot2::ggplot(d, ggplot2::aes(y = category_plot_id)) +
    ggplot2::geom_blank(ggplot2::aes(x = 0, y = category_plot_id)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.35) +
    ggplot2::geom_segment(
      data = segment_df,
      ggplot2::aes(x = display_mean_NES_A, xend = display_mean_NES_B,
        y = category_plot_id, yend = category_plot_id, color = color),
      inherit.aes = FALSE, linewidth = 0.55, alpha = 0.70
    ) +
    ggplot2::geom_point(
      data = supported_points,
      ggplot2::aes(x = mean_NES, y = category_plot_id, shape = contrast,
        size = n_genesets, fill = point_fill),
      inherit.aes = FALSE, color = "grey18", stroke = 0.28, alpha = 0.96
    ) +
    ggplot2::geom_point(
      data = contextual_points,
      ggplot2::aes(x = mean_NES, y = category_plot_id, shape = contrast, size = n_genesets),
      inherit.aes = FALSE, fill = NA, color = "grey35", stroke = 0.75, alpha = 0.38
    )
  if (annotate_gene_sets) {
    p <- p +
      ggplot2::geom_text(
        data = visible_points[visible_points$contrast == contrast_a_label, , drop = FALSE],
        ggplot2::aes(x = label_x, y = category_plot_id, label = point_label, hjust = label_hjust),
        inherit.aes = FALSE, nudge_y = 0.17, size = 3.25, color = "grey25",
        lineheight = 0.9, show.legend = FALSE
      ) +
      ggplot2::geom_text(
        data = visible_points[visible_points$contrast == contrast_b_label, , drop = FALSE],
        ggplot2::aes(x = label_x, y = category_plot_id, label = point_label, hjust = label_hjust),
        inherit.aes = FALSE, nudge_y = -0.17, size = 3.25, color = "grey25",
        lineheight = 0.9, show.legend = FALSE
      )
  }
  contrast_levels <- c(contrast_a_label, contrast_b_label)
  cutoff <- one(d$gsea_padj_cutoff, "0.25")
  p <- p +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_shape_manual(
      values = stats::setNames(c(21, 22), contrast_levels),
      limits = contrast_levels, drop = FALSE
    ) +
    ggplot2::scale_size_continuous(
      range = c(1.7, 6.2), breaks = c(0, 1, 5, 10, 20), name = "n gene sets"
    ) +
    ggplot2::scale_y_discrete(labels = label_values, drop = TRUE) +
    ggplot2::scale_x_continuous(
      limits = x_limits, expand = ggplot2::expansion(mult = c(0.02, 0.02))
    ) +
    ggplot2::labs(
      title = figure_title,
      subtitle = figure_subtitle,
      x = "Mean NES (significant-set mean; contextual pre-threshold mean for an unsupported side)",
      y = NULL, shape = "Contrast",
      caption = sprintf(
        paste0(
          "Filled: significant support at GSEA FDR <= %s.\n",
          "Hollow/translucent: contextual mean without significant support. Blank row: neither side significant."
        ),
        cutoff
      )
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(override.aes = list(
        size = 5.2, fill = "grey75", color = "grey18", stroke = 0.45
      )),
      size = ggplot2::guide_legend(override.aes = list(
        shape = 21, fill = "grey75", color = "grey18"
      ))
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
  if (group_by_supracategory) {
    p <- p +
      ggplot2::facet_grid(
        macrogroup_name ~ ., scales = "free_y", space = "free_y", switch = "y"
      ) +
      ggplot2::theme(
        strip.placement = "outside",
        strip.background.y = ggplot2::element_rect(fill = "grey94", color = "grey80"),
        strip.text.y.left = ggplot2::element_text(
          angle = 0, hjust = 1, size = 8, face = "bold"
        ),
        panel.spacing.y = grid::unit(0.25, "lines")
      )
  }

  panel_rows <- if (group_by_supracategory) {
    as.integer(table(d$macrogroup_name))
  } else {
    nrow(d)
  }
  clean_field <- function(value) gsub("[[:cntrl:]]+", " ", as.character(value))
  cat(paste(c(
    "LISA_DUMBBELL_REPRODUCTION",
    paste0("plot_set=", clean_field(plot_set)),
    paste0("variant=", annotation_variant),
    paste0("grouped=", tolower(as.character(group_by_supracategory))),
    paste0("facets=", length(panel_rows)),
    paste0("panel_rows=", paste(panel_rows, collapse = ",")),
    paste0("categories=", nrow(d)),
    paste0("blank=", sum(!d$plot_has_any_significant_support)),
    paste0("segments=", nrow(segment_df)),
    paste0("supported=", nrow(supported_points)),
    paste0("contextual=", nrow(contextual_points)),
    paste0("labels=", if (annotate_gene_sets) nrow(visible_points) else 0L),
    paste0("contrast_a=", clean_field(contrast_a_label)),
    paste0("contrast_b=", clean_field(contrast_b_label)),
    paste0("colors=", paste(unique(as.character(d$color)), collapse = ",")),
    paste0("title=", clean_field(figure_title)),
    paste0("width=", format(figure_width, trim = TRUE)),
    paste0("height=", format(figure_height, trim = TRUE)),
    paste0("dpi=", format(figure_dpi, trim = TRUE))
  ), collapse = "\t"), "\n", sep = "")
  save_plot(p, figure_width, figure_height, figure_dpi)
}

render_lisa_contrast_heatmap <- function(x) {
  title <- one(x$figure_title, one(x$category_display_name))
  width <- positive_metadata_number(x, "figure_width", 5.5)
  height <- positive_metadata_number(x, "figure_height", 3.8)
  dpi <- positive_metadata_number(x, "figure_dpi", 300)
  support_a <- if ("has_significant_support_A" %in% names(x)) truth(x$has_significant_support_A) else if ("n_genesets_A" %in% names(x)) num(x$n_genesets_A) > 0 else is.finite(num(x$mean_NES_A))
  support_b <- if ("has_significant_support_B" %in% names(x)) truth(x$has_significant_support_B) else if ("n_genesets_B" %in% names(x)) num(x$n_genesets_B) > 0 else is.finite(num(x$mean_NES_B))
  has_signal <- if ("plot_has_any_significant_support" %in% names(x)) truth(x$plot_has_any_significant_support) else support_a | support_b
  if (!any(has_signal)) {
    p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
      ggplot2::geom_blank() +
      ggplot2::annotate("text", x = 1, y = 1,
        label = sprintf("No member gene set met GSEA FDR <= %s on either side", one(x$gsea_padj_cutoff, "0.25")),
        color = "grey35", size = 3.8) +
      ggplot2::labs(title = one(x$figure_title, one(x$category_display_name))) + ggplot2::theme_void()
    return(save_plot(p, width, height, dpi))
  }
  label_a <- one(x$contrast_a_label, "A")
  label_b <- one(x$contrast_b_label, "B")
  display_a <- if ("display_mean_NES_A" %in% names(x)) num(x$display_mean_NES_A) else num(x$mean_NES_A)
  display_b <- if ("display_mean_NES_B" %in% names(x)) num(x$display_mean_NES_B) else num(x$mean_NES_B)
  values <- data.frame(
    side = factor(c(label_a, label_b), levels = c(label_a, label_b)),
    mean_NES = c(display_a[[1]], display_b[[1]]),
    support = c(support_a[[1]], support_b[[1]])
  )
  values$label <- sprintf("%.2f%s", values$mean_NES, ifelse(values$support, "", " (NS)"))
  p <- ggplot2::ggplot(values, ggplot2::aes(x = side, y = title, fill = mean_NES, alpha = support)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 4) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    ggplot2::labs(title = one(x$figure_title, one(x$category_display_name)), x = NULL, y = NULL, fill = "Mean NES") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
  save_plot(p, width, height, dpi)
}

render_network <- function(x) {
  category <- one(x$category_id)
  d <- x[x$category_id == category, , drop = FALSE]
  d <- d[order(-num(d$max_gene_contribution_score), d$symbol), , drop = FALSE]
  d <- head(d[!duplicated(d$symbol), , drop = FALSE], 40)
  d$symbol <- factor(d$symbol, levels = rev(d$symbol))
  long <- rbind(
    data.frame(symbol = d$symbol, side = "A", log2FC = num(d$log2FC_A)),
    data.frame(symbol = d$symbol, side = "B", log2FC = num(d$log2FC_B)))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = side, y = symbol, size = abs(log2FC), color = log2FC)) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0) +
    ggplot2::labs(title = paste("Gene-category network:", one(d$category_display_name, category)), x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 9) + ggplot2::theme(panel.grid = ggplot2::element_blank())
  save_plot(p, 7, max(5, 0.22 * nrow(d) + 2))
}

render_generic <- function(x) {
  numeric_names <- names(x)[vapply(x, function(z) any(is.finite(num(z))), logical(1))]
  x_col <- intersect(c("NES", "mean_NES", "log2FC", "delta", "total_overlap", "category_count",
    "max_gene_contribution_score", "abs_delta_mean_NES", "plot_x"), numeric_names)
  y_col <- intersect(c("category_display_name", "category_id", "pathway", "symbol", "kegg_title"), names(x))
  if (!length(x_col) || !length(y_col)) stop("Source TSV has no supported plot geometry.", call. = FALSE)
  x_col <- x_col[[1]]; y_col <- y_col[[1]]
  p <- ggplot2::ggplot(x, ggplot2::aes(x = num(.data[[x_col]]), y = stats::reorder(.data[[y_col]], num(.data[[x_col]])))) +
    ggplot2::geom_col(fill = "#4C78A8") + ggplot2::labs(x = x_col, y = NULL, title = one(x$figure_title, one(x$figure_id))) +
    ggplot2::theme_minimal(base_size = 9)
  save_plot(p, 9, max(5, 0.24 * nrow(x) + 2))
}

render_lisa_lollipop <- function(x) {
  if (!requireNamespace("lisaR", quietly = TRUE) ||
      !exists("lisa_support_save_lollipop", envir = asNamespace("lisaR"), inherits = FALSE)) {
    stop("Lollipop reproduction requires the corrected lisaR package used to create its source data.", call. = FALSE)
  }
  # Reuse the original saver, including its explicit Cairo device. An automatic
  # ggsave PNG device can select ragg and change pixels despite identical data.
  lisaR:::lisa_support_save_lollipop(x, output_png)
}

render_native_extended <- function(x, type) {
  script <- switch(type, native_recurrent_genes = "build_single_de_recurrent_gene_screen.R",
    native_contrast_network = "build_contrast_gene_category_network.R")
  builder <- system.file("scripts", script, package = "lisaR")
  expected <- constant_metadata(x, "native_builder_sha256", required = TRUE)
  if (!nzchar(builder) || digest::digest(file = builder, algo = "sha256") != expected)
    stop("Install the recorded lisaR builder version to replay this native figure.")
  cfg <- jsonlite::fromJSON(constant_metadata(x, "native_config_json", required = TRUE), simplifyVector = FALSE)
  env <- new.env(parent = globalenv())
  sys.source(builder, env)
  if (type == "native_contrast_network") {
    pdf <- tempfile(fileext = ".pdf")
    on.exit(unlink(pdf), add = TRUE)
    env$plot_network(x, output_png, pdf, cfg)
  } else {
    empty <- "empty_reason" %in% names(x)
    plot <- if (empty) env$build_empty_plot(cfg, one(x$empty_reason)) else env$build_plot(x, cfg)
    ggplot2::ggsave(output_png, plot, width = 8.8, height = if (empty) 4.2 else 7.2,
      dpi = 240, bg = "white", limitsize = FALSE)
  }
}

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.", call. = FALSE)
x <- lisaR:::lisa_read_figure_source_tsv(source_tsv)
if (!nrow(x)) stop("Figure source TSV is empty.", call. = FALSE)
figure_type <- constant_metadata(x, "figure_type", required = TRUE)
if (figure_type %in% c("native_recurrent_genes", "native_contrast_network")) {
  render_native_extended(x, figure_type)
} else if (figure_type == "lisa_category_inference") {
  lisaR:::lisa_category_inference_save_plot(x, output_png)
} else if (figure_type == "lisa_category_gene_sets") {
  render_category_gene_sets(x)
} else if (figure_type == "volcano_overlay") {
  render_volcano(x)
} else if (figure_type == "gene_card") {
  render_saved_native_gene(x,"gene_card")
} else if (figure_type == "contrast_gene_card") {
  render_saved_native_gene(x,"contrast_gene_card")
} else if (figure_type == "contrast_heatmap") {
  render_saved_native_gene(x,"contrast_heatmap")
} else if (figure_type == "kegg_single") {
  render_kegg(x, FALSE)
} else if (figure_type == "kegg_contrast") {
  render_kegg(x, TRUE)
} else if (figure_type == "lisa_dumbbell") {
  render_dumbbell(x)
} else if (figure_type == "lisa_lollipop") {
  render_lisa_lollipop(x)
} else if (figure_type == "lisa_contrast_heatmap") {
  render_lisa_contrast_heatmap(x)
} else if (figure_type == "network") {
  render_network(x)
} else if (grepl("heatmap", figure_type, fixed = TRUE)) {
  render_heatmap(x)
} else {
  render_generic(x)
}

if (!file.exists(output_png) || file.info(output_png)$size <= 0) stop("Figure regeneration produced no PNG.", call. = FALSE)
message(normalizePath(output_png, mustWork = TRUE))
