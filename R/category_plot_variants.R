# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Descriptive NES distributions and optional category-plot annotations.
# The unit is one accepted gene set within one category. Statistics are shared
# by all visual variants; no enrichment, multiplicity adjustment, or category
# inference is performed here. In contrasts the distributions remain A/B.

lisa_nes_number <- function(x, field) {
  # Keep binary numeric inputs numeric. as.character() rounds doubles to about
  # 15 significant digits, changing canonical means/deltas and even inclusive
  # cutoff membership before the exact-source writer has a chance to save them.
  if (is.numeric(x) && !is.complex(x)) {
    out <- as.numeric(x)
    out[is.na(out)] <- NA_real_
    return(out)
  }
  raw <- as.character(x)
  missing <- is.na(raw) | raw %in% c("", "NA", "NaN")
  out <- suppressWarnings(as.numeric(raw))
  if (any(!missing & is.na(out)))
    stop("Invalid numeric NES-variant field: ", field, call. = FALSE)
  out[missing] <- NA_real_
  out
}

lisa_category_nes_statistics <- function(gsea_all, gsea_padj_cutoff = .25,
    category_ids = NULL) {
  if (!is.data.frame(gsea_all)) stop("GSEA members must be a data frame.", call. = FALSE)
  if (!is.numeric(gsea_padj_cutoff) || length(gsea_padj_cutoff) != 1L ||
      !is.finite(gsea_padj_cutoff) || gsea_padj_cutoff < 0 || gsea_padj_cutoff > 1)
    stop("gsea_padj_cutoff must be a number in [0,1].", call. = FALSE)
  set_col <- intersect(c("pathway", "gene_set_id"), names(gsea_all))
  padj_col <- intersect(c("padj", "gsea_fdr"), names(gsea_all))
  if (!all(c("category_id", "NES") %in% names(gsea_all)) ||
      !length(set_col) || !length(padj_col))
    stop("NES statistics require category_id, pathway/gene_set_id, NES and padj/gsea_fdr.", call. = FALSE)
  x <- data.frame(category_id = as.character(gsea_all$category_id),
    pathway = as.character(gsea_all[[set_col[[1L]]]]),
    NES = lisa_nes_number(gsea_all$NES, "NES"),
    padj = lisa_nes_number(gsea_all[[padj_col[[1L]]]], "padj"),
    stringsAsFactors = FALSE)
  keep <- !is.na(x$category_id) & nzchar(x$category_id) &
    x$category_id != "OTHER_UNCLASSIFIED"
  if ("classification_status" %in% names(gsea_all))
    keep <- keep & !is.na(gsea_all$classification_status) &
      gsea_all$classification_status == "classified"
  x <- x[keep, , drop = FALSE]
  if (any(is.na(x$pathway) | !nzchar(x$pathway)))
    stop("Classified NES members require a gene-set identifier.", call. = FALSE)
  if (any(is.finite(x$padj) & (x$padj < 0 | x$padj > 1)))
    stop("Finite GSEA adjusted p-values must lie in [0,1].", call. = FALSE)
  ids <- if (is.null(category_ids)) unique(x$category_id) else as.character(category_ids)
  if (anyNA(ids) || any(!nzchar(ids)) || any(ids == "OTHER_UNCLASSIFIED") || anyDuplicated(ids))
    stop("category_ids must contain unique classified category identifiers.", call. = FALSE)
  empty <- data.frame(category_id = character(), n_genesets_mapped = integer(),
    n_genesets_evaluable = integer(), n_genesets_unavailable = integer(),
    n_genesets_significant = integer(), n_pos_genesets = integer(),
    n_neg_genesets = integer(), n_zero_genesets = integer(),
    mean_NES = numeric(), median_NES = numeric(), p25_NES = numeric(), p75_NES = numeric(),
    positive_pct = numeric(), negative_pct = numeric(), zero_pct = numeric(),
    same_direction_pct = numeric(), mean_NES_direction = character(),
    gsea_padj_cutoff = numeric(), quantile_type = integer(), stringsAsFactors = FALSE)
  if (!length(ids)) return(empty)
  groups <- split(seq_len(nrow(x)), x$category_id)
  rows <- lapply(ids, function(id) {
    z <- x[groups[[id]], , drop = FALSE]
    if (nrow(z)) {
      # Identical repeated assignments do not add support. Conflicting results
      # cannot be resolved by selecting whichever row happened to come first.
      for (ix in split(seq_len(nrow(z)), z$pathway)) {
        if (nrow(unique(z[ix, c("NES", "padj"), drop = FALSE])) > 1L)
          stop("Conflicting duplicate GSEA member in category ", id, ": ",
            z$pathway[ix[[1L]]], call. = FALSE)
      }
      z <- z[!duplicated(z$pathway), , drop = FALSE]
    }
    evaluable <- is.finite(z$NES) & is.finite(z$padj)
    significant <- evaluable & z$padj <= gsea_padj_cutoff
    values <- z$NES[significant]
    n <- length(values); npos <- sum(values > 0); nneg <- sum(values < 0); nzero <- sum(values == 0)
    avg <- if (n) mean(values) else NA_real_
    qs <- if (n) stats::quantile(values, c(.25, .75), type = 7, names = FALSE) else c(NA_real_, NA_real_)
    pct <- function(count) if (n) 100 * count / n else NA_real_
    concordant <- if (!n || avg == 0) NA_real_ else if (avg > 0) pct(npos) else pct(nneg)
    data.frame(category_id = id, n_genesets_mapped = nrow(z),
      n_genesets_evaluable = sum(evaluable), n_genesets_unavailable = sum(!evaluable),
      n_genesets_significant = n, n_pos_genesets = npos, n_neg_genesets = nneg,
      n_zero_genesets = nzero, mean_NES = avg,
      median_NES = if (n) stats::median(values) else NA_real_,
      p25_NES = qs[[1L]], p75_NES = qs[[2L]], positive_pct = pct(npos),
      negative_pct = pct(nneg), zero_pct = pct(nzero), same_direction_pct = concordant,
      mean_NES_direction = if (!n) NA_character_ else if (avg > 0) "Positive NES" else
        if (avg < 0) "Negative NES" else "Zero mean NES",
      gsea_padj_cutoff = gsea_padj_cutoff, quantile_type = 7L, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

lisa_nes_variants_validate <- function(variants = c("clean", "percentages", "direction", "dispersion")) {
  allowed <- c("clean", "percentages", "direction", "dispersion")
  if (!is.character(variants) || !length(variants) || anyNA(variants) ||
      any(!variants %in% allowed) || anyDuplicated(variants))
    stop("NES variants must be a nonempty unique selection of clean, percentages, direction, dispersion.", call. = FALSE)
  variants
}

lisa_nes_variant_summary <- function(summary, stats) {
  if (!is.data.frame(summary) || !all(c("category_id", "mean_NES") %in% names(summary)))
    stop("NES variants require the existing category summary and mean_NES.", call. = FALSE)
  ids <- as.character(summary$category_id)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) || any(ids == "OTHER_UNCLASSIFIED"))
    stop("NES-variant summary must contain unique classified categories.", call. = FALSE)
  if ("classification_status" %in% names(summary) &&
      any(is.na(summary$classification_status) | summary$classification_status != "classified"))
    stop("Unclassified rows cannot be plotted as NES categories.", call. = FALSE)
  idx <- match(ids, stats$category_id)
  if (anyNA(idx)) stop("Category distribution is missing a summary category.", call. = FALSE)
  s <- stats[idx, , drop = FALSE]
  expected <- lisa_nes_number(summary$mean_NES, "mean_NES")
  supported <- s$n_genesets_significant > 0
  mismatch <- supported & (!is.finite(expected) |
    abs(expected - s$mean_NES) > 1e-12 * pmax(1, abs(expected), abs(s$mean_NES)))
  if (any(mismatch) || any(!supported & is.finite(expected)))
    stop("Existing category mean_NES disagrees with unique significant member support.", call. = FALSE)
  for (key in intersect(c("n_genesets", "n_genesets_significant"), names(summary))) {
    n <- lisa_nes_number(summary[[key]], key)
    if (any(!is.finite(n) | n != s$n_genesets_significant))
      stop("Existing category support count disagrees with NES distribution.", call. = FALSE)
  }
  # Canonical mean stays byte-for-byte numeric where supplied; the statistical
  # validation above does not replace an existing rounding/serialization choice.
  s$mean_NES <- expected
  s$category_display_name <- if ("category_display_name" %in% names(summary))
    as.character(summary$category_display_name) else ids
  s$category_display_name[is.na(s$category_display_name) | !nzchar(s$category_display_name)] <-
    ids[is.na(s$category_display_name) | !nzchar(s$category_display_name)]
  s$macrogroup_name <- if ("macrogroup_name" %in% names(summary)) as.character(summary$macrogroup_name) else rep("", nrow(s))
  s$category_color <- if ("color" %in% names(summary)) as.character(summary$color) else rep("#737373", nrow(s))
  s$category_color[is.na(s$category_color) | !nzchar(s$category_color)] <- "#737373"
  s$source_row_order <- seq_len(nrow(s))
  s
}

lisa_nes_variant_long <- function(s, side = "single", display_mean = s$mean_NES) {
  s$side <- rep(side, nrow(s))
  s$display_mean_NES <- display_mean
  s$endpoint_state <- ifelse(s$n_genesets_significant > 0, "significant_support",
    ifelse(is.finite(display_mean), "contextual_not_significant", "unavailable"))
  s$display_color <- ifelse(!is.finite(display_mean) | display_mean == 0, "#687787",
    ifelse(display_mean > 0, "#b2182b", "#2166ac"))
  s
}

lisa_nes_variant_object <- function(source, scope, labels, title, subtitle) {
  limits <- unlist(source[c("display_mean_NES", "median_NES", "p25_NES", "p75_NES")], use.names = FALSE)
  # All variants reserve the same room for near-point labels, including clean.
  span <- 1.52 * max(c(1, abs(limits[is.finite(limits)])))
  metadata <- list(schema_version = "1.0", scope = scope, side_labels = labels,
    presentation_version = "2.0", plot_set = "all",
    distribution_layout = "separate_lanes_v1",
    render_context = list(profile = "fresh_r_process_v1", png_device = "cairo",
      isolation = "One synchronous fresh Rscript --vanilla process for each output figure."),
    title = as.character(title), subtitle = as.character(subtitle), x_limits = c(-span, span),
    mean_definition = "Existing mean NES; significant unique classified member sets.",
    percentages = "Positive, negative and zero counts divided by significant finite-NES member count; absent members excluded.",
    quantiles = "P25-P75: stats::quantile type 7, between-gene-set dispersion, not a confidence interval.",
    contrast_rule = "Separate A/B distributions; existing descriptive A-B delta is unchanged; no delta percentiles.",
    row_order = "Supplied canonical category-summary row order; identical in every variant.")
  structure(list(source = source, metadata = metadata), class = "lisa_category_nes_variants")
}

lisa_prepare_category_nes_variants <- function(summary_df, gsea_all,
    gsea_padj_cutoff = .25, title = "LISA category NES", subtitle = "") {
  if (any(c("mean_NES_A", "mean_NES_B", "delta_mean_NES") %in% names(summary_df)))
    stop("Use the contrast NES-variant preparer for A/B summaries.", call. = FALSE)
  stats <- lisa_category_nes_statistics(gsea_all, gsea_padj_cutoff,
    as.character(summary_df$category_id))
  s <- lisa_nes_variant_summary(summary_df, stats)
  lisa_nes_variant_object(lisa_nes_variant_long(s), "single", c(single = "Analysis"), title, subtitle)
}

lisa_prepare_contrast_nes_variants <- function(contrast_summary, summary_a, summary_b,
    gsea_a, gsea_b, gsea_padj_cutoff = .25, label_a = "A", label_b = "B",
    title = "LISA category contrast", subtitle = "") {
  required <- c("category_id", "display_mean_NES_A", "display_mean_NES_B", "delta_mean_NES")
  if (!is.data.frame(contrast_summary) || !all(required %in% names(contrast_summary)))
    stop("Contrast variants need the existing display endpoints and A-B delta.", call. = FALSE)
  ids <- as.character(contrast_summary$category_id)
  if (anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids)) || any(ids == "OTHER_UNCLASSIFIED"))
    stop("Contrast NES categories must be unique and classified.", call. = FALSE)
  sides <- list(A = list(summary = summary_a, members = gsea_a), B = list(summary = summary_b, members = gsea_b))
  tables <- lapply(names(sides), function(side) {
    input <- sides[[side]]
    idx <- match(ids, as.character(input$summary$category_id))
    # An absent side must be represented explicitly by the canonical summary;
    # this layer does not fabricate a zero endpoint or missing category map.
    if (anyNA(idx)) stop("Contrast side is missing canonical category rows: ", side, call. = FALSE)
    summary <- input$summary[idx, , drop = FALSE]
    stats <- lisa_category_nes_statistics(input$members, gsea_padj_cutoff, ids)
    s <- lisa_nes_variant_summary(summary, stats)
    display <- lisa_nes_number(contrast_summary[[paste0("display_mean_NES_", side)]], "display_mean_NES")
    has <- s$n_genesets_significant > 0
    if (any(has & (!is.finite(display) | abs(display - s$mean_NES) >
        1e-12 * pmax(1, abs(display), abs(s$mean_NES)))))
      stop("Contrast significant endpoint differs from the side mean.", call. = FALSE)
    s <- lisa_nes_variant_long(s, side, display)
    s$delta_mean_NES <- lisa_nes_number(contrast_summary$delta_mean_NES, "delta_mean_NES")
    for (key in intersect(c("is_same_direction", "is_opposite_direction"), names(contrast_summary)))
      s[[key]] <- as.logical(contrast_summary[[key]])
    s
  })
  source <- do.call(rbind, tables); rownames(source) <- NULL
  lisa_nes_variant_object(source, "contrast", c(A = label_a, B = label_b), title, subtitle)
}

lisa_nes_direction_label <- function(source) {
  # paste0() otherwise recycles its fixed text to one row for zero-length
  # columns, which cannot be assigned to a zero-row plotting frame.
  if (!nrow(source)) return(character())
  pct <- function(x) ifelse(is.finite(x), paste0(sub("[.]0$", "",
    formatC(x, digits = 1, format = "f", decimal.mark = ".")), "%"), "NA")
  paste0("n=", source$n_genesets_significant, " | Pos ", pct(source$positive_pct),
    " / Neg ", pct(source$negative_pct),
    ifelse(source$n_zero_genesets > 0, paste0(" / Zero ", pct(source$zero_pct)), ""),
    " | Agree ", pct(source$same_direction_pct))
}

plot_lisa_category_nes_variant_legacy <- function(prepared, variant = "clean") {
  if (!inherits(prepared, "lisa_category_nes_variants"))
    stop("Expected prepared lisa_category_nes_variants.", call. = FALSE)
  lisa_nes_variants_validate(variant)
  if (length(variant) != 1L) stop("Select one NES variant to plot.", call. = FALSE)
  x <- prepared$source; m <- prepared$metadata
  # Local bindings avoid package-check global-variable notes from ggplot NSE.
  plot_y <- display_mean_NES <- median_NES <- p25_NES <- p75_NES <- NULL
  display_color <- n_genesets_significant <- direction_label <- side <- NULL
  ids <- unique(x$category_id)
  ncat <- length(ids)
  x$plot_y <- ncat - match(x$category_id, ids) + 1
  paired <- identical(m$scope, "contrast")
  if (paired) x$plot_y <- x$plot_y + ifelse(x$side == "A", .15, -.15)
  labels <- x$category_display_name[match(ids, x$category_id)]
  visible <- x[is.finite(x$display_mean_NES), , drop = FALSE]
  supported <- visible[visible$n_genesets_significant > 0, , drop = FALSE]
  context <- visible[visible$n_genesets_significant == 0, , drop = FALSE]
  median_rows <- x[is.finite(x$median_NES) & x$n_genesets_significant > 0, , drop = FALSE]
  median_rows$direction_label <- lisa_nes_direction_label(median_rows)
  caption <- if (variant == "clean") "Mean NES is a descriptive category summary, not a category significance test." else
    "Filled symbols: mean NES. Black diamonds: median. Pos/Neg: positive/negative NES. Percentages use significant finite-NES member sets; Agree: direction matches the mean's sign."
  if (variant == "dispersion") caption <- paste(caption,
    "Thick segments: P25-P75 between-gene-set dispersion (type 7), not confidence intervals.")
  if (paired) caption <- paste(caption, "A/B supports remain separate; hollow endpoints are nonsignificant contextual means.")
  p <- ggplot2::ggplot(x, ggplot2::aes(x = display_mean_NES, y = plot_y)) +
    ggplot2::geom_vline(xintercept = 0, color = "#a6afb8", linewidth = .35, linetype = "dashed")
  if (paired) {
    a <- x[x$side == "A", , drop = FALSE]; b <- x[x$side == "B", , drop = FALSE]
    connectors <- data.frame(from = a$display_mean_NES, to = b$display_mean_NES,
      y_a = a$plot_y, y_b = b$plot_y)
    connectors <- connectors[is.finite(connectors$from) & is.finite(connectors$to), , drop = FALSE]
    from <- to <- y_a <- y_b <- NULL
    p <- p + ggplot2::geom_segment(data = connectors,
      ggplot2::aes(x = from, xend = to, y = y_a, yend = y_b), inherit.aes = FALSE,
      linewidth = .5, color = "#c1c8ce")
  } else p <- p + ggplot2::geom_segment(data = visible,
    ggplot2::aes(x = 0, xend = display_mean_NES, yend = plot_y), linewidth = .45, color = "#c1c8ce")
  if (variant == "dispersion") p <- p + ggplot2::geom_segment(data = median_rows,
    ggplot2::aes(x = p25_NES, xend = p75_NES, yend = plot_y, color = display_color),
    linewidth = 3, alpha = .3)
  if (paired) {
    # No endpoints means no shape data or A/B legend. A named manual scale on
    # two empty layers warns about missing shared levels in recent ggplot2.
    if (nrow(visible)) p <- p + ggplot2::geom_point(data = supported,
      ggplot2::aes(fill = display_color, shape = side, size = n_genesets_significant), color = "#293946", stroke = .35) +
      ggplot2::geom_point(data = context, ggplot2::aes(color = display_color, shape = side),
        size = 2.5, fill = "white", stroke = .65) +
      ggplot2::scale_shape_manual(values = c(A = 21, B = 24), labels = m$side_labels,
        name = "Analysis", guide = ggplot2::guide_legend(order = 2L))
  } else p <- p + ggplot2::geom_point(data = supported,
    ggplot2::aes(fill = display_color, size = n_genesets_significant), shape = 21,
    color = "#293946", stroke = .35)
  if (variant != "clean") {
    p <- p + ggplot2::geom_point(data = median_rows, ggplot2::aes(x = median_NES),
      shape = 18, size = 2.3, color = "#18242e") +
      ggplot2::geom_text(data = median_rows,
        ggplot2::aes(x = m$x_limits[[2L]] * 1.04, label = direction_label),
        hjust = 0, size = if (paired) 2.5 else 2.8, color = "#293946")
  }
  if (!ncat) p <- p + ggplot2::annotate("text", x = 0, y = 1,
    label = "No classified categories", color = "#687787")
  p <- p + ggplot2::scale_color_identity() + ggplot2::scale_fill_identity() +
    ggplot2::scale_size_continuous(range = c(2.2, 6.2), name = "Significant gene sets",
      guide = ggplot2::guide_legend(order = 1L),
      breaks = function(limits) {
        if (length(limits) != 2L || any(!is.finite(limits))) return(numeric())
        ticks <- sort(unique(c(round(pretty(limits, n = 4L)),
          ceiling(limits[[1L]]), floor(limits[[2L]]))))
        ticks[ticks >= max(1, ceiling(limits[[1L]])) & ticks <= floor(limits[[2L]])]
      }, labels = function(ticks) format(ticks, scientific = FALSE, trim = TRUE)) +
    ggplot2::scale_y_continuous(breaks = rev(seq_len(ncat)), labels = labels,
      limits = c(.45, max(1, ncat) + .55), expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_x_continuous(limits = NULL, expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_cartesian(xlim = m$x_limits, clip = "off") +
    ggplot2::labs(title = m$title, subtitle = m$subtitle, x = "NES", y = NULL,
      caption = paste(strwrap(caption, width = 122), collapse = "\n")) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(), legend.position = "bottom",
      plot.caption = ggplot2::element_text(hjust = 0, size = 8),
      plot.margin = ggplot2::margin(10, 230, 10, 10),
      plot.background = ggplot2::element_rect(fill = "white", color = NA))
  attr(p, "lisa_nes_variant") <- variant
  attr(p, "lisa_nes_source") <- x
  p
}

lisa_nes_plot_sets_validate <- function(plot_sets, scope) {
  if (is.null(plot_sets)) plot_sets <- "all"
  allowed <- if (identical(scope, "contrast")) c("all", "same_direction", "opposite_direction") else "all"
  if (!is.character(plot_sets) || !length(plot_sets) || anyNA(plot_sets) ||
      anyDuplicated(plot_sets) || any(!plot_sets %in% allowed))
    stop("Invalid nonempty unique NES plot-set selection for this scope.", call. = FALSE)
  plot_sets
}

lisa_nes_select_plot_set <- function(prepared, plot_set = "all") {
  lisa_nes_plot_sets_validate(plot_set, prepared$metadata$scope)
  if (length(plot_set) != 1L) stop("Select one NES plot set.", call. = FALSE)
  result <- prepared
  if (plot_set != "all") {
    key <- paste0("is_", plot_set)
    if (!key %in% names(prepared$source))
      stop("NES subset requires canonical contrast membership: ", key, call. = FALSE)
    keep <- prepared$source[[key]]
    if (!is.logical(keep) || anyNA(keep)) stop("Invalid canonical NES subset membership.", call. = FALSE)
    result$source <- prepared$source[keep, , drop = FALSE]
  }
  result$metadata$plot_set <- plot_set
  if (plot_set != "all") result$metadata$title <- paste(prepared$metadata$title,
    switch(plot_set, same_direction = "Same-direction categories",
      opposite_direction = "Opposite-direction categories"), sep = " - ")
  # The same common x-limits are kept across variants AND contrast subsets.
  result
}

lisa_nes_compact_percentage_label <- function(source) {
  if (!nrow(source)) return(character())
  label <- paste0(ifelse(is.finite(source$same_direction_pct),
    paste0(round(source$same_direction_pct), "%"), "Direction undefined"),
    " (n=", source$n_genesets_significant, ")")
  label[source$n_genesets_significant == 0L] <- paste0("NS (n_eval=",
    source$n_genesets_evaluable[source$n_genesets_significant == 0L], ")")
  label
}

lisa_nes_distribution_lane_offset <- function(metadata, side = "single") {
  # An explicit profile distinguishes the visual correction from earlier
  # emitted 2.0 settings, which keep their original coincident annotations.
  profile <- metadata$distribution_layout
  if (is.null(profile)) return(if (side == "single") 0 else if (side == "A") .14 else -.14)
  if (!identical(profile, "separate_lanes_v1"))
    stop("Unsupported NES distribution layout.", call. = FALSE)
  if (side == "B") -.4 else .4
}

plot_lisa_category_nes_variant <- function(prepared, variant = "clean") {
  if (!inherits(prepared, "lisa_category_nes_variants"))
    stop("Expected prepared lisa_category_nes_variants.", call. = FALSE)
  lisa_nes_variants_validate(variant)
  if (length(variant) != 1L) stop("Select one NES variant to plot.", call. = FALSE)
  version <- prepared$metadata$presentation_version
  # Archived settings retain their original layout and meaning. In particular,
  # direction never becomes the new percentages-only view.
  if (is.null(version)) {
    if (variant == "percentages") stop("Legacy NES settings have no percentages-only view.", call. = FALSE)
    return(plot_lisa_category_nes_variant_legacy(prepared, variant))
  }
  if (!identical(version, "2.0")) stop("Unsupported NES presentation version.", call. = FALSE)
  x <- prepared$source; m <- prepared$metadata
  category_plot_id <- display_mean_NES <- category_color <- n_genesets_significant <- NULL
  point_fill <- side <- label_x <- point_label <- label_hjust <- NULL
  median_NES <- p25_NES <- p75_NES <- from <- to <- NULL
  macrogroup_name <- y <- NULL
  ids <- unique(x$category_id)
  x$category_plot_id <- factor(x$category_id, levels = rev(ids))
  group <- as.character(x$macrogroup_name)
  group[is.na(group) | !nzchar(group)] <- "Categories"
  x$macrogroup_name <- factor(group, levels = unique(group))
  if (!"category_color" %in% names(x)) x$category_color <- rep("#737373", nrow(x))
  labels <- stats::setNames(x$category_display_name[match(ids, x$category_id)], ids)
  paired <- identical(m$scope, "contrast")
  lisa_nes_distribution_lane_offset(m, if (paired) "A" else "single")
  visible <- x[is.finite(x$display_mean_NES), , drop = FALSE]
  visible$point_fill <- visible$category_color
  if (paired && nrow(visible)) {
    is_a <- visible$side == "A"
    visible$point_fill[is_a] <- vapply(visible$category_color[is_a],
      lighten_color, character(1L), amount = .55)
  }
  supported <- visible[visible$n_genesets_significant > 0L, , drop = FALSE]
  context <- visible[visible$n_genesets_significant == 0L, , drop = FALSE]
  medians <- x[is.finite(x$median_NES) & x$n_genesets_significant > 0L, , drop = FALSE]
  visible$point_label <- lisa_nes_compact_percentage_label(visible)
  span <- diff(m$x_limits)
  if (paired && nrow(visible)) {
    a <- x[x$side == "A", , drop = FALSE]
    b <- x[x$side == "B", , drop = FALSE]
    b <- b[match(a$category_id, b$category_id), , drop = FALSE]
    a_left <- a$display_mean_NES <= b$display_mean_NES
    a_left[is.na(a_left)] <- TRUE
    label_sides <- stats::setNames(ifelse(a_left, -1, 1), a$category_id)
    placement <- unname(label_sides[visible$category_id]) * ifelse(visible$side == "A", 1, -1)
  } else placement <- ifelse(visible$display_mean_NES < 0, -1, 1)
  visible$label_x <- visible$display_mean_NES + placement * span * .018
  visible$label_hjust <- ifelse(placement < 0, 1, 0)
  caption <- "Symbols: mean NES; size: significant gene sets. Category colours and supercategory bands follow the LISA map. Mean NES is descriptive, not a category significance test."
  if (variant != "clean") caption <- paste(caption,
    "Labels: percentage of significant finite-NES member sets agreeing with the mean's sign (n is that denominator).")
  if (variant %in% c("direction", "dispersion")) caption <- paste(caption, "Black diamonds: median NES.")
  if (variant == "dispersion") caption <- paste(caption,
    "Thick segments: P25-P75 between-gene-set dispersion (type 7), not confidence intervals.")
  if (variant %in% c("direction", "dispersion") && identical(m$distribution_layout, "separate_lanes_v1"))
    caption <- paste(caption, if (variant == "dispersion") "Medians and intervals" else "Medians",
      if (paired) "sit above (A) and below (B) the category's mean/support row." else
        "sit above the category's mean/support row.")
  if (paired) caption <- paste(caption,
    "Circle: A; square: B. Hollow endpoints: nonsignificant contextual means. A/B distributions remain separate.")
  p <- ggplot2::ggplot(x, ggplot2::aes(y = category_plot_id)) +
    ggplot2::geom_blank(ggplot2::aes(x = 0)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = .35)
  if (!length(ids)) {
    p <- ggplot2::ggplot(data.frame(x = 0, y = 1), ggplot2::aes(x = x, y = y)) +
      ggplot2::annotate("text", x = 0, y = 1, label = "No classified categories in this selection", color = "grey40") +
      ggplot2::theme_void() + ggplot2::labs(title = m$title, subtitle = m$subtitle,
        caption = paste(strwrap(caption, width = 125), collapse = "\n"))
    attr(p, "lisa_nes_variant") <- variant; attr(p, "lisa_nes_source") <- x
    return(p)
  }
  if (paired) {
    a <- x[x$side == "A", , drop = FALSE]; b <- x[x$side == "B", , drop = FALSE]
    b <- b[match(a$category_id, b$category_id), , drop = FALSE]
    connectors <- a
    connectors$from <- a$display_mean_NES; connectors$to <- b$display_mean_NES
    connectors <- connectors[is.finite(connectors$from) & is.finite(connectors$to), , drop = FALSE]
    p <- p + ggplot2::geom_segment(data = connectors,
      ggplot2::aes(x = from, xend = to, yend = category_plot_id, color = category_color),
      linewidth = .55, alpha = .7)
  } else p <- p + ggplot2::geom_segment(data = visible,
    ggplot2::aes(x = 0, xend = display_mean_NES, yend = category_plot_id), color = "grey70", linewidth = .35)
  if (variant == "dispersion") {
    for (current in if (paired) c("A", "B") else "single") {
      rows <- medians[medians$side == current, , drop = FALSE]
      offset <- lisa_nes_distribution_lane_offset(m, current)
      p <- p + ggplot2::geom_segment(data = rows,
        ggplot2::aes(x = p25_NES, xend = p75_NES, yend = category_plot_id, color = category_color),
        position = ggplot2::position_nudge(y = offset), linewidth = 3, alpha = .3)
    }
  }
  if (paired && nrow(visible)) {
    # A/B means remain original common-row dumbbells. Intervals and medians
    # have separate lanes; endpoint percentages stay beside their own mean.
    p <- p + ggplot2::geom_point(data = supported,
      ggplot2::aes(x = display_mean_NES, fill = point_fill, shape = side, size = n_genesets_significant),
      color = "grey18", stroke = .28, alpha = .96) +
      ggplot2::geom_point(data = context,
        ggplot2::aes(x = display_mean_NES, shape = side, size = n_genesets_significant),
        fill = NA, color = "grey35", stroke = .75, alpha = .6) +
      ggplot2::scale_shape_manual(values = c(A = 21, B = 22), limits = c("A", "B"),
        labels = unname(m$side_labels), name = "Analysis", drop = FALSE,
        guide = ggplot2::guide_legend(order = 2L))
  } else if (!paired) p <- p + ggplot2::geom_point(data = supported,
    ggplot2::aes(x = display_mean_NES, fill = point_fill, size = n_genesets_significant),
    shape = 21, color = "grey25", stroke = .25)
  if (variant %in% c("direction", "dispersion")) {
    for (current in if (paired) c("A", "B") else "single") {
      rows <- medians[medians$side == current, , drop = FALSE]
      offset <- lisa_nes_distribution_lane_offset(m, current)
      p <- p + ggplot2::geom_point(data = rows, ggplot2::aes(x = median_NES),
        position = ggplot2::position_nudge(y = offset), shape = 18, size = 2.3, color = "#18242e")
    }
  }
  if (variant != "clean") {
    for (current in if (paired) c("A", "B") else "single") {
      rows <- visible[visible$side == current, , drop = FALSE]
      offset <- if (!paired) 0 else if (current == "A") .2 else -.2
      p <- p + ggplot2::geom_text(data = rows,
        ggplot2::aes(x = label_x, label = point_label, hjust = label_hjust),
        position = ggplot2::position_nudge(y = offset), size = if (paired) 2.5 else 2.7,
        fontface = "bold", color = "grey20", lineheight = .9, show.legend = FALSE)
    }
  }
  p <- p + ggplot2::scale_color_identity() + ggplot2::scale_fill_identity() +
    ggplot2::scale_size_continuous(range = c(1.7, 6.2), name = "Significant gene sets",
      guide = ggplot2::guide_legend(order = 1L), breaks = function(limits) {
        if (length(limits) != 2L || any(!is.finite(limits))) return(numeric())
        ticks <- sort(unique(c(round(pretty(limits, n = 4L)), ceiling(limits[[1L]]), floor(limits[[2L]]))))
        ticks[ticks >= max(1, ceiling(limits[[1L]])) & ticks <= floor(limits[[2L]])]
      }, labels = function(ticks) format(ticks, scientific = FALSE, trim = TRUE)) +
    ggplot2::scale_y_discrete(labels = labels, drop = TRUE) +
    ggplot2::scale_x_continuous(limits = m$x_limits, expand = ggplot2::expansion(mult = c(.02, .02))) +
    ggplot2::facet_grid(macrogroup_name ~ ., scales = "free_y", space = "free_y", switch = "y", drop = TRUE) +
    ggplot2::labs(title = m$title, subtitle = m$subtitle,
      x = "Mean NES across significant LISA-mapped gene sets", y = NULL,
      caption = paste(strwrap(caption, width = 125), collapse = "\n")) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom", strip.placement = "outside",
      strip.background.y = ggplot2::element_rect(fill = "grey94", color = "grey80"),
      strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 8, face = "bold"),
      panel.spacing.y = grid::unit(.25, "lines"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8),
      plot.margin = ggplot2::margin(10, 12, 10, 10),
      plot.background = ggplot2::element_rect(fill = "white", color = NA))
  attr(p, "lisa_nes_variant") <- variant; attr(p, "lisa_nes_source") <- x
  p
}

lisa_nes_save_figure <- function(plot, output_path, settings) {
  ext <- tolower(tools::file_ext(output_path))
  if (!ext %in% c("png", "pdf", "svg")) stop("NES recipe output must be PNG, PDF or SVG.", call. = FALSE)
  context <- settings$render_context
  if (is.null(settings$presentation_version)) {
    device <- switch(ext, png = "png", pdf = grDevices::pdf, svg = grDevices::svg)
    ggplot2::ggsave(output_path, plot, device = device, width = settings$width,
      height = settings$height, dpi = settings$dpi, limitsize = FALSE, bg = "white")
    return(invisible(output_path))
  }
  if (!identical(context$profile, "fresh_r_process_v1") || !identical(context$png_device, "cairo"))
    stop("Unsupported NES rendering context.", call. = FALSE)
  # Device/font caches can depend on previous drawing resolutions. Every new
  # figure therefore starts in a separate R process, both at origin and replay.
  # No global options, other devices or caller jobs are changed.
  folder <- tempfile("lisa-nes-device-")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  bundle <- file.path(folder, "figure.rds")
  paths <- file.path(folder, "libraries.rds")
  script <- file.path(folder, "render.R")
  log <- file.path(folder, "render.log")
  saveRDS(.libPaths(), paths)
  saveRDS(list(plot = plot, settings = settings,
    output = normalizePath(output_path, winslash = "/", mustWork = FALSE), format = ext), bundle)
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    ".libPaths(readRDS(args[[1L]]))",
    "suppressPackageStartupMessages(library(ggplot2))",
    "item <- readRDS(args[[2L]])",
    "if (item$format == 'png' && !capabilities('cairo')) stop('The declared NES PNG context requires Cairo.')",
    "device <- switch(item$format, png = function(filename, width, height, ...) grDevices::png(filename, width = width, height = height, units = 'in', res = item$settings$dpi, type = 'cairo', ...), pdf = grDevices::pdf, svg = grDevices::svg)",
    "ggplot2::ggsave(item$output, item$plot, device = device, width = item$settings$width, height = item$settings$height, dpi = item$settings$dpi, limitsize = FALSE, bg = 'white')"
  ), script)
  executable <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  if (!file.exists(executable)) stop("The pinned Rscript executable is unavailable.", call. = FALSE)
  status <- system2(executable, args = c("--vanilla", shQuote(script), shQuote(paths), shQuote(bundle)),
    stdout = log, stderr = log)
  if (!identical(status, 0L) || !file.exists(output_path))
    stop("Isolated NES figure rendering failed: ", paste(readLines(log, warn = FALSE), collapse = "\n"), call. = FALSE)
  invisible(output_path)
}

render_lisa_category_nes_variants <- function(prepared, output_dir,
    variants = c("clean", "percentages", "direction", "dispersion"), prefix = "category_nes",
    formats = c("png", "pdf", "svg"), width = 13, height = NULL, dpi = 180,
    plot_sets = NULL) {
  if (!inherits(prepared, "lisa_category_nes_variants")) stop("Expected prepared NES variants.", call. = FALSE)
  variants <- lisa_nes_variants_validate(variants)
  plot_sets <- lisa_nes_plot_sets_validate(plot_sets, prepared$metadata$scope)
  if (!is.character(prefix) || length(prefix) != 1L || !grepl("^[A-Za-z0-9_-]+$", prefix))
    stop("NES-variant prefix must be a safe file stem.", call. = FALSE)
  if (!length(formats) || anyNA(formats) || any(!formats %in% c("png", "pdf", "svg")) || anyDuplicated(formats))
    stop("NES-variant formats must be PNG, PDF or SVG.", call. = FALSE)
  row_height <- if (identical(prepared$metadata$scope, "contrast") &&
    identical(prepared$metadata$distribution_layout, "separate_lanes_v1")) .50 else .34
  if (is.null(height)) height <- max(4.5, 3.2 + row_height * length(unique(prepared$source$category_id)) +
    .1 * length(unique(prepared$source$macrogroup_name)))
  dimensions <- c(width, height, dpi)
  if (length(dimensions) != 3L || any(!is.finite(dimensions) | dimensions <= 0))
    stop("Invalid NES-variant dimensions.", call. = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  results <- lapply(plot_sets, function(plot_set) {
  selected <- lisa_nes_select_plot_set(prepared, plot_set)
  selected_prefix <- if (plot_set == "all") prefix else paste0(prefix, "_", plot_set)
  source_path <- file.path(output_dir, paste0(selected_prefix, "_source.tsv"))
  stored_source <- selected$source
  for (field in names(stored_source)[vapply(stored_source, is.numeric, logical(1))]) {
    value <- stored_source[[field]]
    stored_source[[field]] <- ifelse(is.na(value), "", formatC(value, format = "g", digits = 17))
  }
  utils::write.table(stored_source, source_path, sep = "\t", row.names = FALSE,
    quote = TRUE, na = "", fileEncoding = "UTF-8")
  rows <- lapply(variants, function(variant) {
    settings <- selected$metadata
    settings$variant <- variant; settings$width <- width; settings$height <- height; settings$dpi <- dpi
    settings$source_sha256 <- digest::digest(file = source_path, algo = "sha256")
    settings_path <- file.path(output_dir, paste0(selected_prefix, "_", variant, "_settings.json"))
    jsonlite::write_json(settings, settings_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
    p <- plot_lisa_category_nes_variant(selected, variant)
    files <- vapply(formats, function(format) {
      path <- file.path(output_dir, paste0(selected_prefix, "_", variant, ".", format))
      lisa_nes_save_figure(p, path, settings)
      path
    }, character(1))
    data.frame(variant = variant, plot_set = plot_set, format = formats, path = unname(files),
      source = source_path, settings = settings_path, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
  })
  out <- do.call(rbind, results); rownames(out) <- NULL
  invisible(out)
}

lisa_reproduce_category_nes_variant <- function(source_path, settings_path, output_path) {
  settings <- jsonlite::read_json(settings_path, simplifyVector = TRUE)
  if (!identical(settings$schema_version, "1.0") ||
      !settings$scope %in% c("single", "contrast") ||
      length(settings$x_limits) != 2L || any(!is.finite(settings$x_limits)) ||
      settings$x_limits[[1L]] >= settings$x_limits[[2L]])
    stop("Invalid NES-variant settings.", call. = FALSE)
  if (!identical(digest::digest(file = source_path, algo = "sha256"), settings$source_sha256))
    stop("NES-variant source hash mismatch.", call. = FALSE)
  source <- utils::read.delim(source_path, sep = "\t", quote = '"', comment.char = "",
    colClasses = "character", check.names = FALSE, na.strings = NULL, stringsAsFactors = FALSE)
  numeric <- c("n_genesets_mapped", "n_genesets_evaluable", "n_genesets_unavailable",
    "n_genesets_significant", "n_pos_genesets", "n_neg_genesets", "n_zero_genesets",
    "mean_NES", "median_NES", "p25_NES", "p75_NES", "positive_pct", "negative_pct",
    "zero_pct", "same_direction_pct", "gsea_padj_cutoff", "quantile_type", "source_row_order",
    "display_mean_NES", "delta_mean_NES")
  required <- setdiff(numeric, "delta_mean_NES")
  if (!all(c(required, "category_id", "category_display_name", "side", "display_color") %in% names(source)))
    stop("Incomplete NES-variant source.", call. = FALSE)
  for (field in intersect(numeric, names(source))) source[[field]] <- lisa_nes_number(source[[field]], field)
  if (!is.null(settings$presentation_version) &&
      !all(c("macrogroup_name", "category_color") %in% names(source)))
    stop("Grouped NES settings require explicit supercategory and colour source fields.", call. = FALSE)
  for (field in intersect(c("is_same_direction", "is_opposite_direction"), names(source))) {
    if (any(!source[[field]] %in% c("TRUE", "FALSE"))) stop("Invalid canonical NES subset flag.", call. = FALSE)
    source[[field]] <- source[[field]] == "TRUE"
  }
  prepared <- structure(list(source = source, metadata = settings), class = "lisa_category_nes_variants")
  plot <- plot_lisa_category_nes_variant(prepared, settings$variant)
  lisa_nes_save_figure(plot, output_path, settings)
  invisible(output_path)
}
