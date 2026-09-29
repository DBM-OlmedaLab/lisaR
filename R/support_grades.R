# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé. GPL-3 or later.

# One editorial rule for tables, figures and portable recipes. Counts are
# bounded R integers, so every cross-product below is an exact double integer.
# Decimal long division truncates the assertion; no epsilon or round() is used.
lisa_support_grade <- function(d, n) {
  if (!is.numeric(d) || !is.numeric(n) || length(d) != length(n) ||
      anyNA(n) || any(!is.finite(n) | n < 0 | n != floor(n) | n > .Machine$integer.max) ||
      any(!is.na(d) & (!is.finite(d) | d < 0 | d != floor(d) | d > n)) ||
      any(n > 0 & is.na(d)))
    stop("Support grades require integer counts 0 <= d <= N <= .Machine$integer.max; missing d only for N=0.", call. = FALSE)
  positive <- n > 0 & !is.na(d) & d > 0
  grade <- rep("", length(n))
  grade[positive] <- paste0("S", 1L + vapply(which(positive), function(i)
    sum(100 * d[[i]] >= c(10, 25, 50, 75) * n[[i]]), integer(1)))
  percent <- vapply(seq_along(n), function(i) {
    if (!n[[i]]) return(NA_character_)
    numerator <- 100 * d[[i]]
    whole <- floor(numerator / n[[i]])
    remainder <- numerator - whole * n[[i]]
    decimals <- character()
    while (remainder > 0 && (length(decimals) < 2L || (whole == 0 && all(decimals == "0")))) {
      numerator <- remainder * 10
      digit <- floor(numerator / n[[i]])
      remainder <- numerator - digit * n[[i]]
      decimals <- c(decimals, as.character(digit))
    }
    fraction <- sub("0+$", "", paste(decimals, collapse = ""))
    paste0(format(whole, scientific = FALSE, trim = TRUE), if (nzchar(fraction)) paste0(".", fraction), "%")
  }, character(1))
  labels <- c(S1 = "Enrichment supported; minimum below 10%",
    S2 = "Supported in at least one tenth", S3 = "Supported in at least one quarter",
    S4 = "Supported in at least one half", S5 = "Supported in at least three quarters")
  label <- ifelse(n == 0, "Not evaluable: no eligible member with a usable raw P value",
    ifelse(positive, unname(labels[grade]), "No minimum greater than zero could be supported"))
  count <- ifelse(n > 0, paste0(d, "/", n), NA_character_)
  annotation <- ifelse(positive, paste0(grade, " \u00b7 \u2265", percent, " \u00b7 ", count), "")
  data.frame(support_schema_version = rep("support-grades-v1", length(n)), support_grade = grade,
    support_label = label, minimum_support_pct_text = percent,
    support_count = count, support_annotation = annotation, stringsAsFactors = FALSE)
}

lisa_support_key <- function(x, unique = FALSE) {
  fields <- c("analysis_id", "collection", "category_id")
  lisa_require_columns(x, fields, "support keys")
  if (any(vapply(x[fields], function(z) anyNA(z) || any(!nzchar(as.character(z))) ||
      any(grepl("\r", z, fixed = TRUE)), logical(1))))
    stop("Support keys must be explicit nonempty analysis/collection/category identifiers.", call. = FALSE)
  key <- do.call(paste, c(x[fields], sep = "\r"))
  if (unique && anyDuplicated(key)) stop("Duplicate analysis/collection/category support key.", call. = FALSE)
  key
}

lisa_category_support_presentation <- function(x) {
  lisa_support_key(x, unique = TRUE)
  if (any(x$collection == "HALLMARKS")) stop("HALLMARKS is excluded from category support grades.", call. = FALSE)
  if ("confidence_level" %in% names(x) && any(x$confidence_level != .95))
    stop("Support presentation v1 requires simultaneous 95% confidence.", call. = FALSE)
  grades <- lisa_support_grade(x$minimum_enriched_sets, x$n_sets_evaluable)
  positive <- nzchar(grades$support_grade)
  if ("significant" %in% names(x) && any(positive != (x$significant %in% TRUE)))
    stop("Support counts and category significance disagree.", call. = FALSE)
  for (field in names(grades)) x[[field]] <- grades[[field]]
  fields <- c("n_positive_NES", "n_negative_NES", "n_zero_NES", "n_missing_NES")
  lisa_require_columns(x, fields, "descriptive NES counts")
  if (any(rowSums(x[fields]) != x$n_sets_evaluable))
    stop("Observed NES counts must account for all evaluable category members.", call. = FALSE)
  x$observed_NES_pattern <- ifelse(x$n_sets_evaluable == 0, "Not evaluable",
    ifelse(x$n_positive_NES > 0 & x$n_negative_NES > 0, "Observed NES of both signs",
    ifelse(x$n_positive_NES > 0, "Observed positive NES",
    ifelse(x$n_negative_NES > 0, "Observed negative NES",
    ifelse(x$n_zero_NES > 0, "Observed zero NES only", "No observed NES available")))))
  x$observed_NES_counts <- paste0("+", x$n_positive_NES, " / -", x$n_negative_NES,
    " / zero ", x$n_zero_NES, " / missing ", x$n_missing_NES, "; N=", x$n_sets_evaluable)
  x
}

lisa_support_legend <- function() {
  paste0("S1: >0 to <10%; S2: 10 to <25%; S3: 25 to <50%; S4: 50 to <75%; S5: 75 to 100% (bands of the lower bound).\n",
    "All levels use simultaneous 95% confidence within this analysis; they do not indicate effect size, confidence strength or direction.")
}

# Join by scientific identity, never the potentially duplicated display label.
lisa_support_lollipop_source <- function(source, categories, analysis_id, collection) {
  if (length(analysis_id) != 1L || length(collection) != 1L)
    stop("Lollipop support requires one explicit analysis and collection.", call. = FALSE)
  for (field in c("analysis_id", "collection")) {
    expected <- get(field)
    if (field %in% names(source) && any(is.na(source[[field]]) | source[[field]] != expected))
      stop("Lollipop source identity disagrees with the support table.", call. = FALSE)
    source[[field]] <- expected
  }
  source_key <- lisa_support_key(source, unique = TRUE)
  categories <- lisa_category_support_presentation(categories)
  if (any(categories$analysis_id != analysis_id | categories$collection != collection))
    stop("Support table contains a different analysis or collection.", call. = FALSE)
  category_key <- lisa_support_key(categories, unique = TRUE)
  if (any(!source_key %in% category_key)) stop("Lollipop category is absent from the inference catalogue.", call. = FALSE)
  # Native summaries retain empty categories. A legacy source may omit them:
  # retain their rows, but never borrow a contextual mean or manufacture zero.
  ix <- match(category_key, source_key)
  x <- source[ix, , drop = FALSE]
  missing <- is.na(ix)
  for (field in names(source)[vapply(source, function(z) length(unique(z)) == 1L, logical(1))])
    x[[field]][missing] <- source[[field]][[1L]]
  for (field in c("mean_NES", "n_genesets", "same_direction_pct", "consistency"))
    if (field %in% names(x)) x[[field]][missing] <- NA
  for (field in names(categories)) x[[field]] <- categories[[field]]
  if ("macrogroup_order" %in% names(x)) x$macrogroup_order[missing] <- 999
  if ("category_order_within_macrogroup" %in% names(x)) x$category_order_within_macrogroup[missing] <- seq_len(sum(missing))
  if ("macrogroup_name" %in% names(x)) x$macrogroup_name[missing] <- "Without a lollipop summary"
  variant <- unique(as.character(source$annotation_variant))
  if (length(variant) != 1L || !variant %in% c("plain", "direction"))
    stop("Support decoration requires a historical plain or direction lollipop source.", call. = FALSE)
  x$annotation_variant <- if (variant == "plain") "support_v1" else "direction_support_v1"
  x$source_row_order <- seq_len(nrow(x))
  x$selected_for_plot <- TRUE
  x$has_lollipop_geometry <- is.finite(x$mean_NES) & is.finite(x$n_genesets) & x$n_genesets > 0
  x$figure_width <- 17
  x$figure_height <- max(7.5, 4 + .30 * nrow(x))
  rownames(x) <- NULL
  x
}

lisa_plot_support_lollipop <- function(summary_df, group_by_supracategory, plot_order,
    direction = FALSE) {
  lisa_assert_classified_plot_rows(summary_df, "support lollipop")
  # Recompute all displayed assertions with the same versioned helper.
  x <- lisa_category_support_presentation(summary_df)
  x <- prepare_category_plot_df(x, plot_order)
  x$key <- lisa_support_key(x, unique = TRUE)
  x$category_plot_label <- factor(x$key, levels = rev(x$key))
  y_labels <- stats::setNames(x$category_display_name, x$key)
  panels <- c("Minimum support (95%)", "LISA mean NES (individual gene-set threshold)")
  support <- x
  support$panel <- factor(panels[[1]], levels = panels)
  geometry <- x[is.finite(x$mean_NES) & is.finite(x$n_genesets) & x$n_genesets > 0, , drop = FALSE]
  geometry$panel <- factor(rep(panels[[2]], nrow(geometry)), levels = panels)
  rows <- x
  rows$panel <- factor(panels[[2]], levels = panels)
  p <- ggplot2::ggplot() +
    ggplot2::geom_blank(data = support, ggplot2::aes(x = 0, y = category_plot_label)) +
    ggplot2::geom_blank(data = support, ggplot2::aes(x = 100, y = category_plot_label)) +
    ggplot2::geom_blank(data = rows, ggplot2::aes(x = 0, y = category_plot_label)) +
    ggplot2::geom_text(data = support, ggplot2::aes(x = 50, y = category_plot_label, label = support_annotation),
      size = 3.1, color = "#3d4855") +
    ggplot2::geom_vline(data = rows, ggplot2::aes(xintercept = 0), linetype = "dashed", color = "grey55", linewidth = .35) +
    ggplot2::geom_segment(data = geometry, ggplot2::aes(x = 0, xend = mean_NES, y = category_plot_label, yend = category_plot_label),
      color = "grey70", linewidth = .35) +
    ggplot2::geom_point(data = geometry, ggplot2::aes(x = mean_NES, y = category_plot_label, size = n_genesets, fill = color),
      alpha = .95, shape = 21, color = "grey25", stroke = .25) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_size_continuous(range = c(1.5, 7), breaks = c(0, 1, 5, 10, 20)) +
    ggplot2::scale_y_discrete(labels = y_labels, drop = TRUE) +
    ggplot2::scale_x_continuous(breaks = function(limits) if (limits[[2]] >= 100) NULL else pretty(limits, n = 5),
      expand = ggplot2::expansion(mult = c(.18, .18))) +
    ggplot2::labs(x = "Mean NES (right panel)", y = NULL, size = "n gene sets",
      caption = paste0(lisa_support_legend(), "\n",
        "Support: at least d/N evaluable sets; blank = no positive bound or not evaluable (see complete table). Sets may overlap.\n",
        "Points: historical mean NES, size and dictionary colours among individually significant members only; missing means remain blank.",
        if (direction) "\nRight-hand percentages: same-direction members among the individually significant sets, NOT the support bound." else "")) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "right", panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(), strip.text = ggplot2::element_text(face = "bold", size = 9),
      plot.caption = ggplot2::element_text(hjust = 0, size = 9),
      plot.background = ggplot2::element_rect(fill = "white", color = NA))
  if (direction && nrow(geometry)) {
    if (!"same_direction_pct" %in% names(geometry)) geometry$same_direction_pct <- round(100 * geometry$consistency, 1)
    geometry$direction_stats_label <- paste0(round(geometry$same_direction_pct), "% (n=", geometry$n_genesets, ")")
    geometry$label_hjust <- ifelse(geometry$mean_NES >= 0, -.15, 1.15)
    p <- p + ggplot2::geom_text(data = geometry,
      ggplot2::aes(x = mean_NES, y = category_plot_label, label = direction_stats_label, hjust = label_hjust),
      size = 2.7, fontface = "bold", color = "grey20")
  }
  if (isTRUE(group_by_supracategory) && "macrogroup_name" %in% names(x)) {
    p <- p + ggplot2::facet_grid(macrogroup_name ~ panel, scales = "free", space = "free_y", switch = "y") +
      ggplot2::theme(strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 0, hjust = 1, size = 8),
        panel.spacing.y = grid::unit(.25, "lines"))
  } else p <- p + ggplot2::facet_grid(. ~ panel, scales = "free_x")
  p
}
