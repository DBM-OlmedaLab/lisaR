# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé. GPL-3 or later.

# Category inference is separate from the descriptive, FDR-filtered NES layer.
# One closed-testing family contains every eligible gene-set null in one DE
# analysis, including unclassified sets in all applicable collections.
# HALLMARKS is a direct collection, explicitly outside this inference.
lisa_inference_exact_columns <- function(x) {
  # write.table serializes doubles at about 15 digits. Freeze 17 significant
  # digits explicitly so inclusive P <= alpha and figure recipes survive TSV.
  for (field in names(x)[vapply(x, is.numeric, logical(1))]) {
    value <- x[[field]]
    x[[field]] <- ifelse(is.na(value), "", sprintf("%.17g", value))
  }
  x
}

lisa_write_inference_tsv <- function(x, path) {
  lisa_write_figure_source_tsv(lisa_inference_exact_columns(x), path)
}

lisa_inference_closed_bounds <- function(fit, p, ix, alpha) {
  n <- length(ix)
  if (!n) return(list(p = NA_real_, discoveries = NA_integer_))
  zero_count <- sum(p[ix] == 0)
  q <- function(k) {
    # Exact zero shortcut avoids hommel 1.8 out-of-bounds localtest warning.
    if (k < zero_count) return(0)
    as.numeric(hommel::localtest(fit, ix, tdp = if (k == 0L) 0 else (k + .5) / n))
  }
  p_category <- q(0L)
  if (p_category > alpha) return(list(p = p_category, discoveries = 0L))
  low <- 1L; high <- n
  # Invert the same public closed partial-conjunction tests for the bound.
  # discoveries() can disagree at exact machine thresholds (C++ ceil versus
  # localtest multiplications); do not add tolerances or change alpha.
  while (low < high) {
    middle <- as.integer(floor((low + high + 1L) / 2))
    if (q(middle - 1L) <= alpha) low <- middle else high <- middle - 1L
  }
  list(p = p_category, discoveries = as.integer(low))
}

lisa_category_inference <- function(ledger, assignments, categories, alpha = .05) {
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1)
    stop("Category inference alpha must lie strictly between zero and one.", call. = FALSE)
  lisa_require_columns(ledger, c("collection", "pathway", "hypothesis_id", "eligible_for_gsea", "pval", "NES"), "category inference ledger")
  lisa_require_columns(assignments, c("collection", "pathway", "category_id"), "category assignments")
  lisa_require_columns(categories, c("collection", "category_id", "category_display_name"), "category catalogue")
  if (any(ledger$collection == "HALLMARKS") || any(assignments$collection == "HALLMARKS") || any(categories$collection == "HALLMARKS"))
    stop("HALLMARKS is not a LISA category collection and is excluded from category inference.", call. = FALSE)
  if (anyDuplicated(paste(ledger$collection, ledger$pathway, sep = "\r")))
    stop("Duplicate gene-set result within a collection.", call. = FALSE)
  if (anyNA(ledger$hypothesis_id) || any(!nzchar(ledger$hypothesis_id)) || anyNA(ledger$eligible_for_gsea))
    stop("Gene-set hypothesis identity and eligibility must be explicit.", call. = FALSE)
  ledger$pval <- lisa_nes_number(ledger$pval, "raw p-value")
  ledger$NES <- lisa_nes_number(ledger$NES, "NES")
  if (any(!is.na(ledger$pval) & (!is.finite(ledger$pval) | ledger$pval < 0 | ledger$pval > 1)))
    stop("Raw GSEA p-values must be missing or finite numbers in [0,1].", call. = FALSE)
  ledger$raw_p_available <- is.finite(ledger$pval)
  ledger$p_for_inference <- ifelse(ledger$raw_p_available, ledger$pval, 1)
  ledger$included_in_family <- ledger$eligible_for_gsea %in% TRUE
  x <- ledger[ledger$included_in_family, , drop = FALSE]
  ids <- sort(unique(x$hypothesis_id))
  family <- data.frame(hypothesis_id = ids, p_for_inference = numeric(length(ids)),
    n_collection_results = integer(length(ids)), n_missing_raw_p = integer(length(ids)),
    pathway = character(length(ids)), collections = character(length(ids)), stringsAsFactors = FALSE)
  groups <- split(seq_len(nrow(x)), x$hypothesis_id)
  for (i in seq_along(ids)) {
    z <- x[groups[[ids[[i]]]], , drop = FALSE]
    family$p_for_inference[[i]] <- max(z$p_for_inference)
    family$n_collection_results[[i]] <- nrow(z)
    family$n_missing_raw_p[[i]] <- sum(!z$raw_p_available)
    family$pathway[[i]] <- paste(sort(unique(z$pathway)), collapse = ";")
    family$collections[[i]] <- paste(sort(unique(z$collection)), collapse = ";")
  }
  # pmax across stochastic reruns of the SAME null is valid under arbitrary
  # dependence. Replacing unavailable eligible p-values by 1 avoids selecting
  # a family based on successful/finite results. No p-value is selected by NES.
  fit <- if (nrow(family)) hommel::hommel(family$p_for_inference, simes = FALSE) else NULL
  assignments <- unique(assignments[c("collection", "pathway", "category_id")])
  if (anyDuplicated(paste(categories$collection, categories$category_id, sep = "\r")))
    stop("Category catalogue identifiers must be unique within a collection.", call. = FALSE)
  rows <- lapply(seq_len(nrow(categories)), function(i) {
    cat <- categories[i, , drop = FALSE]
    a <- assignments[assignments$collection == cat$collection & assignments$category_id == cat$category_id, , drop = FALSE]
    z <- ledger[ledger$collection == cat$collection & ledger$pathway %in% a$pathway, , drop = FALSE]
    evaluable <- z$included_in_family & z$raw_p_available
    ix <- unique(match(z$hypothesis_id[evaluable], family$hypothesis_id))
    n <- length(ix)
    bounds <- lisa_inference_closed_bounds(fit, family$p_for_inference, ix, alpha)
    p <- bounds$p
    d <- bounds$discoveries
    if (n && !identical(p <= alpha, d >= 1L))
      stop("Category adjusted p-value and simultaneous discovery bound disagree.", call. = FALSE)
    cat$n_sets_total <- nrow(a)
    cat$n_sets_eligible <- sum(z$included_in_family)
    cat$n_sets_evaluable <- n
    cat$n_sets_missing_p <- sum(z$included_in_family & !z$raw_p_available)
    cat$n_sets_not_eligible_or_absent <- nrow(a) - sum(z$included_in_family)
    cat$coverage_pct <- if (nrow(a)) 100 * n / nrow(a) else NA_real_
    cat$category_p_adjusted <- p
    cat$minimum_enriched_sets <- d
    cat$minimum_enriched_pct <- if (n) 100 * d / n else NA_real_
    cat$significant <- n > 0L && p <= alpha
    cat$status <- if (!n) "not_evaluable" else if (p <= alpha) "significant" else "not_significant"
    cat$n_positive_NES <- sum(evaluable & is.finite(z$NES) & z$NES > 0)
    cat$n_negative_NES <- sum(evaluable & is.finite(z$NES) & z$NES < 0)
    cat$n_zero_NES <- sum(evaluable & is.finite(z$NES) & z$NES == 0)
    cat$n_missing_NES <- sum(evaluable & !is.finite(z$NES))
    cat$alpha <- alpha
    cat$confidence_level <- 1 - alpha
    cat$family_size <- nrow(family)
    cat$family_scope <- "All eligible unique gene-set nulls across configured LISA category collections within one DE analysis; HALLMARKS excluded"
    cat$method <- "hommel::hommel(simes=FALSE); invert localtest partial-conjunction tests"
    cat
  })
  table <- if (length(rows)) do.call(rbind, rows) else categories
  members <- merge(assignments, ledger, by = c("collection", "pathway"), all.x = TRUE, sort = FALSE)
  list(categories = table, family = family, raw_results = ledger, members = members)
}

# Readable display only; machine-readable canonical columns remain unchanged.
lisa_category_inference_display <- function(x) {
  x <- lisa_category_support_presentation(x)
  fmt <- function(z) ifelse(is.na(z), "Not available", ifelse(abs(z - .05) <= 1e-6,
    sprintf("%.17g", z), formatC(z, digits = 5, format = "g")))
  data.frame(Category = x$category_display_name,
    `Analysed / total gene sets` = paste(x$n_sets_evaluable, x$n_sets_total, sep = " / "),
    `Adjusted category P` = fmt(x$category_p_adjusted),
    Result = c(significant = "Significant", not_significant = "Not significant", not_evaluable = "Not evaluable")[x$status],
    `Minimum support level (95%)` = x$support_grade,
    `Minimum support: interpretation` = x$support_label,
    `Minimum enriched / evaluable sets (d/N)` = ifelse(is.na(x$support_count), "Not evaluable", x$support_count),
    `Minimum supported percentage (95%)` = ifelse(is.na(x$minimum_support_pct_text), "Not evaluable", paste0("\u2265", x$minimum_support_pct_text)),
    `Eligible sets missing raw P` = x$n_sets_missing_p,
    `Not eligible or absent sets` = x$n_sets_not_eligible_or_absent,
    `Observed NES pattern (descriptive)` = x$observed_NES_pattern,
    `Observed NES signs among evaluable sets` = x$observed_NES_counts,
    check.names = FALSE, stringsAsFactors = FALSE)
}

lisa_category_inference_source <- function(categories, members, title = "LISA category enrichment") {
  categories <- lisa_category_support_presentation(categories)
  if (!nrow(categories)) return(data.frame())
  categories$category_order <- seq_len(nrow(categories))
  if (!"analysis_id" %in% names(members)) {
    if (length(unique(categories$analysis_id)) != 1L) stop("Members need explicit analysis_id for a multi-analysis source.", call. = FALSE)
    members$analysis_id <- categories$analysis_id[[1L]]
  }
  member_key <- lisa_support_key(members)
  members <- members[member_key %in% lisa_support_key(categories) & members$included_in_family %in% TRUE & members$raw_p_available %in% TRUE, , drop = FALSE]
  source <- merge(categories, members[c("analysis_id", "collection", "category_id", "pathway", "NES", "pval", "hypothesis_id")],
    by = c("analysis_id", "collection", "category_id"), all.x = TRUE, sort = FALSE)
  source <- source[order(source$category_order, source$pathway), , drop = FALSE]
  source$figure_type <- "lisa_category_inference"
  source$figure_title <- title
  source$figure_width <- 17
  source$figure_height <- max(7.5, 4 + .30 * nrow(categories))
  source$figure_dpi <- 180
  source$selected_for_plot <- TRUE
  source
}

lisa_plot_category_inference_v1 <- function(source) {
  if (!nrow(source)) stop("No evaluable categories to plot.", call. = FALSE)
  s <- source[!duplicated(source$category_id), , drop = FALSE]
  s <- s[order(s$category_order), , drop = FALSE]
  label <- paste0(s$category_display_name, ifelse(s$significant, " *", ""))
  # Explicit positions avoid using possibly duplicated display names as IDs.
  s$y <- nrow(s) - seq_len(nrow(s)) + 1L
  s$panel <- "Minimum percentage supported"
  s$x <- s$minimum_enriched_pct
  s$count_label <- paste0(s$minimum_enriched_sets, " / ", s$n_sets_evaluable)
  m <- source[is.finite(source$NES), , drop = FALSE]
  m$y <- s$y[match(m$category_id, s$category_id)]
  # Deterministic small offsets: neither random jitter nor another test.
  offsets <- ave(seq_len(nrow(m)), m$category_id, FUN = function(i) ((seq_along(i) * .61803398875) %% 1 - .5) * .48)
  m$y <- m$y + offsets
  m$x <- m$NES
  m$panel <- rep("All analysed member-set NES (descriptive)", nrow(m))
  m$direction <- ifelse(m$NES > 0, "Positive NES", ifelse(m$NES < 0, "Negative NES", "Zero NES"))
  panels <- c("Minimum percentage supported", "All analysed member-set NES (descriptive)")
  s$panel <- factor(s$panel, levels = panels); m$panel <- factor(m$panel, levels = panels)
  limits <- data.frame(panel = factor(rep(panels, each = 2), levels = panels),
    x = c(0, 112, min(c(-1, m$x)), max(c(1, m$x))), y = 1)
  ggplot2::ggplot() +
    ggplot2::geom_blank(data = limits, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_segment(data = s, ggplot2::aes(x = 0, xend = x, y = y, yend = y), color = "#596778", linewidth = 1.1) +
    ggplot2::geom_point(data = s, ggplot2::aes(x = x, y = y), color = "#3d4855", size = 2) +
    ggplot2::geom_text(data = s, ggplot2::aes(x = x + 2, y = y, label = count_label), hjust = 0, size = 2.8, color = "#3d4855") +
    ggplot2::geom_vline(data = data.frame(panel = factor(panels[[2]], levels = panels)), ggplot2::aes(xintercept = 0), color = "grey65", linewidth = .35) +
    ggplot2::geom_point(data = m, ggplot2::aes(x = x, y = y, color = direction), size = .9, alpha = .55) +
    ggplot2::facet_grid(. ~ panel, scales = "free_x") +
    ggplot2::scale_x_continuous(breaks = function(limits) {
      # The fixed support panel reserves 0..112 for count labels; include the
      # biologically readable 100% endpoint rather than automatic 30% steps.
      if (limits[[1]] > -10 && limits[[2]] >= 112) c(0, 25, 50, 75, 100) else pretty(limits, n = 5)
    }) +
    ggplot2::scale_y_continuous(breaks = s$y, labels = label, limits = c(.3, nrow(s) + .7), expand = c(0, 0)) +
    ggplot2::scale_color_manual(values = c(`Positive NES` = "#b2182b", `Negative NES` = "#2166ac", `Zero NES` = "#687787")) +
    ggplot2::labs(title = source$figure_title[[1]],
      subtitle = "* Significant category: adjusted P <= 0.05. Unmarked categories remain visible.",
      x = "Left: minimum percentage (%)     |     Right: normalised enrichment score (NES)", y = NULL, color = NULL,
      caption = paste0("Left: at least the indicated number of analysed sets is enriched (simultaneous 95% confidence).\n",
        "The remainder is uncertain, not proven unenriched. A star does not imply that all members change together or in one direction.\n",
        "Gene sets may overlap; counts are not independent biological mechanisms. Full coverage and exact results are in the table.")) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#e8ebef"),
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      axis.text.y = ggplot2::element_text(size = 9), plot.title = ggplot2::element_text(face = "bold"),
      plot.caption = ggplot2::element_text(hjust = 0, size = 8), plot.margin = ggplot2::margin(12, 14, 12, 12))
}

# Historical sources without a support schema retain the original star plot.
# Newly derived sources show corrected category-P stars; d/N remains in the
# evidence card and machine-readable table, not as a second figure scale.
lisa_plot_category_inference <- function(source) {
  if (!"support_schema_version" %in% names(source)) return(lisa_plot_category_inference_v1(source))
  if (!identical(unique(as.character(source$support_schema_version)), "support-grades-v1"))
    stop("Unsupported support presentation schema.", call. = FALSE)
  key <- lisa_support_key(source)
  s <- lisa_category_support_presentation(source[!duplicated(key), , drop = FALSE])
  s <- s[order(s$category_order), , drop = FALSE]
  s$y <- nrow(s) - seq_len(nrow(s)) + 1L
  panels <- c("Category significance", "All evaluable member-set NES (descriptive)")
  s$panel <- factor(panels[[1]], levels = panels)
  s$significance_mark <- lisa_category_significance_stars(s$category_p_adjusted)
  m <- source[is.finite(source$NES), , drop = FALSE]
  m$y <- s$y[match(lisa_support_key(m), lisa_support_key(s))]
  if (nrow(m)) m$y <- m$y + ave(seq_len(nrow(m)), lisa_support_key(m),
    FUN = function(i) ((seq_along(i) * .61803398875) %% 1 - .5) * .48)
  m$panel <- factor(rep(panels[[2]], nrow(m)), levels = panels)
  m$direction <- ifelse(m$NES > 0, "Positive NES", ifelse(m$NES < 0, "Negative NES", "Zero NES"))
  limits <- data.frame(panel = factor(rep(panels, each = 2), levels = panels),
    x = c(0, 30, min(c(-1, m$NES)), max(c(1, m$NES))), y = 1)
  p <- ggplot2::ggplot() +
    ggplot2::geom_blank(data = limits, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_text(data = s, ggplot2::aes(x = 4, y = y, label = significance_mark), hjust = 0, size = 5.0, color = "#3d4855") +
    ggplot2::geom_vline(data = data.frame(panel = factor(panels[[2]], levels = panels)), ggplot2::aes(xintercept = 0), color = "grey65", linewidth = .35) +
    ggplot2::geom_point(data = m, ggplot2::aes(x = NES, y = y, color = direction), size = .9, alpha = .55) +
    ggplot2::facet_grid(. ~ panel, scales = "free_x") +
    ggplot2::scale_x_continuous(breaks = function(limits) pretty(limits, n = 5)) +
    ggplot2::scale_y_continuous(breaks = s$y, labels = s$category_display_name, limits = c(.3, nrow(s) + .7), expand = c(0, 0)) +
    ggplot2::scale_color_manual(values = c(`Positive NES` = "#b2182b", `Negative NES` = "#2166ac", `Zero NES` = "#687787")) +
    ggplot2::labs(title = source$figure_title[[1]], x = "Left: adjusted category P stars     |     Right: normalised enrichment score (NES)", y = NULL, color = NULL,
      caption = paste0("Stars: * P \u2264 0.05; ** P \u2264 0.01; *** P \u2264 0.001; **** P \u2264 0.0001 (adjusted category P). Blank: not significant or not evaluable.\n",
        "Robust Hommel multiple-testing correction; a significant category has evidence for at least one enriched set.\n",
        "Minimum supported enrichment d/N, with 95% simultaneous confidence, is in Category evidence.\n",
        "Gene sets may overlap: neither independent mechanisms nor additive category counts. Check coverage in the complete table.")) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 9), axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(face = "bold"), plot.caption = ggplot2::element_text(hjust = 0, size = 9),
      plot.margin = ggplot2::margin(12, 14, 12, 12))
  # Suppress numerical ticks only for the layout-only significance panel.
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

lisa_category_inference_save_plot <- function(source, path) {
  plot <- lisa_plot_category_inference(source)
  ext <- tolower(tools::file_ext(path))
  width <- source$figure_width[[1]]; height <- source$figure_height[[1]]
  if (ext == "png") grDevices::png(path, width = width, height = height,
    units = "in", res = source$figure_dpi[[1]], type = "cairo", bg = "white") else
  if (ext == "svg") grDevices::svg(path, width = width, height = height, bg = "white") else
  if (ext == "pdf") grDevices::cairo_pdf(path, width = width, height = height, bg = "white") else
    stop("Category inference figures support png, svg and pdf.", call. = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
  invisible(path)
}
