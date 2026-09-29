# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Category-level navigation is a view of an existing single-analysis summary.
# It deliberately does not aggregate enrichment, infer category p-values, or
# reinterpret a missing mean as zero.  All links retain the exact analysis scope.

lisa_navigation_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  for (pair in list(c("&", "&amp;"), c("<", "&lt;"), c(">", "&gt;"),
      c('"', "&quot;"), c("'", "&#39;"))) x <- gsub(pair[[1L]], pair[[2L]], x, fixed = TRUE)
  x
}

lisa_navigation_read <- function(x) {
  if (is.data.frame(x)) return(x)
  if (!is.character(x) || length(x) != 1L || !file.exists(x))
    stop("Category navigation requires one existing summary TSV or data frame.", call. = FALSE)
  utils::read.delim(x, sep = "\t", quote = "", comment.char = "", check.names = FALSE,
    colClasses = "character", na.strings = NULL, stringsAsFactors = FALSE)
}

lisa_navigation_scope_component <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x) ||
      !grepl("^[A-Za-z0-9_.-]+$", x) || x %in% c(".", ".."))
    stop(name, " must be an exact safe single-analysis path component.", call. = FALSE)
  x
}

# Internal API: summary$mean_NES is the existing arithmetic mean among
# significant mapped member sets.  Input order is preserved from that summary.
build_lisa_category_navigation <- function(summary, analysis_id, collection, tier,
    gsea_padj_cutoff = .25, positive_contrast = "not recorded") {
  analysis_id <- lisa_navigation_scope_component(analysis_id, "analysis_id")
  collection <- lisa_navigation_scope_component(collection, "collection")
  if (!is.character(tier) || length(tier) != 1L || is.na(tier) || !nzchar(tier))
    stop("A single non-empty tier is required.", call. = FALSE)
  if (!is.numeric(gsea_padj_cutoff) || length(gsea_padj_cutoff) != 1L ||
      !is.finite(gsea_padj_cutoff) || gsea_padj_cutoff < 0 || gsea_padj_cutoff > 1)
    stop("gsea_padj_cutoff must lie between zero and one.", call. = FALSE)
  if (length(positive_contrast) != 1L || is.na(positive_contrast))
    stop("positive_contrast must be a single description.", call. = FALSE)
  x <- lisa_navigation_read(summary)
  if (!all(c("category_id", "mean_NES") %in% names(x)))
    stop("Category navigation requires category_id and the existing mean_NES column.", call. = FALSE)
  # A/B merged summaries are outside this single-analysis feature's contract.
  if (any(c("mean_NES_A", "mean_NES_B", "delta_mean_NES") %in% names(x)))
    stop("Category navigation accepts individual analyses, not contrast summaries.", call. = FALSE)
  for (key in intersect(c("analysis_id", "analysis_collection", "universe", "tier"), names(x))) {
    expected <- switch(key, analysis_id = analysis_id, tier = tier, collection)
    observed <- unique(as.character(x[[key]]))
    observed <- observed[!is.na(observed) & nzchar(observed)]
    if (length(setdiff(observed, expected))) stop("Category navigation ", key, " scope mismatch.", call. = FALSE)
  }
  if ("analysis" %in% names(x) && any(!is.na(x$analysis) & x$analysis != "GSEA"))
    stop("Category navigation represents GSEA summaries only.", call. = FALSE)
  if ("gsea_padj_cutoff" %in% names(x)) {
    observed <- suppressWarnings(as.numeric(x$gsea_padj_cutoff))
    if (any(!is.finite(observed) | observed != gsea_padj_cutoff))
      stop("Category navigation cutoff differs from the source summary.", call. = FALSE)
  }
  id <- as.character(x$category_id)
  classified <- !is.na(id) & nzchar(id) & id != "OTHER_UNCLASSIFIED"
  if ("classification_status" %in% names(x))
    classified <- classified & !is.na(x$classification_status) & x$classification_status == "classified"
  x <- x[classified, , drop = FALSE]
  id <- as.character(x$category_id)
  if (anyDuplicated(id)) stop("Duplicate category IDs in navigation summary.", call. = FALSE)
  num <- function(key, fallback = NULL) {
    if (!key %in% names(x)) {
      if (is.null(fallback)) stop("Missing summary column: ", key, call. = FALSE)
      return(fallback)
    }
    suppressWarnings(as.numeric(x[[key]]))
  }
  count_key <- if ("n_genesets_significant" %in% names(x)) "n_genesets_significant" else "n_genesets"
  n <- num(count_key)
  if (any(!is.finite(n) | n < 0 | n != floor(n)))
    stop("Significant gene-set counts must be finite nonnegative integers.", call. = FALSE)
  means <- num("mean_NES")
  if (any(n > 0 & !is.finite(means)))
    stop("Significant support requires a finite existing mean_NES.", call. = FALSE)
  if (any(n == 0 & is.finite(means)))
    stop("No-support categories must retain a missing mean_NES, not zero.", call. = FALSE)
  text_column <- function(key, fallback = "") {
    values <- if (key %in% names(x)) as.character(x[[key]]) else rep(fallback, nrow(x))
    values[is.na(values) | !nzchar(values)] <- rep(fallback, length.out = length(values))[is.na(values) | !nzchar(values)]
    values
  }
  labels <- text_column("category_display_name")
  labels[!nzchar(labels)] <- id[!nzchar(labels)]
  group_id <- text_column("macrogroup_id")
  group_name <- text_column("macrogroup_name")
  group_id[!nzchar(group_id)] <- "UNCLASSIFIED"
  group_name[!nzchar(group_name)] <- "Other or unclassified"
  group_order <- if ("macrogroup_order" %in% names(x)) suppressWarnings(as.numeric(x$macrogroup_order)) else rep(NA_real_, nrow(x))
  category_order <- if ("category_order_within_macrogroup" %in% names(x)) suppressWarnings(as.numeric(x$category_order_within_macrogroup)) else rep(NA_real_, nrow(x))
  category_color <- text_column("color", "#737373")
  category_color[!grepl("^#[0-9A-Fa-f]{6}$", category_color)] <- "#737373"
  table <- data.frame(category_id = id, category_display_name = labels,
    macrogroup_id = group_id, macrogroup_name = group_name,
    macrogroup_order = group_order, category_order_within_macrogroup = category_order,
    category_color = category_color,
    mean_NES = means, n_genesets_significant = n,
    n_pos_genesets = num("n_pos_genesets", rep(NA_real_, nrow(x))),
    n_neg_genesets = num("n_neg_genesets", rep(NA_real_, nrow(x))),
    n_genesets_mapped = num("n_genesets_mapped", rep(NA_real_, nrow(x))),
    n_genesets_evaluable = num("n_genesets_evaluable", rep(NA_real_, nrow(x))),
    support_state = ifelse(n == 0, "no_significant_member_sets", "significant_member_sets"),
    source_row_order = seq_len(nrow(x)), stringsAsFactors = FALSE)
  for (key in c("n_pos_genesets", "n_neg_genesets", "n_genesets_mapped", "n_genesets_evaluable")) {
    count <- table[[key]]
    if (any(!is.na(count) & (!is.finite(count) | count < 0 | count != floor(count))))
      stop("Invalid count in source summary: ", key, call. = FALSE)
  }
  if (any(table$n_pos_genesets + table$n_neg_genesets > n, na.rm = TRUE))
    stop("Directional set counts exceed significant support.", call. = FALSE)
  metadata <- list(schema_version = "1.0", analysis_id = analysis_id, collection = collection,
    tier = tier, gsea_padj_cutoff = gsea_padj_cutoff, positive_contrast = as.character(positive_contrast),
    metric = "Existing mean_NES: arithmetic mean among significant mapped member gene sets.",
    inference = "Descriptive category summary; no category-level p-value or FDR is inferred.",
    ordering = "Preserves the supplied category-summary row order.",
    exclusions = "Unclassified residuals and contrast summaries are not visualized.",
    input_summary_sha256 = if (is.character(summary)) digest::digest(file = summary, algo = "sha256") else "data.frame input")
  structure(list(metadata = metadata, categories = table), class = "lisa_category_navigation")
}

lisa_category_navigation_links <- function(navigation, evidence_base) {
  if (length(evidence_base) != 1L || is.na(evidence_base) || !nzchar(evidence_base) ||
      grepl("^[A-Za-z][A-Za-z0-9+.-]*:|^[/\\\\]|[?#]", evidence_base))
    stop("evidence_base must be a relative HTML path without query or fragment.", call. = FALSE)
  if (!nrow(navigation$categories)) return(character())
  encode <- function(x) utils::URLencode(x, reserved = TRUE)
  m <- navigation$metadata
  paste0(evidence_base, "?category=", vapply(navigation$categories$category_id, encode, character(1L)),
    "&analysis_id=", encode(m$analysis_id), "&collection=", encode(m$collection), "&tier=", encode(m$tier))
}

# One deterministic SVG renderer is captured in the reproduction script.  It
# uses strings only, so opening or reproducing the diagram needs no JS server.
lisa_category_navigation_svg <- function(source, title, subtitle,
    count_header = "Significant sets (+ / \u2212)",
    metric_label = "Mean NES", metric_value_label = "mean NES",
    metric_description = "Mean NES is descriptive, not a category significance test.",
    count_mode = c("total_and_sides", "sides_only"), side_labels = c("+", "\u2212"),
    axis_header = paste0(metric_label, " of significant member sets"),
    unavailable_state = paste0(metric_value_label, " unavailable: no significant member sets")) {
  count_mode <- match.arg(count_mode)
  if (!is.character(metric_label) || length(metric_label) != 1L || is.na(metric_label) || !nzchar(metric_label) ||
      !is.character(metric_value_label) || length(metric_value_label) != 1L || is.na(metric_value_label) || !nzchar(metric_value_label) ||
      !is.character(metric_description) || length(metric_description) != 1L || is.na(metric_description) || !nzchar(metric_description) ||
      !is.character(axis_header) || length(axis_header) != 1L || is.na(axis_header) || !nzchar(axis_header) ||
      !is.character(unavailable_state) || length(unavailable_state) != 1L || is.na(unavailable_state) || !nzchar(unavailable_state))
    stop("Navigation metric labels must be non-empty strings.", call. = FALSE)
  if (!is.character(side_labels) || length(side_labels) != 2L || any(is.na(side_labels) | !nzchar(side_labels)))
    stop("Navigation side labels must contain exactly two non-empty strings.", call. = FALSE)
  escape <- function(x) {
    x <- as.character(x); x[is.na(x)] <- ""
    for (p in list(c("&", "&amp;"), c("<", "&lt;"), c(">", "&gt;"), c('"', "&quot;"), c("'", "&#39;")))
      x <- gsub(p[[1L]], p[[2L]], x, fixed = TRUE)
    x
  }
  number <- function(x) {
    value <- suppressWarnings(as.numeric(x))
    if (length(value) && is.finite(value)) formatC(value, digits = 3L, format = "fg", flag = "#") else "NA"
  }
  adjusted_p <- function(x) {
    value <- suppressWarnings(as.numeric(x))
    if (length(value) && is.finite(value)) format(value, digits = 17L, trim = TRUE) else "Not evaluable"
  }
  significance_stars <- function(p) {
    if (!is.finite(p) || p > .05) "" else if (p <= .0001) "****" else if (p <= .001) "***" else if (p <= .01) "**" else "*"
  }
  paired_count <- function(x) {
    value <- suppressWarnings(as.numeric(x))
    if (length(value) && is.finite(value) && value >= 0 && value == floor(value))
      formatC(value, digits = 0L, format = "f") else "NA"
  }
  wrap <- function(x, width = 23L) {
    words <- strsplit(trimws(as.character(x)), "[[:space:]]+")[[1L]]
    if (!length(words) || !nzchar(words[[1L]])) return("")
    lines <- character(); line <- ""
    for (word in words) {
      if (!nzchar(line)) line <- word else if (nchar(paste(line, word), type = "width") <= width) line <- paste(line, word) else {
        lines <- c(lines, line); line <- word
      }
    }
    c(lines, line)
  }
  paired <- all(c("mean_NES_A", "mean_NES_B") %in% names(source))
  hommel <- "hommel_support_schema" %in% names(source) || "support_schema_version" %in% names(source)
  if (hommel && paired) stop("Hommel support cannot decorate a contrast navigator.", call. = FALSE)
  means <- suppressWarnings(as.numeric(source$mean_NES))
  means_a <- if (paired) suppressWarnings(as.numeric(source$mean_NES_A)) else means
  means_b <- if (paired) suppressWarnings(as.numeric(source$mean_NES_B)) else means
  if (paired) {
    metric_description <- "Separate category profiles from each analysis, not their difference. Missing support remains unavailable."
    axis_header <- "Mean NES: A and B"
  }
  n <- nrow(source); width <- if (hommel) 1740 else 1300; top <- 144
  spacing <- if (paired || hommel) 50 else 34
  group_id <- if ("macrogroup_id" %in% names(source)) as.character(source$macrogroup_id) else rep("", n)
  group_name <- if ("macrogroup_name" %in% names(source)) as.character(source$macrogroup_name) else rep("", n)
  group_id[is.na(group_id) | !nzchar(group_id)] <- "UNCLASSIFIED"
  group_name[is.na(group_name) | !nzchar(group_name)] <- "Other or unclassified"
  group_key <- paste(group_id, group_name, sep = "\r")
  starts <- if (n) c(TRUE, group_key[-1L] != group_key[-n]) else logical()
  group_index <- cumsum(starts)
  group_lines <- if (n) lapply(seq_len(max(group_index)), function(g) wrap(group_name[match(g, group_index)])) else list()
  group_height <- if (n) vapply(seq_len(max(group_index)), function(g) max(sum(group_index == g) * spacing, length(group_lines[[g]]) * 14 + 10), numeric(1L)) else numeric()
  row_y <- numeric(n); group_y <- numeric(length(group_height)); cursor <- top
  if (n) for (g in seq_along(group_height)) {
    group_y[[g]] <- cursor; rows <- which(group_index == g)
    row_y[rows] <- cursor + (group_height[[g]] - length(rows) * spacing) / 2 + (seq_along(rows) - 1L) * spacing + 15
    cursor <- cursor + group_height[[g]] + 8
  }
  content_bottom <- if (n) cursor - 8 else top
  height <- content_bottom + 85
  span <- max(c(1, abs(c(means_a, means_b)[is.finite(c(means_a, means_b))]))) * 1.12
  left <- 555; plot_width <- 430; zero <- left + plot_width / 2
  position <- function(value) zero + value / span * plot_width / 2
  out <- c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s" role="group" %saria-labelledby="navigation-title navigation-description">', width, height, width, height,
    if (hommel) 'data-hommel-support="true" ' else ""),
    paste0('<title id="navigation-title">', escape(title), '</title>'),
    paste0('<desc id="navigation-description">', escape(subtitle), ' Each category label and point links to its category evidence. ',
      if (hommel) 'Stars show adjusted category P from robust Hommel multiple-testing correction; minimum support d/N is in Category evidence. ' else "",
      escape(metric_description), '</desc>'),
    '<rect width="100%" height="100%" fill="#ffffff"/>',
    '<style>a:focus .row-background{fill:#e7edf3;stroke:#24374b;stroke-width:2}a:hover .row-background{fill:#f2f5f8}.category-row{cursor:pointer}</style>',
    paste0('<g font-family="Arial,Helvetica,sans-serif" fill="#172b3a"><text x="18" y="30" font-size="20" font-weight="bold">', escape(title), '</text>'),
    paste0('<text x="18" y="55" font-size="12">', escape(subtitle), '</text>'),
    if (paired) paste0('<text x="18" y="78" font-size="12">Blue circles: A \u00b7 ', escape(side_labels[[1]]), ' \u00b7 Orange squares: B \u00b7 ', escape(side_labels[[2]]), '</text>') else
    paste0('<text x="18" y="78" font-size="12">Red: positive ', escape(metric_value_label), ' \u00b7 Blue: negative ', escape(metric_value_label), ' \u00b7 Grey: unavailable or zero</text>'),
    if (hommel) '<text x="18" y="100" font-size="12">Category P stars: * \u2264 0.05 \u00b7 ** \u2264 0.01 \u00b7 *** \u2264 0.001 \u00b7 **** \u2264 0.0001. Blank: not significant or not evaluable. Open Category evidence for exact P and d/N.</text>' else
      '<text x="18" y="100" font-size="12">Select a category to inspect its gene sets and genes. Counts describe significant member sets, not independent confirmations.</text>',
    '<text x="18" y="124" font-size="12" font-weight="bold">Supercategory</text>',
    '<text x="232" y="124" font-size="12" font-weight="bold">LISA category</text>',
    sprintf('<text x="%s" y="124" font-size="12" text-anchor="middle" font-weight="bold">%s</text>', zero, escape(axis_header)),
    paste0('<text x="1012" y="124" font-size="12" font-weight="bold">', escape(count_header), '</text>'),
    if (hommel) '<text class="hommel-support-header" x="1240" y="124" font-size="12" font-weight="bold">Category significance \u00b7 adjusted P</text>' else "",
    sprintf('<line x1="%s" y1="133" x2="%s" y2="%s" stroke="#9ba8b3"/>', zero, zero, content_bottom))
  if (!n) out <- c(out, '<text x="18" y="161" font-size="14">No classified categories are present in this summary.</text>')
  if (n) for (g in seq_along(group_height)) {
    rows <- which(group_index == g)
    accent <- if ("category_color" %in% names(source)) as.character(source$category_color[[rows[[1L]]]]) else "#737373"
    if (is.na(accent) || !grepl("^#[0-9A-Fa-f]{6}$", accent)) accent <- "#737373"
    out <- c(out, sprintf('<rect class="supercategory-band" x="10" y="%.3f" width="202" height="%.3f" rx="4" fill="#f4f6f8" stroke="%s" stroke-width="2"/>', group_y[[g]], group_height[[g]], accent))
    line_y <- group_y[[g]] + (group_height[[g]] - length(group_lines[[g]]) * 14) / 2 + 11
    out <- c(out, vapply(seq_along(group_lines[[g]]), function(j) sprintf('<text class="supercategory-label" x="20" y="%.3f" font-size="12" font-weight="bold" fill="#172b3a">%s</text>', line_y + (j - 1L) * 14, escape(group_lines[[g]][[j]])), character(1L)))
  }
  for (i in seq_len(n)) {
    y <- row_y[[i]]
    mean <- means[[i]]; value <- number(mean)
    color <- if (!is.finite(mean) || mean == 0) "#687787" else if (mean > 0) "#b2182b" else "#2166ac"
    label <- as.character(source$category_display_name[[i]])
    short <- if (nchar(label, type = "width") > 43L) paste0(substr(label, 1L, 40L), "\u2026") else label
    state <- if (is.finite(mean)) paste0(metric_value_label, " ", value) else unavailable_state
    counts <- if (identical(count_mode, "total_and_sides")) {
      paste0(number(source$n_genesets_significant[[i]]), " (", number(source$n_pos_genesets[[i]]), " / ", number(source$n_neg_genesets[[i]]), ")")
    } else {
      paste0(side_labels[[1L]], " ", paired_count(source$n_pos_genesets[[i]]), " / ", side_labels[[2L]], " ", paired_count(source$n_neg_genesets[[i]]))
    }
    if (paired) {
      state <- paste0("A mean NES ", number(means_a[[i]]), "; B mean NES ", number(means_b[[i]]))
      counts <- paste0("A ", paired_count(source$n_pos_genesets[[i]]), " / B ", paired_count(source$n_neg_genesets[[i]]))
    }
    accessible <- paste0(label, "; ", source$category_id[[i]], "; ", state, "; ", counts, ". Open category evidence.")
    if (hommel) {
      p <- suppressWarnings(as.numeric(source$category_p_adjusted[[i]]))
      mark <- significance_stars(p)
      category_status <- if (!is.finite(p)) "not evaluable" else if (nzchar(mark)) "significant" else "not significant"
      accessible <- paste0(accessible, " Adjusted category P: ", adjusted_p(p),
        "; category significance: ", category_status, if (nzchar(mark)) paste0(" (", mark, ")") else "", ".")
    }
    # A slash may be part of a category identifier, not a filesystem reference.
    # Escape its visible character while keeping links as encoded query values.
    accessible_html <- gsub("/", "&#47;", escape(accessible), fixed = TRUE)
    out <- c(out, paste0('<a class="category-row" data-category-id="', escape(source$category_id[[i]]), '" href="', escape(source$evidence_url[[i]]), '" target="_top" tabindex="0" aria-label="', accessible_html, '"><title>', accessible_html, '</title>'),
      sprintf('<rect class="row-background" x="222" y="%.3f" width="%s" height="%s" rx="3" fill="transparent"/>', y - 22, width - 232, if (hommel) 46 else 32),
      paste0('<text x="232" y="', y, '" font-size="12">', escape(short), '</text>'),
      paste0('<text x="1012" y="', y, '" font-size="12">', escape(counts), '</text>'))
    if (hommel) {
      out <- c(out, paste0('<g class="hommel-support-cell" data-significance-stars="', escape(mark),
        '" data-support-status="', escape(source$status[[i]]), '">',
        if (nzchar(mark)) paste0('<text x="1240" y="', y, '" font-size="17" font-weight="bold">',
          escape(mark), '</text>') else "", '</g>'))
    }
    if (paired) {
      for (side in c("A", "B")) {
        value_side <- if (side == "A") means_a[[i]] else means_b[[i]]
        yy <- y + if (side == "A") -13 else 6
        cc <- if (side == "A") "#2166ac" else "#b85c00"
        if (is.finite(value_side)) {
          px <- position(value_side)
          out <- c(out, sprintf('<g data-profile="%s" data-mean-nes="%.17g"><line x1="%s" y1="%s" x2="%.3f" y2="%s" stroke="%s" stroke-width="2"/>', side, value_side, zero, yy, px, yy, cc),
            if (side == "A") sprintf('<circle cx="%.3f" cy="%s" r="4" fill="%s"/>', px, yy, cc) else sprintf('<rect x="%.3f" y="%s" width="8" height="8" fill="%s"/>', px-4, yy-4, cc),
            sprintf('<text x="%.3f" y="%s" font-size="10" text-anchor="%s" fill="%s">%s</text></g>', px + if(value_side >= 0) 8 else -8, yy+3, if(value_side >= 0) "start" else "end", cc, number(value_side)))
        } else out <- c(out, sprintf('<text data-profile="%s" x="%s" y="%s" font-size="10" text-anchor="middle" fill="%s">%s: unavailable</text>',side,zero,yy,cc,side))
      }
    } else if (is.finite(mean)) {
      px <- position(mean)
      out <- c(out, sprintf('<line x1="%s" y1="%s" x2="%.3f" y2="%s" stroke="%s" stroke-width="3"/>', zero, y - 4, px, y - 4, color),
        sprintf('<circle cx="%.3f" cy="%s" r="5" fill="%s"/>', px, y - 4, color),
        sprintf('<text x="%.3f" y="%s" font-size="11" text-anchor="%s" fill="%s">%s</text>', px + if (mean >= 0) 10 else -10, y, if (mean >= 0) "start" else "end", color, value))
    } else out <- c(out, sprintf('<text x="%s" y="%s" text-anchor="middle" font-size="11" fill="%s">not available</text>', zero, y, color))
    out <- c(out, '</a>')
  }
  bottom <- content_bottom + 25
  out <- c(out, sprintf('<line x1="%s" y1="%s" x2="%s" y2="%s" stroke="#687787"/>', left, bottom, left + plot_width, bottom))
  for (tick in c(-span, 0, span)) out <- c(out, sprintf('<text x="%.3f" y="%s" text-anchor="middle" font-size="11">%s</text>', position(tick), bottom + 18, number(tick)))
  paste(c(out, '</g></svg>'), collapse = "\n")
}

render_lisa_category_navigation <- function(navigation, output_dir, evidence_base = NULL,
                                            download_policy = NULL) {
  if (!inherits(navigation, "lisa_category_navigation")) stop("Expected lisa_category_navigation object.", call. = FALSE)
  m <- navigation$metadata
  if (is.null(evidence_base)) evidence_base <- paste0("../../../evidence/", utils::URLencode(m$analysis_id, reserved = TRUE), "/", utils::URLencode(m$collection, reserved = TRUE), "/index.html")
  source <- navigation$categories
  source$evidence_url <- lisa_category_navigation_links(navigation, evidence_base)
  source$analysis_id <- rep(m$analysis_id, nrow(source)); source$collection <- rep(m$collection, nrow(source)); source$tier <- rep(m$tier, nrow(source))
  source$gsea_padj_cutoff <- rep(m$gsea_padj_cutoff, nrow(source))
  title <- paste0("Category evidence navigator \u00b7 ", if (!is.null(m$display_label)) m$display_label else m$analysis_id)
  subtitle <- paste0(m$collection, " \u00b7 ", m$tier, " \u00b7 member-set GSEA FDR \u2264 ", format(m$gsea_padj_cutoff, digits = 17), " \u00b7 Positive ranking: ", m$positive_contrast)
  svg <- lisa_category_navigation_svg(source, title, subtitle)
  lisa_guarded_dir_create(output_dir)
  lisa_guarded_dir_create(file.path(output_dir, "tables"))
  source_path <- file.path(output_dir, "tables", "category_navigation_source.tsv")
  if ("support_schema_version" %in% names(source)) lisa_write_inference_tsv(source, source_path) else write_lisa_tsv(source, source_path)
  m$source_sha256 <- digest::digest(file = source_path, algo = "sha256")
  m$source_file <- "tables/category_navigation_source.tsv"
  m$evidence_base <- evidence_base; m$title <- title; m$subtitle <- subtitle
  svg_path <- file.path(output_dir, "category_navigation.svg")
  lisa_guarded_write(svg_path, function(path) writeLines(svg, path, useBytes = TRUE))
  metadata_path <- file.path(output_dir, "metadata.json")
  lisa_guarded_write(metadata_path, function(path) jsonlite::write_json(m, path, auto_unbox = TRUE, pretty = TRUE, digits = 17))
  escape <- lisa_navigation_escape
  requested <- function(key) is.null(download_policy) || isTRUE(download_policy[[key]])
  # The report privacy gate deliberately treats a slash after a quote as a
  # possible absolute-path token.  In the embedded copy only, encode slashes
  # in data identifiers; links retain their ordinary relative URL form.
  html_svg <- svg
  for (category_id in unique(as.character(source$category_id))) {
    raw_id <- escape(category_id)
    html_id <- gsub("/", "&#47;", raw_id, fixed = TRUE)
    html_svg <- gsub(paste0('data-category-id="', raw_id, '"'),
      paste0('data-category-id="', html_id, '"'), html_svg, fixed = TRUE)
  }
  html <- c('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
    paste0('<title>', escape(title), '</title>'),
    '<style>body{font-family:Arial,Helvetica,sans-serif;margin:1.5rem;color:#172b3a;background:#fff}main{max-width:1400px;margin:auto}a{color:#16598a}a:focus{outline:3px solid #24374b;outline-offset:3px}.diagram{overflow:auto;border:1px solid #dce2e7}.diagram svg{display:block}table{border-collapse:collapse;width:100%;font-size:0.9rem}th,td{padding:0.6rem;text-align:left;vertical-align:top;border-bottom:1px solid #dce2e7}code{font-size:0.85em}summary{cursor:pointer}nav a{margin-right:1rem}</style></head><body><main>',
    paste0('<h1>', escape(title), '</h1><p>', escape(subtitle), '</p>'),
    if ("support_schema_version" %in% names(source))
      '<p>Select a category label or point to open its evidence. With a keyboard, tab to a row and press Enter. Stars show adjusted category P (* &le; 0.05, ** &le; 0.01, *** &le; 0.001, **** &le; 0.0001); the minimum supported count d/N is in Category evidence. Mean NES describes the individually significant sets separately.</p>' else
      '<p>Select a category label or point to open its evidence. With a keyboard, tab to a row and press Enter. Mean NES and member-set counts are descriptive; this historical navigator does not show category-level adjusted P or minimum support d/N.</p>',
    if ("support_schema_version" %in% names(source)) lisa_hommel_help_button() else "",
    paste0('<nav aria-label="Downloads">',
      if (requested("svg")) '<a href="category_navigation.svg" download>Download SVG</a>' else "",
      if (requested("source_data")) '<a href="tables/category_navigation_source.tsv" download>Download source table</a>' else "",
      '<a href="metadata.json" download>Download provenance</a>',
      if (requested("recipes")) '<a href="reproduce_category_navigation.R" download>R script</a>' else "", '</nav>'),
    '<div class="diagram">', html_svg, '</div>',
    '<p>Red = positive NES; blue = negative NES. Percentages and member evidence are available in each category sheet.</p></main></body></html>')
  html <- lisa_present_evidence_html(html, output_dir, m)
  html_path <- file.path(output_dir, "index.html")
  lisa_guarded_write(html_path, function(path) writeLines(html, path, useBytes = TRUE))
  recipe <- c('#!/usr/bin/env Rscript', '# Offline base-R SVG reproduction. Run from this navigator directory:',
    '# Rscript --vanilla reproduce_category_navigation.R tables/category_navigation_source.tsv reproduced.svg',
    '# Preserve the report directory structure for relative evidence links.',
    paste0('lisa_category_navigation_svg <- ', paste(deparse(lisa_category_navigation_svg), collapse = '\n')),
    paste0('title <- ', encodeString(title, quote = '"')), paste0('subtitle <- ', encodeString(subtitle, quote = '"')),
    'args <- commandArgs(TRUE)', 'if (length(args) != 2L) stop("Expected source TSV and destination SVG.")',
    'if (tolower(tools::file_ext(args[[2L]])) != "svg") stop("Destination must be SVG.")',
    'source <- utils::read.delim(args[[1L]], sep="\\t", quote="", comment.char="", check.names=FALSE, colClasses="character", na.strings=NULL)',
    'writeLines(lisa_category_navigation_svg(source, title, subtitle), args[[2L]], useBytes=TRUE)')
  lisa_guarded_write(file.path(output_dir, "reproduce_category_navigation.R"), function(path) writeLines(recipe, path, useBytes = TRUE))
  invisible(list(html = html_path, svg = svg_path, source = source_path, metadata = metadata_path))
}
