# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Contrast-level navigation is the A/B counterpart of category_navigation.R.
# It reads the already-computed GSEA category-contrast summary (delta_mean_NES
# and the side-specific significant-set counts) and never derives a new
# aggregate, delta or category-level test. Categories link to the exact same
# contrast_id/analysis_a/analysis_b/collection/tier scope already rendered by
# build_contrast_evidence.R, so a click opens that exact category there.

# Internal API: delta_mean_NES is the existing "A minus B" descriptive value
# from the canonical category-contrast table. Input order is preserved.
build_lisa_contrast_navigation <- function(summary, contrast_id, scope, analysis_a, analysis_b,
    collection, tier, gsea_padj_cutoff = .25, label = contrast_id, label_a = analysis_a, label_b = analysis_b) {
  contrast_id <- lisa_navigation_scope_component(contrast_id, "contrast_id")
  scope <- lisa_navigation_scope_component(scope, "scope")
  collection <- lisa_navigation_scope_component(collection, "collection")
  for (id in list(analysis_a = analysis_a, analysis_b = analysis_b))
    lisa_navigation_scope_component(id, "analysis_a/analysis_b")
  if (!is.character(tier) || length(tier) != 1L || is.na(tier) || !nzchar(tier))
    stop("A single non-empty tier is required.", call. = FALSE)
  if (!is.numeric(gsea_padj_cutoff) || length(gsea_padj_cutoff) != 1L ||
      !is.finite(gsea_padj_cutoff) || gsea_padj_cutoff < 0 || gsea_padj_cutoff > 1)
    stop("gsea_padj_cutoff must lie between zero and one.", call. = FALSE)
  if (length(label) != 1L || is.na(label)) stop("label must be a single description.", call. = FALSE)
  x <- lisa_navigation_read(summary)
  if (!"category_id" %in% names(x)) stop("Contrast navigation requires category_id.", call. = FALSE)
  if (!all(c("delta_mean_NES", "n_genesets_significant_A", "n_genesets_significant_B") %in% names(x)))
    stop("Contrast navigation requires the existing category-contrast table (delta_mean_NES and side counts).", call. = FALSE)
  if (!"mean_NES_A" %in% names(x) || !"mean_NES_B" %in% names(x))
    stop("Contrast navigation accepts contrast summaries only, not a single-analysis table.", call. = FALSE)
  for (key in intersect(c("universe"), names(x))) {
    observed <- unique(as.character(x[[key]])); observed <- observed[!is.na(observed) & nzchar(observed)]
    if (length(setdiff(observed, collection))) stop("Contrast navigation ", key, " scope mismatch.", call. = FALSE)
  }
  id <- as.character(x$category_id)
  classified <- !is.na(id) & nzchar(id) & id != "OTHER_UNCLASSIFIED"
  if ("classification_status" %in% names(x))
    classified <- classified & !is.na(x$classification_status) & x$classification_status == "classified"
  x <- x[classified, , drop = FALSE]
  id <- as.character(x$category_id)
  if (anyDuplicated(id)) stop("Duplicate category IDs in contrast navigation summary.", call. = FALSE)
  num <- function(key, fallback = NULL) {
    if (!key %in% names(x)) {
      if (is.null(fallback)) stop("Missing summary column: ", key, call. = FALSE)
      return(fallback)
    }
    suppressWarnings(as.numeric(x[[key]]))
  }
  text_column <- function(key, fallback = "") {
    values <- if (key %in% names(x)) as.character(x[[key]]) else rep(fallback, nrow(x))
    values[is.na(values) | !nzchar(values)] <- rep(fallback, length.out = length(values))[is.na(values) | !nzchar(values)]
    values
  }
  labels <- text_column("category_display_name")
  labels[!nzchar(labels)] <- id[!nzchar(labels)]
  group_id <- text_column("macrogroup_id"); group_name <- text_column("macrogroup_name")
  group_id[!nzchar(group_id)] <- "UNCLASSIFIED"
  group_name[!nzchar(group_name)] <- "Other or unclassified"
  group_order <- if ("macrogroup_order" %in% names(x)) suppressWarnings(as.numeric(x$macrogroup_order)) else rep(NA_real_, nrow(x))
  category_order <- if ("category_order_within_macrogroup" %in% names(x)) suppressWarnings(as.numeric(x$category_order_within_macrogroup)) else rep(NA_real_, nrow(x))
  category_color <- text_column("color", "#737373")
  category_color[!grepl("^#[0-9A-Fa-f]{6}$", category_color)] <- "#737373"
  delta <- num("delta_mean_NES")
  # The leading count has no single non-invented total for a two-sided
  # comparison; only the exact per-side significant counts are shown.
  table <- data.frame(category_id = id, category_display_name = labels,
    macrogroup_id = group_id, macrogroup_name = group_name,
    macrogroup_order = group_order, category_order_within_macrogroup = category_order,
    category_color = category_color,
    mean_NES = delta, mean_NES_A = num("mean_NES_A"), mean_NES_B = num("mean_NES_B"), n_genesets_significant = rep(NA_real_, nrow(x)),
    n_pos_genesets = num("n_genesets_significant_A"), n_neg_genesets = num("n_genesets_significant_B"),
    support_state = ifelse(is.finite(delta), "descriptive_difference_available", "not_available"),
    source_row_order = seq_len(nrow(x)), stringsAsFactors = FALSE)
  metadata <- list(schema_version = "1.0", contrast_id = contrast_id, scope = scope,
    analysis_a = analysis_a, analysis_b = analysis_b, collection = collection, tier = tier,
    gsea_padj_cutoff = gsea_padj_cutoff, label = as.character(label), label_a=as.character(label_a), label_b=as.character(label_b),
    metric = "Existing mean_NES_A and mean_NES_B: separate profiles, without subtracting or recomputing enrichment.",
    inference = "Descriptive category comparison; no category-level p-value, FDR or new test is computed here.",
    ordering = "Preserves the supplied category-contrast row order.",
    exclusions = "Unclassified residuals are not visualized.",
    input_summary_sha256 = if (is.character(summary)) digest::digest(file = summary, algo = "sha256") else "data.frame input")
  structure(list(metadata = metadata, categories = table), class = "lisa_contrast_navigation")
}

lisa_contrast_navigation_links <- function(navigation, evidence_base) {
  if (length(evidence_base) != 1L || is.na(evidence_base) || !nzchar(evidence_base) ||
      grepl("^[A-Za-z][A-Za-z0-9+.-]*:|^[/\\\\]|[?#]", evidence_base))
    stop("evidence_base must be a relative HTML path without query or fragment.", call. = FALSE)
  if (!nrow(navigation$categories)) return(character())
  encode <- function(x) utils::URLencode(x, reserved = TRUE)
  m <- navigation$metadata
  paste0(evidence_base, "?category_id=", vapply(navigation$categories$category_id, encode, character(1L)),
    "&contrast_id=", encode(m$contrast_id), "&analysis_a=", encode(m$analysis_a),
    "&analysis_b=", encode(m$analysis_b), "&collection=", encode(m$collection), "&tier=", encode(m$tier))
}

render_lisa_contrast_navigation <- function(navigation, output_dir, evidence_base = NULL,
                                            download_policy = NULL) {
  if (!inherits(navigation, "lisa_contrast_navigation")) stop("Expected lisa_contrast_navigation object.", call. = FALSE)
  m <- navigation$metadata
  if (is.null(evidence_base)) evidence_base <- paste0("../../../contrast_evidence/",
    utils::URLencode(m$scope, reserved = TRUE), "/", utils::URLencode(m$collection, reserved = TRUE), "/index.html")
  source <- navigation$categories
  source$evidence_url <- lisa_contrast_navigation_links(navigation, evidence_base)
  source$contrast_id <- rep(m$contrast_id, nrow(source)); source$collection <- rep(m$collection, nrow(source))
  source$tier <- rep(m$tier, nrow(source)); source$gsea_padj_cutoff <- rep(m$gsea_padj_cutoff, nrow(source))
  title <- paste0("Category profiles \u00b7 ", m$label)
  subtitle <- paste0(m$collection, " \u00b7 ", m$tier, " \u00b7 member-set GSEA FDR \u2264 ",
    format(m$gsea_padj_cutoff, digits = 17), " \u00b7 Separate profiles; positive NES denotes enrichment toward the positive group in each analysis.")
  svg <- lisa_category_navigation_svg(source, title, subtitle,
    count_header = "Significant sets (A / B)",
    metric_label = "A minus B delta mean NES",
    metric_value_label = "A minus B delta mean NES",
    metric_description = "A minus B delta mean NES is descriptive, not a category significance test.",
    count_mode = "sides_only", side_labels = c(m$label_a, m$label_b),
    axis_header = "Mean NES: A and B",
    unavailable_state = "A minus B delta mean NES unavailable: one or both canonical display endpoints are unavailable")
  lisa_guarded_dir_create(output_dir)
  lisa_guarded_dir_create(file.path(output_dir, "tables"))
  source_path <- file.path(output_dir, "tables", "contrast_navigation_source.tsv")
  write_lisa_tsv(source, source_path)
  m$source_sha256 <- digest::digest(file = source_path, algo = "sha256")
  m$source_file <- "tables/contrast_navigation_source.tsv"
  m$evidence_base <- evidence_base; m$title <- title; m$subtitle <- subtitle
  svg_path <- file.path(output_dir, "category_navigation.svg")
  lisa_guarded_write(svg_path, function(path) writeLines(svg, path, useBytes = TRUE))
  metadata_path <- file.path(output_dir, "metadata.json")
  lisa_guarded_write(metadata_path, function(path) jsonlite::write_json(m, path, auto_unbox = TRUE, pretty = TRUE, digits = 17))
  escape <- lisa_navigation_escape
  requested <- function(key) is.null(download_policy) || isTRUE(download_policy[[key]])
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
    '<p>Click a category label or point to open its contrast evidence. Keyboard users can tab to a row and press Enter. Both analysis profiles are shown separately; the category means are descriptive, not category-level significance tests.</p>',
    paste0('<nav aria-label="Downloads">',
      if (requested("svg")) '<a href="category_navigation.svg" download>Download SVG</a>' else "",
      if (requested("source_data")) '<a href="tables/contrast_navigation_source.tsv" download>Download source table</a>' else "",
      '<a href="metadata.json" download>Download provenance</a>',
      if (requested("recipes")) '<a href="reproduce_contrast_navigation.R" download>R script</a>' else "", '</nav>'),
    '<div class="diagram">', html_svg, '</div>',
    '<p>Blue circles show A; orange squares show B. Both use the same NES axis. Missing support is not plotted as zero. Click a category for paired gene-set and gene evidence.</p></main></body></html>')
  html <- lisa_present_evidence_html(html, output_dir, m, active = "contrasts")
  html_path <- file.path(output_dir, "index.html")
  lisa_guarded_write(html_path, function(path) writeLines(html, path, useBytes = TRUE))
  recipe <- c('#!/usr/bin/env Rscript', '# Offline base-R SVG reproduction. Run from this navigator directory:',
    '# Rscript --vanilla reproduce_contrast_navigation.R tables/contrast_navigation_source.tsv reproduced.svg',
    '# Preserve the report directory structure for relative evidence links.',
    paste0('lisa_category_navigation_svg <- ', paste(deparse(lisa_category_navigation_svg), collapse = '\n')),
    paste0('side_labels <- ', paste(deparse(c(m$label_a,m$label_b)),collapse='')),
    paste0('title <- ', encodeString(title, quote = '"')), paste0('subtitle <- ', encodeString(subtitle, quote = '"')),
    'args <- commandArgs(TRUE)', 'if (length(args) != 2L) stop("Expected source TSV and destination SVG.")',
    'if (tolower(tools::file_ext(args[[2L]])) != "svg") stop("Destination must be SVG.")',
    'source <- utils::read.delim(args[[1L]], sep="\\t", quote="", comment.char="", check.names=FALSE, colClasses="character", na.strings=NULL)',
    'writeLines(lisa_category_navigation_svg(source, title, subtitle, count_header = "Significant sets (A / B)", metric_label = "A minus B delta mean NES", metric_value_label = "A minus B delta mean NES", metric_description = "A minus B delta mean NES is descriptive, not a category significance test.", count_mode = "sides_only", side_labels = side_labels, axis_header = "Mean NES: A and B", unavailable_state = "A minus B delta mean NES unavailable: one or both canonical display endpoints are unavailable"), args[[2L]], useBytes=TRUE)')
  lisa_guarded_write(file.path(output_dir, "reproduce_contrast_navigation.R"), function(path) writeLines(recipe, path, useBytes = TRUE))
  invisible(list(html = html_path, svg = svg_path, source = source_path, metadata = metadata_path))
}
