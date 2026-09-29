# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Shared native leading-edge heatmap drawing and self-contained recipe contract.
# Loading this helper performs no plotting, transformation, or data selection.

lisa_heatmap_num <- function(x) suppressWarnings(as.numeric(x))

lisa_heatmap_subtitle <- function(metadata, prefix = "", suffix = "") {
  values <- c(
    prefix,
    metadata$comparison,
    metadata$positive_direction,
    suffix
  )
  values <- trimws(as.character(values))
  values <- values[!is.na(values) & nzchar(values)]
  lisa_heatmap_device_text(paste(values, collapse = "\n"))
}

lisa_heatmap_device_text <- function(x) {
  # Base PDF and some bitmap devices cannot encode typographic dash glyphs
  # with their default font. Keep the source metadata unchanged and normalize
  # only the human-readable text passed to graphics devices.
  gsub("[\u2013\u2014\u2212]", "-", as.character(x), perl = TRUE)
}

lisa_heatmap_fdr <- function(x) {
  x <- lisa_heatmap_num(x)
  out <- rep("NA", length(x))
  ok <- is.finite(x)
  out[ok & x < 0.001] <- formatC(x[ok & x < 0.001], format = "e", digits = 1)
  out[ok & x >= 0.001] <- formatC(x[ok & x >= 0.001], format = "f", digits = 3)
  out
}

lisa_heatmap_plot <- function(df, cfg, analysis_id, universe, category_label, count_source) {
  ordered_genes <- unique(df$symbol[order(-abs(lisa_heatmap_num(df$log2FC)), df$symbol)])
  df$symbol <- factor(df$symbol, levels = rev(ordered_genes))
  sample_levels <- levels(df$sample)
  x_levels <- c(sample_levels, "log2FC", "FDR")
  df$plot_col <- factor(as.character(df$sample), levels = x_levels)
  gene_stats <- df[!duplicated(as.character(df$symbol)), c("symbol", "log2FC", "padj"), drop = FALSE]
  gene_stats$symbol <- factor(as.character(gene_stats$symbol), levels = levels(df$symbol))
  stat_tiles <- do.call(rbind, lapply(c("log2FC", "FDR"), function(col) {
    data.frame(symbol = gene_stats$symbol, plot_col = factor(col, levels = x_levels), stringsAsFactors = FALSE)
  }))
  stat_text <- rbind(
    data.frame(
      symbol = gene_stats$symbol,
      plot_col = factor("log2FC", levels = x_levels),
      label = formatC(lisa_heatmap_num(gene_stats$log2FC), format = "f", digits = 2),
      stringsAsFactors = FALSE
    ),
    data.frame(
      symbol = gene_stats$symbol,
      plot_col = factor("FDR", levels = x_levels),
      label = lisa_heatmap_fdr(gene_stats$padj),
      stringsAsFactors = FALSE
    )
  )
  fill_name <- if (cfg$scale == "zscore") "per-gene z-score" else count_source
  ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = stat_tiles,
      ggplot2::aes(x = plot_col, y = symbol),
      fill = "grey95",
      color = "white",
      linewidth = 0.12,
      width = 0.86,
      height = 0.86
    ) +
    ggplot2::geom_tile(
      data = df,
      ggplot2::aes(x = plot_col, y = symbol, fill = plot_value),
      color = "white",
      linewidth = 0.12,
      width = 0.86,
      height = 0.86
    ) +
    ggplot2::geom_text(
      data = stat_text,
      ggplot2::aes(x = plot_col, y = symbol, label = label),
      size = 1.75,
      color = "grey20",
      na.rm = TRUE
    ) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, na.value = "grey88", name = fill_name) +
    ggplot2::scale_x_discrete(limits = x_levels, drop = FALSE, expand = ggplot2::expansion(add = 0.08)) +
    ggplot2::labs(
      title = "Leading-edge genes by LISA category",
      subtitle = lisa_heatmap_subtitle(
        cfg$plot_metadata,
        prefix = category_label,
        suffix = sprintf(
          "%s | rows ordered by abs(log2FC); sample values: %s",
          universe, fill_name
        )
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 8.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 6.5),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 6.5),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 10),
      plot.subtitle = ggplot2::element_text(size = 7.2, color = "grey35", lineheight = 1.15),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

lisa_heatmap_size <- function(n_genes, n_samples, n_stat_cols = 2L) {
  n_cols <- n_samples + n_stat_cols
  # Rows and columns need independent space. The former square 7.1-inch cap
  # compressed 52 labelled sample columns into roughly half the width needed
  # for the real dense Riaz matrices. Retain bounded vector/bitmap dimensions
  # while allowing each sample column and gene row to remain readable.
  width <- max(6.4, min(14.5, 4.2 + n_cols * 0.19))
  height <- max(5.6, min(10.5, 4.2 + n_genes * 0.16))
  c(width = width, height = height)
}

# Keep device operations shared as well: PDF must not inherit bitmap dpi and
# SVG uses the native Cairo device rather than an optional svglite dependency.
lisa_heatmap_save_native <- function(path, plot, width, height) {
  format <- tolower(tools::file_ext(path))
  if (format == "png") {
    ggplot2::ggsave(path, plot, width = width, height = height, dpi = 220, bg = "white", limitsize = FALSE)
  } else if (format == "pdf") {
    ggplot2::ggsave(path, plot, width = width, height = height, bg = "white", limitsize = FALSE)
  } else if (format == "svg") {
    if (!isTRUE(capabilities("cairo")))
      stop("LISA-HEATMAP-RECIPE-006 native SVG requires this R runtime's Cairo capability.", call. = FALSE)
    grDevices::svg(filename = path, width = width, height = height, bg = "white")
    on.exit(grDevices::dev.off(), add = TRUE)
    print(plot)
  } else {
    stop("LISA-HEATMAP-RECIPE-004 output extension must be PNG, PDF or SVG.", call. = FALSE)
  }
  invisible(path)
}

# Drawing context is serialized into each emitted recipe, not added to the
# scientific matrix. The recipe can be copied with that matrix and run offline.
lisa_heatmap_context <- function(df, cfg, analysis_id, universe, category_label,
    count_source, source_path) {
  required <- c("symbol", "sample", "category_id", "analysis_id", "collection",
    "category_display_name", "scale", "log2FC", "padj", "plot_value")
  if (!nrow(df) || !all(required %in% names(df)))
    stop("LISA-HEATMAP-RECIPE-001 native source lacks exact scope or saved geometry.", call. = FALSE)
  scalar <- function(key) {
    value <- unique(as.character(df[[key]]))
    if (length(value) != 1L || is.na(value) || !nzchar(value))
      stop("LISA-HEATMAP-RECIPE-001 source scope is ambiguous: ", key, call. = FALSE)
    value
  }
  scope <- lapply(c("category_id", "analysis_id", "collection", "category_display_name", "scale"), scalar)
  names(scope) <- c("category_id", "analysis_id", "collection", "category_display_name", "scale")
  if (!identical(scope$analysis_id, analysis_id) || !identical(scope$collection, universe) ||
      !identical(scope$category_display_name, category_label) || !identical(scope$scale, cfg$scale))
    stop("LISA-HEATMAP-RECIPE-001 plot context differs from its saved category matrix.", call. = FALSE)
  genes <- unique(as.character(df$symbol[order(-abs(lisa_heatmap_num(df$log2FC)), df$symbol)]))
  samples <- if (is.factor(df$sample)) levels(df$sample) else unique(as.character(df$sample))
  dims <- lisa_heatmap_size(length(genes), length(samples), 2L)
  fill_name <- if (cfg$scale == "zscore") "per-gene z-score" else count_source
  if (!is.character(fill_name) || length(fill_name) != 1L || is.na(fill_name) || !nzchar(fill_name))
    stop("LISA-HEATMAP-RECIPE-001 explicit display-scale label is required.", call. = FALSE)
  list(schema_version = "lisa-native-leading-edge-heatmap-v1", scope = scope,
    source_sha256 = digest::digest(file = source_path, algo = "sha256"),
    source_bytes = unname(file.info(source_path)$size), source_rows = nrow(df),
    gene_order = genes, sample_order = samples, plot_metadata = cfg$plot_metadata,
    count_source = count_source, fill_name = fill_name,
    title = "Leading-edge genes by LISA category",
    subtitle = lisa_heatmap_subtitle(cfg$plot_metadata, prefix = category_label,
      suffix = sprintf("%s | rows ordered by abs(log2FC); sample values: %s", universe, fill_name)),
    width = unname(dims[["width"]]), height = unname(dims[["height"]]), dpi = 220,
    background = "white", value_policy = "saved plot_value; no new transformations or selection")
}

lisa_heatmap_read_bound_source <- function(source_path, context) {
  if (!identical(context$schema_version, "lisa-native-leading-edge-heatmap-v1") ||
      !identical(digest::digest(file = source_path, algo = "sha256"), context$source_sha256) ||
      !identical(as.numeric(file.info(source_path)$size), as.numeric(context$source_bytes)))
    stop("LISA-HEATMAP-RECIPE-002 source matrix does not match the recipe's recorded SHA-256/contract.", call. = FALSE)
  df <- utils::read.delim(source_path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = FALSE, colClasses = "character", na.strings = NULL)
  required <- c("symbol", "sample", "category_id", "analysis_id", "collection",
    "category_display_name", "scale", "log2FC", "padj", "plot_value")
  if (anyDuplicated(names(df)) || !all(required %in% names(df)) || nrow(df) != context$source_rows)
    stop("LISA-HEATMAP-RECIPE-003 saved matrix geometry is incomplete.", call. = FALSE)
  for (key in names(context$scope)) {
    if (!identical(unique(df[[key]]), context$scope[[key]]))
      stop("LISA-HEATMAP-RECIPE-003 exact source scope mismatch: ", key, call. = FALSE)
  }
  genes <- context$gene_order; samples <- context$sample_order
  if (!length(genes) || length(samples) < 2L || anyDuplicated(genes) || anyDuplicated(samples) ||
      anyNA(genes) || anyNA(samples) || any(!nzchar(genes)) || any(!nzchar(samples)) ||
      any(grepl("[[:cntrl:]]", c(genes, samples))) || any(samples %in% c("log2FC", "FDR")) ||
      !identical(unique(df$symbol), genes) || !identical(unique(df$sample), samples) ||
      nrow(df) != length(genes) * length(samples) || anyDuplicated(paste(df$symbol, df$sample, sep = "\r")))
    stop("LISA-HEATMAP-RECIPE-003 source has no unique exact gene/sample grid and orders.", call. = FALSE)
  for (gene in genes) {
    rows <- df[df$symbol == gene, , drop = FALSE]
    if (!identical(rows$sample, samples) || length(unique(rows$log2FC)) != 1L || length(unique(rows$padj)) != 1L)
      stop("LISA-HEATMAP-RECIPE-003 sample order or DE annotation varies within a gene.", call. = FALSE)
  }
  for (key in c("plot_value", "log2FC", "padj")) {
    value <- lisa_heatmap_num(df[[key]])
    missing <- df[[key]] %in% c("NA", "NaN", "")
    if (any(!missing & !is.finite(value)))
      stop("LISA-HEATMAP-RECIPE-003 nonnumeric/infinite saved value: ", key, call. = FALSE)
    df[[key]] <- value
  }
  if (!any(is.finite(df$plot_value)))
    stop("LISA-HEATMAP-RECIPE-003 no finite saved display values.", call. = FALSE)
  # Preserve explicit serialized orders, rather than silently changing tied
  # labels under a different locale. Generation and reproduction share plotting.
  native_order <- unique(df$symbol[order(-abs(lisa_heatmap_num(df$log2FC)), df$symbol)])
  if (!identical(native_order, genes))
    stop("LISA-HEATMAP-RECIPE-003 runtime gene ordering differs from recorded native order.", call. = FALSE)
  df$sample <- factor(df$sample, levels = samples)
  df
}

lisa_heatmap_render_saved <- function(source_path, output_path, context) {
  for (package in c("ggplot2", "digest")) if (!requireNamespace(package, quietly = TRUE))
    stop("Required offline recipe package: ", package, call. = FALSE)
  source_path <- normalizePath(source_path, winslash = "/", mustWork = TRUE)
  output_parent <- normalizePath(dirname(output_path), winslash = "/", mustWork = TRUE)
  output_path <- file.path(output_parent, basename(output_path))
  link <- Sys.readlink(output_path)
  if (!tolower(tools::file_ext(output_path)) %in% c("png", "pdf", "svg") || file.exists(output_path) ||
      dir.exists(output_path) || (!is.na(link) && nzchar(link)) || identical(source_path, output_path))
    stop("LISA-HEATMAP-RECIPE-004 use a fresh PNG, PDF or SVG output; original files are never overwritten.", call. = FALSE)
  df <- lisa_heatmap_read_bound_source(source_path, context)
  cfg <- list(scale = context$scope$scale, plot_metadata = context$plot_metadata)
  fill_name <- if (cfg$scale == "zscore") "per-gene z-score" else context$count_source
  subtitle <- lisa_heatmap_subtitle(cfg$plot_metadata, prefix = context$scope$category_display_name,
    suffix = sprintf("%s | rows ordered by abs(log2FC); sample values: %s", context$scope$collection, fill_name))
  dims <- lisa_heatmap_size(length(context$gene_order), length(context$sample_order), 2L)
  if (!cfg$scale %in% c("zscore", "raw") || !identical(fill_name, context$fill_name) ||
      !identical(subtitle, context$subtitle) || !identical(context$title, "Leading-edge genes by LISA category") ||
      !isTRUE(all.equal(unname(dims), c(context$width, context$height), tolerance = 1e-12)) ||
      context$dpi != 220 || !identical(context$background, "white"))
    stop("LISA-HEATMAP-RECIPE-003 recorded native presentation context is inconsistent.", call. = FALSE)
  plot <- lisa_heatmap_plot(df, cfg, context$scope$analysis_id, context$scope$collection,
    context$scope$category_display_name, context$count_source)
  lisa_heatmap_save_native(output_path, plot, width = context$width, height = context$height)
  if (!file.exists(output_path) || file.info(output_path)$size <= 0 ||
      !identical(digest::digest(file = source_path, algo = "sha256"), context$source_sha256))
    stop("LISA-HEATMAP-RECIPE-005 Figure generation failed or source changed.", call. = FALSE)
  invisible(output_path)
}

lisa_heatmap_write_recipe <- function(recipe_path, context, run_root = NULL) {
  # Embed the same native functions used by the builder, with safely serialized
  # data. The emitted file has no source(), project/config or lisaR dependency.
  functions <- c("lisa_heatmap_num", "lisa_heatmap_fdr", "lisa_heatmap_subtitle",
    "lisa_heatmap_device_text", "lisa_heatmap_plot", "lisa_heatmap_size", "lisa_heatmap_save_native",
    "lisa_heatmap_read_bound_source", "lisa_heatmap_render_saved")
  env <- environment(lisa_heatmap_write_recipe)
  definitions <- unlist(lapply(functions, function(name) c(
    paste0(name, " <- ", paste(deparse(get(name, envir = env, inherits = FALSE), width.cutoff = 500L), collapse = "\n")), "")), use.names = FALSE)
  serialized <- paste(capture.output(dput(context)), collapse = "\n")
  code <- c("#!/usr/bin/env Rscript",
    "# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones Científicas (CSIC).",
    "# Author: David Olmeda Casadomé. GNU General Public License version 3 (GPL-3).",
    "# Native leading-edge heatmap recipe. Offline ggplot2 and digest are required.",
    "# Usage: Rscript --vanilla recipe.R SOURCE_MATRIX.tsv NEW_OUTPUT.png|pdf|svg",
    "# Saved plot_value only; no new transformation, selection, DE or enrichment.",
    "options(stringsAsFactors = FALSE)", definitions,
    paste0("lisa_heatmap_saved_context <- ", serialized),
    "args <- commandArgs(trailingOnly = TRUE)",
    "if (length(args) != 2L) stop('Usage: Rscript --vanilla recipe.R SOURCE_MATRIX.tsv NEW_OUTPUT.png|pdf|svg', call. = FALSE)",
    "lisa_heatmap_render_saved(args[[1L]], args[[2L]], lisa_heatmap_saved_context)",
    "message(normalizePath(args[[2L]], winslash = '/', mustWork = TRUE))")
  lisaR:::lisa_guarded_write(recipe_path, function(target) writeLines(code, target, useBytes = TRUE),
    run_root = run_root)
  invisible(recipe_path)
}
