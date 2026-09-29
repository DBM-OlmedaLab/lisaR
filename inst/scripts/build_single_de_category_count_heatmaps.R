#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

default_project_dir <- NA_character_

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
}

write_tsv <- function(df, path) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

save_plot_svg <- function(path, plot, width, height) {
  grDevices::svg(filename = path, width = width, height = height, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot)
}

parse_args <- function(args) {
  out <- list(
    project_dir = default_project_dir,
    universe = "GOBP-C2",
    max_categories = 0
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    if (i == length(args)) stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    val <- args[[i + 1]]
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(out)) stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    if (key %in% c("max_categories")) {
      out[[key]] <- as.integer(val)
    } else {
      out[[key]] <- val
    }
    i <- i + 2
  }
  if (is.na(out$project_dir) || out$project_dir == "") stop("Required argument: --project-dir", call. = FALSE)
  if (is.na(out$universe) || out$universe == "") stop("Required argument: --universe", call. = FALSE)
  out
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_file_component <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

short_label <- function(x, max_chars = 46) {
  x <- as.character(x)
  too_long <- nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1, max_chars - 1), "...")
  x
}

find_summary_files <- function(project_dir, universe, suffix) {
  output_root <- file.path(project_dir, "outputs")
  patterns <- c(
    file.path(output_root, "single_de", "*", paste0("collection_", universe), "lisa_tables", paste0("*_", suffix, "_category_summary.tsv")),
    file.path(output_root, "lisa", "single_de", "*", paste0("collection_", universe), "lisa_tables", paste0("*_", suffix, "_category_summary.tsv"))
  )
  unique(unlist(lapply(patterns, Sys.glob), use.names = FALSE))
}

parse_analysis_id <- function(path) {
  parts <- strsplit(normalizePath(path, winslash = "/", mustWork = FALSE), "/", fixed = TRUE)[[1]]
  idx <- match("single_de", parts)
  if (is.na(idx) || idx + 1 > length(parts)) return("")
  parts[[idx + 1]]
}

load_gsea_counts <- function(project_dir, universe) {
  files <- find_summary_files(project_dir, universe, "GSEA")
  rows <- lapply(files, function(path) {
    df <- read_tsv(path)
    if (nrow(df) == 0) return(data.frame())
    df$analysis_id <- parse_analysis_id(path)
    df$count_type <- "GSEA_n_genesets"
    df$count_value <- safe_num(df$n_genesets)
    df$direction_value <- safe_num(df$mean_NES)
    df
  })
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

load_ora_counts <- function(project_dir, universe) {
  files <- find_summary_files(project_dir, universe, "ORA")
  rows <- lapply(files, function(path) {
    df <- read_tsv(path)
    if (nrow(df) == 0) return(data.frame())
    df$analysis_id <- parse_analysis_id(path)
    df$count_type <- paste0("ORA_", df$direction, "_n_genesets")
    df$count_value <- safe_num(df$n_genesets)
    df$direction_value <- ifelse(df$direction == "UP", 1, -1)
    df
  })
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

select_categories <- function(df, cfg) {
  score <- stats::aggregate(
    count_value ~ category_id + category_display_name + macrogroup_id + macrogroup_name +
      macrogroup_order + category_order_within_macrogroup + color,
    data = df,
    FUN = max
  )
  score <- score[order(
    -safe_num(score$count_value),
    safe_num(score$macrogroup_order),
    safe_num(score$category_order_within_macrogroup),
    score$category_id
  ), , drop = FALSE]
  if (cfg$max_categories > 0) {
    score <- head(score, cfg$max_categories)
  }
  score$category_id
}

plot_count_heatmap <- function(df, cfg, title, count_label) {
  categories <- select_categories(df, cfg)
  df <- df[df$category_id %in% categories, , drop = FALSE]
  category_order <- unique(df[order(
    safe_num(df$macrogroup_order),
    safe_num(df$category_order_within_macrogroup),
    df$category_id
  ), "category_id"])
  df$category_label <- factor(short_label(df$category_display_name), levels = rev(short_label(df$category_display_name[match(category_order, df$category_id)])))
  df$analysis_label <- factor(df$analysis_id, levels = unique(df$analysis_id))
  df$count_value[is.na(df$count_value)] <- 0
  df$cell_label <- ifelse(df$count_value > 0, as.character(df$count_value), "")
  max_count <- max(df$count_value, na.rm = TRUE)
  if (!is.finite(max_count) || max_count <= 0) max_count <- 1

  ggplot2::ggplot(df, ggplot2::aes(x = analysis_label, y = category_label, fill = count_value)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = cell_label), size = 2.7, color = "grey15") +
    ggplot2::scale_fill_gradient(low = "#F7F8F9", high = "#9D332F", limits = c(0, max_count), name = count_label) +
    ggplot2::labs(
      title = title,
      subtitle = sprintf("%s | category-count heatmap; cell values are source gene-set counts", cfg$universe),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1, color = "grey15"),
      axis.text.y = ggplot2::element_text(color = "grey15"),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey35"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "right"
    )
}

write_layer <- function(df, paths, cfg, layer_name, title, count_label) {
  if (nrow(df) == 0) return(data.frame())
  max_count <- max(safe_num(df$count_value), na.rm = TRUE)
  if (!is.finite(max_count) || max_count <= 0) return(data.frame())
  categories <- select_categories(df, cfg)
  selected <- df[df$category_id %in% categories, , drop = FALSE]
  summary_path <- file.path(paths$out_dir, paste0(paths$prefix, "_", layer_name, "_category_count_matrix.tsv"))
  write_tsv(selected, summary_path)
  plot <- plot_count_heatmap(selected, cfg, title, count_label)
  n_categories <- length(unique(selected$category_id))
  n_columns <- length(unique(selected$analysis_id))
  width <- max(7.2, 2.2 + n_columns * 1.35)
  height <- max(6.4, 2.2 + n_categories * 0.18)
  png_path <- file.path(paths$out_dir, paste0(paths$prefix, "_", layer_name, "_category_count_heatmap.png"))
  pdf_path <- file.path(paths$out_dir, paste0(paths$prefix, "_", layer_name, "_category_count_heatmap.pdf"))
  svg_path <- file.path(paths$out_dir, paste0(paths$prefix, "_", layer_name, "_category_count_heatmap.svg"))
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 230, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(pdf_path, plot, width = width, height = height, bg = "white", limitsize = FALSE)
  save_plot_svg(svg_path, plot, width = width, height = height)
  data.frame(
    universe = cfg$universe,
    layer = layer_name,
    n_categories = n_categories,
    n_columns = n_columns,
    max_categories = cfg$max_categories,
    matrix_tsv = summary_path,
    plot_png = png_path,
    plot_pdf = pdf_path,
    plot_svg = svg_path,
    stringsAsFactors = FALSE
  )
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Required package: ggplot2", call. = FALSE)
  paths <- list(
    out_dir = file.path(cfg$project_dir, "outputs", "category_count_heatmaps", paste0("collection_", cfg$universe)),
    prefix = safe_file_component(cfg$universe)
  )
  dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)
  index <- list()
  gsea <- load_gsea_counts(cfg$project_dir, cfg$universe)
  if (nrow(gsea) > 0) {
    index[[length(index) + 1L]] <- write_layer(gsea, paths, cfg, "GSEA", "GSEA category count heatmap", "GSEA gene sets")
  }
  ora <- load_ora_counts(cfg$project_dir, cfg$universe)
  if (nrow(ora) > 0) {
    index[[length(index) + 1L]] <- write_layer(ora, paths, cfg, "ORA", "ORA category count heatmap", "ORA gene sets")
  }
  index_df <- if (length(index) == 0) data.frame() else do.call(rbind, index)
  index_path <- file.path(paths$out_dir, paste0(paths$prefix, "_category_count_heatmaps_index.tsv"))
  write_tsv(index_df, index_path)
  if (nrow(index_df) == 0) {
    message("No category count heatmaps generated; wrote empty index.")
  } else {
    print(index_df[, c("universe", "layer", "n_categories", "n_columns")])
  }
}

if (sys.nframe() == 0L) main()
