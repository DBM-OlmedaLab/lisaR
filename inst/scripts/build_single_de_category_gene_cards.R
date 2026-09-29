#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

default_project_dir <- NA_character_
script_args <- commandArgs(FALSE)
script_candidates <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
script_dir <- if (length(script_candidates) > 0L) {
  dirname(normalizePath(script_candidates[[1]], mustWork = TRUE))
} else {
  getwd()
}

read_tsv <- function(path) {
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
}

write_tsv <- function(df, path) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

parse_args <- function(args) {
  out <- list(
    project_dir = default_project_dir,
    analysis_id = NA_character_,
    universe = "GOBP-C2",
    max_categories = 0,
    top_genes = 15,
    categories = NA_character_,
    formats = "png",
    lisa_internal_renderer_sha256 = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    }
    val <- args[[i + 1]]
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(out)) {
      stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    }
    if (key %in% c("max_categories", "top_genes")) {
      out[[key]] <- as.integer(val)
    } else {
      out[[key]] <- val
    }
    i <- i + 2
  }
  if (is.na(out$project_dir) || out$project_dir == "") {
    stop("Required argument: --project-dir", call. = FALSE)
  }
  if (is.na(out$analysis_id) || out$analysis_id == "") {
    stop("Required argument: --analysis-id", call. = FALSE)
  }
  out$formats <- unique(trimws(strsplit(tolower(out$formats), ",", fixed = TRUE)[[1]]))
  if (any(!out$formats %in% c("png", "svg", "pdf"))) stop("--formats must contain png, svg and/or pdf", call. = FALSE)
  out
}

gene_level_paths <- function(project_dir, analysis_id, universe) {
  collection_dir <- file.path(project_dir, "outputs", "single_de", analysis_id, paste0("collection_", universe))
  out_dir <- file.path(project_dir, "outputs", "gene_level", "single_de", analysis_id, paste0("collection_", universe))
  prefix <- paste(analysis_id, universe, "gene_level", sep = "_")
  summary_candidates <- list.files(
    file.path(collection_dir, "lisa_tables"),
    pattern = "_GSEA_category_summary\\.tsv$",
    full.names = TRUE
  )
  summary_candidates <- summary_candidates[
    startsWith(basename(summary_candidates), paste0(analysis_id, "_"))
  ]
  list(
    out_dir = out_dir,
    gsea_summary = if (length(summary_candidates) == 1) summary_candidates[[1]] else "",
    evidence = file.path(out_dir, paste0(prefix, "_gene_category_contributions.tsv")),
    summary = file.path(out_dir, paste0(prefix, "_category_gene_support_summary.tsv")),
    cards_dir = file.path(out_dir, "category_gene_cards"),
    prefix = prefix
  )
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

gene_card_lfc_breaks <- function(x) {
  values <- safe_num(x)
  values <- values[is.finite(values)]
  if (length(values) == 0L) return(numeric())

  value_range <- range(values)
  if (diff(value_range) == 0) return(value_range[[1L]])
  if (value_range[[1]] < 0 && value_range[[2]] > 0) {
    return(c(value_range[[1]], 0, value_range[[2]]))
  }
  c(value_range[[1]], mean(value_range), value_range[[2]])
}

safe_label <- function(x, max_chars = 42) {
  x <- as.character(x)
  too_long <- nchar(x) > max_chars
  x[too_long] <- paste0(substr(x[too_long], 1, max_chars - 1), "...")
  x
}

safe_file_component <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

format_padj <- function(x) {
  x <- safe_num(x)
  ifelse(is.na(x), "NA", ifelse(x < 0.001, formatC(x, format = "e", digits = 1), sprintf("%.3f", x)))
}

select_categories <- function(evidence, summary, cfg) {
  requested <- character()
  if (!is.na(cfg$categories) && cfg$categories != "") {
    requested <- trimws(strsplit(cfg$categories, ",", fixed = TRUE)[[1]])
  }
  if (length(requested) > 0) {
    found <- requested[requested %in% unique(evidence$category_id)]
    missing <- setdiff(requested, found)
    if (length(missing) > 0) {
      warning(sprintf("Requested categories not found and skipped: %s", paste(missing, collapse = ", ")), call. = FALSE)
    }
    return(found)
  }

  paths <- gene_level_paths(cfg$project_dir, cfg$analysis_id, cfg$universe)
  if (file.exists(paths$gsea_summary)) {
    gsea <- read_tsv(paths$gsea_summary)
    if ("n_genesets" %in% names(gsea)) {
      gsea <- gsea[safe_num(gsea$n_genesets) > 0, , drop = FALSE]
    }
    if (nrow(gsea) == 0) {
      return(character())
    }
    if (nrow(gsea) > 0) {
      if (all(c("macrogroup_order", "category_order_within_macrogroup") %in% names(gsea))) {
        gsea <- gsea[order(
          safe_num(gsea$macrogroup_order),
          safe_num(gsea$category_order_within_macrogroup),
          gsea$category_id
        ), , drop = FALSE]
      }
      category_ids <- unique(as.character(gsea$category_id))
      category_ids <- category_ids[category_ids %in% unique(evidence$category_id)]
      if (cfg$max_categories > 0) {
        category_ids <- head(category_ids, cfg$max_categories)
      }
      return(category_ids)
    }
  }

  summary <- merge(
    summary,
    unique(evidence[, c("category_id", "macrogroup_id", "macrogroup_name")]),
    by = "category_id",
    all.x = TRUE,
    sort = FALSE
  )
  summary <- summary[order(-safe_num(summary$max_gene_contribution_score),
                           -safe_num(summary$n_supporting_genes),
                           summary$category_id), , drop = FALSE]
  if (cfg$max_categories > 0) {
    return(head(summary$category_id, cfg$max_categories))
  }
  summary$category_id
}

build_card_plot <- function(evidence, category_id, top_genes) {
  category_rows <- evidence[evidence$category_id == category_id, , drop = FALSE]
  category_rows <- category_rows[order(-safe_num(category_rows$gene_contribution_score),
                                       safe_num(category_rows$padj),
                                       -abs(safe_num(category_rows$log2FC)),
                                       category_rows$symbol), , drop = FALSE]
  if (top_genes > 0) {
    category_rows <- head(category_rows, top_genes)
  }
  category_rows$gene_label <- factor(category_rows$symbol, levels = rev(category_rows$symbol))
  category_rows$score <- safe_num(category_rows$gene_contribution_score)
  category_rows$abs_rank <- abs(safe_num(category_rows$rank_value))
  category_rows$minus_log10_padj <- -log10(pmax(safe_num(category_rows$padj), 1e-300))
  category_rows$minus_log10_padj[is.infinite(category_rows$minus_log10_padj)] <- NA_real_
  category_rows$evidence_label <- ifelse(
    category_rows$in_gsea_leading_edge & category_rows$in_ora_overlap, "GSEA+ORA",
    ifelse(category_rows$in_gsea_leading_edge, "GSEA",
           ifelse(category_rows$in_ora_overlap, "ORA", "other"))
  )

  title <- unique(category_rows$category_display_name)[1]
  macrogroup <- unique(category_rows$macrogroup_name)[1]
  category_color <- unique(category_rows$color)[1]
  if (is.na(category_color) || category_color == "") category_color <- "#525252"

  header <- sprintf(
    "%s\n%s | %s supporting genes in full table | showing top %s",
    title,
    macrogroup,
    length(unique(evidence$symbol[evidence$category_id == category_id])),
    nrow(category_rows)
  )
  lfc_breaks <- gene_card_lfc_breaks(category_rows$log2FC)

  header_grob <- grid::grobTree(
    grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)),
    grid::textGrob(
      header,
      x = grid::unit(8, "pt"),
      y = grid::unit(0.95, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(
        col = "grey10",
        fontsize = 11,
        fontface = "bold",
        lineheight = 1.05
      )
    )
  )
  p_header <- patchwork::wrap_elements(full = header_grob, clip = FALSE)

  p_lfc <- ggplot2::ggplot(category_rows, ggplot2::aes(x = gene_label, y = 1, fill = log2FC)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35, height = 0.82, width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", log2FC)), size = 2.55, color = "grey10") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient2(
      low = "#3B66B0",
      mid = "white",
      high = "#B73A3A",
      midpoint = 0,
      name = "log2FC",
      breaks = lfc_breaks,
      guide = ggplot2::guide_colorbar(
        barwidth = grid::unit(30, "mm"),
        barheight = grid::unit(3, "mm"),
        title.position = "top",
        label.position = "bottom"
      )
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", color = "grey15"),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(2, 2, 2, 6)
    )

  p_score <- ggplot2::ggplot(category_rows, ggplot2::aes(x = gene_label, y = score)) +
    ggplot2::geom_col(fill = category_color, width = 0.72) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", score)), hjust = -0.1, size = 2.45, color = "grey20") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.20))) +
    ggplot2::labs(x = NULL, y = "DE-weighted contribution") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.margin = ggplot2::margin(2, 8, 2, 8)
    )

  marks <- rbind(
    data.frame(symbol = category_rows$symbol, gene_label = category_rows$gene_label, mark = "DE FDR <= 0.05", mark_short = "FDR<0.05", value = category_rows$de_is_significant),
    data.frame(symbol = category_rows$symbol, gene_label = category_rows$gene_label, mark = "GSEA leading edge", mark_short = "LE", value = category_rows$in_gsea_leading_edge),
    data.frame(symbol = category_rows$symbol, gene_label = category_rows$gene_label, mark = "ORA overlap", mark_short = "ORA", value = category_rows$in_ora_overlap)
  )
  marks$mark_short <- factor(marks$mark_short, levels = c("FDR<0.05", "LE", "ORA"))
  marks$value <- as.logical(marks$value)
  p_marks <- ggplot2::ggplot(marks, ggplot2::aes(x = mark_short, y = gene_label)) +
    ggplot2::geom_point(
      ggplot2::aes(fill = value, color = value),
      shape = 21,
      size = 2.9,
      stroke = 0.85
    ) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#0AA36E", `FALSE` = "white")) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#0A7F58", `FALSE` = "#C8CDD3")) +
    ggplot2::labs(x = NULL, y = NULL, title = "Evidence") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 7.3, color = "grey25"),
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 8.5, face = "bold", color = "grey20", hjust = 0.5),
      plot.margin = ggplot2::margin(2, 6, 2, 2)
    )

  recurrence_df <- category_rows
  recurrence_df$category_recurrence <- safe_num(recurrence_df$category_recurrence)
  recurrence_df$recurrence_label <- sprintf("%s", recurrence_df$category_recurrence)
  max_recurrence <- max(recurrence_df$category_recurrence, na.rm = TRUE)
  if (!is.finite(max_recurrence)) max_recurrence <- 1
  p_recur <- ggplot2::ggplot(recurrence_df, ggplot2::aes(y = gene_label)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = category_recurrence, yend = gene_label),
      color = "#D2D7DD",
      linewidth = 0.55,
      lineend = "round",
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = category_recurrence),
      color = "#4C6A92",
      fill = "#DDE8F7",
      shape = 21,
      size = 2.35,
      stroke = 0.55,
      na.rm = TRUE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = max_recurrence * 1.08, label = recurrence_label),
      hjust = 0,
      size = 2.35,
      color = "grey20",
      na.rm = TRUE
    ) +
    ggplot2::scale_x_continuous(limits = c(0, max_recurrence * 1.55), expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL, title = "LISA categories\ncontaining this gene") +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.title = ggplot2::element_text(size = 7.6, face = "bold", color = "grey20", hjust = 0.5, lineheight = 0.95),
      plot.margin = ggplot2::margin(2, 6, 2, 2)
    )

  table_cols <- c("symbol", "log2FC", "padj", "evidence_label", "source_geneset_count", "category_recurrence")
  if ("lisa_support_score" %in% names(category_rows)) {
    table_cols <- c("symbol", "log2FC", "padj", "gene_contribution_score", "lisa_support_score", "evidence_label", "source_geneset_count", "category_recurrence")
  }
  table_rows <- category_rows[, table_cols, drop = FALSE]
  table_rows$log2FC <- sprintf("%.2f", safe_num(table_rows$log2FC))
  table_rows$padj <- format_padj(table_rows$padj)
  if ("lisa_support_score" %in% names(category_rows)) {
    table_rows$gene_contribution_score <- sprintf("%.2f", safe_num(table_rows$gene_contribution_score))
    table_rows$lisa_support_score <- sprintf("%.2f", safe_num(table_rows$lisa_support_score))
    names(table_rows) <- c("Gene", "log2FC", "DE FDR", "DE score", "LISA support", "Evidence", "sets", "LISA categories")
  } else {
    names(table_rows) <- c("Gene", "log2FC", "DE FDR", "Evidence", "sets", "LISA categories")
  }
  table_grob <- gridExtra::tableGrob(table_rows, rows = NULL, theme = gridExtra::ttheme_minimal(base_size = 7))
  p_table <- patchwork::wrap_elements(full = table_grob, clip = FALSE)

  p_panels <- patchwork::wrap_plots(
    list(p_lfc, p_score, p_marks, p_recur),
    nrow = 1,
    widths = c(0.42, 0.60, 0.28, 0.38),
    guides = "keep"
  )

  # Reserve a device-independent white footer. Without an explicit outer row,
  # Windows raster rounding can place the table on the final pixel and make a
  # complete-looking card indistinguishable from a clipped one.
  p_header / p_panels / p_table / patchwork::plot_spacer() +
    patchwork::plot_layout(heights = c(0.28, 1.88, 1.44, 0.03))
}

write_cards <- function(evidence, category_ids, paths, cfg) {
  dir.create(paths$cards_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(paths$cards_dir, "*_category_gene_card.*"))
  index_path <- file.path(paths$cards_dir, paste0(paths$prefix, "_category_gene_cards_index.tsv"))
  if (length(category_ids) == 0) {
    empty_index <- data.frame(
      analysis_id = character(),
      universe = character(),
      card_rank = integer(),
      category_id = character(),
      category_display_name = character(),
      macrogroup_id = character(),
      macrogroup_name = character(),
      n_supporting_genes = integer(),
      n_genes_shown = integer(),
      max_gene_contribution_score = numeric(),
      png_path = character(),
      pdf_path = character(),
      stringsAsFactors = FALSE
    )
    write_tsv(empty_index, index_path)
    return(empty_index)
  }
  card_width <- 8.9
  card_height <- max(8.2, 5.3 + (cfg$top_genes * 0.34))
  index <- vector("list", length(category_ids))
  for (i in seq_along(category_ids)) {
    category_id <- category_ids[[i]]
    card <- build_card_plot(evidence, category_id, cfg$top_genes)
    png_name <- sprintf("%s_category_gene_card.png", safe_file_component(category_id))
    png_path <- file.path(paths$cards_dir, png_name)
    pdf_path <- sub("[.]png$", ".pdf", png_path)
    svg_path <- sub("[.]png$", ".svg", png_path)
    if ("png" %in% cfg$formats) ggplot2::ggsave(png_path, card, width = card_width, height = card_height, dpi = 220, bg = "white", limitsize = FALSE)
    if ("pdf" %in% cfg$formats) ggplot2::ggsave(pdf_path, card, width = card_width, height = card_height, bg = "white", limitsize = FALSE)
    if ("svg" %in% cfg$formats) ggplot2::ggsave(svg_path, card, width = card_width, height = card_height, bg = "white", limitsize = FALSE)
    category_rows <- evidence[evidence$category_id == category_id, , drop = FALSE]
    category_rows <- category_rows[order(-safe_num(category_rows$gene_contribution_score),
      safe_num(category_rows$padj), -abs(safe_num(category_rows$log2FC)), category_rows$symbol), , drop = FALSE]
    category_rows$figure_id <- paste0("gene_card__", tools::file_path_sans_ext(png_name))
    category_rows$figure_type <- "gene_card"
    category_rows$source_row_order <- seq_len(nrow(category_rows))
    category_rows$category_color <- unique(category_rows$color)[1]
    category_rows$selected_for_plot <- seq_len(nrow(category_rows)) <= if (cfg$top_genes > 0) cfg$top_genes else nrow(category_rows)
    category_rows$highlighted <- category_rows$selected_for_plot
    category_rows$labelled <- category_rows$selected_for_plot
    category_rows$plotted_order <- ifelse(category_rows$selected_for_plot, seq_len(nrow(category_rows)), NA_integer_)
    category_rows$plot_log2FC <- safe_num(category_rows$log2FC)
    category_rows$plot_minus_log10_fdr <- -log10(pmax(safe_num(category_rows$padj), 1e-300))
    category_rows$selection_top_genes <- cfg$top_genes
    source_path <- paste0(tools::file_path_sans_ext(png_path), "_source.tsv")
    recipe_path <- paste0(tools::file_path_sans_ext(png_path), "_recipe.R")
    write_tsv(category_rows, source_path)
    renderer <- file.path(script_dir, "reproduce_lisa_figure.R")
    if (!file.exists(renderer)) {
      stop("LISA-FIGURE-SOURCE-002 cannot install GeneCard reproduction recipe.", call. = FALSE)
    }
    lisaR:::lisa_copy_verified_figure_recipe(
      renderer, recipe_path, cfg$lisa_internal_renderer_sha256,
      run_root = cfg$project_dir
    )
    n_genes_total <- length(unique(category_rows$symbol))
    n_genes_visible <- ifelse(cfg$top_genes > 0, min(cfg$top_genes, n_genes_total), n_genes_total)
    index[[i]] <- data.frame(
      analysis_id = cfg$analysis_id,
      universe = cfg$universe,
      card_rank = i,
      category_id = category_id,
      category_display_name = unique(category_rows$category_display_name)[1],
      macrogroup_id = unique(category_rows$macrogroup_id)[1],
      macrogroup_name = unique(category_rows$macrogroup_name)[1],
      n_supporting_genes = n_genes_total,
      n_genes_shown = n_genes_visible,
      n_genes_visible = n_genes_visible,
      n_genes_total = n_genes_total,
      gene_visibility = sprintf("%s/%s", n_genes_visible, n_genes_total),
      max_gene_contribution_score = max(safe_num(category_rows$gene_contribution_score), na.rm = TRUE),
      png_path = png_path,
      pdf_path = pdf_path,
      svg_path = svg_path,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, index)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE) ||
      !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Required packages: ggplot2, patchwork, gridExtra", call. = FALSE)
  }

  paths <- gene_level_paths(cfg$project_dir, cfg$analysis_id, cfg$universe)
  if (!file.exists(paths$evidence)) stop(sprintf("Missing gene-level evidence table: %s", paths$evidence), call. = FALSE)
  if (!file.exists(paths$summary)) stop(sprintf("Missing category summary table: %s", paths$summary), call. = FALSE)

  evidence <- read_tsv(paths$evidence)
  summary <- read_tsv(paths$summary)
  category_ids <- select_categories(evidence, summary, cfg)

  index <- write_cards(evidence, category_ids, paths, cfg)
  index_path <- file.path(paths$cards_dir, paste0(paths$prefix, "_category_gene_cards_index.tsv"))
  write_tsv(index, index_path)
  if (nrow(index) == 0) {
    message("No GSEA/lollipop categories with n_genesets > 0; wrote empty category gene-card index.")
  } else {
    print(index[, c("card_rank", "category_id", "n_supporting_genes", "n_genes_shown")])
  }
}

if (sys.nframe() == 0L) main()
