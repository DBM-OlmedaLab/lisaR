#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

default_project_dir <- NA_character_
script_file <- sub(
  "^--file=", "",
  commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1]
)
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "lisa_plot_metadata.R"))

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
    label_genes = 16,
    categories = NA_character_,
    formats = "png",
    source_only = "false",
    de_padj_cutoff = 0.05,
    lfc_cutoff = 0.58,
    lisa_internal_renderer_sha256 = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    if (i == length(args)) stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    val <- args[[i + 1]]
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(out)) stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    if (key %in% c("max_categories", "label_genes")) {
      out[[key]] <- as.integer(val)
    } else if (key %in% c("de_padj_cutoff", "lfc_cutoff")) {
      out[[key]] <- as.numeric(val)
    } else {
      out[[key]] <- val
    }
    i <- i + 2
  }
  if (is.na(out$project_dir) || out$project_dir == "") stop("Required argument: --project-dir", call. = FALSE)
  if (is.na(out$analysis_id) || out$analysis_id == "") stop("Required argument: --analysis-id", call. = FALSE)
  if (!out$source_only %in% c("true", "false")) stop("--source-only must be true or false", call. = FALSE)
  out$formats <- unique(trimws(strsplit(tolower(out$formats), ",", fixed = TRUE)[[1]]))
  if (any(!out$formats %in% c("png", "svg", "pdf"))) stop("--formats must contain png, svg and/or pdf", call. = FALSE)
  out
}

paths_for <- function(project_dir, analysis_id, universe) {
  collection_dir <- file.path(project_dir, "outputs", "single_de", analysis_id, paste0("collection_", universe))
  gene_dir <- file.path(project_dir, "outputs", "gene_level", "single_de", analysis_id, paste0("collection_", universe))
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
    standardized_de = file.path(collection_dir, "inputs", paste0(analysis_id, "_standardized_DE.tsv")),
    gsea_summary = if (length(summary_candidates) == 1) summary_candidates[[1]] else "",
    evidence = file.path(gene_dir, paste0(prefix, "_gene_category_contributions.tsv")),
    summary = file.path(gene_dir, paste0(prefix, "_category_gene_support_summary.tsv")),
    volcano_dir = file.path(gene_dir, "category_volcano_overlays"),
    prefix = prefix
  )
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_file_component <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

select_categories <- function(evidence, summary, cfg) {
  requested <- character()
  if (!is.na(cfg$categories) && cfg$categories != "") {
    requested <- trimws(strsplit(cfg$categories, ",", fixed = TRUE)[[1]])
  }
  if (length(requested) > 0) {
    found <- requested[requested %in% unique(evidence$category_id)]
    missing <- setdiff(requested, found)
    if (length(missing) > 0) warning(sprintf("Requested categories not found and skipped: %s", paste(missing, collapse = ", ")), call. = FALSE)
    return(found)
  }

  paths <- paths_for(cfg$project_dir, cfg$analysis_id, cfg$universe)
  if (file.exists(paths$gsea_summary)) {
    gsea <- read_tsv(paths$gsea_summary)
    if ("n_genesets" %in% names(gsea)) {
      gsea <- gsea[safe_num(gsea$n_genesets) > 0, , drop = FALSE]
    }
    if (nrow(gsea) == 0) {
      return(character())
    }
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

  summary <- summary[order(-safe_num(summary$max_gene_contribution_score),
                           -safe_num(summary$n_supporting_genes),
                           summary$category_id), , drop = FALSE]
  if (cfg$max_categories > 0) {
    return(head(summary$category_id, cfg$max_categories))
  }
  summary$category_id
}

prepare_de <- function(de, cfg) {
  de$symbol <- toupper(as.character(de$symbol))
  de$log2FC <- safe_num(de$log2FoldChange)
  de$padj_num <- safe_num(de$padj)
  de$pvalue_num <- safe_num(de$pvalue)
  de$plot_p <- ifelse(!is.na(de$padj_num), de$padj_num, de$pvalue_num)
  de$plot_p <- pmax(de$plot_p, 1e-300)
  de$neg_log10_fdr <- -log10(de$plot_p)
  de$de_class <- ifelse(
    !is.na(de$padj_num) & de$padj_num <= cfg$de_padj_cutoff & de$log2FC >= cfg$lfc_cutoff, "DE up",
    ifelse(!is.na(de$padj_num) & de$padj_num <= cfg$de_padj_cutoff & de$log2FC <= -cfg$lfc_cutoff, "DE down", "not DE")
  )
  de
}

build_volcano <- function(de, evidence, category_id, cfg) {
  category_rows <- evidence[evidence$category_id == category_id, , drop = FALSE]
  category_rows <- category_rows[order(-safe_num(category_rows$gene_contribution_score),
                                       safe_num(category_rows$padj),
                                       category_rows$symbol), , drop = FALSE]
  category_rows <- category_rows[!duplicated(category_rows$symbol), , drop = FALSE]
  highlight <- unique(category_rows$symbol)
  plot_df <- de
  plot_df$category_support <- plot_df$symbol %in% highlight
  label_symbols <- if (cfg$label_genes > 0) head(category_rows$symbol, cfg$label_genes) else category_rows$symbol
  label_df <- plot_df[plot_df$symbol %in% label_symbols, , drop = FALSE]

  title <- unique(category_rows$category_display_name)[1]
  macrogroup <- unique(category_rows$macrogroup_name)[1]
  category_color <- unique(category_rows$color)[1]
  if (is.na(category_color) || category_color == "") category_color <- "#009E73"

  ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FC, y = neg_log10_fdr)) +
    ggplot2::geom_point(ggplot2::aes(color = de_class), alpha = 0.15, size = 0.55, na.rm = TRUE) +
    ggplot2::geom_vline(xintercept = c(-cfg$lfc_cutoff, cfg$lfc_cutoff), linetype = "dashed", color = "grey70", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(cfg$de_padj_cutoff), linetype = "dashed", color = "grey70", linewidth = 0.3) +
    ggplot2::geom_point(
      data = plot_df[plot_df$category_support, , drop = FALSE],
      color = "grey10",
      fill = category_color,
      shape = 21,
      size = 1.55,
      alpha = 0.9,
      stroke = 0.22,
      na.rm = TRUE
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      ggplot2::aes(label = symbol),
      size = 3.25,
      max.overlaps = Inf,
      box.padding = 0.42,
      point.padding = 0.16,
      seed = 1,
      force = 2,
      force_pull = 0.15,
      min.segment.length = 0,
      segment.size = 0.22,
      segment.color = "grey55",
      na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(values = c(`DE down` = "#4A72B8", `DE up` = "#C44E52", `not DE` = "#B8B8B8")) +
    ggplot2::labs(
      title = title,
      subtitle = lisa_plot_subtitle(
        cfg$plot_metadata,
        prefix = macrogroup,
        suffix = sprintf(
          "%s | %s supporting genes highlighted", cfg$universe, length(highlight)
        )
      ),
      x = "log2FC",
      y = "-log10(FDR or p-value)",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey35"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "bottom"
    )
}

write_volcanoes <- function(de, evidence, category_ids, paths, cfg) {
  dir.create(paths$volcano_dir, recursive = TRUE, showWarnings = FALSE)
  source_only <- identical(cfg$source_only, "true")
  if (!source_only) unlink(file.path(paths$volcano_dir, "*_category_volcano_overlay.*"))
  index_path <- file.path(paths$volcano_dir, paste0(paths$prefix, "_category_volcano_overlays_index.tsv"))
  if (length(category_ids) == 0) {
    empty_index <- data.frame(
      analysis_id = character(),
      universe = character(),
      volcano_rank = integer(),
      category_id = character(),
      category_display_name = character(),
      macrogroup_id = character(),
      macrogroup_name = character(),
      n_supporting_genes = integer(),
      n_labelled_genes = integer(),
      n_genes_visible = integer(),
      n_genes_total = integer(),
      gene_visibility = character(),
      png_path = character(),
      pdf_path = character(),
      stringsAsFactors = FALSE
    )
    write_tsv(empty_index, index_path)
    return(empty_index)
  }
  index <- vector("list", length(category_ids))
  for (i in seq_along(category_ids)) {
    category_id <- category_ids[[i]]
    if (!source_only) plot <- build_volcano(de, evidence, category_id, cfg)
    png_name <- sprintf("%s_category_volcano_overlay.png", safe_file_component(category_id))
    png_path <- file.path(paths$volcano_dir, png_name)
    pdf_path <- sub("[.]png$", ".pdf", png_path)
    svg_path <- sub("[.]png$", ".svg", png_path)
    if (source_only && !all(file.exists(c(if ("png" %in% cfg$formats) png_path,
        if ("pdf" %in% cfg$formats) pdf_path, if ("svg" %in% cfg$formats) svg_path))))
      stop("Source-only repair requires all existing requested figure formats.", call. = FALSE)
    if (!source_only && "png" %in% cfg$formats) ggplot2::ggsave(png_path, plot, width = 8.5, height = 7.2, dpi = 240, bg = "white", limitsize = FALSE)
    if (!source_only && "pdf" %in% cfg$formats) ggplot2::ggsave(pdf_path, plot, width = 8.5, height = 7.2, bg = "white", limitsize = FALSE)
    if (!source_only && "svg" %in% cfg$formats) ggplot2::ggsave(svg_path, plot, width = 8.5, height = 7.2, bg = "white", limitsize = FALSE)
    rows <- evidence[evidence$category_id == category_id, , drop = FALSE]
    rows <- rows[order(-safe_num(rows$gene_contribution_score), safe_num(rows$padj), rows$symbol), , drop = FALSE]
    rows <- rows[!duplicated(rows$symbol), , drop = FALSE]
    labels <- if (cfg$label_genes > 0) head(rows$symbol, cfg$label_genes) else rows$symbol
    source_rows <- de
    source_rows$figure_id <- paste0("volcano_overlay__", tools::file_path_sans_ext(png_name))
    source_rows$figure_type <- "volcano_overlay"
    source_rows$native_volcano_style <- TRUE
    source_rows$label_layout_seed <- if (source_only) NA_integer_ else 1L
    source_rows$label_layout_provenance <- if (source_only)
      "Original native label seed was not recorded; replay positions may differ."
      else "Native labels rendered with fixed seed 1."
    source_rows$category_id <- category_id
    source_rows$category_display_name <- unique(rows$category_display_name)[1]
    source_rows$macrogroup_name <- unique(rows$macrogroup_name)[1]
    source_rows$category_color <- unique(rows$color)[1]
    source_rows$plot_subtitle <- lisa_plot_subtitle(
      cfg$plot_metadata,
      prefix = unique(rows$macrogroup_name)[1],
      suffix = sprintf("%s | %s supporting genes highlighted", cfg$universe, length(unique(rows$symbol)))
    )
    source_rows$selected_for_plot <- TRUE
    source_rows$highlighted <- source_rows$symbol %in% rows$symbol
    source_rows$labelled <- source_rows$symbol %in% labels
    source_rows$source_row_order <- seq_len(nrow(source_rows))
    source_rows$plot_x_log2FC <- source_rows$log2FC
    source_rows$plot_y_neg_log10_fdr <- source_rows$neg_log10_fdr
    source_rows$threshold_abs_log2FC <- cfg$lfc_cutoff
    source_rows$threshold_de_fdr <- cfg$de_padj_cutoff
    source_rows$label_rank <- match(source_rows$symbol, labels)
    source_path <- paste0(tools::file_path_sans_ext(png_path), "_source.tsv")
    recipe_path <- paste0(tools::file_path_sans_ext(png_path), "_recipe.R")
    lisaR:::lisa_write_figure_source_tsv(source_rows, source_path)
    renderer <- file.path(script_dir, "reproduce_lisa_figure.R")
    if (!file.exists(renderer)) {
      stop("LISA-FIGURE-SOURCE-003 cannot install volcano reproduction recipe.", call. = FALSE)
    }
    lisaR:::lisa_copy_verified_figure_recipe(
      renderer, recipe_path, cfg$lisa_internal_renderer_sha256,
      run_root = cfg$project_dir
    )
    n_genes_total <- length(unique(rows$symbol))
    n_genes_visible <- ifelse(cfg$label_genes > 0, min(cfg$label_genes, n_genes_total), n_genes_total)
    index[[i]] <- data.frame(
      analysis_id = cfg$analysis_id,
      universe = cfg$universe,
      volcano_rank = i,
      category_id = category_id,
      category_display_name = unique(rows$category_display_name)[1],
      macrogroup_id = unique(rows$macrogroup_id)[1],
      macrogroup_name = unique(rows$macrogroup_name)[1],
      n_supporting_genes = n_genes_total,
      n_labelled_genes = n_genes_visible,
      n_genes_visible = n_genes_visible,
      n_genes_total = n_genes_total,
      gene_visibility = sprintf("%s/%s", n_genes_visible, n_genes_total),
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
  cfg$plot_metadata <- lisa_plot_metadata(cfg$project_dir, cfg$analysis_id)
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Required packages: ggplot2, ggrepel", call. = FALSE)
  }
  paths <- paths_for(cfg$project_dir, cfg$analysis_id, cfg$universe)
  if (!file.exists(paths$standardized_de)) stop(sprintf("Missing standardized DE table: %s", paths$standardized_de), call. = FALSE)
  if (!file.exists(paths$evidence)) stop(sprintf("Missing gene-level evidence table: %s", paths$evidence), call. = FALSE)
  if (!file.exists(paths$summary)) stop(sprintf("Missing category summary table: %s", paths$summary), call. = FALSE)
  de <- prepare_de(read_tsv(paths$standardized_de), cfg)
  evidence <- read_tsv(paths$evidence)
  summary <- read_tsv(paths$summary)
  category_ids <- select_categories(evidence, summary, cfg)
  index <- write_volcanoes(de, evidence, category_ids, paths, cfg)
  index_path <- file.path(paths$volcano_dir, paste0(paths$prefix, "_category_volcano_overlays_index.tsv"))
  write_tsv(index, index_path)
  if (nrow(index) == 0) {
    message("No GSEA/lollipop categories with n_genesets > 0; wrote empty category volcano index.")
  } else {
    print(index[, c("volcano_rank", "category_id", "n_supporting_genes", "n_labelled_genes", "gene_visibility")])
  }
}

if (sys.nframe() == 0L) main()
