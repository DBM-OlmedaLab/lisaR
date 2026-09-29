#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

default_project_dir <- NA_character_

read_tsv <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing TSV: %s", path), call. = FALSE)
  if (file.info(path)$size == 0L) return(data.frame())
  first_lines <- readLines(path, n = 2L, warn = FALSE)
  if (length(first_lines) == 0L || !any(nzchar(trimws(first_lines)))) return(data.frame())
  utils::read.delim(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE)
}

write_tsv <- function(df, path) {
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

parse_args <- function(args) {
  out <- list(
    project_dir = default_project_dir,
    contrast_id = NA_character_,
    universe = "GOBP-C2",
    top_pathways = 0,
    top_genes = 45,
    lfc_cap = 1.5
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    if (i == length(args)) stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    val <- args[[i + 1]]
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(out)) stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    if (key %in% c("top_pathways", "top_genes")) {
      out[[key]] <- as.integer(val)
    } else if (key %in% c("lfc_cap")) {
      out[[key]] <- as.numeric(val)
    } else {
      out[[key]] <- val
    }
    i <- i + 2
  }
  if (is.na(out$contrast_id) || out$contrast_id == "") stop("Required argument: --contrast-id", call. = FALSE)
  out
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

first_nonempty <- function(x) {
  x <- safe_chr(x)
  x <- x[x != ""]
  if (length(x) == 0) "" else x[[1]]
}

min_or_na <- function(x) {
  x <- safe_num(x)
  out <- suppressWarnings(min(x, na.rm = TRUE))
  if (!is.finite(out)) NA_real_ else out
}

max_or_na <- function(x) {
  x <- safe_num(x)
  out <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(out)) NA_real_ else out
}

collapse_unique <- function(x, sep = ";") {
  x <- sort(unique(safe_chr(x)))
  x <- x[x != ""]
  paste(x, collapse = sep)
}

short_label <- function(x, max_chars = 70) {
  x <- safe_chr(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

pathway_label <- function(x) {
  x <- gsub("^KEGG_MEDICUS_(REFERENCE|VARIANT)_", "KEGG_MEDICUS_", safe_chr(x))
  x <- gsub("^KEGG_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

contrast_row <- function(project_dir, contrast_id) {
  idx <- read_tsv(file.path(project_dir, "config", "contrast_index.tsv"))
  row <- idx[idx$contrast_id == contrast_id | paste(idx$contrast_id, idx$output_id, sep = "_") == contrast_id, , drop = FALSE]
  if (nrow(row) != 1) stop(sprintf("Contrast id did not resolve to one row: %s", contrast_id), call. = FALSE)
  row
}

paths_for <- function(project_dir, row, universe) {
  contrast_name <- paste(row$contrast_id, row$output_id, sep = "_")
  out_dir <- file.path(project_dir, "outputs", "gene_level", "category_contrasts", contrast_name, paste0("collection_", universe))
  prefix <- paste(row$output_id, universe, "contrast_gene_level", sep = "_")
  list(
    contrast_name = contrast_name,
    out_dir = out_dir,
    paired = file.path(out_dir, paste0(prefix, "_paired_gene_evidence.tsv")),
    summary = file.path(out_dir, paste0(prefix, "_contrast_category_gene_summary.tsv")),
    kegg_dir = file.path(out_dir, "contrast_kegg_map_layer"),
    prefix = prefix
  )
}

split_source_pathways <- function(x) {
  x <- safe_chr(x)
  if (length(x) == 0 || x[[1]] == "") return(character())
  vals <- unique(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
  vals <- trimws(vals)
  vals[grepl("^KEGG_", vals) & !grepl("^KEGG_MEDICUS", vals)]
}

extract_side_rows <- function(paired, side) {
  path_col <- paste0("source_pathways_", side)
  lfc_col <- paste0("log2FC_", side)
  padj_col <- paste0("padj_", side)
  present_col <- paste0("present_", side)
  score_col <- paste0("gene_contribution_score_", side)
  if (!all(c(path_col, lfc_col, padj_col, present_col, score_col) %in% names(paired))) return(data.frame())
  rows <- vector("list", nrow(paired))
  for (i in seq_len(nrow(paired))) {
    if (!isTRUE(paired[[present_col]][[i]])) next
    pathways <- split_source_pathways(paired[[path_col]][[i]])
    if (length(pathways) == 0) next
    rows[[i]] <- data.frame(
      side = side,
      pathway = pathways,
      pathway_display_label = pathway_label(pathways),
      category_id = paired$category_id[[i]],
      category_display_name = paired$category_display_name[[i]],
      macrogroup_id = paired$macrogroup_id[[i]],
      macrogroup_name = paired$macrogroup_name[[i]],
      direction_class = paired$direction_class[[i]],
      contrast_gene_class = paired$contrast_gene_class[[i]],
      symbol = paired$symbol[[i]],
      log2FC = safe_num(paired[[lfc_col]][[i]]),
      padj = safe_num(paired[[padj_col]][[i]]),
      gene_contribution_score = safe_num(paired[[score_col]][[i]]),
      abs_delta_mean_NES = safe_num(paired$abs_delta_mean_NES[[i]]),
      mean_abs_NES = safe_num(paired$mean_abs_NES[[i]]),
      stringsAsFactors = FALSE
    )
  }
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0) data.frame() else do.call(rbind, rows)
}

build_kegg_rows <- function(paired) {
  rows <- rbind(extract_side_rows(paired, "A"), extract_side_rows(paired, "B"))
  if (nrow(rows) == 0) return(data.frame())
  rows$symbol <- toupper(safe_chr(rows$symbol))
  rows <- rows[rows$pathway != "" & rows$symbol != "", , drop = FALSE]
  rows
}

summarize_pathways <- function(kegg_rows, cfg, row) {
  if (nrow(kegg_rows) == 0) return(data.frame())
  split_rows <- split(kegg_rows, paste(kegg_rows$pathway, kegg_rows$category_id, sep = "\r"))
  out <- lapply(split_rows, function(df) {
    a <- df[df$side == "A", , drop = FALSE]
    b <- df[df$side == "B", , drop = FALSE]
    shared <- intersect(unique(a$symbol), unique(b$symbol))
    shared_classes <- unique(df$contrast_gene_class[df$symbol %in% shared])
    data.frame(
      contrast_id = row$contrast_id,
      output_id = row$output_id,
      universe = cfg$universe,
      pathway = first_nonempty(df$pathway),
      pathway_display_label = first_nonempty(df$pathway_display_label),
      category_id = first_nonempty(df$category_id),
      category_display_name = first_nonempty(df$category_display_name),
      macrogroup_id = first_nonempty(df$macrogroup_id),
      macrogroup_name = first_nonempty(df$macrogroup_name),
      direction_class = first_nonempty(df$direction_class),
      n_genes_A = length(unique(a$symbol)),
      n_genes_B = length(unique(b$symbol)),
      n_shared_genes = length(shared),
      n_shared_opposite_genes = sum(shared_classes == "shared_opposite_direction"),
      n_shared_same_genes = sum(shared_classes == "shared_same_direction"),
      min_padj_A = min_or_na(a$padj),
      min_padj_B = min_or_na(b$padj),
      max_abs_log2FC_A = max_or_na(abs(a$log2FC)),
      max_abs_log2FC_B = max_or_na(abs(b$log2FC)),
      max_gene_contribution_score = max_or_na(df$gene_contribution_score),
      abs_delta_mean_NES = max_or_na(df$abs_delta_mean_NES),
      mean_abs_NES = max_or_na(df$mean_abs_NES),
      genes_A = collapse_unique(a$symbol),
      genes_B = collapse_unique(b$symbol),
      shared_genes = paste(sort(shared), collapse = ";"),
      interpretation_note = "KEGG-derived source pathway inside paired contrast gene-level evidence; not causal evidence.",
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out$pathway_gene_total <- out$n_genes_A + out$n_genes_B
  out$contrast_kegg_score <- safe_num(out$abs_delta_mean_NES) +
    0.20 * safe_num(out$n_shared_opposite_genes) +
    0.03 * safe_num(out$n_shared_genes) +
    0.01 * safe_num(out$pathway_gene_total) +
    0.05 * safe_num(out$max_gene_contribution_score)
  out <- out[order(
    -safe_num(out$contrast_kegg_score),
    -safe_num(out$n_shared_opposite_genes),
    -safe_num(out$n_shared_genes),
    -safe_num(out$pathway_gene_total),
    out$pathway,
    out$category_id
  ), , drop = FALSE]
  rownames(out) <- NULL
  out$contrast_kegg_rank <- seq_len(nrow(out))
  out[, c("contrast_kegg_rank", setdiff(names(out), "contrast_kegg_rank")), drop = FALSE]
}

limit_rows <- function(df, n) {
  if (n <= 0 || nrow(df) <= n) return(df)
  utils::head(df, n)
}

build_gene_matrix <- function(pathway_summary, kegg_rows, cfg) {
  selected <- limit_rows(pathway_summary, cfg$top_pathways)
  empty <- list(matrix = data.frame(), n_genes_visible = 0L, n_genes_total = 0L)
  if (nrow(selected) == 0 || nrow(kegg_rows) == 0) return(empty)
  keys <- paste(selected$pathway, selected$category_id, sep = "\r")
  kegg_rows$key <- paste(kegg_rows$pathway, kegg_rows$category_id, sep = "\r")
  rows <- kegg_rows[kegg_rows$key %in% keys, , drop = FALSE]
  if (nrow(rows) == 0) return(data.frame())
  gene_priority <- stats::aggregate(
    cbind(pathway_count = rep(1, nrow(rows)), abs_lfc = abs(rows$log2FC)) ~ symbol,
    data = rows,
    FUN = function(x) sum(!is.na(x))
  )
  max_abs <- stats::aggregate(abs(rows$log2FC), list(symbol = rows$symbol), max, na.rm = TRUE)
  names(max_abs)[2] <- "max_abs_lfc"
  gene_priority <- merge(gene_priority[, c("symbol", "pathway_count")], max_abs, by = "symbol", all.x = TRUE)
  gene_priority <- gene_priority[order(-safe_num(gene_priority$pathway_count), -safe_num(gene_priority$max_abs_lfc), gene_priority$symbol), , drop = FALSE]
  total_genes <- length(unique(gene_priority$symbol))
  keep_genes <- if (cfg$top_genes > 0) utils::head(gene_priority$symbol, cfg$top_genes) else gene_priority$symbol
  rows <- rows[rows$symbol %in% keep_genes, , drop = FALSE]
  rows$pathway_label <- short_label(paste0(rows$pathway_display_label, " [", rows$category_display_name, "]"), 68)
  list(
    matrix = rows,
    n_genes_visible = length(unique(rows$symbol)),
    n_genes_total = total_genes
  )
}

build_pathway_plot <- function(pathway_summary, cfg, row, gene_stats) {
  plot_df <- limit_rows(pathway_summary, cfg$top_pathways)
  plot_df$label <- short_label(paste0(plot_df$pathway_display_label, " [", plot_df$category_display_name, "]"), 78)
  plot_df$label <- factor(plot_df$label, levels = rev(plot_df$label))
  plot_df$evidence_class <- ifelse(
    plot_df$n_shared_opposite_genes > 0,
    "shared opposite genes",
    ifelse(plot_df$n_shared_genes > 0, "shared same/mixed genes", "one-sided genes")
  )
  ggplot2::ggplot(plot_df, ggplot2::aes(x = contrast_kegg_score, y = label)) +
    ggplot2::geom_point(
      ggplot2::aes(size = pathway_gene_total, fill = evidence_class),
      shape = 21,
      color = "grey20",
      stroke = 0.35,
      alpha = 0.92
    ) +
    ggplot2::scale_size_continuous(range = c(1.7, 5.5), name = "A+B gene rows") +
    ggplot2::scale_fill_manual(values = c(
      "shared opposite genes" = "#0AA36E",
      "shared same/mixed genes" = "#4C78A8",
      "one-sided genes" = "#C7C7C7"
    )) +
    ggplot2::labs(
      title = "Contrast KEGG-derived pathway layer",
      subtitle = sprintf(
        "%s | A = %s; B = %s | pathways shown %s/%s | heatmap genes %s/%s",
        row$output_id,
        row$label_a,
        row$label_b,
        nrow(plot_df),
        nrow(pathway_summary),
        gene_stats$n_genes_visible,
        gene_stats$n_genes_total
      ),
      x = "contrast KEGG priority score",
      y = NULL,
      fill = "Gene evidence"
    ) +
    ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey35"),
      axis.text.y = ggplot2::element_text(size = 7.2),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "right"
    )
}

build_gene_heatmap <- function(gene_matrix, cfg, row, gene_stats) {
  if (nrow(gene_matrix) == 0) return(NULL)
  gene_matrix$pathway_label <- factor(gene_matrix$pathway_label, levels = rev(unique(gene_matrix$pathway_label)))
  gene_matrix$symbol_side <- paste(gene_matrix$symbol, gene_matrix$side, sep = "\n")
  symbol_order <- unique(gene_matrix$symbol)
  gene_matrix$symbol_side <- factor(
    gene_matrix$symbol_side,
    levels = as.vector(rbind(paste(symbol_order, "A", sep = "\n"), paste(symbol_order, "B", sep = "\n")))
  )
  gene_matrix$plot_lfc <- pmax(pmin(safe_num(gene_matrix$log2FC), cfg$lfc_cap), -cfg$lfc_cap)
  ggplot2::ggplot(gene_matrix, ggplot2::aes(x = symbol_side, y = pathway_label, fill = plot_lfc)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.22) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-cfg$lfc_cap, cfg$lfc_cap), na.value = "grey88", name = "log2FC"
    ) +
    ggplot2::labs(
      title = "Genes mapped to selected contrast KEGG-derived pathways",
      subtitle = sprintf("%s | A = %s; B = %s | genes shown %s/%s", row$output_id, row$label_a, row$label_b, gene_stats$n_genes_visible, gene_stats$n_genes_total),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 8) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 8, color = "grey35"),
      axis.text.y = ggplot2::element_text(size = 7),
      axis.text.x = ggplot2::element_text(angle = 60, hjust = 1, vjust = 1, size = 5.7),
      panel.grid = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

write_readme <- function(path, cfg, row) {
  txt <- c(
    "Contrast KEGG-derived pathway layer",
    "",
    sprintf("Contrast: %s_%s", row$contrast_id, row$output_id),
    sprintf("Universe: %s", cfg$universe),
    sprintf("A: %s (%s)", row$contrast_a, row$label_a),
    sprintf("B: %s (%s)", row$contrast_b, row$label_b),
    "",
    "Purpose",
    "This output extracts canonical KEGG source pathways already present in the paired contrast gene-level evidence.",
    "It summarizes pathway/category/gene support across A and B and provides compact pathway and gene heatmaps.",
    "",
    "Important caution",
    "This is not a rerun of pathway enrichment and does not add causal evidence.",
    "It is an auditable contrast summary over KEGG terms that were already present in LISA source pathways.",
    "",
    "Files",
    "- *_contrast_kegg_pathway_layer.tsv: all summarized KEGG pathway-category rows.",
    "- *_top_contrast_kegg_pathways.tsv: top selected rows used for compact review.",
    "- *_contrast_kegg_gene_matrix.tsv: gene/pathway/side matrix used for the heatmap.",
    "- *_top_contrast_kegg_pathways.png/pdf: pathway-level priority plot.",
    "- *_contrast_kegg_gene_heatmap.png/pdf: gene-by-pathway A/B heatmap.",
    "",
    "This is post-processing only; no DE, GSEA, ORA, or LISA step is rerun."
  )
  writeLines(txt, con = path, useBytes = TRUE)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Required package: ggplot2", call. = FALSE)

  row <- contrast_row(cfg$project_dir, cfg$contrast_id)
  paths <- paths_for(cfg$project_dir, row, cfg$universe)
  missing <- c(paths$paired, paths$summary)[!file.exists(c(paths$paired, paths$summary))]
  if (length(missing) > 0) stop(sprintf("Missing required input files:\n%s", paste(missing, collapse = "\n")), call. = FALSE)
  dir.create(paths$kegg_dir, recursive = TRUE, showWarnings = FALSE)

  paired <- read_tsv(paths$paired)
  kegg_rows <- build_kegg_rows(paired)
  if (nrow(kegg_rows) == 0) stop("No KEGG source pathways found in paired contrast evidence.", call. = FALSE)
  pathway_summary <- summarize_pathways(kegg_rows, cfg, row)
  gene_result <- build_gene_matrix(pathway_summary, kegg_rows, cfg)
  gene_matrix <- gene_result$matrix
  visible_pathways <- limit_rows(pathway_summary, cfg$top_pathways)

  summary_path <- file.path(paths$kegg_dir, paste0(paths$prefix, "_contrast_kegg_pathway_layer.tsv"))
  top_path <- file.path(paths$kegg_dir, paste0(paths$prefix, "_top_contrast_kegg_pathways.tsv"))
  matrix_path <- file.path(paths$kegg_dir, paste0(paths$prefix, "_contrast_kegg_gene_matrix.tsv"))
  write_tsv(pathway_summary, summary_path)
  write_tsv(visible_pathways, top_path)
  write_tsv(gene_matrix, matrix_path)

  pathway_plot <- build_pathway_plot(pathway_summary, cfg, row, gene_result)
  pathway_png <- file.path(paths$kegg_dir, paste0(paths$prefix, "_top_contrast_kegg_pathways.png"))
  pathway_pdf <- file.path(paths$kegg_dir, paste0(paths$prefix, "_top_contrast_kegg_pathways.pdf"))
  pathway_height <- max(7.8, 2.6 + nrow(visible_pathways) * 0.18)
  ggplot2::ggsave(pathway_png, pathway_plot, width = 11.2, height = pathway_height, dpi = 240, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(pathway_pdf, pathway_plot, width = 11.2, height = pathway_height, bg = "white", limitsize = FALSE)

  heatmap_png <- ""
  heatmap_pdf <- ""
  if (nrow(gene_matrix) > 0) {
    heatmap <- build_gene_heatmap(gene_matrix, cfg, row, gene_result)
    heatmap_png <- file.path(paths$kegg_dir, paste0(paths$prefix, "_contrast_kegg_gene_heatmap.png"))
    heatmap_pdf <- file.path(paths$kegg_dir, paste0(paths$prefix, "_contrast_kegg_gene_heatmap.pdf"))
    ggplot2::ggsave(heatmap_png, heatmap, width = 12.5, height = 8.0, dpi = 240, bg = "white", limitsize = FALSE)
    ggplot2::ggsave(heatmap_pdf, heatmap, width = 12.5, height = 8.0, bg = "white", limitsize = FALSE)
  }

  readme_path <- file.path(paths$kegg_dir, "README_contrast_kegg_map_layer.txt")
  write_readme(readme_path, cfg, row)
  index <- data.frame(
    contrast_id = row$contrast_id,
    output_id = row$output_id,
    universe = cfg$universe,
    n_kegg_pathway_category_rows = nrow(pathway_summary),
    n_pathway_rows_visible = nrow(visible_pathways),
    n_pathway_rows_total = nrow(pathway_summary),
    n_gene_matrix_rows = nrow(gene_matrix),
    n_genes_visible = gene_result$n_genes_visible,
    n_genes_total = gene_result$n_genes_total,
    gene_visibility = sprintf("%s/%s", gene_result$n_genes_visible, gene_result$n_genes_total),
    summary_tsv = summary_path,
    top_kegg_tsv = top_path,
    gene_matrix_tsv = matrix_path,
    pathway_plot_png = pathway_png,
    pathway_plot_pdf = pathway_pdf,
    gene_heatmap_png = heatmap_png,
    gene_heatmap_pdf = heatmap_pdf,
    readme = readme_path,
    stringsAsFactors = FALSE
  )
  index_path <- file.path(paths$kegg_dir, paste0(paths$prefix, "_contrast_kegg_map_layer_index.tsv"))
  write_tsv(index, index_path)
  print(index[, c("contrast_id", "output_id", "universe", "n_kegg_pathway_category_rows", "n_pathway_rows_visible", "n_gene_matrix_rows", "gene_visibility")])
}

if (!identical(Sys.getenv("KEGG_MAP_LAYER_LIBRARY_ONLY"), "1")) {
  main()
}
