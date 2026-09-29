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

parse_args <- function(args) {
  out <- list(
    project_dir = default_project_dir,
    contrast_id = NA_character_,
    universe = "GOBP-C2",
    max_macrogroups = 0,
    max_categories_per_macrogroup = 0,
    max_genes = 45,
    lfc_cap = 1.5,
    macrogroups = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    if (i == length(args)) stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    val <- args[[i + 1]]
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(out)) stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    if (key %in% c("max_macrogroups", "max_categories_per_macrogroup", "max_genes")) {
      out[[key]] <- as.integer(val)
    } else if (key %in% c("lfc_cap")) {
      out[[key]] <- as.numeric(val)
    } else {
      out[[key]] <- val
    }
    i <- i + 2
  }
  if (is.na(out$project_dir) || out$project_dir == "") stop("Required argument: --project-dir", call. = FALSE)
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

safe_file_component <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

short_label <- function(x, max_chars = 30) {
  x <- safe_chr(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
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
    heatmap_dir = file.path(out_dir, "contrast_macrogroup_heatmaps"),
    prefix = prefix
  )
}

pick_macrogroups <- function(summary, paired, cfg) {
  requested <- character()
  if (!is.na(cfg$macrogroups) && cfg$macrogroups != "") {
    requested <- trimws(strsplit(cfg$macrogroups, ",", fixed = TRUE)[[1]])
  }
  available <- unique(safe_chr(paired$macrogroup_id))
  available <- available[available != ""]
  if (length(requested) > 0) {
    found <- requested[requested %in% available]
    missing <- setdiff(requested, found)
    if (length(missing) > 0) {
      warning(sprintf("Requested macrogroups not found and skipped: %s", paste(missing, collapse = ", ")), call. = FALSE)
    }
    return(found)
  }

  summary$macrogroup_id <- safe_chr(summary$macrogroup_id)
  summary$rank_score <- safe_num(summary$abs_delta_mean_NES) +
    0.06 * safe_num(summary$n_shared_opposite_direction) +
    0.01 * safe_num(summary$n_genes_total)
  macro_scores <- stats::aggregate(rank_score ~ macrogroup_id + macrogroup_name, data = summary, FUN = max)
  macro_counts <- stats::aggregate(n_genes_total ~ macrogroup_id, data = summary, FUN = sum)
  names(macro_counts)[2] <- "n_gene_category_rows"
  macro_scores <- merge(macro_scores, macro_counts, by = "macrogroup_id", all.x = TRUE, sort = FALSE)
  macro_scores <- macro_scores[macro_scores$macrogroup_id %in% available, , drop = FALSE]
  macro_scores <- macro_scores[order(-safe_num(macro_scores$rank_score), -safe_num(macro_scores$n_gene_category_rows), macro_scores$macrogroup_id), , drop = FALSE]
  if (cfg$max_macrogroups > 0) {
    return(head(macro_scores$macrogroup_id, cfg$max_macrogroups))
  }
  macro_scores$macrogroup_id
}

select_macrogroup_data <- function(paired, summary, macrogroup_id, cfg) {
  macro_summary <- summary[summary$macrogroup_id == macrogroup_id, , drop = FALSE]
  macro_summary <- macro_summary[order(
    -safe_num(macro_summary$abs_delta_mean_NES),
    -safe_num(macro_summary$n_shared_opposite_direction),
    -safe_num(macro_summary$n_genes_total),
    macro_summary$category_id
  ), , drop = FALSE]
  selected_categories <- macro_summary$category_id
  if (cfg$max_categories_per_macrogroup > 0) {
    selected_categories <- head(selected_categories, cfg$max_categories_per_macrogroup)
  }
  df <- paired[paired$macrogroup_id == macrogroup_id & paired$category_id %in% selected_categories, , drop = FALSE]
  if (nrow(df) == 0) return(list(df = data.frame(), gene_scores = data.frame(), categories = macro_summary))

  gene_scores <- stats::aggregate(
    max_gene_contribution_score ~ symbol,
    data = df,
    FUN = function(x) max(safe_num(x), na.rm = TRUE)
  )
  gene_scores$max_gene_contribution_score[!is.finite(gene_scores$max_gene_contribution_score)] <- NA_real_
  gene_recur <- stats::aggregate(category_id ~ symbol, data = unique(df[, c("symbol", "category_id")]), FUN = length)
  names(gene_recur)[2] <- "category_recurrence"
  gene_shared_opp <- stats::aggregate(contrast_gene_class ~ symbol, data = df, FUN = function(x) any(x == "shared_opposite_direction"))
  names(gene_shared_opp)[2] <- "has_shared_opposite"
  gene_shared_same <- stats::aggregate(contrast_gene_class ~ symbol, data = df, FUN = function(x) any(x == "shared_same_direction"))
  names(gene_shared_same)[2] <- "has_shared_same"
  gene_abs <- stats::aggregate(max_abs_log2FC ~ symbol, data = df, FUN = function(x) max(safe_num(x), na.rm = TRUE))
  gene_abs$max_abs_log2FC[!is.finite(gene_abs$max_abs_log2FC)] <- NA_real_
  gene_scores <- Reduce(function(x, y) merge(x, y, by = "symbol", all.x = TRUE, sort = FALSE),
                        list(gene_scores, gene_recur, gene_shared_opp, gene_shared_same, gene_abs))
  gene_scores <- gene_scores[order(
    -as.integer(gene_scores$has_shared_opposite %in% TRUE),
    -safe_num(gene_scores$max_gene_contribution_score),
    -safe_num(gene_scores$category_recurrence),
    -safe_num(gene_scores$max_abs_log2FC),
    gene_scores$symbol
  ), , drop = FALSE]
  selected_genes <- head(gene_scores$symbol, cfg$max_genes)
  df <- df[df$symbol %in% selected_genes, , drop = FALSE]
  gene_scores <- gene_scores[gene_scores$symbol %in% selected_genes, , drop = FALSE]

  category_order <- selected_categories[selected_categories %in% unique(df$category_id)]
  category_labels <- macro_summary$category_display_name[match(category_order, macro_summary$category_id)]
  category_labels <- short_label(category_labels, 26)
  names(category_labels) <- category_order
  gene_order <- rev(gene_scores$symbol)

  long <- rbind(
    data.frame(
      category_id = df$category_id,
      category_display_name = df$category_display_name,
      symbol = df$symbol,
      side = "A",
      log2FC = safe_num(df$log2FC_A),
      present = df$present_A,
      gene_class = df$contrast_gene_class,
      stringsAsFactors = FALSE
    ),
    data.frame(
      category_id = df$category_id,
      category_display_name = df$category_display_name,
      symbol = df$symbol,
      side = "B",
      log2FC = safe_num(df$log2FC_B),
      present = df$present_B,
      gene_class = df$contrast_gene_class,
      stringsAsFactors = FALSE
    )
  )
  long$category_side <- paste0(category_labels[long$category_id], "\n", long$side)
  long$category_side <- factor(long$category_side, levels = as.vector(rbind(paste0(category_labels, "\nA"), paste0(category_labels, "\nB"))))
  long$symbol <- factor(long$symbol, levels = gene_order)
  gene_scores$symbol <- factor(gene_scores$symbol, levels = gene_order)
  list(df = long, gene_scores = gene_scores, categories = macro_summary[macro_summary$category_id %in% category_order, , drop = FALSE])
}

plot_macrogroup_heatmap <- function(paired, summary, macrogroup_id, cfg, row) {
  prepared <- select_macrogroup_data(paired, summary, macrogroup_id, cfg)
  df <- prepared$df
  gene_scores <- prepared$gene_scores
  if (nrow(df) == 0) return(NULL)

  macrogroup_name <- first_nonempty(summary$macrogroup_name[summary$macrogroup_id == macrogroup_id])
  subtitle <- sprintf(
    "%s | A = %s; B = %s | %s categories | %s genes",
    row$output_id,
    row$label_a,
    row$label_b,
    length(unique(safe_chr(df$category_id))),
    length(unique(safe_chr(df$symbol)))
  )
  df$plot_lfc <- pmax(pmin(safe_num(df$log2FC), cfg$lfc_cap), -cfg$lfc_cap)
  df$text_label <- ifelse(is.na(df$log2FC), "", sprintf("%.2f", df$log2FC))
  df$text_color <- ifelse(!is.na(df$log2FC) & abs(df$log2FC) >= cfg$lfc_cap * 0.72, "white", "grey12")
  n_categories <- length(unique(safe_chr(df$category_id)))
  pair_separators <- if (n_categories > 1) {
    data.frame(x = seq(2.5, by = 2, length.out = n_categories - 1))
  } else {
    data.frame(x = numeric())
  }

  p_heatmap <- ggplot2::ggplot(df, ggplot2::aes(x = category_side, y = symbol, fill = plot_lfc)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.20)
  if (nrow(pair_separators) > 0) {
    p_heatmap <- p_heatmap +
      ggplot2::geom_vline(
        xintercept = pair_separators$x,
        color = "grey45",
        linewidth = 0.32
      )
  }
  p_heatmap <- p_heatmap +
    ggplot2::geom_text(ggplot2::aes(label = text_label, color = text_color), size = 1.75) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-cfg$lfc_cap, cfg$lfc_cap), na.value = "grey92", name = "log2FC"
    ) +
    ggplot2::labs(title = macrogroup_name, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 8.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, color = "grey15", size = 7.2),
      axis.text.y = ggplot2::element_text(face = "bold", color = "grey15", size = 6.4),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 8.2, color = "grey35"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "bottom"
    )

  marks <- rbind(
    data.frame(symbol = gene_scores$symbol, mark = "opposite\nsign", value = gene_scores$has_shared_opposite, stringsAsFactors = FALSE),
    data.frame(symbol = gene_scores$symbol, mark = "same\nsign", value = gene_scores$has_shared_same, stringsAsFactors = FALSE),
    data.frame(symbol = gene_scores$symbol, mark = "category\ncount", value = NA, stringsAsFactors = FALSE)
  )
  marks$mark <- factor(marks$mark, levels = c("opposite\nsign", "same\nsign", "category\ncount"))
  marks$count_label <- ""
  marks$count_label[marks$mark == "category\ncount"] <- as.character(gene_scores$category_recurrence)
  p_marks <- ggplot2::ggplot(marks, ggplot2::aes(x = mark, y = symbol)) +
    ggplot2::geom_point(
      data = marks[marks$mark != "category\ncount", , drop = FALSE],
      ggplot2::aes(fill = value, color = value),
      shape = 21,
      size = 2.15,
      stroke = 0.6
    ) +
    ggplot2::geom_text(
      data = marks[marks$mark == "category\ncount", , drop = FALSE],
      ggplot2::aes(label = count_label),
      size = 2.2,
      color = "grey20"
    ) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#0AA36E", `FALSE` = "white")) +
    ggplot2::scale_color_manual(values = c(`TRUE` = "#0A7F58", `FALSE` = "#C8CDD3")) +
    ggplot2::labs(x = NULL, y = NULL, title = "gene-level\nsummary") +
    ggplot2::theme_minimal(base_size = 8.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 0, hjust = 0.5, color = "grey20", size = 6.4, lineheight = 0.92),
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(size = 7.2, face = "bold", color = "grey20", hjust = 0.5, lineheight = 0.95),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )

  p_heatmap | p_marks + patchwork::plot_layout(widths = c(3.3, 0.72))
}

write_heatmaps <- function(paired, summary, macrogroups, paths, cfg, row) {
  dir.create(paths$heatmap_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(Sys.glob(file.path(paths$heatmap_dir, "*_contrast_macrogroup_heatmap.png")), force = TRUE)
  pdf_path <- file.path(paths$heatmap_dir, paste0(paths$prefix, "_contrast_macrogroup_heatmaps.pdf"))
  figure_width <- 13.2
  figure_height <- max(7.5, 3.6 + cfg$max_genes * 0.15)
  grDevices::pdf(pdf_path, width = figure_width, height = figure_height, onefile = TRUE, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)

  index <- vector("list", length(macrogroups))
  for (i in seq_along(macrogroups)) {
    macrogroup_id <- macrogroups[[i]]
    plot <- plot_macrogroup_heatmap(paired, summary, macrogroup_id, cfg, row)
    if (is.null(plot)) next
    print(plot)
    png_name <- sprintf("%02d_%s_contrast_macrogroup_heatmap.png", i, safe_file_component(macrogroup_id))
    png_path <- file.path(paths$heatmap_dir, png_name)
    ggplot2::ggsave(png_path, plot, width = figure_width, height = figure_height, dpi = 220, bg = "white", limitsize = FALSE)
    rows <- paired[paired$macrogroup_id == macrogroup_id, , drop = FALSE]
    index[[i]] <- data.frame(
      contrast_id = row$contrast_id,
      output_id = row$output_id,
      universe = cfg$universe,
      heatmap_rank = i,
      macrogroup_id = macrogroup_id,
      macrogroup_name = first_nonempty(rows$macrogroup_name),
      n_categories_available = length(unique(rows$category_id)),
      n_genes_available = length(unique(rows$symbol)),
      n_genes_visible = ifelse(cfg$max_genes > 0, min(cfg$max_genes, length(unique(rows$symbol))), length(unique(rows$symbol))),
      n_genes_total = length(unique(rows$symbol)),
      gene_visibility = sprintf("%s/%s", ifelse(cfg$max_genes > 0, min(cfg$max_genes, length(unique(rows$symbol))), length(unique(rows$symbol))), length(unique(rows$symbol))),
      max_categories_per_macrogroup = cfg$max_categories_per_macrogroup,
      max_genes = cfg$max_genes,
      lfc_cap = cfg$lfc_cap,
      png_path = png_path,
      pdf_path = pdf_path,
      stringsAsFactors = FALSE
    )
  }
  index <- index[!vapply(index, is.null, logical(1))]
  if (length(index) == 0) data.frame() else do.call(rbind, index)
}

write_readme <- function(paths, cfg, row, n_heatmaps) {
  lines <- c(
    "Contrast macrogroup heatmaps",
    "",
    sprintf("Contrast: %s_%s", row$contrast_id, row$output_id),
    sprintf("Universe: %s", cfg$universe),
    sprintf("A: %s (%s)", row$contrast_a, row$label_a),
    sprintf("B: %s (%s)", row$contrast_b, row$label_b),
    "",
    "Purpose",
    "These plots summarize paired gene-level contrast evidence by LISA macrogroup.",
    "Each selected category is shown as two adjacent columns: A and B.",
    "Tile color encodes log2FC, capped for display only; the printed number is the uncapped value.",
    "",
    "Right-side markers",
    "- opposite sign: gene appears in both analyses with opposite log2FC signs.",
    "- same sign: gene appears in both analyses with the same log2FC sign.",
    "- category count: number of selected categories in this macrogroup containing the gene.",
    "",
    sprintf("Heatmaps generated: %s", n_heatmaps),
    "",
    "Interpretation caution: this is a compact triage view over post-processed gene-category evidence, not a causal network."
  )
  writeLines(lines, file.path(paths$heatmap_dir, "README_contrast_macrogroup_heatmaps.txt"), useBytes = TRUE)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Required packages: ggplot2, patchwork", call. = FALSE)
  }
  row <- contrast_row(cfg$project_dir, cfg$contrast_id)
  paths <- paths_for(cfg$project_dir, row, cfg$universe)
  missing <- c(paths$paired, paths$summary)[!file.exists(c(paths$paired, paths$summary))]
  if (length(missing) > 0) stop(sprintf("Missing required input files:\n%s", paste(missing, collapse = "\n")), call. = FALSE)
  paired <- read_tsv(paths$paired)
  summary <- read_tsv(paths$summary)
  macrogroups <- pick_macrogroups(summary, paired, cfg)
  if (length(macrogroups) == 0) stop("No macrogroups selected.", call. = FALSE)
  index <- write_heatmaps(paired, summary, macrogroups, paths, cfg, row)
  index_path <- file.path(paths$heatmap_dir, paste0(paths$prefix, "_contrast_macrogroup_heatmaps_index.tsv"))
  write_tsv(index, index_path)
  write_readme(paths, cfg, row, nrow(index))
  print(index[, c("heatmap_rank", "macrogroup_id", "n_categories_available", "n_genes_available")])
}

if (sys.nframe() == 0L) main()
