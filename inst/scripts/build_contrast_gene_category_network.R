#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

default_project_dir <- NA_character_

required_pkgs <- "ggplot2"
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(sprintf("Missing required R packages: %s", paste(missing_pkgs, collapse = ", ")), call. = FALSE)
}

read_tsv <- function(path) {
  if (file.exists(path) && isTRUE(file.info(path)$size == 0)) return(data.frame())
  first_line <- readLines(path, n = 1, warn = FALSE)
  if (length(first_line) == 0 || trimws(first_line[[1]]) == "") return(data.frame())
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
    top_categories = 0,
    top_genes_per_category = 7,
    # --- exploration seams (H3) --------------------------------------------
    # Opt-in. Without them this script behaves exactly as it always has.
    #
    #   --variant        the single accepted exploration variant,
    #                    `all_categories_top7`. It MAPS ONTO the accepted
    #                    settings already defaulted above -- every category and
    #                    up to seven genes per category -- so naming it selects
    #                    the accepted figure rather than introducing a new
    #                    scientific control. There is deliberately no second
    #                    variant and no category-local crop.
    #   --paired-input   read an already prepared paired evidence table and its
    #   --summary-input  complete ranked category summary instead of the default
    #                    sibling paths, so this network never has to make the
    #                    card/heatmap builder run to satisfy its dependency.
    variant = NA_character_,
    paired_input = NA_character_,
    summary_input = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    if (i == length(args)) stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    val <- args[[i + 1]]
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(out)) stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    if (key %in% c("top_categories", "top_genes_per_category")) {
      out[[key]] <- as.integer(val)
    } else {
      out[[key]] <- val
    }
    i <- i + 2
  }
  if (is.na(out$project_dir) || out$project_dir == "") stop("Required argument: --project-dir", call. = FALSE)
  if (is.na(out$contrast_id) || out$contrast_id == "") stop("Required argument: --contrast-id", call. = FALSE)
  variant <- if (is.na(out$variant)) "" else out$variant
  if (nzchar(variant)) {
    if (!identical(variant, "all_categories_top7")) {
      stop(sprintf("LISA-CONTRAST-050 --variant accepts only all_categories_top7; received: %s.", variant), call. = FALSE)
    }
    out$top_categories <- 0L
    out$top_genes_per_category <- 7L
  }
  out$variant <- variant
  supplied <- c(!is.na(out$paired_input) && nzchar(out$paired_input),
                !is.na(out$summary_input) && nzchar(out$summary_input))
  if (any(supplied) && !all(supplied)) {
    stop("LISA-CONTRAST-051 --paired-input and --summary-input must be given together; a paired table without its ranked summary is not a preparation.", call. = FALSE)
  }
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

# ggplot2 calls oob functions with the values and the scale limits. Clamp only
# finite values, matching the behaviour previously supplied by scales::squish,
# while keeping this installed script free of an otherwise unused dependency.
squish_finite <- function(x, range) {
  finite <- is.finite(x)
  x[finite] <- pmax(range[[1L]], pmin(range[[2L]], x[finite]))
  x
}

first_nonempty <- function(x) {
  x <- safe_chr(x)
  x <- x[x != ""]
  if (length(x) == 0) "" else x[[1]]
}

short_label <- function(x, max_chars = 34) {
  x <- safe_chr(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

first_non_na_num <- function(x) {
  x <- safe_num(x)
  x <- x[!is.na(x)]
  if (length(x) == 0) NA_real_ else x[[1]]
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
    network_dir = file.path(out_dir, "contrast_gene_category_network"),
    paired = file.path(out_dir, paste0(prefix, "_paired_gene_evidence.tsv")),
    summary = file.path(out_dir, paste0(prefix, "_contrast_category_gene_summary.tsv")),
    prefix = prefix
  )
}

select_edges <- function(paired, summary, cfg) {
  if (nrow(paired) == 0 || nrow(summary) == 0) return(data.frame())
  category_ids <- summary$category_id
  if (cfg$top_categories > 0) {
    category_ids <- head(category_ids, cfg$top_categories)
  }
  rows <- lapply(category_ids, function(cid) {
    df <- paired[paired$category_id == cid, , drop = FALSE]
    if (nrow(df) == 0) return(NULL)
    class_priority <- match(df$contrast_gene_class, c(
      "shared_opposite_direction", "shared_same_direction",
      "A_specific", "B_specific", "shared_mixed_or_missing_direction"
    ))
    class_priority[is.na(class_priority)] <- 99
    df <- df[order(class_priority, -safe_num(df$max_gene_contribution_score), -safe_num(df$max_abs_log2FC), df$symbol), , drop = FALSE]
    df <- df[!duplicated(df$symbol), , drop = FALSE]
    if (cfg$top_genes_per_category > 0) head(df, cfg$top_genes_per_category) else df
  })
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out)) data.frame() else out
}

plot_network <- function(edges, out_png, out_pdf, cfg) {
  if (nrow(edges) == 0) return(FALSE)
  edges$category_label <- short_label(edges$category_display_name, 32)
  edges$gene_label <- edges$symbol
  category_order <- unique(edges$category_id)
  gene_order <- unique(edges$symbol)

  cat_nodes <- data.frame(
    node_id = category_order,
    label = vapply(category_order, function(cid) first_nonempty(edges$category_label[edges$category_id == cid]), character(1)),
    x = 0,
    y = seq(length(gene_order), 1, length.out = length(category_order)),
    type = "category",
    stringsAsFactors = FALSE
  )
  gene_nodes <- data.frame(
    node_id = gene_order,
    label = gene_order,
    x = 1,
    y = seq(length(gene_order), 1),
    type = "gene",
    stringsAsFactors = FALSE
  )
  node_y <- c(setNames(cat_nodes$y, cat_nodes$node_id), setNames(gene_nodes$y, gene_nodes$node_id))
  edges$x <- 0.08
  edges$xend <- 0.92
  edges$y <- node_y[edges$category_id]
  edges$yend <- node_y[edges$symbol]
  edges$class_label <- gsub("_", " ", edges$contrast_gene_class)

  gene_fc <- data.frame(
    symbol = gene_order,
    log2FC_A = vapply(gene_order, function(g) first_non_na_num(edges$log2FC_A_num[edges$symbol == g]), numeric(1)),
    log2FC_B = vapply(gene_order, function(g) first_non_na_num(edges$log2FC_B_num[edges$symbol == g]), numeric(1)),
    stringsAsFactors = FALSE
  )
  heatmap_df <- rbind(
    data.frame(symbol = gene_fc$symbol, contrast = "A", x = 1.24, y = node_y[gene_fc$symbol], log2FC = gene_fc$log2FC_A, stringsAsFactors = FALSE),
    data.frame(symbol = gene_fc$symbol, contrast = "B", x = 1.37, y = node_y[gene_fc$symbol], log2FC = gene_fc$log2FC_B, stringsAsFactors = FALSE)
  )
  heatmap_df$value_label <- ifelse(is.na(heatmap_df$log2FC), "", sprintf("%.1f", heatmap_df$log2FC))
  heatmap_df$text_is_light <- !is.na(heatmap_df$log2FC) & abs(heatmap_df$log2FC) >= 1.1
  max_abs_fc_raw <- max(abs(heatmap_df$log2FC), na.rm = TRUE)
  if (!is.finite(max_abs_fc_raw) || max_abs_fc_raw <= 0) max_abs_fc_raw <- 1
  max_abs_fc <- max(1, min(3, ceiling(max_abs_fc_raw * 2) / 2))
  fill_legend_title <- if (max_abs_fc_raw > max_abs_fc) "Gene log2FC\n(color capped)" else "Gene log2FC"

  class_colors <- c(
    shared_opposite_direction = "#7B3294",
    shared_same_direction = "#008837",
    A_specific = "#2C7FB8",
    B_specific = "#D95F02",
    shared_mixed_or_missing_direction = "grey45"
  )
  class_labels <- c(
    shared_opposite_direction = "shared opposite",
    shared_same_direction = "shared same",
    A_specific = "A only",
    B_specific = "B only",
    shared_mixed_or_missing_direction = "shared mixed"
  )
  h <- max(5.2, 0.16 * length(gene_order) + 0.22 * length(category_order))
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = x, xend = xend, y = y, yend = yend, color = contrast_gene_class),
      linewidth = 0.35,
      alpha = 0.42
    ) +
    ggplot2::geom_point(data = cat_nodes, ggplot2::aes(x = x, y = y), size = 3.0, color = "grey20") +
    ggplot2::geom_point(data = gene_nodes, ggplot2::aes(x = x, y = y), size = 2.0, color = "grey30") +
    ggplot2::geom_text(data = cat_nodes, ggplot2::aes(x = x - 0.025, y = y, label = label), hjust = 1, size = 2.6, color = "grey18") +
    ggplot2::geom_text(data = gene_nodes, ggplot2::aes(x = x + 0.025, y = y, label = label), hjust = 0, size = 2.45, color = "grey18") +
    ggplot2::geom_tile(
      data = heatmap_df,
      ggplot2::aes(x = x, y = y, fill = log2FC),
      width = 0.105,
      height = 0.72,
      color = "white",
      linewidth = 0.2
    ) +
    ggplot2::geom_text(
      data = heatmap_df[!heatmap_df$text_is_light, , drop = FALSE],
      ggplot2::aes(x = x, y = y, label = value_label),
      size = 1.75,
      color = "grey15",
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = heatmap_df[heatmap_df$text_is_light, , drop = FALSE],
      ggplot2::aes(x = x, y = y, label = value_label),
      size = 1.75,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::annotate("text", x = 1.24, y = max(gene_nodes$y) + 0.68, label = "A", size = 2.45, fontface = "bold", color = "grey20") +
    ggplot2::annotate("text", x = 1.37, y = max(gene_nodes$y) + 0.68, label = "B", size = 2.45, fontface = "bold", color = "grey20") +
    ggplot2::annotate("text", x = 1.305, y = max(gene_nodes$y) + 0.98, label = "log2FC", size = 2.2, color = "grey35") +
    ggplot2::scale_color_manual(values = class_colors, labels = class_labels, drop = FALSE, name = "Gene class") +
    ggplot2::scale_fill_gradient2(
      low = "#2C7FB8",
      mid = "white",
      high = "#D95F02",
      midpoint = 0,
      limits = c(-max_abs_fc, max_abs_fc),
      oob = squish_finite,
      na.value = "grey92",
      name = fill_legend_title
    ) +
    ggplot2::coord_cartesian(xlim = c(-0.38, 1.53), ylim = c(0.5, max(cat_nodes$y, gene_nodes$y) + 1.08), clip = "off") +
    ggplot2::labs(
      title = "Category-gene network",
      subtitle = sprintf("A = %s; B = %s\n%s categories; up to %s genes/category", cfg$label_a, cfg$label_b, ifelse(cfg$top_categories > 0, paste("Top", cfg$top_categories), "All"), cfg$top_genes_per_category),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 9) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "grey35"),
      plot.margin = ggplot2::margin(10, 60, 10, 120)
    )
  ggplot2::ggsave(out_png, p, width = 9.25, height = h, dpi = 180, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(out_pdf, p, width = 9.25, height = h, bg = "white", limitsize = FALSE)
  save_svg_twin(out_png, p, 9.25, h)
  TRUE
}

# D4: every PNG product also gets a vector SVG twin drawn with the cairo svg()
# device (no extra dependency). The run's format policy removes it afterwards
# when SVG was not requested, exactly as for the other generators.
save_svg_twin <- function(out_png, plot, width, height) {
  out_svg <- sub("[.]png$", ".svg", out_png)
  tryCatch(ggplot2::ggsave(out_svg, plot, width = width, height = height,
    device = grDevices::svg, bg = "white", limitsize = FALSE),
    error = function(e) warning("SVG twin not written for ", basename(out_png), ": ",
      conditionMessage(e), call. = FALSE))
  invisible(out_svg)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  row <- contrast_row(cfg$project_dir, cfg$contrast_id)
  cfg$label_a <- first_nonempty(row$label_a)
  cfg$label_b <- first_nonempty(row$label_b)
  paths <- paths_for(cfg$project_dir, row, cfg$universe)
  # An explicitly supplied prepared pair takes precedence over the default
  # sibling paths. This is the whole reason the network can be requested on its
  # own: its dependency is satisfied by the shared preparation, so it never has
  # to invoke the card/heatmap builder and never produces a sibling figure.
  paired_path <- if (!is.na(cfg$paired_input) && nzchar(cfg$paired_input)) cfg$paired_input else paths$paired
  summary_path <- if (!is.na(cfg$summary_input) && nzchar(cfg$summary_input)) cfg$summary_input else paths$summary
  if (!file.exists(paired_path) || !file.exists(summary_path)) {
    stop("Run build_contrast_gene_level_product.R before building the network.", call. = FALSE)
  }
  dir.create(paths$network_dir, recursive = TRUE, showWarnings = FALSE)
  paired <- read_tsv(paired_path)
  summary <- read_tsv(summary_path)
  category_ids_total <- summary$category_id
  if (cfg$top_categories > 0) {
    category_ids_total <- head(category_ids_total, cfg$top_categories)
  }
  total_rows <- paired[paired$category_id %in% category_ids_total, , drop = FALSE]
  edges <- select_edges(paired, summary, cfg)
  # Owner and collection already appear in the containing directories.
  network_prefix <- cfg$universe
  write_tsv(edges, file.path(paths$network_dir, paste0(network_prefix, "_network_edges.tsv")))

  out_png <- file.path(paths$network_dir, paste0(network_prefix, "_category_gene_network.png"))
  out_pdf <- file.path(paths$network_dir, paste0(network_prefix, "_category_gene_network.pdf"))
  network_written <- plot_network(edges, out_png, out_pdf, cfg)
  # The recipe records the EFFECTIVE configuration. The H3 seams are dropped
  # when unset -- exactly as `project_dir` and the renderer hash already are --
  # so a default command-line or FULL run writes the same `native_config_json`
  # it wrote before these flags existed.
  recipe_cfg <- cfg
  for (field in c("variant", "paired_input", "summary_input")) {
    value <- recipe_cfg[[field]]
    if (is.null(value) || is.na(value) || !nzchar(value)) recipe_cfg[[field]] <- NULL
  }
  if (network_written) lisaR:::lisa_write_native_extended_recipe(edges, recipe_cfg, out_png,
    "build_contrast_gene_category_network.R", "native_contrast_network")
  if (nzchar(cfg$variant) && !isTRUE(network_written)) {
    stop("LISA-CONTRAST-052 this contrast/collection has no plottable paired gene evidence; no network was drawn.", call. = FALSE)
  }
  manifest <- data.frame(
    contrast_id = row$contrast_id,
    output_id = row$output_id,
    universe = cfg$universe,
    top_categories = cfg$top_categories,
    top_genes_per_category = cfg$top_genes_per_category,
    n_edges = nrow(edges),
    n_categories = length(unique(edges$category_id)),
    n_genes = length(unique(edges$symbol)),
    n_genes_visible = length(unique(edges$symbol)),
    n_genes_total = length(unique(total_rows$symbol)),
    gene_visibility = sprintf("%s/%s", length(unique(edges$symbol)), length(unique(total_rows$symbol))),
    network_written = network_written,
    output_dir = paths$network_dir,
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(paths$network_dir, paste0(network_prefix, "_network_manifest.tsv")))
  # The exploration variant writes its own scoped product inventory naming the
  # files it actually produced, so collection never has to guess from directory
  # contents. It is written only in variant mode, leaving the default output set
  # unchanged.
  if (nzchar(cfg$variant)) {
    write_tsv(data.frame(
      contrast_id = row$contrast_id,
      output_id = row$output_id,
      universe = cfg$universe,
      variant = cfg$variant,
      top_categories = cfg$top_categories,
      top_genes_per_category = cfg$top_genes_per_category,
      n_categories = length(unique(edges$category_id)),
      n_genes = length(unique(edges$symbol)),
      output_png = out_png,
      output_pdf = out_pdf,
      stringsAsFactors = FALSE
    ), file.path(paths$network_dir,
      paste0(network_prefix, "_contrast_gene_category_network_index.tsv")))
  }
  writeLines(c(
    "Contrast category-gene network",
    "",
    sprintf("Contrast: %s_%s", row$contrast_id, row$output_id),
    sprintf("Universe: %s", cfg$universe),
    "This is a deliberately small network view of the top contrast categories and selected supporting genes.",
    "The two-column heatmap next to each gene shows gene-level log2FC in A and B, using the same A/B mapping shown in the figure subtitle.",
    "It is meant for visual triage, not exhaustive gene discovery.",
    "",
    "Use the paired gene evidence table for complete data."
  ), file.path(paths$network_dir, "README_contrast_gene_category_network.txt"))
  print(manifest)
}

if (sys.nframe() == 0L) main()
