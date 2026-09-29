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
source(system.file("scripts/lisa_plot_metadata.R", package = "lisaR"), local = TRUE)

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
    analysis_id = NA_character_,
    universe = "GOBP-C2",
    top_genes = 40,
    top_categories = 0,
    min_recurrent_categories = 3,
    min_high_categories = 5,
    min_high_macrogroups = 2,
    de_padj_cutoff = 0.05,
    lfc_cutoff = 0.58
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    if (i == length(args)) stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    val <- args[[i + 1]]
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(out)) stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    if (key %in% c("top_genes", "top_categories", "min_recurrent_categories", "min_high_categories", "min_high_macrogroups")) {
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
  out
}

first_existing_path <- function(candidates) {
  hits <- candidates[file.exists(candidates)]
  if (length(hits) > 0) return(hits[[1]])
  ""
}

first_glob_path <- function(patterns) {
  for (pattern in patterns) {
    hits <- Sys.glob(pattern)
    hits <- hits[file.exists(hits)]
    if (length(hits) > 0) return(sort(hits)[[1]])
  }
  ""
}

paths_for <- function(project_dir, analysis_id, universe) {
  out_dir <- file.path(project_dir, "outputs", "gene_level", "single_de", analysis_id, paste0("collection_", universe))
  prefix <- paste(analysis_id, universe, "gene_level", sep = "_")
  evidence <- first_existing_path(c(
    file.path(out_dir, paste0(prefix, "_gene_category_contributions.tsv")),
    first_glob_path(c(
      file.path(out_dir, paste0(prefix, "_gene_category*.tsv")),
      file.path(out_dir, paste0(prefix, "_gene_catego*.tsv"))
    ))
  ))
  list(
    out_dir = out_dir,
    evidence = evidence,
    recurrent_dir = file.path(out_dir, "recurrent_gene_screen"),
    prefix = prefix
  )
}

safe_num <- function(x) {
  suppressWarnings(as.numeric(x))
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

collapse_unique <- function(x, sep = ";") {
  x <- sort(unique(safe_chr(x)))
  x <- x[x != ""]
  paste(x, collapse = sep)
}

collapse_top_categories <- function(df, n = 0) {
  df <- df[order(-safe_num(df$gene_contribution_score),
                 safe_num(df$min_source_padj),
                 df$category_display_name), , drop = FALSE]
  df <- df[!duplicated(df$category_id), , drop = FALSE]
  if (nrow(df) == 0) return("")
  vals <- sprintf("%s [%s]", df$category_display_name, df$macrogroup_name)
  if (n > 0) {
    vals <- head(vals, n)
  }
  paste(vals, collapse = " | ")
}

support_tier <- function(category_count, macrogroup_count, de_sig, abs_lfc, le_count, ora_count, cfg) {
  evidence_supported <- (le_count >= 2) | (ora_count >= 1)
  effect_supported <- de_sig & !is.na(abs_lfc) & abs_lfc >= cfg$lfc_cutoff
  if (category_count >= cfg$min_high_categories &&
      macrogroup_count >= cfg$min_high_macrogroups &&
      effect_supported &&
      evidence_supported) {
    return("high_support_recurrent_candidate")
  }
  if (category_count >= cfg$min_recurrent_categories &&
      effect_supported &&
      evidence_supported) {
    return("moderate_support_recurrent_candidate")
  }
  if (category_count >= cfg$min_recurrent_categories) {
    return("recurrent_only")
  }
  "limited_recurrence"
}

build_recurrent_summary <- function(evidence, cfg) {
  evidence$symbol <- toupper(safe_chr(evidence$symbol))
  evidence$category_id <- safe_chr(evidence$category_id)
  evidence$macrogroup_id <- safe_chr(evidence$macrogroup_id)
  evidence$macrogroup_name <- safe_chr(evidence$macrogroup_name)
  evidence$category_display_name <- safe_chr(evidence$category_display_name)
  evidence$log2FC_num <- safe_num(evidence$log2FC)
  evidence$padj_num <- safe_num(evidence$padj)
  evidence$gene_contribution_score_num <- safe_num(evidence$gene_contribution_score)
  evidence$source_geneset_count_num <- safe_num(evidence$source_geneset_count)
  evidence$min_source_padj_num <- safe_num(evidence$min_source_padj)
  evidence$in_gsea_leading_edge_bool <- as.logical(evidence$in_gsea_leading_edge)
  evidence$in_ora_overlap_bool <- as.logical(evidence$in_ora_overlap)
  evidence$de_is_significant_bool <- as.logical(evidence$de_is_significant)

  split_rows <- split(evidence, evidence$symbol)
  rows <- lapply(split_rows, function(df) {
    category_count <- length(unique(df$category_id[df$category_id != ""]))
    macrogroup_count <- length(unique(df$macrogroup_id[df$macrogroup_id != ""]))
    gsea_count <- length(unique(df$category_id[df$in_gsea_leading_edge_bool %in% TRUE]))
    ora_count <- length(unique(df$category_id[df$in_ora_overlap_bool %in% TRUE]))
    de_sig <- any(df$de_is_significant_bool %in% TRUE)
    lfc_values <- df$log2FC_num[!is.na(df$log2FC_num)]
    lfc <- if (length(lfc_values) > 0) lfc_values[1] else NA_real_
    padj <- suppressWarnings(min(df$padj_num, na.rm = TRUE))
    if (!is.finite(padj)) padj <- NA_real_
    contribution_max <- suppressWarnings(max(df$gene_contribution_score_num, na.rm = TRUE))
    if (!is.finite(contribution_max)) contribution_max <- NA_real_
    contribution_mean <- suppressWarnings(mean(df$gene_contribution_score_num, na.rm = TRUE))
    if (!is.finite(contribution_mean)) contribution_mean <- NA_real_
    min_source_padj <- suppressWarnings(min(df$min_source_padj_num, na.rm = TRUE))
    if (!is.finite(min_source_padj)) min_source_padj <- NA_real_
    source_geneset_total <- suppressWarnings(sum(df$source_geneset_count_num, na.rm = TRUE))
    if (!is.finite(source_geneset_total)) source_geneset_total <- NA_real_
    tier <- support_tier(
      category_count = category_count,
      macrogroup_count = macrogroup_count,
      de_sig = de_sig,
      abs_lfc = abs(lfc),
      le_count = gsea_count,
      ora_count = ora_count,
      cfg = cfg
    )
    data.frame(
      analysis_id = df$analysis_id[1],
      universe = cfg$universe,
      symbol = df$symbol[1],
      category_count = category_count,
      macrogroup_count = macrogroup_count,
      gsea_leading_edge_category_count = gsea_count,
      ora_overlap_category_count = ora_count,
      source_geneset_count_total = source_geneset_total,
      max_gene_contribution_score = contribution_max,
      mean_gene_contribution_score = contribution_mean,
      log2FC = lfc,
      abs_log2FC = abs(lfc),
      padj = padj,
      de_is_significant = de_sig,
      min_source_padj = min_source_padj,
      recurrent_support_tier = tier,
      macrogroups = collapse_unique(df$macrogroup_name),
      top_contributing_categories = collapse_top_categories(df, cfg$top_categories),
      interpretation_note = "Recurrent association across LISA categories; this is a prioritization signal, not causal driver evidence.",
      stringsAsFactors = FALSE
    )
  })
  summary <- do.call(rbind, rows)
  summary <- summary[order(
    -safe_num(summary$category_count),
    -safe_num(summary$macrogroup_count),
    -safe_num(summary$de_is_significant),
    -safe_num(summary$gsea_leading_edge_category_count),
    -safe_num(summary$ora_overlap_category_count),
    -safe_num(summary$abs_log2FC),
    safe_num(summary$padj),
    summary$symbol
  ), , drop = FALSE]
  rownames(summary) <- NULL
  summary$recurrent_rank <- seq_len(nrow(summary))
  summary[, c("recurrent_rank", setdiff(names(summary), "recurrent_rank")), drop = FALSE]
}

empty_recurrent_summary <- function() {
  data.frame(
    recurrent_rank = integer(),
    analysis_id = character(),
    universe = character(),
    symbol = character(),
    category_count = integer(),
    macrogroup_count = integer(),
    gsea_leading_edge_category_count = integer(),
    ora_overlap_category_count = integer(),
    source_geneset_count_total = numeric(),
    max_gene_contribution_score = numeric(),
    mean_gene_contribution_score = numeric(),
    log2FC = numeric(),
    abs_log2FC = numeric(),
    padj = numeric(),
    de_is_significant = logical(),
    min_source_padj = numeric(),
    recurrent_support_tier = character(),
    macrogroups = character(),
    top_contributing_categories = character(),
    interpretation_note = character(),
    stringsAsFactors = FALSE
  )
}

build_plot <- function(summary, cfg) {
  plot_df <- summary[summary$category_count >= cfg$min_recurrent_categories, , drop = FALSE]
  if (nrow(plot_df) == 0) {
    plot_df <- head(summary, cfg$top_genes)
  } else {
    plot_df <- head(plot_df, cfg$top_genes)
  }
  plot_df$symbol <- factor(plot_df$symbol, levels = rev(plot_df$symbol))
  plot_df$effect_direction <- ifelse(safe_num(plot_df$log2FC) >= 0, "up", "down")
  plot_df$neg_log10_padj <- -log10(pmax(safe_num(plot_df$padj), 1e-300))
  plot_df$neg_log10_padj[!is.finite(plot_df$neg_log10_padj)] <- NA_real_
  tier_levels <- c(
    "high_support_recurrent_candidate",
    "moderate_support_recurrent_candidate",
    "recurrent_only",
    "limited_recurrence"
  )
  plot_df$recurrent_support_tier <- factor(plot_df$recurrent_support_tier, levels = tier_levels)

  ggplot2::ggplot(plot_df, ggplot2::aes(y = symbol)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = category_count, yend = symbol),
      color = "grey80",
      linewidth = 0.55
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = category_count, size = macrogroup_count, fill = recurrent_support_tier, shape = effect_direction),
      color = "grey15",
      stroke = 0.35,
      alpha = 0.95
    ) +
    ggplot2::scale_shape_manual(values = c(up = 24, down = 25), name = "log2FC direction") +
    ggplot2::scale_fill_manual(
      values = c(
        high_support_recurrent_candidate = "#0072B2",
        moderate_support_recurrent_candidate = "#009E73",
        recurrent_only = "#E69F00",
        limited_recurrence = "#BDBDBD"
      ),
      labels = c(
        high_support_recurrent_candidate = "High-support recurrent candidate",
        moderate_support_recurrent_candidate = "Moderate-support recurrent candidate",
        recurrent_only = "Recurrent only",
        limited_recurrence = "Limited recurrence"
      ),
      name = "Support tier"
    ) +
    ggplot2::scale_size_continuous(range = c(1.6, 4.8), breaks = sort(unique(plot_df$macrogroup_count)), name = "Macrogroup count") +
    ggplot2::labs(
      title = "Cross-category recurrent genes",
      subtitle = lisa_plot_subtitle(
        cfg$plot_metadata,
        suffix = sprintf(
          "%s | recurrence is a prioritization signal, not causal evidence",
          cfg$universe
        )
      ),
      x = "LISA categories containing this gene",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey35"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      legend.position = "right"
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(shape = 21, size = 3.2, color = "grey15")),
      shape = ggplot2::guide_legend(override.aes = list(size = 3.2, fill = "white", color = "grey15"))
    )
}

build_empty_plot <- function(cfg, reason) {
  ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_blank() +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 1.05,
      label = "No recurrent gene screen signal",
      fontface = "bold",
      size = 5,
      color = "grey25"
    ) +
    ggplot2::annotate(
      "text",
      x = 1,
      y = 0.95,
      label = reason,
      size = 3.5,
      color = "grey40"
    ) +
    ggplot2::labs(
      title = "Cross-category recurrent genes",
      subtitle = lisa_plot_subtitle(
        cfg$plot_metadata, suffix = cfg$universe
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle = ggplot2::element_text(size = 9, color = "grey35"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

write_readme <- function(path, cfg) {
  txt <- c(
    "Cross-category recurrent gene screen",
    "",
    "Purpose",
    "This output summarizes genes that appear across multiple altered LISA categories.",
    "The recurrence signal is intended for prioritization and interpretation, not as causal evidence.",
    "",
    "Important caution",
    "Genes in the high or moderate support tiers should not be called drivers based on this output alone.",
    "They are recurrent candidates for follow-up because they combine category breadth with DE/effect and enrichment support.",
    "",
    "Support tiers",
    sprintf("- high_support_recurrent_candidate: category_count >= %s, macrogroup_count >= %s, DE significant at padj <= %.3g, |log2FC| >= %.3g, and enrichment support.",
            cfg$min_high_categories, cfg$min_high_macrogroups, cfg$de_padj_cutoff, cfg$lfc_cutoff),
    sprintf("- moderate_support_recurrent_candidate: category_count >= %s, DE significant at padj <= %.3g, |log2FC| >= %.3g, and enrichment support.",
            cfg$min_recurrent_categories, cfg$de_padj_cutoff, cfg$lfc_cutoff),
    sprintf("- recurrent_only: category_count >= %s but incomplete DE/effect/enrichment support.", cfg$min_recurrent_categories),
    "- limited_recurrence: below the recurrence threshold.",
    "",
    "Enrichment support means at least two GSEA leading-edge category hits or at least one ORA-overlap category hit.",
    "",
    "Key columns",
    "- category_count: number of LISA categories containing the gene.",
    "- macrogroup_count: number of LISA macrogroups containing the gene.",
    "- gsea_leading_edge_category_count: number of categories where the gene is in a GSEA leading edge.",
    "- ora_overlap_category_count: number of categories where the gene appears in ORA overlap.",
    "- max_gene_contribution_score and mean_gene_contribution_score: previous display-ranking score summarized across categories.",
    "- top_contributing_categories: highest-scoring categories supporting the recurrent signal.",
    "",
    "This is post-processing of the existing gene-level evidence table; no DE, GSEA, ORA, or LISA step is rerun."
  )
  writeLines(txt, con = path, useBytes = TRUE)
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  cfg$plot_metadata <- lisa_plot_metadata(cfg$project_dir, cfg$analysis_id)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Required package: ggplot2", call. = FALSE)

  paths <- paths_for(cfg$project_dir, cfg$analysis_id, cfg$universe)
  if (!file.exists(paths$evidence)) stop(sprintf("Missing gene-level evidence table: %s", paths$evidence), call. = FALSE)
  dir.create(paths$recurrent_dir, recursive = TRUE, showWarnings = FALSE)

  evidence <- read_tsv(paths$evidence)
  summary <- if (nrow(evidence) == 0) empty_recurrent_summary() else build_recurrent_summary(evidence, cfg)
  summary_path <- file.path(paths$recurrent_dir, paste0(paths$prefix, "_recurrent_gene_screen.tsv"))
  write_tsv(summary, summary_path)

  recurrent_only <- summary[summary$category_count >= cfg$min_recurrent_categories, , drop = FALSE]
  top_path <- file.path(paths$recurrent_dir, paste0(paths$prefix, "_top_recurrent_genes.tsv"))
  top_recurrent <- if (cfg$top_genes > 0) head(recurrent_only, cfg$top_genes) else recurrent_only
  write_tsv(top_recurrent, top_path)
  n_recurrent_visible <- ifelse(cfg$top_genes > 0, min(cfg$top_genes, nrow(recurrent_only)), nrow(recurrent_only))

  png_path <- file.path(paths$recurrent_dir, paste0(paths$prefix, "_top_recurrent_genes.png"))
  pdf_path <- file.path(paths$recurrent_dir, paste0(paths$prefix, "_top_recurrent_genes.pdf"))
  svg_path <- file.path(paths$recurrent_dir, paste0(paths$prefix, "_top_recurrent_genes.svg"))
  status <- "completed"
  empty_reason <- ""
  if (nrow(summary) > 0) {
    plot <- build_plot(summary, cfg)
    ggplot2::ggsave(png_path, plot, width = 8.8, height = 7.2, dpi = 240, bg = "white", limitsize = FALSE)
    ggplot2::ggsave(pdf_path, plot, width = 8.8, height = 7.2, bg = "white", limitsize = FALSE)
    save_plot_svg(svg_path, plot, width = 8.8, height = 7.2)
  } else {
    status <- "empty"
    empty_reason <- "Gene-category contribution table is empty."
    plot <- build_empty_plot(cfg, empty_reason)
    ggplot2::ggsave(png_path, plot, width = 8.8, height = 4.2, dpi = 240, bg = "white", limitsize = FALSE)
    ggplot2::ggsave(pdf_path, plot, width = 8.8, height = 4.2, bg = "white", limitsize = FALSE)
    save_plot_svg(svg_path, plot, width = 8.8, height = 4.2)
  }

  lisaR:::lisa_write_native_extended_recipe(summary, cfg, png_path,
    "build_single_de_recurrent_gene_screen.R", "native_recurrent_genes")

  readme_path <- file.path(paths$recurrent_dir, "README_recurrent_gene_screen.txt")
  write_readme(readme_path, cfg)

  index <- data.frame(
    analysis_id = cfg$analysis_id,
    universe = cfg$universe,
    n_genes_total = nrow(summary),
    n_recurrent_genes = nrow(recurrent_only),
    n_genes_visible = n_recurrent_visible,
    n_genes_available_for_plot = nrow(recurrent_only),
    gene_visibility = sprintf("%s/%s", n_recurrent_visible, nrow(recurrent_only)),
    status = status,
    empty_reason = empty_reason,
    min_recurrent_categories = cfg$min_recurrent_categories,
    summary_tsv = summary_path,
    top_recurrent_tsv = top_path,
    plot_png = png_path,
    plot_pdf = pdf_path,
    plot_svg = svg_path,
    readme = readme_path,
    stringsAsFactors = FALSE
  )
  index_path <- file.path(paths$recurrent_dir, paste0(paths$prefix, "_recurrent_gene_screen_index.tsv"))
  write_tsv(index, index_path)
  print(index[, c("analysis_id", "universe", "n_genes_total", "n_recurrent_genes", "gene_visibility")])
}

if (sys.nframe() == 0L) main()
