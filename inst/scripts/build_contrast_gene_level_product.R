#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


options(stringsAsFactors = FALSE)

default_project_dir <- NA_character_
script_file <- sub("^--file=", "", commandArgs(FALSE)[grepl("^--file=", commandArgs(FALSE))][1])
script_dir <- dirname(normalizePath(script_file))

required_pkgs <- c("ggplot2", "openxlsx")
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
    top_genes_per_category = 28,
    lfc_cap = 1.5,
    # --- plot-free preparation and exact product selection (H3) ------------
    # All five are opt-in. With none of them this script behaves exactly as it
    # always has: prepare the paired table and the category summary, draw every
    # selected/default card, and draw the native paired heatmap.
    #
    #   --tables-only        prepare the paired evidence and the complete ranked
    #                        category summary and stop. No plot function is
    #                        called and no figure job is submitted. This is the
    #                        shared-preparation mode; enumerating what could be
    #                        drawn must not draw anything.
    #   --exact-product      draw exactly ONE product:
    #                        `contrast_gene_card` (one named category) or
    #                        `contrast_paired_heatmap` (the collection-wide
    #                        heatmap). The other product is not drawn and no
    #                        sibling output is removed.
    #   --category-id        the exact category for `contrast_gene_card`.
    #   --paired-input       reuse an already prepared paired evidence table
    #   --summary-input      and its complete ranked category summary instead of
    #                        recomputing them, so the card, the heatmap, the
    #                        network and the map all read ONE preparation.
    exact_product = NA_character_,
    category_id = NA_character_,
    paired_input = NA_character_,
    summary_input = NA_character_,
    tables_only = "false",
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
    if (key %in% c("top_categories", "top_genes_per_category")) {
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
  if (!out$tables_only %in% c("true", "false")) {
    stop("--tables-only must be true or false.", call. = FALSE)
  }
  out$tables_only <- identical(out$tables_only, "true")
  exact <- if (is.na(out$exact_product)) "" else out$exact_product
  if (nzchar(exact) && !exact %in% c("contrast_gene_card", "contrast_paired_heatmap")) {
    stop(sprintf("LISA-CONTRAST-040 --exact-product must be contrast_gene_card or contrast_paired_heatmap; received: %s.", exact), call. = FALSE)
  }
  if (nzchar(exact) && out$tables_only) {
    stop("LISA-CONTRAST-041 --exact-product and --tables-only true are mutually exclusive: one draws a figure, the other refuses to draw anything.", call. = FALSE)
  }
  has_category <- !is.na(out$category_id) && nzchar(out$category_id)
  if (identical(exact, "contrast_gene_card") && !has_category) {
    stop("LISA-CONTRAST-042 --exact-product contrast_gene_card needs --category-id naming exactly one category.", call. = FALSE)
  }
  if (has_category && !identical(exact, "contrast_gene_card")) {
    stop("LISA-CONTRAST-043 --category-id applies only to --exact-product contrast_gene_card; the paired heatmap is one collection-wide figure.", call. = FALSE)
  }
  supplied <- c(!is.na(out$paired_input) && nzchar(out$paired_input),
                !is.na(out$summary_input) && nzchar(out$summary_input))
  if (any(supplied) && !all(supplied)) {
    stop("LISA-CONTRAST-044 --paired-input and --summary-input must be given together; a paired table without its ranked summary is not a preparation.", call. = FALSE)
  }
  out$exact_product <- exact
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

collapse_unique <- function(x, sep = ";") {
  x <- sort(unique(safe_chr(x)))
  x <- x[x != ""]
  paste(x, collapse = sep)
}

first_nonempty <- function(x) {
  x <- safe_chr(x)
  x <- x[x != ""]
  if (length(x) == 0) "" else x[[1]]
}

max_or_na <- function(x) {
  x <- safe_num(x)
  out <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(out)) NA_real_ else out
}

min_or_na <- function(x) {
  x <- safe_num(x)
  out <- suppressWarnings(min(x, na.rm = TRUE))
  if (!is.finite(out)) NA_real_ else out
}

mean_or_na <- function(x) {
  x <- safe_num(x)
  out <- suppressWarnings(mean(x, na.rm = TRUE))
  if (!is.finite(out)) NA_real_ else out
}

sign_class <- function(x) {
  x <- safe_num(x)
  ifelse(is.na(x) | x == 0, "zero_or_missing", ifelse(x > 0, "up", "down"))
}

short_label <- function(x, max_chars = 58) {
  x <- safe_chr(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

gene_class_display <- function(class, label_a, label_b) {
  class <- safe_chr(class)
  label_a <- first_nonempty(label_a)
  label_b <- first_nonempty(label_b)
  out <- gsub("_", " ", class)
  out[class == "A_specific"] <- sprintf("%s only", label_a)
  out[class == "B_specific"] <- sprintf("%s only", label_b)
  out[class == "shared_opposite_direction"] <- "shared opposite direction"
  out[class == "shared_same_direction"] <- "shared same direction"
  out[class == "shared_mixed_or_missing_direction"] <- "shared mixed/missing direction"
  out
}

contrast_row <- function(project_dir, contrast_id) {
  idx <- read_tsv(file.path(project_dir, "config", "contrast_index.tsv"))
  row <- idx[idx$contrast_id == contrast_id | paste(idx$contrast_id, idx$output_id, sep = "_") == contrast_id, , drop = FALSE]
  if (nrow(row) != 1) stop(sprintf("Contrast id did not resolve to one row: %s", contrast_id), call. = FALSE)
  row
}

find_contrast_table <- function(contrast_dir, output_id) {
  table_dir <- file.path(contrast_dir, "lisa_tables")
  if (!dir.exists(table_dir)) return("")
  candidates <- list.files(
    table_dir,
    pattern = "_GSEA_category_contrast\\.tsv$",
    full.names = TRUE
  )
  candidates <- candidates[startsWith(basename(candidates), paste0(output_id, "_"))]
  if (length(candidates) > 1) {
    stop(sprintf(
      "Multiple GSEA category contrast tables found for %s:\n%s",
      output_id, paste(candidates, collapse = "\n")
    ), call. = FALSE)
  }
  if (length(candidates) == 0) "" else candidates[[1]]
}

paths_for <- function(project_dir, row, universe) {
  contrast_name <- paste(row$contrast_id, row$output_id, sep = "_")
  contrast_dir <- file.path(project_dir, "outputs", "category_contrasts", contrast_name, paste0("collection_", universe))
  contrast_table <- find_contrast_table(contrast_dir, row$output_id)
  single_dir_a <- file.path(project_dir, "outputs", "gene_level", "single_de", row$contrast_a, paste0("collection_", universe))
  single_dir_b <- file.path(project_dir, "outputs", "gene_level", "single_de", row$contrast_b, paste0("collection_", universe))
  prefix_a <- paste(row$contrast_a, universe, "gene_level", sep = "_")
  prefix_b <- paste(row$contrast_b, universe, "gene_level", sep = "_")
  out_dir <- file.path(project_dir, "outputs", "gene_level", "category_contrasts", contrast_name, paste0("collection_", universe))
  prefix <- paste(row$output_id, universe, "contrast_gene_level", sep = "_")
  list(
    contrast_name = contrast_name,
    contrast_table = contrast_table,
    evidence_a = file.path(single_dir_a, paste0(prefix_a, "_gene_category_contributions.tsv")),
    evidence_b = file.path(single_dir_b, paste0(prefix_b, "_gene_category_contributions.tsv")),
    out_dir = out_dir,
    cards_dir = file.path(out_dir, "contrast_category_cards"),
    prefix = prefix
  )
}

require_files <- function(paths) {
  need <- unlist(paths)
  missing <- need[!file.exists(need)]
  if (length(missing) > 0) {
    stop(sprintf("Missing required input files:\n%s", paste(missing, collapse = "\n")), call. = FALSE)
  }
}

aggregate_evidence <- function(df, side_label) {
  common_cols <- c("category_id", "symbol", "category_display_name", "macrogroup_id", "macrogroup_name", "color")
  side_cols <- paste0(c(
    "gene_id", "log2FC", "padj", "de_direction", "de_is_significant",
    "in_gsea_leading_edge", "in_ora_overlap", "source_geneset_count",
    "gene_contribution_score", "min_source_padj", "mean_source_NES",
    "source_evidence", "source_pathways"
  ), "_", side_label)
  if (nrow(df) == 0) {
    return(setNames(data.frame(matrix(nrow = 0, ncol = length(c(common_cols, side_cols)))), c(common_cols, side_cols)))
  }
  df$symbol <- toupper(safe_chr(df$symbol))
  df$category_id <- safe_chr(df$category_id)
  df <- df[df$symbol != "" & df$category_id != "", , drop = FALSE]
  if (nrow(df) == 0) {
    return(setNames(data.frame(matrix(nrow = 0, ncol = length(c(common_cols, side_cols)))), c(common_cols, side_cols)))
  }
  keys <- paste(df$category_id, df$symbol, sep = "\r")
  split_rows <- split(df, keys)
  rows <- lapply(split_rows, function(x) {
    data.frame(
      category_id = first_nonempty(x$category_id),
      symbol = first_nonempty(x$symbol),
      category_display_name = first_nonempty(x$category_display_name),
      macrogroup_id = first_nonempty(x$macrogroup_id),
      macrogroup_name = first_nonempty(x$macrogroup_name),
      color = first_nonempty(x$color),
      gene_id = first_nonempty(x$gene_id),
      log2FC = mean_or_na(x$log2FC),
      padj = min_or_na(x$padj),
      de_direction = first_nonempty(x$de_direction),
      de_is_significant = any(as.logical(x$de_is_significant) %in% TRUE),
      in_gsea_leading_edge = any(as.logical(x$in_gsea_leading_edge) %in% TRUE),
      in_ora_overlap = any(as.logical(x$in_ora_overlap) %in% TRUE),
      source_geneset_count = max_or_na(x$source_geneset_count),
      gene_contribution_score = max_or_na(x$gene_contribution_score),
      min_source_padj = min_or_na(x$min_source_padj),
      mean_source_NES = mean_or_na(x$mean_source_NES),
      source_evidence = collapse_unique(x$source_evidence),
      source_pathways = collapse_unique(x$source_pathways),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  names(out)[!names(out) %in% c("category_id", "symbol", "category_display_name", "macrogroup_id", "macrogroup_name", "color")] <-
    paste0(names(out)[!names(out) %in% c("category_id", "symbol", "category_display_name", "macrogroup_id", "macrogroup_name", "color")], "_", side_label)
  out
}

build_paired_table <- function(a, b, contrast, row) {
  paired <- merge(
    a,
    b,
    by = c("category_id", "symbol", "category_display_name", "macrogroup_id", "macrogroup_name", "color"),
    all = TRUE,
    sort = FALSE
  )
  if (nrow(paired) == 0) {
    out <- data.frame(
      category_id = character(),
      symbol = character(),
      category_display_name = character(),
      macrogroup_id = character(),
      macrogroup_name = character(),
      color = character(),
      contrast_id = character(),
      output_id = character(),
      analysis_A = character(),
      analysis_B = character(),
      label_A = character(),
      label_B = character(),
      present_A = logical(),
      present_B = logical(),
      contrast_gene_class = character(),
      contrast_gene_class_display = character(),
      interpretation_note = character(),
      stringsAsFactors = FALSE
    )
    return(out)
  }
  paired$present_A <- !is.na(paired$gene_contribution_score_A)
  paired$present_B <- !is.na(paired$gene_contribution_score_B)
  paired$log2FC_A_num <- safe_num(paired$log2FC_A)
  paired$log2FC_B_num <- safe_num(paired$log2FC_B)
  paired$direction_A_gene <- sign_class(paired$log2FC_A_num)
  paired$direction_B_gene <- sign_class(paired$log2FC_B_num)
  paired$contrast_gene_class <- ifelse(
    paired$present_A & paired$present_B & paired$direction_A_gene == paired$direction_B_gene & paired$direction_A_gene != "zero_or_missing",
    "shared_same_direction",
    ifelse(
      paired$present_A & paired$present_B & paired$direction_A_gene != paired$direction_B_gene &
        !paired$direction_A_gene %in% "zero_or_missing" & !paired$direction_B_gene %in% "zero_or_missing",
      "shared_opposite_direction",
      ifelse(
        paired$present_A & paired$present_B,
        "shared_mixed_or_missing_direction",
        ifelse(paired$present_A, "A_specific", "B_specific")
      )
    )
  )
  paired$max_gene_contribution_score <- pmax(safe_num(paired$gene_contribution_score_A), safe_num(paired$gene_contribution_score_B), na.rm = TRUE)
  paired$max_gene_contribution_score[!is.finite(paired$max_gene_contribution_score)] <- NA_real_
  paired$abs_log2FC_A <- abs(paired$log2FC_A_num)
  paired$abs_log2FC_B <- abs(paired$log2FC_B_num)
  paired$max_abs_log2FC <- pmax(paired$abs_log2FC_A, paired$abs_log2FC_B, na.rm = TRUE)
  paired$max_abs_log2FC[!is.finite(paired$max_abs_log2FC)] <- NA_real_
  paired <- merge(
    contrast[, c(
      "category_id", "contrast_a_label", "contrast_b_label", "direction_A", "direction_B",
      "direction_class", "is_same_direction", "is_opposite_direction", "mean_NES_A",
      "mean_NES_B", "min_padj_A", "min_padj_B", "abs_delta_mean_NES", "mean_abs_NES"
    )],
    paired,
    by = "category_id",
    all.y = TRUE,
    sort = FALSE
  )
  paired$contrast_id <- row$contrast_id
  paired$output_id <- row$output_id
  paired$analysis_A <- row$contrast_a
  paired$analysis_B <- row$contrast_b
  paired$label_A <- row$label_a
  paired$label_B <- row$label_b
  paired$contrast_gene_class_display <- gene_class_display(paired$contrast_gene_class, paired$label_A, paired$label_B)
  paired$interpretation_note <- "Paired post-processing of single-DE gene-level evidence; gene class is a prioritization label, not causal evidence."
  paired <- paired[order(
    paired$category_id,
    match(paired$contrast_gene_class, c("shared_opposite_direction", "shared_same_direction", "A_specific", "B_specific", "shared_mixed_or_missing_direction")),
    -safe_num(paired$max_gene_contribution_score),
    -safe_num(paired$max_abs_log2FC),
    paired$symbol
  ), , drop = FALSE]
  rownames(paired) <- NULL
  paired
}

top_symbols <- function(df, class_name, n = 8) {
  x <- df[df$contrast_gene_class == class_name, , drop = FALSE]
  if (nrow(x) == 0) return("")
  x <- x[order(-safe_num(x$max_gene_contribution_score), -safe_num(x$max_abs_log2FC), x$symbol), , drop = FALSE]
  paste(head(unique(x$symbol), n), collapse = ", ")
}

build_category_summary <- function(paired) {
  if (nrow(paired) == 0) {
    return(data.frame(
      contrast_category_rank = integer(),
      contrast_id = character(),
      output_id = character(),
      universe = character(),
      category_id = character(),
      category_display_name = character(),
      macrogroup_id = character(),
      macrogroup_name = character(),
      direction_A = character(),
      direction_B = character(),
      direction_class = character(),
      mean_NES_A = numeric(),
      mean_NES_B = numeric(),
      min_padj_A = numeric(),
      min_padj_B = numeric(),
      abs_delta_mean_NES = numeric(),
      mean_abs_NES = numeric(),
      n_genes_total = integer(),
      n_present_A = integer(),
      n_present_B = integer(),
      n_shared_any = integer(),
      n_shared_same_direction = integer(),
      n_shared_opposite_direction = integer(),
      n_A_specific = integer(),
      n_B_specific = integer(),
      n_shared_mixed_or_missing_direction = integer(),
      top_shared_opposite_genes = character(),
      top_shared_same_genes = character(),
      top_A_specific_genes = character(),
      top_B_specific_genes = character(),
      interpretation_tag = character(),
      gene_evidence_balance_tag = character(),
      stringsAsFactors = FALSE
    ))
  }
  split_rows <- split(paired, paired$category_id)
  rows <- lapply(split_rows, function(df) {
    data.frame(
      contrast_id = first_nonempty(df$contrast_id),
      output_id = first_nonempty(df$output_id),
      universe = first_nonempty(df$universe),
      category_id = first_nonempty(df$category_id),
      category_display_name = first_nonempty(df$category_display_name),
      macrogroup_id = first_nonempty(df$macrogroup_id),
      macrogroup_name = first_nonempty(df$macrogroup_name),
      direction_A = first_nonempty(df$direction_A),
      direction_B = first_nonempty(df$direction_B),
      direction_class = first_nonempty(df$direction_class),
      mean_NES_A = mean_or_na(df$mean_NES_A),
      mean_NES_B = mean_or_na(df$mean_NES_B),
      min_padj_A = min_or_na(df$min_padj_A),
      min_padj_B = min_or_na(df$min_padj_B),
      abs_delta_mean_NES = max_or_na(df$abs_delta_mean_NES),
      mean_abs_NES = max_or_na(df$mean_abs_NES),
      n_genes_total = length(unique(df$symbol)),
      n_present_A = sum(df$present_A),
      n_present_B = sum(df$present_B),
      n_shared_any = sum(df$present_A & df$present_B),
      n_shared_same_direction = sum(df$contrast_gene_class == "shared_same_direction"),
      n_shared_opposite_direction = sum(df$contrast_gene_class == "shared_opposite_direction"),
      n_A_specific = sum(df$contrast_gene_class == "A_specific"),
      n_B_specific = sum(df$contrast_gene_class == "B_specific"),
      n_shared_mixed_or_missing_direction = sum(df$contrast_gene_class == "shared_mixed_or_missing_direction"),
      top_shared_opposite_genes = top_symbols(df, "shared_opposite_direction"),
      top_shared_same_genes = top_symbols(df, "shared_same_direction"),
      top_A_specific_genes = top_symbols(df, "A_specific"),
      top_B_specific_genes = top_symbols(df, "B_specific"),
      interpretation_tag = first_nonempty(make_interpretation_tag(df)),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(-safe_num(out$abs_delta_mean_NES), -safe_num(out$n_shared_opposite_direction), -safe_num(out$n_genes_total), out$category_id), , drop = FALSE]
  rownames(out) <- NULL
  out$contrast_category_rank <- seq_len(nrow(out))
  min_side <- pmin(out$n_present_A, out$n_present_B)
  max_side <- pmax(out$n_present_A, out$n_present_B)
  out$gene_evidence_balance_tag <- ifelse(
    max_side == 0,
    "no_gene_level_evidence",
    ifelse(min_side == 0, "one_sided_gene_evidence",
      ifelse(min_side / max_side < 0.15, "strongly_unbalanced_gene_evidence", "balanced_or_moderately_unbalanced_gene_evidence")
    )
  )
  out[, c("contrast_category_rank", setdiff(names(out), "contrast_category_rank")), drop = FALSE]
}

make_interpretation_tag <- function(df) {
  dc <- first_nonempty(df$direction_class)
  n_opp <- sum(df$contrast_gene_class == "shared_opposite_direction")
  n_same <- sum(df$contrast_gene_class == "shared_same_direction")
  n_a <- sum(df$contrast_gene_class == "A_specific")
  n_b <- sum(df$contrast_gene_class == "B_specific")
  if (grepl("opposite", dc) && n_opp > 0) return("reciprocal_category_with_opposite_gene_support")
  if (grepl("opposite", dc)) return("reciprocal_category_without_shared_opposite_gene_support")
  if (grepl("same", dc) && n_same > 0) return("conserved_category_with_shared_gene_support")
  if (n_a > 0 && n_b == 0) return("A_context_specific_gene_support")
  if (n_b > 0 && n_a == 0) return("B_context_specific_gene_support")
  "mixed_or_weak_gene_support"
}

select_top_category_ids <- function(summary, n) {
  if (n <= 0) return(summary$category_id)
  head(summary$category_id, n)
}

select_card_genes <- function(df, n) {
  df <- df[order(
    match(df$contrast_gene_class, c("shared_opposite_direction", "shared_same_direction", "A_specific", "B_specific", "shared_mixed_or_missing_direction")),
    -safe_num(df$max_gene_contribution_score),
    -safe_num(df$max_abs_log2FC),
    df$symbol
  ), , drop = FALSE]
  df <- df[!duplicated(df$symbol), , drop = FALSE]
  if (n <= 0) df else head(df, n)
}

plot_category_card <- function(df, out_png, out_pdf, cfg) {
  df <- select_card_genes(df, cfg$top_genes_per_category)
  if (nrow(df) == 0) return(FALSE)
  df$symbol <- factor(df$symbol, levels = rev(df$symbol))
  long <- rbind(
    data.frame(symbol = df$symbol, analysis = "A", log2FC = safe_num(df$log2FC_A), gene_class = df$contrast_gene_class, stringsAsFactors = FALSE),
    data.frame(symbol = df$symbol, analysis = "B", log2FC = safe_num(df$log2FC_B), gene_class = df$contrast_gene_class, stringsAsFactors = FALSE)
  )
  long$analysis <- factor(long$analysis, levels = unique(long$analysis))
  long$symbol <- factor(long$symbol, levels = levels(df$symbol))
  long$text_color <- ifelse(!is.na(long$log2FC) & abs(long$log2FC) >= cfg$lfc_cap * 0.72, "white", "grey12")
  title <- sprintf("%s | %s", first_nonempty(df$category_display_name), first_nonempty(df$direction_class))
  subtitle <- sprintf("A = %s; B = %s; gene class shown at right", first_nonempty(df$label_A), first_nonempty(df$label_B))
  class_df <- df
  class_df$x <- 2.54
  class_df$class_label <- gene_class_display(class_df$contrast_gene_class, class_df$label_A, class_df$label_B)
  p <- ggplot2::ggplot(long, ggplot2::aes(x = analysis, y = symbol, fill = pmax(pmin(log2FC, cfg$lfc_cap), -cfg$lfc_cap))) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35, width = 0.92, height = 0.92) +
    ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(log2FC), "", sprintf("%.2f", log2FC)), color = text_color), size = 2.5) +
    ggplot2::scale_color_identity() +
    ggplot2::geom_text(data = class_df, ggplot2::aes(x = x, y = symbol, label = class_label), inherit.aes = FALSE, hjust = 0, size = 2.45, color = "grey20") +
    ggplot2::coord_cartesian(xlim = c(0.5, 3.55), clip = "off") +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-cfg$lfc_cap, cfg$lfc_cap), breaks = c(-cfg$lfc_cap, 0, cfg$lfc_cap),
      na.value = "grey92", name = "log2FC"
    ) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 9),
      axis.text.y = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 8.5, color = "grey30"),
      plot.margin = ggplot2::margin(8, 100, 8, 8)
  )
  h <- max(4.5, 1.0 + 0.22 * nrow(df))
  ggplot2::ggsave(out_png, p, width = 7.1, height = h, dpi = 180, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(out_pdf, p, width = 7.1, height = h, bg = "white", limitsize = FALSE)
  save_svg_twin(out_png, p, 7.1, h)
  TRUE
}

plot_combined_heatmap <- function(paired, category_ids, out_png, out_pdf, cfg) {
  rows <- lapply(category_ids, function(cid) {
    df <- paired[paired$category_id == cid, , drop = FALSE]
    df <- select_card_genes(df, max(8, floor(cfg$top_genes_per_category / 2)))
    if (nrow(df) == 0) return(NULL)
    df$row_label <- paste(short_label(first_nonempty(df$category_display_name), 22), df$symbol, sep = " | ")
    df
  })
  df <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(df) || nrow(df) == 0) return(FALSE)
  df$row_label <- factor(df$row_label, levels = rev(unique(df$row_label)))
  label_a <- first_nonempty(df$label_A)
  label_b <- first_nonempty(df$label_B)
  long <- rbind(
    data.frame(row_label = df$row_label, analysis = "A", log2FC = safe_num(df$log2FC_A), gene_class = df$contrast_gene_class, stringsAsFactors = FALSE),
    data.frame(row_label = df$row_label, analysis = "B", log2FC = safe_num(df$log2FC_B), gene_class = df$contrast_gene_class, stringsAsFactors = FALSE)
  )
  long$analysis <- factor(long$analysis, levels = unique(long$analysis))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = analysis, y = row_label, fill = pmax(pmin(log2FC, cfg$lfc_cap), -cfg$lfc_cap))) +
    ggplot2::geom_tile(color = "white", linewidth = 0.18, width = 0.92, height = 0.92) +
    ggplot2::coord_fixed(ratio = 1, clip = "off") +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-cfg$lfc_cap, cfg$lfc_cap), na.value = "grey92", name = "log2FC") +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(add = 0.08)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = 0.08)) +
    ggplot2::labs(
      title = "Paired heatmap",
      subtitle = sprintf("A = %s\nB = %s", label_a, label_b),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 8.5) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 8, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 5.5),
      plot.title = ggplot2::element_text(face = "bold", size = 9.5),
      plot.subtitle = ggplot2::element_text(size = 7.5, color = "grey30"),
      legend.position = "none",
      plot.margin = ggplot2::margin(6, 6, 6, 6)
    )
  h <- max(5.8, 0.125 * nrow(df))
  ggplot2::ggsave(out_png, p, width = 3.9, height = h, dpi = 180, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(out_pdf, p, width = 3.9, height = h, bg = "white", limitsize = FALSE)
  save_svg_twin(out_png, p, 3.9, h)
  TRUE
}

write_readme <- function(paths, cfg, row, summary, card_index) {
  n_a <- if ("n_present_A" %in% names(summary)) sum(summary$n_present_A, na.rm = TRUE) else NA_integer_
  n_b <- if ("n_present_B" %in% names(summary)) sum(summary$n_present_B, na.rm = TRUE) else NA_integer_
  n_shared <- if ("n_shared_any" %in% names(summary)) sum(summary$n_shared_any, na.rm = TRUE) else NA_integer_
  balance_note <- c(
    "",
    "Gene-evidence balance:",
    sprintf("- Category-level present rows in A: %s", n_a),
    sprintf("- Category-level present rows in B: %s", n_b),
    sprintf("- Shared gene-category rows: %s", n_shared)
  )
  if (isTRUE(n_shared == 0) && (isTRUE(n_a > 0) || isTRUE(n_b > 0))) {
    balance_note <- c(
      balance_note,
      "- Warning: no shared gene-category rows were detected. Treat figures as one-sided diagnostic summaries, not as paired shared-gene support."
    )
  }
  cross_species_note <- if (is_cross_species_contrast(row)) {
    c(
      "",
      "Cross-species caution:",
      "This contrast includes B16 murine and external human material. Shared gene symbols are useful for triage but are not an explicit orthology validation."
    )
  } else {
    character()
  }
  lines <- c(
    "Contrast gene-level product",
    "",
    sprintf("Contrast: %s_%s", row$contrast_id, row$output_id),
    sprintf("Universe: %s", cfg$universe),
    sprintf("Analysis A: %s (%s)", row$contrast_a, row$label_a),
    sprintf("Analysis B: %s (%s)", row$contrast_b, row$label_b),
    "",
    "This output pairs the single-DE gene-category evidence from both sides of an existing LISA category contrast.",
    "It does not rerun DE, GSEA, ORA, or LISA.",
    "",
    "Gene classes:",
    "- shared_opposite_direction: present in both analyses with opposite log2FC signs.",
    "- shared_same_direction: present in both analyses with the same log2FC sign.",
    sprintf("- A_specific / B_specific: present in only one side's gene-level evidence; A = %s, B = %s.", row$label_a, row$label_b),
    "- shared_mixed_or_missing_direction: present in both but direction is zero/missing.",
    "",
    "Files:",
    sprintf("- %s_paired_gene_evidence.tsv: all paired gene-category rows.", paths$prefix),
    sprintf("- %s_contrast_category_gene_summary.tsv: one row per category, report-ready.", paths$prefix),
    sprintf("- %s_paired_gene_heatmap.png/.pdf: compact heatmap for representative contrast categories.", paths$prefix),
    "- contrast_category_cards/: one card per contrast category.",
    "- contrast_category_cards/contrast_category_cards_index.tsv: card labels and figure paths.",
    "",
    sprintf("Categories summarized: %s", nrow(summary)),
    sprintf("Cards generated: %s", nrow(card_index)),
    balance_note,
    cross_species_note,
    "",
    "Interpretation caution: paired gene classes are prioritization labels for explaining semantic category contrasts; they are not causal driver calls."
  )
  writeLines(lines, file.path(paths$out_dir, "README_contrast_gene_level.txt"))
}

is_cross_species_contrast <- function(row) {
  a <- paste(row$contrast_a, row$label_a)
  b <- paste(row$contrast_b, row$label_b)
  has_b16 <- grepl("B16", a, ignore.case = TRUE) || grepl("B16", b, ignore.case = TRUE)
  has_external_human <- grepl("Huh7|MS751|Qiu|human", a, ignore.case = TRUE) || grepl("Huh7|MS751|Qiu|human", b, ignore.case = TRUE)
  has_b16 && has_external_human
}

# --- preparation, split from plotting (H3) ----------------------------------
#
# H3_SCOPE section 4: "split the current paired-table/category-summary
# preparation from plotting". This is the split. The scientific logic is the
# accepted code, moved and not rewritten: the same `aggregate_evidence()`,
# `build_paired_table()` and `build_category_summary()` calls in the same order
# with the same arguments. Nothing here draws.
#
# When a verified prepared pair is supplied the tables are read instead of
# recomputed, which is what lets one preparation serve the card, the heatmap,
# the network and the map without any of them recomputing the evidence.
prepare_contrast_tables <- function(cfg, row, paths) {
  if (!is.na(cfg$paired_input) && nzchar(cfg$paired_input)) {
    require_files(list(paired = cfg$paired_input, summary = cfg$summary_input))
    return(list(paired = read_tsv(cfg$paired_input),
                summary = read_tsv(cfg$summary_input),
                reused = TRUE))
  }
  require_files(paths[c("contrast_table", "evidence_a", "evidence_b")])
  contrast <- read_tsv(paths$contrast_table)
  evidence_a <- aggregate_evidence(read_tsv(paths$evidence_a), "A")
  evidence_b <- aggregate_evidence(read_tsv(paths$evidence_b), "B")
  paired <- build_paired_table(evidence_a, evidence_b, contrast, row)
  paired$universe <- rep(cfg$universe, nrow(paired))
  list(paired = paired, summary = build_category_summary(paired), reused = FALSE)
}

# The heatmap's category selection, stated once so the exact and the default
# paths cannot drift. It is taken from the COMPLETE ranked summary -- never from
# a requested or browsed category -- which is the accepted independent top-six
# selection H3_SCOPE requires the exact mode to preserve.
heatmap_category_ids_for <- function(summary) {
  head(summary$category_id, min(6, length(summary$category_id)))
}

category_slug <- function(cid) gsub("[^A-Za-z0-9]+", "_", cid)

install_figure_recipe <- function(png, cfg, code) {
  renderer <- file.path(script_dir, "reproduce_lisa_figure.R")
  if (!file.exists(renderer)) {
    stop(sprintf("%s cannot install contrast reproduction recipe.", code), call. = FALSE)
  }
  lisaR:::lisa_copy_verified_figure_recipe(
    renderer, paste0(tools::file_path_sans_ext(png), "_recipe.R"),
    cfg$lisa_internal_renderer_sha256, run_root = cfg$project_dir
  )
}

# The card's own source contract, extracted verbatim from the default loop so
# the exact single card and the default loop write byte-identical contracts.
write_card_contract <- function(df, png, cfg) {
  selected_df <- select_card_genes(df, cfg$top_genes_per_category)
  df$figure_id <- paste0("contrast_gene_card__", tools::file_path_sans_ext(basename(png)))
  df$figure_type <- "contrast_gene_card"
  df$source_row_order <- seq_len(nrow(df))
  df$selected_for_plot <- df$symbol %in% selected_df$symbol
  df$highlighted <- df$selected_for_plot
  df$labelled <- df$selected_for_plot
  df$plotted_order <- match(df$symbol, selected_df$symbol)
  df$lfc_cap <- cfg$lfc_cap
  write_tsv(df, paste0(tools::file_path_sans_ext(png), "_source.tsv"))
  install_figure_recipe(png, cfg, "LISA-FIGURE-SOURCE-006")
  invisible(df)
}

# The paired heatmap's own source contract, likewise extracted verbatim.
write_heatmap_contract <- function(paired, heatmap_category_ids, png, cfg) {
  heatmap_rows <- lapply(heatmap_category_ids, function(cid) {
    df <- paired[paired$category_id == cid, , drop = FALSE]
    select_card_genes(df, max(8, floor(cfg$top_genes_per_category / 2)))
  })
  heatmap_source <- do.call(rbind, heatmap_rows[vapply(heatmap_rows, nrow, integer(1)) > 0])
  heatmap_source$figure_id <- paste0("contrast_heatmap__", tools::file_path_sans_ext(basename(png)))
  heatmap_source$figure_type <- "contrast_heatmap"
  heatmap_source$source_row_order <- seq_len(nrow(heatmap_source))
  heatmap_source$selected_for_plot <- TRUE
  heatmap_source$highlighted <- TRUE
  heatmap_source$labelled <- TRUE
  heatmap_source$lfc_cap <- cfg$lfc_cap
  write_tsv(heatmap_source, paste0(tools::file_path_sans_ext(png), "_source.tsv"))
  install_figure_recipe(png, cfg, "LISA-FIGURE-SOURCE-007")
  invisible(heatmap_source)
}

# --- exact modes ------------------------------------------------------------
#
# Each draws ONE figure and writes its own scoped inventory naming exactly the
# files it produced. Neither removes any sibling output: in exact mode the
# default sweep would be destructive rather than hygienic, because the sibling
# it would delete is a perfectly valid product of a different request.

exact_gene_card <- function(cfg, paired, summary, paths) {
  cid <- cfg$category_id
  known <- as.character(summary$category_id)
  if (!cid %in% known) {
    stop(sprintf("LISA-CONTRAST-045 category %s is not among the %d categories of this contrast/collection summary.", cid, length(known)), call. = FALSE)
  }
  # The card's rank is its rank in the COMPLETE ranked summary, so the exact
  # card carries the same rank and file stem the default loop would give it.
  rank <- match(cid, known)
  dir.create(paths$cards_dir, recursive = TRUE, showWarnings = FALSE)
  slug <- category_slug(cid)
  png <- file.path(paths$cards_dir, sprintf("%02d_%s_contrast_category_card.png", rank, slug))
  pdf <- file.path(paths$cards_dir, sprintf("%02d_%s_contrast_category_card.pdf", rank, slug))
  df <- paired[paired$category_id == cid, , drop = FALSE]
  n_genes_total <- length(unique(df$symbol))
  if (!isTRUE(plot_category_card(df, png, pdf, cfg))) {
    stop(sprintf("LISA-CONTRAST-046 category %s has no plottable paired gene evidence; no card was drawn.", cid), call. = FALSE)
  }
  write_card_contract(df, png, cfg)
  n_genes_visible <- ifelse(cfg$top_genes_per_category > 0,
    min(cfg$top_genes_per_category, n_genes_total), n_genes_total)
  index <- data.frame(
    rank = rank, category_id = cid,
    category_display_name = first_nonempty(df$category_display_name),
    n_genes_visible = n_genes_visible, n_genes_total = n_genes_total,
    gene_visibility = sprintf("%s/%s", n_genes_visible, n_genes_total),
    output_png = png, output_pdf = pdf, stringsAsFactors = FALSE)
  write_tsv(index, file.path(paths$cards_dir,
    sprintf("%s_contrast_gene_card_index.tsv", slug)))
  message(sprintf("Wrote exactly one contrast gene card: %s", basename(png)))
  invisible(index)
}

exact_paired_heatmap <- function(cfg, paired, summary, paths) {
  heatmap_category_ids <- heatmap_category_ids_for(summary)
  png <- file.path(paths$out_dir, paste0(paths$prefix, "_paired_gene_heatmap.png"))
  pdf <- file.path(paths$out_dir, paste0(paths$prefix, "_paired_gene_heatmap.pdf"))
  if (!isTRUE(plot_combined_heatmap(paired, heatmap_category_ids, png, pdf, cfg))) {
    stop("LISA-CONTRAST-047 this contrast/collection has no plottable paired gene evidence; no heatmap was drawn.", call. = FALSE)
  }
  write_heatmap_contract(paired, heatmap_category_ids, png, cfg)
  index <- data.frame(
    contrast_id = cfg$contrast_id, universe = cfg$universe,
    n_categories_total = nrow(summary),
    n_heatmap_categories = length(heatmap_category_ids),
    heatmap_categories = paste(heatmap_category_ids, collapse = ";"),
    output_png = png, output_pdf = pdf, stringsAsFactors = FALSE)
  write_tsv(index, file.path(paths$out_dir,
    paste0(paths$prefix, "_paired_gene_heatmap_index.tsv")))
  message(sprintf("Wrote exactly one contrast paired heatmap: %s", basename(png)))
  invisible(index)
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
  paths <- paths_for(cfg$project_dir, row, cfg$universe)
  dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- prepare_contrast_tables(cfg, row, paths)
  paired <- tables$paired
  summary <- tables$summary

  # Plot-free preparation. Writes the two tables the card, the heatmap, the
  # network and the map all read, and returns before any plot call.
  if (isTRUE(cfg$tables_only)) {
    write_tsv(paired, file.path(paths$out_dir, paste0(paths$prefix, "_paired_gene_evidence.tsv")))
    write_tsv(summary, file.path(paths$out_dir, paste0(paths$prefix, "_contrast_category_gene_summary.tsv")))
    message(sprintf("Prepared %d paired gene-category rows and %d ranked categories; no figure was drawn.",
      nrow(paired), nrow(summary)))
    return(invisible(NULL))
  }

  # Exactly one product, selected before any plot call. The sibling product's
  # output is neither drawn nor removed.
  if (nzchar(cfg$exact_product)) {
    if (identical(cfg$exact_product, "contrast_gene_card")) {
      return(invisible(exact_gene_card(cfg, paired, summary, paths)))
    }
    return(invisible(exact_paired_heatmap(cfg, paired, summary, paths)))
  }

  # --- default FULL/CLI behaviour, unchanged --------------------------------
  dir.create(paths$cards_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(file.path(paths$cards_dir, "*_contrast_category_card.png"))
  unlink(file.path(paths$cards_dir, "*_contrast_category_card.pdf"))

  category_ids <- select_top_category_ids(summary, cfg$top_categories)
  heatmap_category_ids <- heatmap_category_ids_for(summary)

  write_tsv(paired, file.path(paths$out_dir, paste0(paths$prefix, "_paired_gene_evidence.tsv")))
  write_tsv(summary, file.path(paths$out_dir, paste0(paths$prefix, "_contrast_category_gene_summary.tsv")))
  openxlsx::write.xlsx(
    list(paired_gene_evidence = paired, category_gene_summary = summary),
    file.path(paths$out_dir, paste0(paths$prefix, "_contrast_gene_level_tables.xlsx")),
    overwrite = TRUE
  )

  card_rows <- list()
  for (i in seq_along(category_ids)) {
    cid <- category_ids[[i]]
    df <- paired[paired$category_id == cid, , drop = FALSE]
    n_genes_total <- length(unique(df$symbol))
    n_genes_visible <- ifelse(cfg$top_genes_per_category > 0, min(cfg$top_genes_per_category, n_genes_total), n_genes_total)
    slug <- category_slug(cid)
    png <- file.path(paths$cards_dir, sprintf("%02d_%s_contrast_category_card.png", i, slug))
    pdf <- file.path(paths$cards_dir, sprintf("%02d_%s_contrast_category_card.pdf", i, slug))
    ok <- plot_category_card(df, png, pdf, cfg)
    if (ok) {
      write_card_contract(df, png, cfg)
      card_rows[[length(card_rows) + 1]] <- data.frame(
        rank = i,
        category_id = cid,
        category_display_name = first_nonempty(df$category_display_name),
        n_genes_visible = n_genes_visible,
        n_genes_total = n_genes_total,
        gene_visibility = sprintf("%s/%s", n_genes_visible, n_genes_total),
        png = png,
        pdf = pdf,
        stringsAsFactors = FALSE
      )
    }
  }
  card_index <- if (length(card_rows) > 0) {
    do.call(rbind, card_rows)
  } else {
    data.frame(
      rank = integer(),
      category_id = character(),
      category_display_name = character(),
      png = character(),
      pdf = character(),
      stringsAsFactors = FALSE
    )
  }
  # The directory already identifies the contrast and collection. Repeating
  # those IDs in the index filename exceeds the Windows path budget.
  write_tsv(card_index, file.path(paths$cards_dir, "contrast_category_cards_index.tsv"))

  heatmap_png <- file.path(paths$out_dir, paste0(paths$prefix, "_paired_gene_heatmap.png"))
  heatmap_pdf <- file.path(paths$out_dir, paste0(paths$prefix, "_paired_gene_heatmap.pdf"))
  heatmap_written <- plot_combined_heatmap(paired, heatmap_category_ids, heatmap_png, heatmap_pdf, cfg)
  if (isTRUE(heatmap_written)) {
    write_heatmap_contract(paired, heatmap_category_ids, heatmap_png, cfg)
  }

  manifest <- data.frame(
    contrast_id = row$contrast_id,
    output_id = row$output_id,
    universe = cfg$universe,
    analysis_A = row$contrast_a,
    analysis_B = row$contrast_b,
    n_paired_gene_category_rows = nrow(paired),
    n_categories = nrow(summary),
    n_card_categories = nrow(card_index),
    n_heatmap_categories = length(heatmap_category_ids),
    n_card_genes_visible_total = ifelse(nrow(card_index) > 0, sum(card_index$n_genes_visible), 0),
    n_card_genes_available_total = ifelse(nrow(card_index) > 0, sum(card_index$n_genes_total), 0),
    card_gene_visibility = ifelse(nrow(card_index) > 0, sprintf("%s/%s", sum(card_index$n_genes_visible), sum(card_index$n_genes_total)), "0/0"),
    heatmap_written = heatmap_written,
    n_gene_category_rows_present_A = sum(paired$present_A, na.rm = TRUE),
    n_gene_category_rows_present_B = sum(paired$present_B, na.rm = TRUE),
    n_shared_gene_category_rows = sum(paired$present_A & paired$present_B, na.rm = TRUE),
    cross_species_symbol_caution = is_cross_species_contrast(row),
    output_dir = paths$out_dir,
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(paths$out_dir, paste0(paths$prefix, "_manifest.tsv")))
  write_readme(paths, cfg, row, summary, card_index)
  print(manifest)
}

if (sys.nframe() == 0L) main()
