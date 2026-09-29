# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

lisa_split_leading_edge <- function(value) {
  if (length(value) == 0 || is.na(value) || !nzchar(value)) {
    return(character())
  }
  value <- gsub("^c\\(|\\)$", "", value)
  value <- gsub("[\"']", "", value)
  parts <- unlist(strsplit(value, "[,;|/[:space:]]+"), use.names = FALSE)
  unique(parts[nzchar(parts)])
}

build_lisa_gene_level_tables <- function(
  de_table,
  gsea_table,
  ora_table = NULL,
  output_dir,
  symbol_col = "symbol",
  logfc_col = NULL,
  padj_col = "padj",
  analysis_id = NULL,
  collection = NULL,
  de_padj_cutoff = 0.05,
  lfc_cutoff = 0.58
) {
  if (is.character(de_table) && length(de_table) == 1) {
    de_table <- read_lisa_tsv(de_table)
  }
  if (is.character(gsea_table) && length(gsea_table) == 1) {
    gsea_table <- read_lisa_tsv(gsea_table)
  }
  if (is.character(ora_table) && length(ora_table) == 1) {
    ora_table <- read_lisa_tsv(ora_table, required = FALSE)
  }
  if (is.null(ora_table)) {
    ora_table <- data.frame()
  }

  if (nrow(ora_table) > 0 && "pathway" %in% names(ora_table)) {
    ora_table <- ora_table[as.character(ora_table$pathway) != "pathway", , drop = FALSE]
  }

  lisa_require_columns(de_table, symbol_col, "de_table")
  gsea_required <- c("category_id", "pathway", "padj", "NES", "leadingEdge")
  if (nrow(gsea_table) == 0 && length(names(gsea_table)) == 0) {
    gsea_table <- lisa_empty_df(gsea_required)
  }
  lisa_require_columns(gsea_table, gsea_required, "gsea_table")
  gsea_table <- lisa_classified_gsea_rows(gsea_table)
  lisa_assert_classified_plot_rows(gsea_table, "gene-level GSEA evidence")
  ora_required <- c("category_id", "pathway", "padj", "genes")
  if (nrow(ora_table) == 0 && length(names(ora_table)) == 0) {
    ora_table <- lisa_empty_df(ora_required)
  }
  if (!"category_id" %in% names(ora_table) && "pathway" %in% names(ora_table)) {
    ora_table$category_id <- as.character(ora_table$pathway)
  }
  lisa_require_columns(ora_table, ora_required, "ora_table")

  logfc_col <- logfc_col %||% if ("log2FoldChange" %in% names(de_table)) "log2FoldChange" else "log2FC"
  if (!logfc_col %in% names(de_table)) {
    lisa_structural_error("LISA-EVIDENCE-002", "logfc_col", logfc_col, "required effect column is absent", "a source effect column", "Add the effect column or use a reduced-evidence product.")
  }
  if (!padj_col %in% names(de_table)) {
    lisa_structural_error("LISA-EVIDENCE-002", "padj_col", padj_col, "required significance column is absent", "a source adjusted-p-value column", "Add the significance column or use a reduced-evidence product.")
  }
  if (!"gene_id" %in% names(de_table)) {
    de_table$gene_id <- de_table[[symbol_col]]
  }
  if (!"category_display_name" %in% names(gsea_table)) {
    gsea_table$category_display_name <- gsea_table$category_id
  }
  if (!"macrogroup_id" %in% names(gsea_table)) {
    gsea_table$macrogroup_id <- rep("", nrow(gsea_table))
  }
  if (!"macrogroup_name" %in% names(gsea_table)) {
    gsea_table$macrogroup_name <- rep("", nrow(gsea_table))
  }
  if (!"color" %in% names(gsea_table)) {
    gsea_table$color <- rep("#999999", nrow(gsea_table))
  }
  if (!"source_family" %in% names(gsea_table)) {
    gsea_table$source_family <- rep("", nrow(gsea_table))
  }
  for (col in c("category_display_name", "macrogroup_id", "macrogroup_name", "color", "source_family")) {
    if (!col %in% names(ora_table)) {
      ora_table[[col]] <- rep("", nrow(ora_table))
    }
  }
  if (!"direction" %in% names(ora_table)) {
    ora_table$direction <- rep("", nrow(ora_table))
  }

  de_table[[symbol_col]] <- toupper(as.character(de_table[[symbol_col]]))
  de_symbols <- as.character(de_table[[symbol_col]])
  de_lookup <- split(de_table, de_symbols)

  pathway_rows <- list()
  row_i <- 0L

  add_gene_source_row <- function(gene, source, row) {
    gene_symbol <- toupper(gene)
    de_row <- de_lookup[[gene_symbol]]
    if (is.null(de_row)) {
      logfc <- NA_real_
      gene_padj <- NA_real_
      gene_id <- gene_symbol
    } else {
      logfc <- if (logfc_col %in% names(de_row)) suppressWarnings(as.numeric(de_row[[logfc_col]][1])) else NA_real_
      gene_padj <- if (padj_col %in% names(de_row)) suppressWarnings(as.numeric(de_row[[padj_col]][1])) else NA_real_
      gene_id <- if ("gene_id" %in% names(de_row)) as.character(de_row$gene_id[[1]]) else gene_symbol
    }
    data.frame(
      analysis_id = analysis_id %||% "",
      collection = collection %||% "",
      symbol = gene_symbol,
      gene_id = gene_id,
      category_id = row$category_id[[1]],
      category_display_name = row$category_display_name[[1]],
      macrogroup_id = row$macrogroup_id[[1]],
      macrogroup_name = row$macrogroup_name[[1]],
      color = row$color[[1]],
      pathway = row$pathway[[1]],
      source_family = row$source_family[[1]],
      source_evidence = source$evidence,
      source_pathway = row$pathway[[1]],
      source_padj = suppressWarnings(as.numeric(row$padj[[1]])),
      source_NES = source$nes,
      log2FC = logfc,
      rank_value = if (!is.null(de_row) && "rank_value" %in% names(de_row)) suppressWarnings(as.numeric(de_row$rank_value[[1]])) else logfc,
      padj = gene_padj,
      de_direction = ifelse(is.na(logfc) | logfc == 0, "neutral", ifelse(logfc > 0, "up", "down")),
      de_is_significant = !is.na(gene_padj) & gene_padj <= de_padj_cutoff,
      in_gsea_leading_edge = source$gsea,
      in_ora_overlap = source$ora,
      stringsAsFactors = FALSE
    )
  }

  for (i in seq_len(nrow(gsea_table))) {
    genes <- lisa_split_leading_edge(gsea_table$leadingEdge[[i]])
    if (length(genes) == 0) {
      next
    }
    for (gene in genes) {
      row_i <- row_i + 1L
      pathway_rows[[row_i]] <- add_gene_source_row(
        gene,
        list(evidence = "GSEA_LE", nes = suppressWarnings(as.numeric(gsea_table$NES[[i]])), gsea = TRUE, ora = FALSE),
        gsea_table[i, , drop = FALSE]
      )
    }
  }

  for (i in seq_len(nrow(ora_table))) {
    genes <- lisa_split_leading_edge(ora_table$genes[[i]])
    if (length(genes) == 0) {
      next
    }
    evidence <- if (nzchar(as.character(ora_table$direction[[i]]))) {
      paste0("ORA_", toupper(as.character(ora_table$direction[[i]])))
    } else {
      "ORA_OVERLAP"
    }
    for (gene in genes) {
      row_i <- row_i + 1L
      pathway_rows[[row_i]] <- add_gene_source_row(
        gene,
        list(evidence = evidence, nes = NA_real_, gsea = FALSE, ora = TRUE),
        ora_table[i, , drop = FALSE]
      )
    }
  }

  leading_edge_cols <- c(
    "analysis_id", "collection", "symbol", "gene_id", "category_id",
    "category_display_name", "macrogroup_id", "macrogroup_name", "color",
    "pathway", "source_family", "source_evidence", "source_pathway",
    "source_padj", "source_NES", "log2FC", "rank_value", "padj", "de_direction",
    "de_is_significant", "in_gsea_leading_edge", "in_ora_overlap"
  )
  leading_edge <- if (length(pathway_rows) == 0) {
    lisa_empty_df(leading_edge_cols)
  } else {
    do.call(rbind, pathway_rows)
  }

  contribution_cols <- c(
    "analysis_id", "collection", "symbol", "gene_id", "category_id",
    "category_display_name", "macrogroup_id", "macrogroup_name", "color",
    "log2FC", "rank_value", "padj", "de_direction", "de_is_significant",
    "in_gsea_leading_edge", "in_ora_overlap", "source_geneset_count",
    "gene_contribution_score", "lisa_support_score", "de_magnitude",
    "de_fdr_score", "min_source_padj", "mean_source_NES",
    "source_evidence", "source_pathways", "category_recurrence"
  )
  if (nrow(leading_edge) == 0) {
    gene_category <- lisa_empty_df(contribution_cols)
  } else {
    keys <- paste(leading_edge$symbol, leading_edge$category_id, sep = "\r")
    grouped <- split(leading_edge, keys)
    category_rows <- lapply(grouped, function(x) {
      min_padj <- lisa_min_or_na(x$source_padj)
      mean_nes <- lisa_mean_or_na(x$source_NES)
      mean_nes_factor <- ifelse(is.na(mean_nes), 0, abs(mean_nes))
      lisa_support_score <- length(unique(x$source_pathway)) * (1 + mean_nes_factor) *
        ifelse(is.na(min_padj), 1, -log10(pmax(min_padj, 1e-300)))
      gene_logfc <- if (logfc_col %in% names(x)) suppressWarnings(as.numeric(x[[logfc_col]][1])) else suppressWarnings(as.numeric(x$log2FC[1]))
      gene_padj <- suppressWarnings(as.numeric(x$padj[1]))
      de_magnitude <- ifelse(is.na(gene_logfc), 0, abs(gene_logfc))
      de_fdr_score <- ifelse(is.na(gene_padj), 0, -log10(pmax(gene_padj, 1e-300)))
      score <- sqrt(pmax(lisa_support_score, 0)) * de_magnitude * (1 + de_fdr_score)
      data.frame(
        analysis_id = x$analysis_id[[1]],
        collection = x$collection[[1]],
        symbol = x$symbol[[1]],
        gene_id = x$gene_id[[1]],
        category_id = x$category_id[[1]],
        category_display_name = x$category_display_name[[1]],
        macrogroup_id = x$macrogroup_id[[1]],
        macrogroup_name = x$macrogroup_name[[1]],
        color = x$color[[1]],
        log2FC = x$log2FC[[1]],
        rank_value = x$rank_value[[1]],
        padj = x$padj[[1]],
        de_direction = x$de_direction[[1]],
        de_is_significant = any(as.logical(x$de_is_significant) %in% TRUE),
        in_gsea_leading_edge = any(as.logical(x$in_gsea_leading_edge) %in% TRUE),
        in_ora_overlap = any(as.logical(x$in_ora_overlap) %in% TRUE),
        source_geneset_count = length(unique(x$source_pathway)),
        gene_contribution_score = score,
        lisa_support_score = lisa_support_score,
        de_magnitude = de_magnitude,
        de_fdr_score = de_fdr_score,
        min_source_padj = min_padj,
        mean_source_NES = mean_nes,
        source_evidence = paste(sort(unique(x$source_evidence)), collapse = ";"),
        source_pathways = paste(sort(unique(x$source_pathway)), collapse = ";"),
        category_recurrence = NA_integer_,
        stringsAsFactors = FALSE
      )
    })
    gene_category <- do.call(rbind, category_rows)
    recurrence <- stats::aggregate(category_id ~ symbol, data = unique(gene_category[, c("symbol", "category_id")]), FUN = length)
    names(recurrence)[names(recurrence) == "category_id"] <- "category_recurrence"
    gene_category <- merge(gene_category, recurrence, by = "symbol", all.x = TRUE, sort = FALSE, suffixes = c("", ".calc"))
    gene_category$category_recurrence <- gene_category$category_recurrence.calc
    gene_category$category_recurrence.calc <- NULL
    gene_category <- gene_category[order(gene_category$category_id, gene_category$min_source_padj, -gene_category$source_geneset_count), , drop = FALSE]
  }

  summary_cols <- c(
    "analysis_id", "collection", "category_id", "category_display_name",
    "macrogroup_id", "macrogroup_name", "color", "n_supporting_genes",
    "max_gene_contribution_score", "mean_gene_contribution_score",
    "min_source_padj", "mean_source_NES", "top_supporting_genes"
  )
  if (nrow(gene_category) == 0) {
    category_summary <- lisa_empty_df(summary_cols)
  } else {
    grouped_categories <- split(gene_category, gene_category$category_id)
    category_summary <- do.call(rbind, lapply(grouped_categories, function(x) {
      x <- x[order(-suppressWarnings(as.numeric(x$gene_contribution_score)), x$min_source_padj, x$symbol), , drop = FALSE]
      data.frame(
        analysis_id = x$analysis_id[[1]],
        collection = x$collection[[1]],
        category_id = x$category_id[[1]],
        category_display_name = x$category_display_name[[1]],
        macrogroup_id = x$macrogroup_id[[1]],
        macrogroup_name = x$macrogroup_name[[1]],
        color = x$color[[1]],
        n_supporting_genes = length(unique(x$symbol)),
        max_gene_contribution_score = lisa_max_or_na(x$gene_contribution_score),
        mean_gene_contribution_score = lisa_mean_or_na(x$gene_contribution_score),
        min_source_padj = lisa_min_or_na(x$min_source_padj),
        mean_source_NES = lisa_mean_or_na(x$mean_source_NES),
        top_supporting_genes = paste(head(unique(x$symbol), 20), collapse = ";"),
        stringsAsFactors = FALSE
      )
    }))
    category_summary <- category_summary[order(-category_summary$max_gene_contribution_score, category_summary$min_source_padj, category_summary$category_id), , drop = FALSE]
  }

  lisa_guarded_dir_create(output_dir, output_dir)
  leading_edge_path <- file.path(output_dir, "leading_edge_gene_pathways.tsv")
  gene_category_path <- file.path(output_dir, "gene_category_contributions.tsv")
  category_summary_path <- file.path(output_dir, "category_gene_support_summary.tsv")
  write_lisa_tsv(leading_edge, leading_edge_path)
  write_lisa_tsv(gene_category, gene_category_path)
  write_lisa_tsv(category_summary, category_summary_path)

  if (!is.null(analysis_id) && !is.null(collection)) {
    prefix <- paste(analysis_id, collection, "gene_level", sep = "_")
    write_lisa_tsv(leading_edge, file.path(output_dir, paste0(prefix, "_leading_edge_gene_pathways.tsv")))
    write_lisa_tsv(gene_category, file.path(output_dir, paste0(prefix, "_gene_category_contributions.tsv")))
    write_lisa_tsv(category_summary, file.path(output_dir, paste0(prefix, "_category_gene_support_summary.tsv")))
  }

  invisible(list(
    leading_edge = leading_edge,
    gene_category = gene_category,
    category_summary = category_summary,
    files = c(
      leading_edge_gene_pathways = leading_edge_path,
      gene_category_contributions = gene_category_path,
      category_gene_support_summary = category_summary_path
    )
  ))
}

lisa_min_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- suppressWarnings(min(x, na.rm = TRUE))
  if (!is.finite(out)) NA_real_ else out
}

lisa_max_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(out)) NA_real_ else out
}

lisa_mean_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  out <- suppressWarnings(mean(x, na.rm = TRUE))
  if (!is.finite(out)) NA_real_ else out
}
