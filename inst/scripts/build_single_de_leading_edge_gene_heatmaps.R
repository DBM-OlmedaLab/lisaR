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
source(file.path(dirname(normalizePath(script_file)), "lisa_plot_metadata.R"))
source(file.path(dirname(normalizePath(script_file)), "lisa_leading_edge_heatmap_native.R"))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || (length(x) == 1 && is.na(x))) y else x
}

read_tsv <- function(path) {
  if (!file.exists(path)) stop(sprintf("Missing TSV: %s", path), call. = FALSE)
  if (file.info(path)$size == 0L) return(data.frame())
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
    expression_matrix = "",
    expression_matrix_column = "",
    scale = "zscore",
    top_genes = 30,
    categories = NA_character_,
    formats = "png"
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    if (i == length(args)) stop(sprintf("Missing value for argument: %s", key), call. = FALSE)
    val <- args[[i + 1]]
    key <- gsub("-", "_", sub("^--", "", key))
    if (!key %in% names(out)) stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    if (key %in% c("top_genes")) {
      out[[key]] <- as.integer(val)
    } else {
      out[[key]] <- val
    }
    i <- i + 2
  }
  if (is.na(out$project_dir) || out$project_dir == "") stop("Required argument: --project-dir", call. = FALSE)
  if (is.na(out$analysis_id) || out$analysis_id == "") stop("Required argument: --analysis-id", call. = FALSE)
  if (!out$scale %in% c("zscore", "raw")) stop("--scale must be zscore or raw", call. = FALSE)
  out$formats <- unique(trimws(strsplit(tolower(out$formats), ",", fixed = TRUE)[[1]]))
  if (any(!out$formats %in% c("png", "svg", "pdf"))) stop("--formats must contain png, svg and/or pdf", call. = FALSE)
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
  single_dir <- file.path(project_dir, "outputs", "single_de", analysis_id, paste0("collection_", universe))
  gene_dir <- file.path(project_dir, "outputs", "gene_level", "single_de", analysis_id, paste0("collection_", universe))
  prefix <- paste(analysis_id, universe, "gene_level", sep = "_")
  standardized_de <- first_existing_path(c(
    file.path(single_dir, "inputs", paste0(analysis_id, "_standardized_DE.tsv")),
    file.path(single_dir, paste0(analysis_id, "_standardized_DE.tsv"))
  ))
  evidence <- first_existing_path(c(
    file.path(gene_dir, paste0(prefix, "_gene_category_contributions.tsv")),
    first_glob_path(c(
      file.path(gene_dir, paste0(prefix, "_gene_category*.tsv")),
      file.path(gene_dir, paste0(prefix, "_gene_catego*.tsv"))
    ))
  ))
  list(
    single_dir = single_dir,
    gene_dir = gene_dir,
    standardized_de = standardized_de,
    original_de = file.path(project_dir, "input", paste0(analysis_id, ".tsv")),
    de_index = file.path(project_dir, "config", "de_index.tsv"),
    evidence = evidence,
    out_dir = file.path(gene_dir, "leading_edge_gene_heatmaps"),
    prefix = prefix
  )
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x
}

safe_file_component <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

first_existing_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) "" else hit[[1]]
}

first_nonempty_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  for (col in hit) {
    vals <- trimws(safe_chr(df[[col]]))
    vals <- vals[!is.na(vals) & vals != ""]
    if (length(vals) > 0) return(col)
  }
  ""
}

norm_path <- function(path, base_dir) {
  path <- trimws(as.character(path %||% ""))
  if (path == "") return("")
  lisaR:::lisa_norm_path(path, base_dir)
}

expression_matrix_source <- function(path, source_column = "") {
  row <- data.frame()
  if (nzchar(source_column)) row <- stats::setNames(data.frame(path), source_column)
  lisaR:::lisa_expression_matrix_source(path, row)
}

indexed_expression_matrix_source <- function(path, row, index_path) {
  lisaR:::lisa_expression_matrix_source(path, row, dirname(index_path))
}

discover_matrix_source <- function(paths, cfg) {
  explicit <- trimws(as.character(cfg$expression_matrix %||% ""))
  matrix_columns <- c(
    "expression_matrix", "expression_matrix_path", "counts_matrix", "counts_matrix_path",
    "vst_matrix", "vst_matrix_path", "normalized_matrix", "normalized_matrix_path",
    "tpm_matrix", "tpm_matrix_path", "sample_counts_path", "sample_counts"
  )
  row <- analysis_index_row(paths, cfg$analysis_id)
  if (explicit != "") {
    if (!file.exists(explicit)) stop(sprintf("Expression/count matrix does not exist: %s", explicit), call. = FALSE)
    source_column <- trimws(as.character(
      cfg$expression_matrix_column %||% ""
    ))
    if (nzchar(source_column)) {
      if (nrow(row) != 1L || !source_column %in% matrix_columns ||
          !source_column %in% names(row) ||
          is.na(row[[source_column]][[1L]]) ||
          !nzchar(trimws(as.character(row[[source_column]][[1L]])))) {
        stop("--expression-matrix-column does not identify one saved de_index matrix declaration.",
             call. = FALSE)
      }
      # The orchestrator has staged a manifest-bound byte-identical effective
      # matrix. Rebind every alias of the originally selected declaration to the
      # staged path only for scale/provenance resolution; direct CLI behavior is
      # unchanged when --expression-matrix-column is absent.
      selected_path <- norm_path(row[[source_column]][[1L]],
                                 dirname(paths$de_index))
      rebound <- row
      for (column in intersect(matrix_columns, names(rebound))) {
        declared <- as.character(rebound[[column]][[1L]])
        if (is.na(declared) || !nzchar(trimws(declared))) next
        declared <- norm_path(declared, dirname(paths$de_index))
        if (identical(lisaR:::lisa_path_key(declared),
                      lisaR:::lisa_path_key(selected_path))) {
          rebound[[column]][[1L]] <- explicit
        }
      }
      return(indexed_expression_matrix_source(explicit, rebound,
                                               paths$de_index))
    }
    return(indexed_expression_matrix_source(explicit, row, paths$de_index))
  }
  if (nrow(row) > 0) {
    col <- first_nonempty_col(row, matrix_columns)
    if (col != "") {
      path <- norm_path(row[[col]][[1]], dirname(paths$de_index))
      if (!file.exists(path)) stop(sprintf("de_index %s points to missing expression/count matrix: %s", col, path), call. = FALSE)
      # Preserve file-selection precedence, but inspect every declaration of
      # that file. A generic alias must not hide an explicit VST/count scale.
      return(indexed_expression_matrix_source(path, row, paths$de_index))
    }
  }
  candidates <- c(
    file.path(dirname(paths$original_de), paste0(cfg$analysis_id, ".vst.tsv")),
    file.path(dirname(paths$original_de), paste0(cfg$analysis_id, ".normalized.tsv")),
    file.path(dirname(paths$original_de), paste0(cfg$analysis_id, ".counts.tsv")),
    file.path(dirname(paths$original_de), paste0(cfg$analysis_id, "_vst.tsv")),
    file.path(dirname(paths$original_de), paste0(cfg$analysis_id, "_normalized.tsv")),
    file.path(dirname(paths$original_de), paste0(cfg$analysis_id, "_counts.tsv"))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) return(expression_matrix_source(hit[[1]]))
  expression_matrix_source("")
}

discover_matrix_path <- function(paths, cfg) {
  discover_matrix_source(paths, cfg)$path
}

read_expression_matrix <- function(path, source = expression_matrix_source(path)) {
  mat <- read_tsv(path)
  if (nrow(mat) == 0) stop(sprintf("Expression/count matrix is empty: %s", path), call. = FALSE)
  gene_col <- first_existing_col(mat, c("symbol", "Symbol", "Genes", "gene_name", "gene", "Gene", "external_gene_name", "gene_id", "Geneid"))
  if (gene_col == "") stop(sprintf("Expression/count matrix has no recognizable gene column: %s", path), call. = FALSE)
  mat$symbol <- toupper(safe_chr(mat[[gene_col]]))
  gene_aliases <- c(gene_col, "symbol", "Symbol", "gene_name", "gene", "Gene", "Genes", "gene_id", "Geneid", "external_gene_name")
  de_stat_cols <- c(
    "baseMean", "log2FoldChange", "log2FC", "lfcSE", "stat", "pvalue", "p_value",
    "padj", "FDR", "qvalue", "q_value", "P.Value", "adj.P.Val", "rank_value",
    "score", "NES", "ES", "leadingEdge", "core_enrichment"
  )
  sample_cols <- setdiff(names(mat), unique(c(gene_aliases, de_stat_cols)))
  sample_cols <- sample_cols[vapply(mat[sample_cols], function(x) any(is.finite(safe_num(x))), logical(1))]
  if (length(sample_cols) < 2) {
    stop(
      sprintf(
        "Expression/count matrix must contain at least two numeric sample columns after excluding DE statistics: %s",
        path
      ),
      call. = FALSE
    )
  }
  mat <- mat[mat$symbol != "", c("symbol", sample_cols), drop = FALSE]
  mat <- mat[!duplicated(mat$symbol), , drop = FALSE]
  for (col in sample_cols) mat[[col]] <- safe_num(mat[[col]])
  attr(mat, "sample_cols") <- sample_cols
  attr(mat, "source_path") <- path
  attr(mat, "expression_source") <- source
  mat
}

analysis_index_row <- function(paths, analysis_id) {
  if (!file.exists(paths$de_index)) return(data.frame())
  idx <- read_tsv(paths$de_index)
  if (nrow(idx) == 0 || !"analysis_id" %in% names(idx)) return(data.frame())
  idx[idx$analysis_id == analysis_id, , drop = FALSE]
}

sample_token_regex <- function(token) {
  token <- tolower(as.character(token %||% ""))
  token <- gsub("[^a-z0-9]+", "_", token)
  token <- gsub("_+", "_", token)
  token <- gsub("^_|_$", "", token)
  if (token == "") return(character())
  aliases <- token
  if (grepl("lof_sh2", token, fixed = TRUE)) aliases <- c(aliases, "shfa5_2")
  if (grepl("lof_sh3", token, fixed = TRUE)) aliases <- c(aliases, "shfa5_3")
  if (grepl("shctrl", token, fixed = TRUE)) aliases <- c(aliases, "shctrl")
  if (grepl("gof", token, fixed = TRUE)) aliases <- c(aliases, "f5")
  if (grepl("neg", token, fixed = TRUE)) aliases <- c(aliases, "neg")
  unique(aliases[nzchar(aliases)])
}

infer_analysis_sample_regex <- function(analysis_id) {
  parts <- strsplit(as.character(analysis_id), "_vs_", fixed = TRUE)[[1]]
  if (length(parts) != 2) return("")
  tokens <- unique(c(sample_token_regex(parts[[1]]), sample_token_regex(parts[[2]])))
  tokens <- tokens[nzchar(tokens)]
  if (length(tokens) == 0) return("")
  paste(tokens, collapse = "|")
}

filter_expression_matrix_samples <- function(expr, paths, cfg) {
  sample_cols <- attr(expr, "sample_cols")
  if (length(sample_cols) < 2) return(expr)

  idx_row <- analysis_index_row(paths, cfg$analysis_id)
  include_regex <- ""
  exclude_regex <- ""
  method <- "all_matrix_samples"
  if (nrow(idx_row) > 0) {
    include_col <- first_nonempty_col(idx_row, c("sample_include_regex", "heatmap_sample_include_regex"))
    exclude_col <- first_nonempty_col(idx_row, c("sample_exclude_regex", "heatmap_sample_exclude_regex"))
    if (include_col != "") {
      include_regex <- as.character(idx_row[[include_col]][[1]])
      method <- include_col
    }
    if (exclude_col != "") exclude_regex <- as.character(idx_row[[exclude_col]][[1]])
  }
  if (include_regex == "") {
    include_regex <- infer_analysis_sample_regex(cfg$analysis_id)
    if (include_regex != "") method <- "analysis_id_inference"
  }
  if (include_regex == "") {
    attr(expr, "sample_selection_method") <- method
    attr(expr, "available_sample_cols") <- sample_cols
    return(expr)
  }

  source_path <- attr(expr, "source_path") %||% ""
  expression_source <- attr(expr, "expression_source")
  keep <- grepl(include_regex, sample_cols, ignore.case = TRUE)
  if (exclude_regex != "") keep <- keep & !grepl(exclude_regex, sample_cols, ignore.case = TRUE)
  selected <- sample_cols[keep]
  if (length(selected) < 2) {
    stop(
      sprintf(
        "Sample selection for %s kept fewer than two columns. regex=%s; available=%s",
        cfg$analysis_id,
        include_regex,
        paste(sample_cols, collapse = ",")
      ),
      call. = FALSE
    )
  }
  expr <- expr[, c("symbol", selected), drop = FALSE]
  attr(expr, "sample_cols") <- selected
  attr(expr, "source_path") <- source_path
  attr(expr, "expression_source") <- expression_source
  attr(expr, "available_sample_cols") <- sample_cols
  attr(expr, "sample_selection_method") <- method
  attr(expr, "sample_include_regex") <- include_regex
  attr(expr, "sample_exclude_regex") <- exclude_regex
  expr
}

load_de_context <- function(paths) {
  de <- read_tsv(paths$standardized_de)
  de$symbol <- toupper(safe_chr(de$symbol))
  de$log2FC <- safe_num(de$log2FoldChange)
  de$padj <- safe_num(de$padj)
  de
}

zscore_vec <- function(x) {
  x <- safe_num(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

scale_by_gene <- function(df, value_col = "plot_value_raw") {
  df$plot_value <- NA_real_
  for (gene in unique(df$symbol)) {
    idx <- which(df$symbol == gene)
    df$plot_value[idx] <- zscore_vec(df[[value_col]][idx])
  }
  df
}

build_category_matrix <- function(evidence, de, expr, category_id, cfg) {
  rows <- evidence[evidence$category_id == category_id & evidence$in_gsea_leading_edge_bool %in% TRUE, , drop = FALSE]
  if (nrow(rows) == 0) return(data.frame())
  rows <- rows[, setdiff(names(rows), c("log2FC", "padj")), drop = FALSE]
  rows <- merge(rows, de[, c("symbol", "log2FC", "padj"), drop = FALSE], by = "symbol", all.x = TRUE, sort = FALSE)
  rows <- rows[order(-abs(safe_num(rows$log2FC)), safe_num(rows$min_source_padj), rows$symbol), , drop = FALSE]
  rows <- rows[!duplicated(rows$symbol), , drop = FALSE]
  if (cfg$top_genes > 0) rows <- head(rows, cfg$top_genes)
  sample_cols <- attr(expr, "sample_cols")
  wide <- merge(rows, expr, by = "symbol", all.x = FALSE, sort = FALSE)
  if (nrow(wide) == 0 || length(sample_cols) < 2) return(data.frame())
  long <- do.call(rbind, lapply(sample_cols, function(sample) {
    data.frame(
      wide[, setdiff(names(wide), sample_cols), drop = FALSE],
      sample = sample,
      value_raw = wide[[sample]],
      stringsAsFactors = FALSE
    )
  }))
  source_path <- attr(expr, "source_path") %||% ""
  expression_source <- attr(expr, "expression_source") %||% expression_matrix_source(source_path)
  use_raw_values <- identical(expression_source$input_scale, "transformed")
  long$plot_value_raw <- if (use_raw_values) {
    safe_num(long$value_raw)
  } else {
    log2(pmax(safe_num(long$value_raw), 0) + 1)
  }
  long$input_scale <- expression_source$input_scale
  long$input_scale_source <- expression_source$input_scale_source
  long$expression_transform <- expression_source$expression_transform
  long <- long[order(-abs(safe_num(long$log2FC)), long$symbol, match(long$sample, sample_cols)), , drop = FALSE]
  long$sample <- factor(long$sample, levels = sample_cols)
  long$scale <- cfg$scale
  if (cfg$scale == "zscore") {
    long <- scale_by_gene(long, "plot_value_raw")
  } else {
    long$plot_value <- long$plot_value_raw
  }
  long
}

main <- function() {
  cfg <- parse_args(commandArgs(trailingOnly = TRUE))
  cfg$plot_metadata <- lisa_plot_metadata(cfg$project_dir, cfg$analysis_id)
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Required package: ggplot2", call. = FALSE)
  paths <- paths_for(cfg$project_dir, cfg$analysis_id, cfg$universe)
  if (!file.exists(paths$evidence)) stop(sprintf("Missing gene-category contribution table: %s", paths$evidence), call. = FALSE)
  dir.create(paths$out_dir, recursive = TRUE, showWarnings = FALSE)

  evidence <- read_tsv(paths$evidence)
  if (nrow(evidence) == 0) stop("Gene-category contribution table is empty.", call. = FALSE)
  evidence$symbol <- toupper(safe_chr(evidence$symbol))
  evidence$category_id <- safe_chr(evidence$category_id)
  evidence$category_display_name <- safe_chr(evidence$category_display_name)
  evidence$macrogroup_id <- safe_chr(evidence$macrogroup_id)
  evidence$macrogroup_name <- safe_chr(evidence$macrogroup_name)
  evidence$in_gsea_leading_edge_bool <- as.logical(evidence$in_gsea_leading_edge)
  de <- load_de_context(paths)
  expression_source <- discover_matrix_source(paths, cfg)
  expression_matrix_path <- expression_source$path
  if (expression_matrix_path == "") {
    stop(
      paste(
        "No expression/count matrix found for leading-edge heatmaps.",
        "Provide --expression-matrix or add one of expression_matrix_path/counts_matrix_path/vst_matrix_path/normalized_matrix_path/tpm_matrix_path to config/de_index.tsv.",
        "Refusing to fall back to baseMean/log2FC single-column plots."
      ),
      call. = FALSE
    )
  }
  expr <- filter_expression_matrix_samples(read_expression_matrix(expression_matrix_path, expression_source), paths, cfg)
  sample_cols <- attr(expr, "sample_cols")
  count_source <- if (expression_source$input_scale == "transformed") {
    basename(expression_matrix_path)
  } else if (grepl("tpm", basename(expression_matrix_path), ignore.case = TRUE)) {
    paste0("log2(", basename(expression_matrix_path), " TPM+1)")
  } else {
    paste0("log2(", basename(expression_matrix_path), "+1)")
  }

  categories <- unique(evidence$category_id[evidence$in_gsea_leading_edge_bool %in% TRUE & evidence$category_id != ""])
  if (!is.na(cfg$categories) && nzchar(cfg$categories)) {
    categories <- intersect(categories, trimws(strsplit(cfg$categories, ",", fixed = TRUE)[[1]]))
  }
  index <- list()
  for (category_id in categories) {
    matrix <- build_category_matrix(evidence, de, expr, category_id, cfg)
    if (nrow(matrix) == 0) next
    meta <- evidence[evidence$category_id == category_id, , drop = FALSE][1, , drop = FALSE]
    category_label <- meta$category_display_name[[1]]
    safe_cat <- safe_file_component(category_id)
    prefix <- paste(paths$prefix, safe_cat, "leading_edge_gene_heatmap", sep = "_")
    matrix_path <- file.path(paths$out_dir, paste0(prefix, "_matrix.tsv"))
    recipe_path <- file.path(paths$out_dir, paste0(prefix, "_recipe.R"))
    png_path <- file.path(paths$out_dir, paste0(prefix, "_", cfg$scale, ".png"))
    pdf_path <- file.path(paths$out_dir, paste0(prefix, "_", cfg$scale, ".pdf"))
    svg_path <- file.path(paths$out_dir, paste0(prefix, "_", cfg$scale, ".svg"))
    write_tsv(matrix, matrix_path)
    context <- lisa_heatmap_context(matrix, cfg, cfg$analysis_id, cfg$universe,
      category_label, count_source, matrix_path)
    lisa_heatmap_write_recipe(recipe_path, context, run_root = cfg$project_dir)
    plot <- lisa_heatmap_plot(matrix, cfg, cfg$analysis_id, cfg$universe, category_label, count_source)
    n_genes <- length(unique(matrix$symbol))
    n_samples <- length(unique(as.character(matrix$sample)))
    dims <- lisa_heatmap_size(n_genes, n_samples, 2L)
    w <- dims[["width"]]
    h <- dims[["height"]]
    for (format in cfg$formats) {
      figure_path <- switch(format, png = png_path, pdf = pdf_path, svg = svg_path)
      lisa_heatmap_save_native(figure_path, plot, width = w, height = h)
    }
    index[[length(index) + 1L]] <- data.frame(
      analysis_id = cfg$analysis_id,
      universe = cfg$universe,
      category_id = category_id,
      category_display_name = category_label,
      macrogroup_id = meta$macrogroup_id[[1]],
      macrogroup_name = meta$macrogroup_name[[1]],
      n_leading_edge_genes = n_genes,
      n_samples = n_samples,
      scale = cfg$scale,
      count_source = count_source,
      expression_matrix = if (nzchar(cfg$expression_matrix_column)) {
        basename(expression_matrix_path)
      } else {
        expression_matrix_path
      },
      input_scale = expression_source$input_scale,
      input_scale_source = expression_source$input_scale_source,
      expression_transform = expression_source$expression_transform,
      n_matrix_sample_columns = length(attr(expr, "available_sample_cols") %||% sample_cols),
      selected_samples = paste(sample_cols, collapse = ";"),
      sample_selection_method = attr(expr, "sample_selection_method") %||% "",
      sample_include_regex = attr(expr, "sample_include_regex") %||% "",
      sample_exclude_regex = attr(expr, "sample_exclude_regex") %||% "",
      matrix_tsv = matrix_path,
      recipe_r = recipe_path,
      plot_png = png_path,
      plot_pdf = pdf_path,
      plot_svg = svg_path,
      stringsAsFactors = FALSE
    )
  }
  index_df <- if (length(index) == 0) data.frame() else do.call(rbind, index)
  index_path <- file.path(paths$out_dir, paste0(paths$prefix, "_leading_edge_gene_heatmaps_index.tsv"))
  write_tsv(index_df, index_path)
  readme <- file.path(paths$out_dir, "README_leading_edge_gene_heatmaps.txt")
  writeLines(c(
    "Leading-edge gene heatmaps by LISA category",
    "",
    "Each plot shows unique genes found in the GSEA leading edge for one LISA category.",
    "Rows are ordered by absolute log2 fold change from the DE table.",
    "Columns are biological samples from the provided expression/count matrix.",
    "Sample columns are displayed in their input-matrix order after filtering; the renderer does not cluster, alphabetize or otherwise reorder samples.",
    "When the matrix contains a broader panel, sample columns are selected for the current analysis using sample_include_regex/heatmap_sample_include_regex from config/de_index.tsv when present, otherwise by conservative analysis_id inference.",
    "The two right-hand annotation columns show DE log2FC and DE FDR for each row gene; they are not expression samples.",
    "The de_index column declares the input scale: vst_matrix/vst_matrix_path preserves stored transformed values; counts_matrix/tpm_matrix/sample_counts paths use log2(value + 1).",
    "Only generic or unindexed paths retain the legacy VST/rlog/RLD/log2 basename heuristic. Declare a VST matrix via vst_matrix_path to make its scale independent of its filename.",
    "Matrix TSV files and the heatmap index record input_scale, input_scale_source and expression_transform.",
    "The default plot scale is per-gene z-score across samples. Matrix TSV files preserve raw and plotted values.",
    "Each matrix has a self-contained native R recipe: Rscript --vanilla *_recipe.R SOURCE_MATRIX.tsv NEW_OUTPUT.png|pdf|svg. Offline ggplot2 and digest are required.",
    "The two-argument recipe binds the exact matrix SHA-256 and embeds native ordering, labels and dimensions; it draws saved plot_value without new transformations or analysis. Use a fresh PNG, PDF or native Cairo SVG destination.",
    "The script refuses to fall back to baseMean/log2FC single-column plots."
  ), readme, useBytes = TRUE)
  if (nrow(index_df) == 0) {
    message("No leading-edge category heatmaps generated.")
  } else {
    print(index_df[, c("analysis_id", "universe", "n_leading_edge_genes", "scale", "count_source")])
  }
}

if (sys.nframe() == 0L) main()
