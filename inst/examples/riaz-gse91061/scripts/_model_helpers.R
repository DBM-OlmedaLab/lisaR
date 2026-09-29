# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Statistical-model helpers shared by the Riaz worked-example scripts.
#
# These helpers preserve two distinct identifiers:
#   * Entrez ID is the native row key of the public count matrix.
#   * HGNC symbol is added only at the lisaR export boundary.
#
# When historical Entrez IDs map to the same current symbol, the representative
# row is selected by mean normalized abundance. This rule is independent of
# effect size and significance, preventing significance-driven deduplication.

riaz_read_counts <- function(sample_ids) {
  path <- file.path(
    riaz_input_cache(),
    "GSE91061_BMS038109Sample.hg19KnownGene.raw.csv.gz"
  )
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  table <- read.csv(
    connection,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  entrez_id <- as.character(table[[1]])
  table[[1]] <- NULL
  if (anyDuplicated(entrez_id)) {
    stop("The public count matrix contains duplicated Entrez IDs.")
  }
  rownames(table) <- entrez_id
  missing <- setdiff(sample_ids, colnames(table))
  if (length(missing)) {
    stop("Samples absent from count matrix: ", paste(missing, collapse = ", "))
  }
  matrix <- as.matrix(table[, sample_ids, drop = FALSE])
  storage.mode(matrix) <- "integer"
  if (anyNA(matrix) || any(matrix < 0L)) {
    stop("The selected raw-count matrix contains missing or negative values.")
  }
  matrix
}

riaz_fit_deseq2_matrix <- function(count_matrix, design, model_matrix) {
  if (nrow(design) != ncol(count_matrix)) {
    stop("Design rows and count-matrix columns differ.")
  }
  if (!identical(design$sample_id, colnames(count_matrix))) {
    stop("Design and count matrix are not in the same sample order.")
  }
  rank <- qr(model_matrix)$rank
  if (rank != ncol(model_matrix)) {
    stop(
      "Model matrix is not full rank: rank=", rank,
      ", columns=", ncol(model_matrix), "."
    )
  }
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = design,
    design = model_matrix
  )
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds, quiet = TRUE)
  dds <- DESeq2::nbinomWaldTest(
    dds,
    modelMatrix = model_matrix,
    # A generous iteration ceiling is important for this heterogeneous bulk
    # tumour cohort. It lets difficult coefficients converge instead of
    # silently contributing an unstable statistic to a downstream rank list.
    maxit = 10000,
    quiet = TRUE
  )
  list(dds = dds, rank = rank, columns = colnames(model_matrix))
}

riaz_numeric_contrast <- function(dds, weights, alpha = 0.05) {
  names_available <- DESeq2::resultsNames(dds)
  unknown <- setdiff(names(weights), names_available)
  if (length(unknown)) {
    stop(
      "Contrast references unknown coefficients: ",
      paste(unknown, collapse = ", "),
      "; available: ", paste(names_available, collapse = ", ")
    )
  }
  vector <- setNames(rep(0, length(names_available)), names_available)
  vector[names(weights)] <- weights
  DESeq2::results(dds, contrast = vector, alpha = alpha)
}

riaz_annotate_result <- function(result, dds, analysis_id) {
  table <- as.data.frame(result)
  table$entrez_id <- rownames(table)
  table$symbol <- AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = table$entrez_id,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  table$beta_converged <- as.logical(
    SummarizedExperiment::mcols(dds)$betaConv
  )
  table$analysis_id <- analysis_id
  table[
    ,
    c(
      "analysis_id", "entrez_id", "symbol", "baseMean",
      "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
      "beta_converged"
    ),
    drop = FALSE
  ]
}

riaz_collapse_symbols <- function(table) {
  usable <- table[!is.na(table$symbol) & nzchar(table$symbol), , drop = FALSE]
  usable <- usable[
    order(
      usable$symbol,
      -usable$baseMean,
      suppressWarnings(as.numeric(usable$entrez_id)),
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  keep <- !duplicated(usable$symbol)
  duplicated_symbol <- duplicated(usable$symbol) |
    duplicated(usable$symbol, fromLast = TRUE)
  audit <- usable[duplicated_symbol, , drop = FALSE]
  if (nrow(audit) > 0L) {
    # Recompute duplication status inside the selected audit rows. Using the
    # full-table logical vector here would have the wrong length when there are
    # no duplicate symbols, which is precisely the empty-edge case guarded by
    # this branch.
    audit$resolution <- ifelse(
      duplicated(audit$symbol),
      "removed_duplicate",
      "kept_highest_baseMean"
    )
  } else {
    audit$resolution <- character()
  }
  list(table = usable[keep, , drop = FALSE], audit = audit)
}

riaz_export_model_result <- function(result, dds, analysis_id, output_dir) {
  annotated <- riaz_annotate_result(result, dds, analysis_id)
  collapsed <- riaz_collapse_symbols(annotated)
  riaz_write_tsv(
    annotated,
    file.path(output_dir, paste0(analysis_id, "_entrez_full.tsv"))
  )
  riaz_write_tsv(
    collapsed$table,
    file.path(output_dir, paste0(analysis_id, "_symbol_full.tsv"))
  )
  riaz_write_tsv(
    collapsed$audit,
    file.path(output_dir, paste0(analysis_id, "_symbol_resolution.tsv"))
  )
  invisible(collapsed$table)
}
