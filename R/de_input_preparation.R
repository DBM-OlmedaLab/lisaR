# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

#' Prepare an auditable differential-expression table for LISA
#'
#' Adapts an already computed data frame without loading DESeq2, edgeR or limma.
#' Every tested row, including non-significant rows and permitted missing
#' statistics, remains in `full_data`. No adjusted-P threshold is applied.
#' `data` differs only if an explicitly requested duplicate policy selects or
#' remaps rows. This function does not fit models or infer sample metadata.
#'
#' @param data A data frame containing the complete tested DE results.
#' @param method Declared table convention: `"DESeq2"`, `"edgeR"` or `"limma"`.
#'   This is a caller declaration, not verification of the producing software.
#' @param columns Named list mapping `symbol`, `logfc`, `rank`, `pvalue` and
#'   `padj` to source columns. Use `symbol = ".rownames"` explicitly to take
#'   identifiers from non-default row names. Standard effect/P-value names and
#'   DESeq2 signed Wald `stat` or limma `t` are defaults. Unsigned DESeq2 LRT
#'   results require a separately computed signed rank. edgeR requires an explicit signed
#'   `rank` column: multi-df `F` or unsigned `LR` must never be arbitrarily signed.
#' @param positive_direction One non-empty description of what positive effects
#'   mean, for example `"treated relative to control"`. No contrast is inferred.
#' @param rank_direction Whether the supplied signed ranking has the same or
#'   opposite orientation as the effect. Opposite ranks are negated once; finite
#'   nonzero ranks must then agree in sign with finite nonzero effects.
#' @param duplicate_policy An existing LISA duplicate policy: `"error"`, a
#'   deterministic `select` policy, or a row-level `mapping_file` policy.
#'   DE statistics are never aggregated. Identifier collisions are checked
#'   case-insensitively, consistently with LISA ranking normalization.
#' @param na_policy `"retain"` keeps missing statistics and records rank
#'   eligibility; `"error"` rejects missing statistics. Missing/empty identifiers,
#'   malformed numbers and infinite statistics are always errors.
#' @param matrix_scale Optional separately declared expression-matrix scale:
#'   `"transformed"` or `"count_like"`. Recorded only; no expression matrix is
#'   inferred, generated, transformed or validated by this adapter.
#'
#' @return A list with canonical `data` and unfiltered `full_data` columns
#'   `symbol`, `log2FoldChange`, `rank_value`, `pvalue`, `padj`; an input-declaration
#'   `receipt`; explicit `column_mapping`; `row_audit` including rank eligibility;
#'   and `duplicate_resolution`. Row order is preserved except explicit duplicate
#'   resolution. The ranking and effect both describe `positive_direction`.
#' @export
#'
#' @examples
#' results <- data.frame(gene = c("A", "B", "C"),
#'   log2FoldChange = c(1, -0.5, 0.1), stat = c(3, -1, 0.2),
#'   pvalue = c(0.001, 0.3, 0.8), padj = c(0.01, 0.6, NA))
#' prepared <- prepare_lisa_de_input(results, "DESeq2",
#'   columns = list(symbol = "gene"),
#'   positive_direction = "treated relative to control")
#' prepared$data
#' prepared$receipt
prepare_lisa_de_input <- function(data, method = c("DESeq2", "edgeR", "limma"),
                                  columns = list(), positive_direction,
                                  rank_direction = c("same_as_effect", "opposite_to_effect"),
                                  duplicate_policy = "error",
                                  na_policy = c("retain", "error"), matrix_scale = NULL) {
  method <- match.arg(method)
  rank_direction <- match.arg(rank_direction)
  na_policy <- match.arg(na_policy)
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data))) {
    stop("LISA-DE-INPUT-001 data must be a non-empty data frame with unique column names.", call. = FALSE)
  }
  if (missing(positive_direction) || !is.character(positive_direction) ||
      length(positive_direction) != 1L || is.na(positive_direction) ||
      !nzchar(trimws(positive_direction))) {
    stop("LISA-DE-INPUT-002 declare one non-empty positive_direction; no contrast orientation is inferred.", call. = FALSE)
  }
  if (!is.null(matrix_scale) && (!is.character(matrix_scale) ||
      length(matrix_scale) != 1L || is.na(matrix_scale) ||
      !matrix_scale %in% c("transformed", "count_like"))) {
    stop("LISA-DE-INPUT-003 matrix_scale must be NULL, transformed or count_like.", call. = FALSE)
  }
  allowed <- c("symbol", "logfc", "rank", "pvalue", "padj")
  if (!is.list(columns) || (length(columns) && (is.null(names(columns)) ||
      anyDuplicated(names(columns)) || any(!names(columns) %in% allowed)))) {
    stop("LISA-DE-INPUT-004 columns must be a named mapping of symbol, logfc, rank, pvalue and padj.", call. = FALSE)
  }
  defaults <- switch(method,
    DESeq2 = list(symbol = "symbol", logfc = "log2FoldChange", rank = "stat", pvalue = "pvalue", padj = "padj"),
    edgeR = list(symbol = "symbol", logfc = "logFC", rank = NULL, pvalue = "PValue", padj = "FDR"),
    limma = list(symbol = "symbol", logfc = "logFC", rank = "t", pvalue = "P.Value", padj = "adj.P.Val"))
  mapping <- defaults
  for (key in names(columns)) mapping[key] <- columns[key]
  if (is.null(mapping$rank)) {
    stop("LISA-DE-INPUT-005 edgeR requires an explicit signed rank column; an F or LR statistic has no defensible sign for an unspecified contrast.", call. = FALSE)
  }
  for (key in allowed) {
    value <- mapping[[key]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("LISA-DE-INPUT-004 mapping for ", key, " must name exactly one source column.", call. = FALSE)
    }
    if (!(key == "symbol" && value == ".rownames") && !value %in% names(data)) {
      stop("LISA-DE-INPUT-004 source column is absent: ", key, " = ", value, ".", call. = FALSE)
    }
  }
  if (mapping$rank %in% c("F", "LR", "F.statistic", "Chisq", "chisq") ||
      mapping$rank %in% unlist(mapping[c("pvalue", "padj")], use.names = FALSE)) {
    stop("LISA-DE-INPUT-005 rank must name an explicitly signed statistic, not unsigned F/LR or significance.", call. = FALSE)
  }
  if (mapping$symbol == ".rownames") {
    ids <- rownames(data)
    if (is.null(ids) || identical(ids, as.character(seq_len(nrow(data))))) {
      stop("LISA-DE-INPUT-006 default sequential row names are ambiguous gene identifiers; supply an explicit ID column.", call. = FALSE)
    }
  } else {
    ids <- as.character(data[[mapping$symbol]])
  }
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || any(ids != trimws(ids))) {
    stop("LISA-DE-INPUT-006 identifiers must be non-missing, non-empty and free of surrounding whitespace.", call. = FALSE)
  }
  numeric_column <- function(key) {
    raw <- data[[mapping[[key]]]]
    if (is.factor(raw)) raw <- as.character(raw)
    if (!is.numeric(raw) && !is.character(raw)) {
      stop("LISA-DE-INPUT-007 ", key, " must be numeric, not logical/list data.", call. = FALSE)
    }
    values <- suppressWarnings(as.numeric(raw))
    if (any(!is.na(raw) & is.na(values)) || any(is.infinite(values))) {
      stop("LISA-DE-INPUT-007 ", key, " contains malformed or infinite numeric statistics.", call. = FALSE)
    }
    if (na_policy == "error" && anyNA(values)) {
      stop("LISA-DE-INPUT-008 missing statistics violate na_policy = error (", key, ").", call. = FALSE)
    }
    if (key %in% c("pvalue", "padj") && any(values < 0 | values > 1, na.rm = TRUE)) {
      stop("LISA-DE-INPUT-009 pvalue and padj must lie in [0, 1].", call. = FALSE)
    }
    values
  }
  full <- data.frame(symbol = ids, log2FoldChange = numeric_column("logfc"),
    rank_value = numeric_column("rank"), pvalue = numeric_column("pvalue"),
    padj = numeric_column("padj"), stringsAsFactors = FALSE)
  if (rank_direction == "opposite_to_effect") full$rank_value <- -full$rank_value
  comparable <- is.finite(full$rank_value) & is.finite(full$log2FoldChange) &
    full$rank_value != 0 & full$log2FoldChange != 0
  if (any(sign(full$rank_value[comparable]) != sign(full$log2FoldChange[comparable]))) {
    stop("LISA-DE-INPUT-010 rank/effect signs conflict with the declared rank_direction; supply a compatible signed contrast statistic.", call. = FALSE)
  }
  if (!any(is.finite(full$rank_value))) {
    stop("LISA-DE-INPUT-011 at least one finite signed rank is required.", call. = FALSE)
  }
  type <- lisa_validate_duplicate_policy(duplicate_policy, "duplicate_policy")
  if (type == "aggregate") {
    stop("LISA-DE-INPUT-012 DE statistics must not be aggregated; use error, select or mapping_file.", call. = FALSE)
  }
  # Include caller metadata while applying an explicit selection/tie-break rule,
  # then return only canonical columns. Never overwrite a metadata column.
  work <- data
  reserved <- c(".lisa_source_row", ".lisa_symbol", ".lisa_logfc", ".lisa_rank", ".lisa_pvalue", ".lisa_padj")
  if (any(reserved %in% names(work))) stop("LISA-DE-INPUT-004 source uses reserved .lisa_* adapter column names.", call. = FALSE)
  work$.lisa_source_row <- seq_len(nrow(full))
  work$.lisa_symbol <- full$symbol
  work$.lisa_logfc <- full$log2FoldChange
  work$.lisa_rank <- full$rank_value
  work$.lisa_pvalue <- full$pvalue
  work$.lisa_padj <- full$padj
  resolution <- if (type == "mapping_file") {
    lisa_de_input_mapping(work, duplicate_policy)
  } else {
    lisa_resolve_duplicates(work, ".lisa_symbol", duplicate_policy,
      input_class = "de_table", columns = list(padj = ".lisa_padj", rank = ".lisa_rank"))
  }
  retained <- resolution$data$.lisa_source_row
  effective <- full[retained, , drop = FALSE]
  effective$symbol <- as.character(resolution$data$.lisa_symbol)
  if (anyNA(effective$symbol) || any(!nzchar(trimws(effective$symbol))) ||
      any(effective$symbol != trimws(effective$symbol))) {
    stop("LISA-DE-INPUT-006 mapping produced missing, empty or whitespace-padded identifiers.", call. = FALSE)
  }
  rownames(effective) <- NULL
  audit <- data.frame(source_row = seq_len(nrow(full)), source_id = ids,
    effective_id = effective$symbol[match(seq_len(nrow(full)), retained)],
    retained = seq_len(nrow(full)) %in% retained,
    rank_eligible = is.finite(full$rank_value),
    missing_effect = is.na(full$log2FoldChange), missing_rank = is.na(full$rank_value),
    missing_pvalue = is.na(full$pvalue), missing_padj = is.na(full$padj),
    stringsAsFactors = FALSE)
  receipt <- data.frame(adapter_version = "1", declared_method = method,
    positive_direction = positive_direction, supplied_rank_direction = rank_direction,
    effective_rank_direction = "same_as_effect", rank_source = mapping$rank,
    rank_source_declaration = if ("rank" %in% names(columns)) "explicit" else "method_convention",
    duplicate_policy = type, na_policy = na_policy, source_rows = nrow(full),
    effective_rows = nrow(effective), full_rows_retained = nrow(full),
    significance_filter = "none", finite_rank_rows = sum(is.finite(effective$rank_value)),
    matrix_scale = if (is.null(matrix_scale)) "not_declared" else matrix_scale,
    producing_software_verified = FALSE, stringsAsFactors = FALSE)
  list(data = effective, full_data = full, receipt = receipt,
    column_mapping = data.frame(canonical = allowed, source = unlist(mapping[allowed], use.names = FALSE),
      declaration = ifelse(allowed %in% names(columns), "explicit", "method_convention"), stringsAsFactors = FALSE),
    row_audit = audit, duplicate_resolution = resolution$evidence)
}

# This adapter reads identifier mappings as text, so e.g. "001" is never
# silently recoded as "1". Fractional source-row numbers also fail closed.
lisa_de_input_mapping <- function(data, policy) {
  mapping <- utils::read.delim(policy$mapping_file, sep = "\t", quote = "",
    comment.char = "", colClasses = "character", check.names = FALSE,
    stringsAsFactors = FALSE, na.strings = character())
  row_col <- policy$source_row_col %||% "source_row"
  id_col <- policy$output_id_col %||% "output_id"
  if (!all(c(row_col, id_col) %in% names(mapping))) {
    stop("LISA-DE-INPUT-013 mapping_file requires source-row and output-ID columns.", call. = FALSE)
  }
  rows <- suppressWarnings(as.numeric(mapping[[row_col]]))
  if (anyNA(rows) || any(!is.finite(rows)) || any(rows != floor(rows)) ||
      any(rows < 1 | rows > nrow(data)) || anyDuplicated(rows)) {
    stop("LISA-DE-INPUT-013 mapping source rows must be unique in-range integers.", call. = FALSE)
  }
  canonical <- lisa_canonical_duplicate_id(data$.lisa_symbol)
  affected <- which(duplicated(canonical) | duplicated(canonical, fromLast = TRUE))
  if (!all(affected %in% rows)) {
    stop("LISA-DE-INPUT-013 mapping must resolve every duplicate source row.", call. = FALSE)
  }
  output <- data
  output$.lisa_symbol[affected] <- mapping[[id_col]][match(affected, rows)]
  if (anyDuplicated(lisa_canonical_duplicate_id(output$.lisa_symbol))) {
    stop("LISA-DE-INPUT-013 mapping must produce unique identifiers.", call. = FALSE)
  }
  list(data = output, evidence = lisa_duplicate_audit(data, ".lisa_symbol", output,
    type = "mapping_file", parameters = paste0("input_class=de_table;mapping_file=",
      normalizePath(policy$mapping_file, mustWork = TRUE))))
}
