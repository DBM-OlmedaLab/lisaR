# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# calibration read-only sensitivity helpers for the provisional lisa-gps contract.
# These functions never select weights, pool biological contrasts, or create a
# biological gold standard.  They only summarize rank sensitivity.

#' Build the compact calibration provisional-weight sensitivity grid
#'
#' @param evidence pathway-compatible gene-set evidence.
#' @param weights Optional data frame of explicitly documented weights. Zero
#'   weights exclude that component, including its missing values, from the
#'   scenario score. All positive-weight components must be finite; a zero in
#'   any positive-weight component yields score zero. Equal-weight Pareto
#'   metadata still describes all three components, not the custom scenario.
#' @return A list containing the weights, per-gene scores, rank stability, and
#'   the equal-weight first-Pareto-front metadata.
#' @keywords internal
lisa_gps_sensitivity_grid <- function(evidence, weights = NULL) {
  if (is.null(weights)) {
    weights <- data.frame(
      scenario_id = c("equal", "effect_emphasis", "confidence_emphasis", "breadth_emphasis"),
      effect_weight = c(1/3, .50, .25, .25),
      confidence_weight = c(1/3, .25, .50, .25),
      breadth_weight = c(1/3, .25, .25, .50),
      stringsAsFactors = FALSE
    )
  }
  weights <- as.data.frame(weights, stringsAsFactors = FALSE)
  required <- c("scenario_id", "effect_weight", "confidence_weight", "breadth_weight")
  missing <- setdiff(required, names(weights))
  if (length(missing)) stop("LISA-GPS-6B-001 weights are missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  if (anyDuplicated(weights$scenario_id) || any(!nzchar(as.character(weights$scenario_id)))) stop("LISA-GPS-6B-002 scenario_id values must be unique and non-empty.", call. = FALSE)
  numeric_weights <- as.matrix(data.frame(lapply(weights[required[-1]], as.numeric)))
  if (any(!is.finite(numeric_weights)) || any(numeric_weights < 0) || any(abs(rowSums(numeric_weights) - 1) > 1e-10)) stop("LISA-GPS-6B-003 each non-negative weight row must sum to one.", call. = FALSE)
  base <- lisa_gps_prioritization(evidence)
  component_names <- c("effect_strength_normalized", "de_confidence_normalized", "independent_category_breadth_normalized")
  rows <- lapply(seq_len(nrow(weights)), function(i) {
    # Select components before checking completeness or taking logarithms.
    # Otherwise log(0) * 0 is NaN, and a missing ignored component can wrongly
    # suppress a valid score. Positive-weight missing components stay NA.
    active <- numeric_weights[i, ] > 0
    components <- as.matrix(base[, component_names[active], drop = FALSE])
    complete <- base$applicability == "applicable" &
      rowSums(!is.finite(components)) == 0L
    score <- rep(NA_real_, nrow(base))
    zero <- complete & rowSums(components == 0, na.rm = TRUE) > 0L
    positive <- complete & !zero
    score[positive] <- exp(rowSums(sweep(
      log(components[positive, , drop = FALSE]), 2L,
      numeric_weights[i, active], `*`
    )))
    score[zero] <- 0
    rank <- rep(NA_integer_, nrow(base)); rank[complete] <- rank(-score[complete], ties.method = "min")
    audit_columns <- c(
      "priority_id", "raw_gene_set_count", "independent_category_breadth",
      "gene_set_ids", "category_ids", "effect_signed_raw",
      "effect_strength_raw", "padj_raw", "de_confidence_raw",
      component_names, "pareto_front", "applicability"
    )
    cbind(weights[i, , drop = FALSE][rep(1, nrow(base)), , drop = FALSE], base[, audit_columns, drop = FALSE], lisa_gps = score, rank = rank)
  })
  scores <- do.call(rbind, rows); rownames(scores) <- NULL
  list(weights = weights, scores = scores, rank_stability = lisa_gps_rank_stability(scores), equal_weight_pareto = scores[scores$scenario_id == "equal", c("priority_id", "pareto_front"), drop = FALSE])
}

#' Summarize rank stability across explicitly documented sensitivity scenarios
#' @param scores Result component `scores` from [lisa_gps_sensitivity_grid()].
#' @return One row per priority ID; this is not a pooled biological ranking.
#' @keywords internal
lisa_gps_rank_stability <- function(scores) {
  x <- as.data.frame(scores, stringsAsFactors = FALSE)
  if (!all(c("priority_id", "rank", "lisa_gps") %in% names(x))) stop("LISA-GPS-6B-004 scores require priority_id, rank, and lisa_gps.", call. = FALSE)
  ids <- sort(unique(as.character(x$priority_id)))
  out <- lapply(ids, function(id) {
    z <- x[x$priority_id == id & is.finite(x$rank), , drop = FALSE]
    data.frame(priority_id = id, scenarios_ranked = nrow(z), rank_min = if (nrow(z)) min(z$rank) else NA_integer_, rank_max = if (nrow(z)) max(z$rank) else NA_integer_, rank_sd = if (nrow(z) > 1) stats::sd(z$rank) else 0, lisa_gps_min = if (nrow(z)) min(z$lisa_gps) else NA_real_, lisa_gps_max = if (nrow(z)) max(z$lisa_gps) else NA_real_, stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

#' Report concordance for two kept-separate biological contrasts
#' @param first First per-gene score table.
#' @param second Second per-gene score table.
#' @param top_n Number of leading ranks used for overlap metadata.
#' @return One-row non-pooled concordance metadata.
#' @keywords internal
lisa_gps_concordance <- function(first, second, top_n = 25L) {
  required <- c("priority_id", "rank", "lisa_gps")
  if (!all(required %in% names(first)) || !all(required %in% names(second))) stop("LISA-GPS-6B-005 concordance inputs require priority_id, rank, and lisa_gps.", call. = FALSE)
  x <- merge(first[, required], second[, required], by = "priority_id", suffixes = c("_first", "_second"))
  x <- x[stats::complete.cases(x), , drop = FALSE]
  top_n <- as.integer(top_n)
  first_top <- x$priority_id[x$rank_first <= top_n]; second_top <- x$priority_id[x$rank_second <= top_n]
  data.frame(comparison = "separate_inputs_no_pooling", n_shared = nrow(x), rank_spearman = if (nrow(x) > 1) stats::cor(x$rank_first, x$rank_second, method = "spearman") else NA_real_, score_spearman = if (nrow(x) > 1) stats::cor(x$lisa_gps_first, x$lisa_gps_second, method = "spearman") else NA_real_, top_n = top_n, top_n_overlap = length(intersect(first_top, second_top)), stringsAsFactors = FALSE)
}
