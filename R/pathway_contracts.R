# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# pathway prioritization, PATHWAYS-universe, and offline KEGG map contracts.
#
# These helpers are development contracts.  They preserve raw inputs and make
# no inferential or causal claim.

lisa_gps_development_contract <- function() {
  list(
    development_contract_version = "lisa-gps-development@0.1.0-dev",
    status = "provisional; this is not a frozen final score",
    label = "prioritization heuristic; not inferential or causal",
    redundancy_unit = "versioned LISA category_id",
    transformations = list(
      effect_strength = "abs(effect); finite values min-max normalized; a finite constant column normalizes to 1",
      de_confidence = "-log10(padj); finite values min-max normalized; a finite constant column normalizes to 1",
      independent_category_breadth = "distinct category_id; min-max normalized; a finite constant column normalizes to 1",
      score = "geometric mean of the three normalized components with equal weights (1/3 each)"
    ),
    na_behavior = "missing effect or padj remains NA; padj=0 is deterministically clamped to .Machine$double.xmin before -log10; a score is NA unless all three normalized components are finite",
    tie_behavior = "input rows are ordered by stable lexical priority_id for display; equal component values retain the same score and Pareto front",
    evidence_mode_applicability = "full_de only; all other modes are not_applicable",
    pareto = "exact first-front membership only; robustness metadata, never a competing score or ranking"
  )
}

lisa_pathway_minmax <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out)
  lo <- min(x[ok]); hi <- max(x[ok])
  if (identical(lo, hi)) out[ok] <- 1 else out[ok] <- (x[ok] - lo) / (hi - lo)
  out
}

#' Construct the provisional pathway lisa-gps prioritization table
#'
#' @param evidence Gene-set evidence with `priority_id`, `gene_set_id`,
#'   `category_id`, `effect`, and `padj` columns.
#' @param evidence_mode Declared input evidence mode.
#' @return A reproducible prioritization table with raw and normalized inputs.
#' @keywords internal
lisa_gps_prioritization <- function(evidence, evidence_mode = "full_de") {
  required <- c("priority_id", "gene_set_id", "category_id", "effect", "padj")
  evidence <- as.data.frame(evidence, stringsAsFactors = FALSE)
  missing <- setdiff(required, names(evidence))
  if (length(missing)) stop("LISA-GPS-001 evidence is missing required column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  id_columns <- c("priority_id", "gene_set_id", "category_id")
  if (any(vapply(id_columns, function(column) any(is.na(evidence[[column]]) | !nzchar(as.character(evidence[[column]]))), logical(1)))) {
    stop("LISA-GPS-002 priority_id, gene_set_id, and category_id must be non-empty.", call. = FALSE)
  }
  evidence$priority_id <- as.character(evidence$priority_id)
  evidence$gene_set_id <- as.character(evidence$gene_set_id)
  evidence$category_id <- as.character(evidence$category_id)
  evidence$effect <- suppressWarnings(as.numeric(evidence$effect))
  evidence$padj <- suppressWarnings(as.numeric(evidence$padj))
  if (any(!is.na(evidence$padj) & (evidence$padj < 0 | evidence$padj > 1))) {
    stop("LISA-GPS-003 padj values must be within [0, 1]. Repair: provide adjusted p-values, retaining NA only when the value is unavailable.", call. = FALSE)
  }
  # LISA dictionaries are legitimately many-to-many: one gene set may support
  # more than one versioned LISA category. Count the source gene set once in
  # raw breadth while retaining every distinct category assignment in the
  # independent-category component.
  same_value <- function(x) length(unique(x)) == 1L
  groups <- split(evidence, evidence$priority_id, drop = TRUE)
  inconsistent <- vapply(groups, function(x) !same_value(x$effect) || !same_value(x$padj), logical(1))
  if (any(inconsistent)) stop("LISA-GPS-005 repeated gene-set evidence has inconsistent effect or padj for priority_id(s): ", paste(names(groups)[inconsistent], collapse = ", "), ". Repair: supply one consistent gene-level value per priority_id; lisaR will not select or aggregate DE statistics.", call. = FALSE)
  out <- do.call(rbind, lapply(groups, function(x) {
    effect_signed <- x$effect[[1]]
    padj_raw <- x$padj[[1]]
    padj_for_confidence <- if (isTRUE(padj_raw == 0)) .Machine$double.xmin else padj_raw
    data.frame(
      priority_id = x$priority_id[[1]],
      raw_gene_set_count = length(unique(x$gene_set_id)),
      independent_category_breadth = length(unique(x$category_id)),
      gene_set_ids = paste(sort(unique(x$gene_set_id)), collapse = ";"),
      category_ids = paste(sort(unique(x$category_id)), collapse = ";"),
      effect_signed_raw = effect_signed,
      effect_strength_raw = abs(effect_signed),
      padj_raw = padj_raw,
      padj_for_confidence = padj_for_confidence,
      zero_padj_clamped = !is.na(padj_raw) && padj_raw == 0,
      de_confidence_raw = if (is.na(padj_for_confidence)) NA_real_ else -log10(padj_for_confidence),
      stringsAsFactors = FALSE
    )
  }))
  out <- out[order(out$priority_id), , drop = FALSE]
  rownames(out) <- NULL
  effect_component <- lisa_pathway_minmax(out$effect_strength_raw)
  confidence_component <- lisa_pathway_minmax(out$de_confidence_raw)
  breadth_component_rows <- lisa_pathway_minmax(out$independent_category_breadth)
  mode <- tolower(as.character(evidence_mode[[1]]))
  applicable <- identical(mode, "full_de")
  component_complete <- is.finite(effect_component) & is.finite(confidence_component) & is.finite(breadth_component_rows)
  score <- rep(NA_real_, nrow(out))
  if (applicable) score[component_complete] <- (effect_component[component_complete] * confidence_component[component_complete] * breadth_component_rows[component_complete])^(1 / 3)
  out$effect_strength_normalized <- effect_component
  out$de_confidence_normalized <- confidence_component
  out$independent_category_breadth_normalized <- breadth_component_rows
  out$lisa_gps <- score
  out$evidence_mode <- mode
  out$applicability <- if (applicable) "applicable" else "not_applicable"
  out$interpretation <- "prioritization heuristic; not inferential or causal"
  out$pareto_front <- lisa_gps_pareto_front(out)
  out$pareto_role <- ifelse(is.na(out$pareto_front), "not_applicable", "robustness_metadata_only")
  out
}

#' Calculate Pareto-front metadata for a lisa-gps table
#'
#' Pareto fronts are metadata only and must not replace `lisa_gps` ordering.
#' @param score_table Output of [lisa_gps_prioritization()].
#' @return Integer Pareto fronts, or `NA` for incomplete/not-applicable rows.
#' @keywords internal
lisa_gps_pareto_front <- function(score_table) {
  x <- as.data.frame(score_table, stringsAsFactors = FALSE)
  cols <- c("effect_strength_normalized", "de_confidence_normalized", "independent_category_breadth_normalized")
  if (!all(cols %in% names(x))) stop("LISA-GPS-004 score table lacks normalized components.", call. = FALSE)
  valid <- stats::complete.cases(x[, cols, drop = FALSE])
  if ("applicability" %in% names(x)) valid <- valid & x$applicability == "applicable"
  front <- rep(NA_integer_, nrow(x))
  idx <- which(valid)
  if (!length(idx)) return(front)

  # Exact three-dimensional skyline in O(n log n). The earlier development
  # implementation calculated every Pareto layer by repeated all-pairs
  # comparisons, which is not usable for real transcriptomic datasets. For
  # the approved robustness role only first-front membership is required.
  values <- as.matrix(x[idx, cols, drop = FALSE])
  ord <- order(-values[, 1], -values[, 2], -values[, 3], idx)
  values <- values[ord, , drop = FALSE]
  original <- idx[ord]
  y_levels <- sort(unique(values[, 2]), decreasing = TRUE)
  y_rank <- match(values[, 2], y_levels)
  bit <- rep(-Inf, length(y_levels))
  bit_query <- function(position) {
    answer <- -Inf
    while (position > 0L) {
      answer <- max(answer, bit[[position]])
      position <- position - bitwAnd(position, -position)
    }
    answer
  }
  bit_update <- function(position, value) {
    while (position <= length(bit)) {
      bit[[position]] <<- max(bit[[position]], value)
      position <- position + bitwAnd(position, -position)
    }
  }

  group <- interaction(values[, 1], values[, 2], values[, 3], drop = TRUE, lex.order = TRUE)
  starts <- c(1L, which(group[-1L] != group[-length(group)]) + 1L)
  ends <- c(starts[-1L] - 1L, length(group))
  for (g in seq_along(starts)) {
    rows <- starts[[g]]:ends[[g]]
    dominated <- bit_query(y_rank[[rows[[1L]]]]) >= values[rows[[1L]], 3]
    if (!dominated) front[original[rows]] <- 1L
    bit_update(y_rank[[rows[[1L]]]], values[rows[[1L]], 3])
  }
  front
}

#' Create a complete PATHWAYS category-status universe
#'
#' @param dictionary_categories Versioned PATHWAYS dictionary categories.
#' @param tested_results Optional tested result rows keyed by `category_id`.
#' @param padj_cutoff Existing analysis cutoff used only to label tested rows.
#' @param category_display Display policy; testing default is `"all"`.
#' @param expected_categories Expected active dictionary breadth (69 after
#'   removing zero-assignment categories from the approved map).
#' @return One status row for every dictionary category.
#' @keywords internal
lisa_pathways_status_table <- function(dictionary_categories, tested_results = NULL,
                                       padj_cutoff = 0.05, category_display = "all",
                                       expected_categories = 69L) {
  dict <- as.data.frame(dictionary_categories, stringsAsFactors = FALSE)
  if (!"category_id" %in% names(dict) || any(!nzchar(as.character(dict$category_id))) || anyDuplicated(dict$category_id)) {
    stop("LISA-PATHWAYS-001 dictionary categories require unique non-empty category_id values.", call. = FALSE)
  }
  if (!identical(tolower(as.character(category_display[[1]])), "all")) stop("LISA-PATHWAYS-002 category_display must be 'all' for the pathway testing default.", call. = FALSE)
  if (!is.null(expected_categories) && nrow(dict) != as.integer(expected_categories)) {
    stop("LISA-PATHWAYS-003 dictionary category count is ", nrow(dict), "; expected ", expected_categories, ". Do not modify the approved dictionary contract.", call. = FALSE)
  }
  out <- dict
  if (!"category_display_name" %in% names(out)) out$category_display_name <- out$category_id
  out$category_display <- "all"
  out$status <- "not_tested"; out$padj <- NA_real_; out$tested <- FALSE
  if (is.null(tested_results) || nrow(as.data.frame(tested_results)) == 0) return(out)
  result <- as.data.frame(tested_results, stringsAsFactors = FALSE)
  if (!"category_id" %in% names(result) || anyDuplicated(result$category_id)) stop("LISA-PATHWAYS-004 tested results require unique category_id values.", call. = FALSE)
  unknown <- setdiff(result$category_id, out$category_id)
  if (length(unknown)) stop("LISA-PATHWAYS-005 tested results contain category IDs absent from the versioned dictionary: ", paste(unknown, collapse = ", "), call. = FALSE)
  hit <- match(out$category_id, result$category_id)
  present <- !is.na(hit)
  if ("available" %in% names(result)) {
    available_result <- as.logical(result$available[hit[present]])
    out$status[which(present)[!available_result %in% TRUE]] <- "not_available"
  }
  if ("tested" %in% names(result)) out$tested[present] <- as.logical(result$tested[hit[present]]) else out$tested[present] <- TRUE
  if ("padj" %in% names(result)) out$padj[present] <- suppressWarnings(as.numeric(result$padj[hit[present]]))
  available <- out$status != "not_available"
  tested <- out$tested & available
  out$status[tested & !is.na(out$padj) & out$padj <= padj_cutoff] <- "significant"
  out$status[tested & (is.na(out$padj) | out$padj > padj_cutoff)] <- "tested_not_significant"
  out
}

#' Select KEGG maps from existing enrichment results
#'
#' This performs no independent KEGG enrichment.
#' @param enrichment Existing enrichment result rows with `padj`.
#' @param species Homo sapiens or Mus musculus.
#' @param padj_cutoff Existing enrichment cutoff.
#' @return Significant KEGG map candidates only.
#' @keywords internal
lisa_kegg_map_selection <- function(enrichment, species, padj_cutoff = 0.05) {
  lisa_species_contract(species)
  x <- as.data.frame(enrichment, stringsAsFactors = FALSE)
  if (!"padj" %in% names(x)) stop("LISA-KEGG-011 enrichment requires padj for significant-only KEGG selection.", call. = FALSE)
  x$padj <- suppressWarnings(as.numeric(x$padj))
  x$selected_for_kegg_map <- !is.na(x$padj) & x$padj <= padj_cutoff
  x[x$selected_for_kegg_map, , drop = FALSE]
}

lisa_assert_selected_kegg_painter_outputs <- function(index, node_table) {
  index <- as.data.frame(index, stringsAsFactors = FALSE)
  node_table <- as.data.frame(node_table, stringsAsFactors = FALSE)
  required <- c("kegg_id", "kegg_title", "output_png", "output_pdf")
  if (!all(required %in% names(index))) {
    stop("LISA-KEGG-026 selected KEGG painter index is incomplete.", call. = FALSE)
  }
  safe <- function(x) {
    x <- as.character(x)
    x[is.na(x)] <- ""
    x
  }
  missing <- which(!nzchar(safe(index$output_png)) | !nzchar(safe(index$output_pdf)))
  if (!length(missing)) return(invisible(TRUE))
  labels <- paste0(safe(index$kegg_id[missing]), " (", safe(index$kegg_title[missing]), ")")
  failures <- if ("error" %in% names(node_table)) {
    unique(safe(node_table$error[nzchar(safe(node_table$error))]))
  } else character()
  stop(
    "LISA-KEGG-026 required selected KEGG pathway renders are absent: ",
    paste(labels, collapse = ", "),
    if (length(failures)) paste0(". Cache/render diagnostics: ", paste(failures, collapse = " | ")) else "",
    ". Supply every selected KGML/PNG resource in the declared immutable snapshot or change the upstream significant evidence; lisaR does not silently omit required selected maps.",
    call. = FALSE
  )
}

#' Build a dual A/B KEGG node table with one common display scale
#'
#' @param side_a Side A gene values (`gene_id`, `log2fc`).
#' @param side_b Side B gene values (`gene_id`, `log2fc`).
#' @param mappings Offline mapping rows (`gene_id`, `node_id`); one-to-many is retained.
#' @param species Homo sapiens or Mus musculus.
#' @return Node-level dual-map data. Delta is data only and is never a paint value.
#' @keywords internal
lisa_kegg_dual_node_table <- function(side_a, side_b, mappings, species) {
  contract <- lisa_species_contract(species)
  normalize_side <- function(x, side) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    if (!all(c("gene_id", "log2fc") %in% names(x))) stop("LISA-KEGG-012 ", side, " requires gene_id and log2fc.", call. = FALSE)
    x$gene_id <- as.character(x$gene_id)
    if (any(is.na(x$gene_id) | !nzchar(x$gene_id))) stop("LISA-KEGG-012 ", side, " requires non-empty gene_id values.", call. = FALSE)
    if (anyDuplicated(x$gene_id)) stop("LISA-KEGG-015 ", side, " contains duplicate gene_id values. Repair: resolve gene-level values before KEGG node summarization.", call. = FALSE)
    x$log2fc <- suppressWarnings(as.numeric(x$log2fc)); x$side <- side; x
  }
  maps <- as.data.frame(mappings, stringsAsFactors = FALSE)
  if (!all(c("gene_id", "node_id") %in% names(maps)) || any(!nzchar(as.character(maps$node_id)))) stop("LISA-KEGG-013 mappings require non-empty gene_id and node_id.", call. = FALSE)
  maps$gene_id <- as.character(maps$gene_id); maps$node_id <- as.character(maps$node_id)
  if (anyDuplicated(maps[, c("gene_id", "node_id")])) stop("LISA-KEGG-014 duplicate gene-to-node mappings are ambiguous; retain one explicit mapping row.", call. = FALSE)
  values <- rbind(normalize_side(side_a, "A"), normalize_side(side_b, "B"))
  joined <- merge(maps, values, by = "gene_id", all.x = FALSE, all.y = FALSE, sort = FALSE)
  if (!nrow(joined)) return(data.frame())
  node_side <- aggregate(log2fc ~ node_id + side, joined, function(x) mean(x[is.finite(x)]))
  wide <- reshape(node_side, idvar = "node_id", timevar = "side", direction = "wide")
  names(wide) <- sub("log2fc\\.", "log2fc_", names(wide))
  for (name in c("log2fc_A", "log2fc_B")) if (!name %in% names(wide)) wide[[name]] <- NA_real_
  counts <- aggregate(gene_id ~ node_id + side, joined, function(x) length(unique(x)))
  counts <- reshape(counts, idvar = "node_id", timevar = "side", direction = "wide")
  names(counts) <- sub("gene_id\\.", "mapped_gene_count_", names(counts))
  out <- merge(wide, counts, by = "node_id", all = TRUE, sort = TRUE)
  for (name in c("mapped_gene_count_A", "mapped_gene_count_B")) if (!name %in% names(out)) out[[name]] <- 0L
  out$mapping_rows <- vapply(out$node_id, function(id) sum(joined$node_id == id), integer(1))
  out$mapped_gene_ids <- vapply(out$node_id, function(id) paste(sort(unique(joined$gene_id[joined$node_id == id])), collapse = ";"), character(1))
  out$multi_id_mapping_present <- vapply(out$node_id, function(id) length(unique(joined$gene_id[joined$node_id == id])) > 1L, logical(1))
  out$compound_node <- grepl("^(cpd:|C[0-9]{5}$)", out$node_id)
  out$delta_log2fc_data_only <- out$log2fc_A - out$log2fc_B
  out$common_scale_max_abs_log2fc <- max(abs(c(out$log2fc_A, out$log2fc_B)), na.rm = TRUE)
  if (!is.finite(out$common_scale_max_abs_log2fc[[1]])) out$common_scale_max_abs_log2fc <- NA_real_
  out$paint_encoding <- "dual A/B values on common scale; delta_log2fc_data_only is not painted"
  out$node_summarization <- "arithmetic mean of finite, unique gene-level log2fc values per node and side"
  out$species <- contract$scientific_name
  out
}
