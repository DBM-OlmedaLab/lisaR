# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Category evidence is descriptive: it never combines enrichment and DE into a
# prioritisation score. Identifiers are joined exactly, without case conversion.
lisa_evidence_read <- function(x) {
  if (is.null(x)) return(data.frame())
  if (is.character(x) && length(x) == 1L) {
    return(utils::read.delim(x, sep = "\t", quote = "", comment.char = "",
      check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character"))
  }
  as.data.frame(x, stringsAsFactors = FALSE)
}

lisa_evidence_column <- function(x, choices, default = "") {
  found <- intersect(choices, names(x))
  if (!length(found)) return(rep(default, nrow(x)))
  x[[found[[1L]]]]
}

lisa_evidence_num <- function(x) suppressWarnings(as.numeric(x))

# Preserve the exact IEEE double when a heterogeneous figure source must store
# numeric fields as text (ordinary as.character rounds to 15 digits).
lisa_evidence_text <- function(x) if (is.numeric(x)) ifelse(is.na(x), NA_character_, sprintf("%.17g", x)) else as.character(x)

lisa_evidence_genes <- function(x) {
  if (is.list(x)) x <- unlist(x, use.names = FALSE)
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(character())
  x <- gsub('^c\\(|\\)$|["\']', "", x)
  out <- trimws(unlist(strsplit(x, "[,;|/[:space:]]+"), use.names = FALSE))
  sort(unique(out[nzchar(out)]), method = "radix")
}

lisa_evidence_order <- function(...) order(..., na.last = TRUE, method = "radix")

#' Build descriptive category evidence from complete GSEA results
#'
#' The full annotated GSEA table is authoritative. An optional universe ledger
#' preserves eligibility states, and dictionary assignments can restore mapped
#' sets absent from older, filtered result files. No ORA evidence is consumed.
#' The return value contains normalised tables, not a gene-prioritisation score.
#' @param gsea_table Full annotated GSEA data frame or TSV path.
#' @param de_table Standardised DE table or TSV path; may be absent for ranks-only input.
#' @param universe_ledger Optional complete GSEA universe ledger.
#' @param category_dictionary Optional dictionary/category metadata with category_id.
#' @param memberships Optional complete TERM2GENE table (pathway/symbol or gs_name/gene_symbol).
#' @param analysis_id,collection,tier Exact scope labels; one analysis and tier only.
#' @param gsea_padj_cutoff,de_padj_cutoff Adjusted-p-value thresholds.
#' @param positive_contrast Meaning of positive NES and positive DE effect, if recorded.
#' @param max_sets,max_genes Positive display caps, independent of table completeness.
#' @return A lisa_category_evidence object containing normalised evidence tables.
build_lisa_category_evidence <- function(gsea_table, de_table = NULL,
    universe_ledger = NULL, category_dictionary = NULL, memberships = NULL,
    analysis_id = "", collection = "", tier = "not recorded",
    gsea_padj_cutoff = 0.25, de_padj_cutoff = 0.05,
    positive_contrast = "not recorded", max_sets = 25L, max_genes = 40L) {
  for (value in list(gsea_padj_cutoff, de_padj_cutoff)) {
    if (length(value) != 1L || !is.finite(value) || value < 0 || value > 1)
      stop("Evidence FDR cutoffs must be finite numbers in [0, 1].", call. = FALSE)
  }
  for (value in list(max_sets, max_genes)) {
    if (length(value) != 1L || !is.finite(value) || value < 1 || value != floor(value))
      stop("Evidence display caps must be positive integers.", call. = FALSE)
  }
  for (value in list(analysis_id, collection, tier, positive_contrast)) {
    if (length(value) != 1L || is.na(value)) stop("Evidence scope labels must be single strings.", call. = FALSE)
  }
  g <- lisa_evidence_read(gsea_table)
  ledger <- lisa_evidence_read(universe_ledger)
  dictionary <- lisa_evidence_read(category_dictionary)
  de <- lisa_evidence_read(de_table)
  member <- lisa_evidence_read(memberships)
  if (!all(c("pathway", "category_id") %in% names(g)))
    stop("Full annotated GSEA requires pathway and category_id columns.", call. = FALSE)
  if (nrow(dictionary) && !"category_id" %in% names(dictionary))
    stop("Category dictionary requires category_id.", call. = FALSE)
  # Scope metadata is evidence, not a label to overwrite. Check every supplied
  # source before assignment recovery or rendering. An expanded dictionary may
  # contain core assignments; the reverse would silently widen a core request.
  scope_tables <- list(gsea_table = g, universe_ledger = ledger,
    category_dictionary = dictionary)
  scope_values <- function(columns) {
    values <- unlist(lapply(scope_tables, function(x) {
      unlist(lapply(intersect(columns, names(x)), function(column) as.character(x[[column]])), use.names = FALSE)
    }), use.names = FALSE)
    unique(values[!is.na(values) & nzchar(trimws(values))])
  }
  observed_tiers <- scope_values("tier")
  if (identical(tier, "not recorded") || !nzchar(tier)) {
    if (length(observed_tiers) > 1L) stop("Specify the analysis tier when dictionary row tiers differ.", call. = FALSE)
    if (length(observed_tiers) == 1L) tier <- observed_tiers[[1L]]
  }
  allowed_tiers <- if (identical(tier, "expanded")) c("core", "expanded") else tier
  incompatible_tiers <- setdiff(observed_tiers, allowed_tiers)
  if (length(incompatible_tiers)) {
    stop("Evidence tier scope '", tier, "' conflicts with input row tiers: ",
      paste(incompatible_tiers, collapse = ", "), ". Supply inputs from the requested tier.", call. = FALSE)
  }
  observed_collections <- scope_values(c("analysis_collection", "universe"))
  if (!nzchar(collection) && length(observed_collections) == 1L) collection <- observed_collections[[1L]]
  if (length(setdiff(observed_collections, collection))) {
    stop("Evidence collection scope '", collection, "' conflicts with input analysis_collection/universe: ",
      paste(observed_collections, collapse = ", "), ". Supply one matching collection.", call. = FALSE)
  }
  meta_cols <- c("category_display_name", "category_description", "macrogroup_id", "macrogroup_name", "color")
  aliases <- list(category_display_name = c("category_display_name", "display_name"),
    category_description = c("category_description", "description", "definition"),
    macrogroup_id = "macrogroup_id", macrogroup_name = "macrogroup_name", color = "color")
  assignments <- unique(data.frame(category_id = as.character(g$category_id),
    pathway = as.character(g$pathway), stringsAsFactors = FALSE))
  dict_set <- intersect(c("gene_set_id", "pathway", "gs_name"), names(dictionary))
  if (length(dict_set)) {
    extra <- data.frame(category_id = as.character(dictionary$category_id), pathway = as.character(dictionary[[dict_set[[1L]]]]))
    if (nrow(ledger) && "pathway" %in% names(ledger)) extra <- extra[extra$pathway %in% ledger$pathway, , drop = FALSE]
    assignments <- unique(rbind(assignments, extra))
  }
  valid_category <- function(x) !is.na(x) & nzchar(x) & x != "OTHER_UNCLASSIFIED"
  excluded_rows <- sum(!valid_category(as.character(g$category_id)))
  assignments <- assignments[valid_category(assignments$category_id) & !is.na(assignments$pathway) & nzchar(assignments$pathway), , drop = FALSE]
  assignments <- assignments[lisa_evidence_order(assignments$category_id, assignments$pathway), , drop = FALSE]
  category_ids <- sort(unique(c(assignments$category_id, as.character(dictionary$category_id))), method = "radix")
  category_ids <- category_ids[valid_category(category_ids)]
  categories <- data.frame(category_id = category_ids, stringsAsFactors = FALSE)
  for (column in meta_cols) {
    categories[[column]] <- vapply(category_ids, function(id) {
      vals <- c(as.character(lisa_evidence_column(dictionary[dictionary$category_id %in% id, , drop = FALSE], aliases[[column]])),
        as.character(lisa_evidence_column(g[g$category_id %in% id, , drop = FALSE], aliases[[column]])))
      vals <- unique(vals[!is.na(vals) & nzchar(vals)])
      if (length(vals)) vals[[1L]] else if (column == "category_display_name") id else ""
    }, character(1L))
  }
  # Enrichment statistics must agree across all category annotations of a set.
  # This detects accidentally concatenated analyses instead of silently mixing.
  g$pathway <- as.character(g$pathway)
  for (column in c("NES", "padj", "leadingEdge")) {
    if (!column %in% names(g)) next
    vals <- if (is.list(g[[column]])) vapply(g[[column]], function(x) paste(lisa_evidence_genes(x), collapse = ";"), character(1L)) else lisa_evidence_text(g[[column]])
    if (any(vapply(split(vals, g$pathway), function(x) length(unique(x[!is.na(x)])) > 1L, logical(1L))))
      stop("Conflicting GSEA ", column, " values for an exact pathway ID; provide one analysis/tier.", call. = FALSE)
  }
  if (nrow(ledger) && (!"pathway" %in% names(ledger) || anyDuplicated(ledger$pathway)))
    stop("Universe ledger requires unique pathway IDs.", call. = FALSE)
  sets <- assignments
  gi <- match(sets$pathway, g$pathway)
  li <- match(sets$pathway, ledger$pathway)
  sets$NES <- lisa_evidence_num(lisa_evidence_column(g, "NES", NA_real_))[gi]
  sets$gsea_fdr <- lisa_evidence_num(lisa_evidence_column(g, "padj", NA_real_))[gi]
  for (column in c("gsea_eligibility_status", "gsea_result_status", "source_family", "source_gene_count", "ranked_gene_count")) {
    fallback <- if (column == "source_family") c("source_family", "term2gene_source_family") else column
    value <- lisa_evidence_column(g, fallback, NA_character_)[gi]
    lv <- lisa_evidence_column(ledger, fallback, NA_character_)[li]
    missing <- is.na(value) | !nzchar(as.character(value))
    value[missing] <- lv[missing]
    sets[[column]] <- value
  }
  eligible <- lisa_evidence_column(g, "eligible_for_gsea", NA)[gi]
  ledger_eligible <- lisa_evidence_column(ledger, "eligible_for_gsea", NA)[li]
  eligible[is.na(eligible)] <- ledger_eligible[is.na(eligible)]
  sets$eligible_for_gsea <- as.logical(eligible)
  sets$evaluable <- is.finite(sets$NES) & is.finite(sets$gsea_fdr) & sets$gsea_fdr >= 0 & sets$gsea_fdr <= 1 & !(sets$eligible_for_gsea %in% FALSE)
  sets$significant <- sets$evaluable & sets$gsea_fdr <= gsea_padj_cutoff
  sets$state <- ifelse(!sets$evaluable, "non_evaluable", ifelse(sets$significant, "significant", "non_significant"))
  sets$non_evaluable_reason <- ifelse(sets$evaluable, "", ifelse(!is.na(sets$gsea_eligibility_status) & sets$gsea_eligibility_status != "eligible", sets$gsea_eligibility_status,
    ifelse(is.na(gi), "missing_result", "missing_or_invalid_statistics")))
  sets$direction <- ifelse(!sets$evaluable, "not_evaluable", ifelse(sets$NES > 0, "positive", ifelse(sets$NES < 0, "negative", "zero")))
  edges <- data.frame(category_id = character(), pathway = character(), symbol = character(), stringsAsFactors = FALSE)
  edge_list <- vector("list", nrow(sets))
  sets$leading_edge_state <- rep("not_recorded", nrow(sets))
  sets$n_leading_edge_genes <- integer(nrow(sets))
  for (i in seq_len(nrow(sets))) {
    present <- !is.na(gi[[i]]) && "leadingEdge" %in% names(g)
    raw <- if (present) g$leadingEdge[[gi[[i]]]] else NA_character_
    genes <- lisa_evidence_genes(raw)
    sets$leading_edge_state[[i]] <- if (!present || all(is.na(raw))) "not_recorded" else if (!length(genes)) "empty" else "recorded"
    sets$n_leading_edge_genes[[i]] <- length(genes)
    if (length(genes)) edge_list[[i]] <- data.frame(category_id = sets$category_id[[i]], pathway = sets$pathway[[i]], symbol = genes, stringsAsFactors = FALSE)
  }
  edge_list <- Filter(Negate(is.null), edge_list)
  if (length(edge_list)) edges <- do.call(rbind, edge_list)
  if (nrow(de) && !"symbol" %in% names(de)) stop("Standardised DE requires symbol for exact joins.", call. = FALSE)
  if (nrow(de) && anyDuplicated(as.character(de$symbol))) stop("Standardised DE has duplicate symbols; resolve them upstream.", call. = FALSE)
  if (nrow(de) && any(is.na(de$symbol) | !nzchar(as.character(de$symbol)))) stop("Standardised DE contains empty symbol IDs.", call. = FALSE)
  symbols <- sort(unique(c(as.character(de$symbol), edges$symbol)), method = "radix")
  di <- match(symbols, as.character(de$symbol))
  genes <- data.frame(symbol = symbols,
    gene_id = as.character(lisa_evidence_column(de, "gene_id", NA_character_))[di],
    log2FC = lisa_evidence_num(lisa_evidence_column(de, c("log2FC", "log2FoldChange", "logFC"), NA_real_))[di],
    de_fdr = lisa_evidence_num(lisa_evidence_column(de, "padj", NA_real_))[di],
    stringsAsFactors = FALSE)
  genes$de_state <- ifelse(is.na(di), "not_in_standardized_DE", ifelse(!is.finite(genes$log2FC) | !is.finite(genes$de_fdr) | genes$de_fdr < 0 | genes$de_fdr > 1, "not_evaluable", ifelse(genes$de_fdr <= de_padj_cutoff, "significant", "non_significant")))
  genes$in_leading_edge <- genes$symbol %in% edges$symbol
  # Raw complete membership remains separate from GSEA leading-edge membership.
  complete <- data.frame(pathway = character(), symbol = character(), stringsAsFactors = FALSE)
  if (nrow(member)) {
    pc <- intersect(c("pathway", "gs_name", "gene_set_id"), names(member))
    gc <- intersect(c("symbol", "gene_symbol"), names(member))
    if (!length(pc) || !length(gc)) stop("Complete memberships require pathway and symbol columns.", call. = FALSE)
    complete <- unique(data.frame(pathway = as.character(member[[pc[[1L]]]]), symbol = as.character(member[[gc[[1L]]]]), stringsAsFactors = FALSE))
    complete <- complete[complete$pathway %in% sets$pathway & !is.na(complete$symbol) & nzchar(complete$symbol), , drop = FALSE]
    complete <- complete[lisa_evidence_order(complete$pathway, complete$symbol), , drop = FALSE]
  }
  sets$display_order <- rep(NA_integer_, nrow(sets))
  selection <- data.frame(category_id = character(), symbol = character(), display_order = integer(), n_significant_set_connections = integer(), n_significant_sets_denominator = integer(), stringsAsFactors = FALSE)
  selection_list <- list()
  categories$mean_significant_set_NES <- rep(NA_real_, nrow(categories))
  categories$category_mean_source <- rep("derived_all_significant_unique_sets", nrow(categories))
  counts <- c("n_mapped_sets", "n_evaluable_sets", "n_significant_sets", "n_positive_sets", "n_negative_sets", "n_significant_positive_sets", "n_significant_negative_sets", "n_unique_le_genes", "n_le_assignments", "n_significant_unique_le_genes", "n_significant_le_assignments", "n_plotted_sets", "n_available_plot_genes", "n_plotted_genes")
  for (column in counts) categories[[column]] <- integer(nrow(categories))
  set_indices <- split(seq_len(nrow(sets)), sets$category_id)
  edge_indices <- split(seq_len(nrow(edges)), edges$category_id)
  for (i in seq_len(nrow(categories))) {
    id <- categories$category_id[[i]]
    si <- set_indices[[id]]
    if (is.null(si)) si <- integer()
    ei <- edge_indices[[id]]
    if (is.null(ei)) ei <- integer()
    ss <- sets[si, , drop = FALSE]
    ee <- edges[ei, , drop = FALSE]
    sig_edges <- ee[ee$pathway %in% ss$pathway[ss$significant], , drop = FALSE]
    plot_rows <- si[ss$significant]
    plot_rows <- head(plot_rows[lisa_evidence_order(sets$gsea_fdr[plot_rows], -abs(sets$NES[plot_rows]), sets$pathway[plot_rows])], max_sets)
    sets$display_order[plot_rows] <- seq_along(plot_rows)
    available <- sort(unique(ee$symbol[ee$pathway %in% sets$pathway[plot_rows]]), method = "radix")
    # Recurrence is counted over ALL significant unique sets in the category,
    # never only the capped visible rows. It is descriptive, not a gene score.
    recurrence <- table(sig_edges$symbol)
    available_counts <- as.integer(recurrence[available]); available_counts[is.na(available_counts)] <- 0L
    available <- available[lisa_evidence_order(-available_counts, available)]
    selected <- head(available, max_genes)
    if (any(ss$significant)) categories$mean_significant_set_NES[[i]] <- mean(ss$NES[ss$significant])
    # The pipeline supplies the canonical category-summary mean. Verify it,
    # rather than silently using a different visible-set mean for alignment.
    if ("mean_NES" %in% names(dictionary)) {
      canonical <- unique(lisa_evidence_num(dictionary$mean_NES[dictionary$category_id %in% id]))
      canonical <- canonical[is.finite(canonical)]
      if (length(canonical)) {
        derived <- categories$mean_significant_set_NES[[i]]
        if (!is.finite(derived) || any(abs(canonical - derived) > 1e-10 * max(1, abs(derived))))
          stop("Category mean_NES conflicts with all significant unique member sets for ", id, ". Check the analysis and FDR cutoff.", call. = FALSE)
        categories$mean_significant_set_NES[[i]] <- canonical[[1L]]
        categories$category_mean_source[[i]] <- "verified_category_summary_mean_NES"
      }
    }
    if (length(selected)) selection_list[[id]] <- data.frame(category_id = id, symbol = selected,
      display_order = seq_along(selected), n_significant_set_connections = as.integer(recurrence[selected]),
      n_significant_sets_denominator = sum(ss$significant), stringsAsFactors = FALSE)
    categories[i, counts] <- as.list(c(nrow(ss), sum(ss$evaluable), sum(ss$significant), sum(ss$direction == "positive"), sum(ss$direction == "negative"), sum(ss$significant & ss$direction == "positive"), sum(ss$significant & ss$direction == "negative"), length(unique(ee$symbol)), nrow(ee), length(unique(sig_edges$symbol)), nrow(sig_edges), length(plot_rows), length(available), length(selected)))
  }
  if (length(selection_list)) selection <- do.call(rbind, selection_list)
  categories$sets_truncated <- categories$n_plotted_sets < categories$n_significant_sets
  categories$genes_truncated <- categories$n_plotted_genes < categories$n_available_plot_genes
  categories$support_state <- ifelse(categories$n_mapped_sets == 0L, "no_mapped_sets", ifelse(categories$n_evaluable_sets == 0L, "no_evaluable_sets", ifelse(categories$n_significant_sets == 0L, "no_significant_sets", ifelse(categories$n_significant_unique_le_genes == 0L, "no_significant_leading_edge", ifelse(categories$n_significant_positive_sets > 0L & categories$n_significant_negative_sets > 0L, "mixed_direction_support", "directional_support")))))
  for (name in c("categories", "sets", "genes", "edges", "selection", "complete")) {
    value <- get(name); rownames(value) <- NULL; assign(name, value)
  }
  structure(list(metadata = list(analysis_id = analysis_id, collection = collection, tier = tier,
    gsea_padj_cutoff = gsea_padj_cutoff, de_padj_cutoff = de_padj_cutoff,
    positive_contrast = positive_contrast, max_sets = as.integer(max_sets), max_genes = as.integer(max_genes),
    set_selection = "Significant sets by GSEA FDR ascending, then absolute NES descending, then exact pathway ID; no combined score.",
    gene_selection = "Leading-edge recurrence across ALL significant category sets descending, then exact symbol; candidates are leading-edge genes of displayed sets.",
    excluded_unclassified_rows = excluded_rows, universe_ledger_supplied = nrow(ledger) > 0L,
    complete_membership_supplied = !is.null(memberships), schema_version = "1.1"),
    categories = categories, sets = sets, genes = genes, leading_edges = edges,
    gene_selection = selection, complete_memberships = complete), class = "lisa_category_evidence")
}

# Indexed gene postings calculate exact intersections without constructing a
# dense set-by-gene matrix. An emit callback streams at most n_sets rows at a
# time, even for explicitly requested all-pair exports of a large category.
lisa_category_evidence_overlaps <- function(evidence, category_id,
    membership = c("leading_edge", "complete"), displayed = FALSE, emit = NULL) {
  membership <- match.arg(membership)
  sets <- evidence$sets[evidence$sets$category_id %in% category_id, , drop = FALSE]
  if (isTRUE(displayed)) sets <- sets[!is.na(sets$display_order), , drop = FALSE]
  if (length(unique(sets$category_id)) > 1L) stop("Overlap comparisons are within one category only.", call. = FALSE)
  empty <- data.frame(category_id = character(), pathway_a = character(), pathway_b = character(),
    direction_a = character(), direction_b = character(), direction_relation = character(),
    membership_type = character(), pair_scope = character(), n_a = integer(), n_b = integer(), intersection_n = integer(),
    union_n = integer(), jaccard = numeric(), overlap_state = character(), stringsAsFactors = FALSE)
  if (nrow(sets) < 2L) return(empty)
  mm <- if (membership == "leading_edge") evidence$leading_edges[evidence$leading_edges$category_id %in% category_id, , drop = FALSE] else evidence$complete_memberships
  mm <- unique(mm[mm$pathway %in% sets$pathway, c("pathway", "symbol"), drop = FALSE])
  ids <- match(mm$pathway, sets$pathway)
  postings <- split(ids, mm$symbol)
  groups <- split(mm$symbol, ids)
  sizes <- tabulate(ids, nbins = nrow(sets))
  chunks <- if (is.null(emit)) vector("list", nrow(sets) - 1L) else NULL
  for (a in seq_len(nrow(sets) - 1L)) {
    aa <- groups[[as.character(a)]]
    touched <- unlist(postings[aa], use.names = FALSE)
    intersection <- tabulate(as.integer(touched), nbins = nrow(sets))
    b <- seq.int(a + 1L, nrow(sets))
    union <- sizes[[a]] + sizes[b] - intersection[b]
    da <- sets$direction[[a]]; db <- sets$direction[b]
    relation <- ifelse(!da %in% c("positive", "negative") | !db %in% c("positive", "negative"), "not_directional", ifelse(da == db, "same", "opposite"))
    defined <- sizes[[a]] > 0L & sizes[b] > 0L
    chunk <- data.frame(category_id = category_id, pathway_a = sets$pathway[[a]], pathway_b = sets$pathway[b],
      direction_a = da, direction_b = db, direction_relation = relation, membership_type = membership,
      pair_scope = if (isTRUE(displayed)) "default_displayed_sets" else "all_mapped_sets",
      n_a = sizes[[a]], n_b = sizes[b], intersection_n = intersection[b], union_n = union,
      jaccard = ifelse(defined, intersection[b] / union, NA_real_),
      overlap_state = ifelse(defined, "defined", "empty_or_unrecorded_membership"), stringsAsFactors = FALSE)
    if (is.null(emit)) chunks[[a]] <- chunk else emit(chunk)
  }
  if (!is.null(emit)) return(invisible(NULL))
  out <- do.call(rbind, chunks); rownames(out) <- NULL; out
}

lisa_category_evidence_figure_source <- function(evidence, category_id) {
  ss <- evidence$sets[evidence$sets$category_id == category_id & !is.na(evidence$sets$display_order), , drop = FALSE]
  ss <- ss[order(ss$display_order), , drop = FALSE]
  pick <- evidence$gene_selection[evidence$gene_selection$category_id == category_id, , drop = FALSE]
  gg <- evidence$genes[match(pick$symbol, evidence$genes$symbol), , drop = FALSE]
  category_edges <- evidence$leading_edges[evidence$leading_edges$category_id == category_id, , drop = FALSE]
  cat <- evidence$categories[evidence$categories$category_id == category_id, , drop = FALSE]
  columns <- c("figure_id", "figure_type", "row_type", "category_id", "category_display_name", "pathway", "symbol", "set_order", "gene_order", "NES", "gsea_fdr", "log2FC", "de_fdr", "de_state", "leading_edge_state", "n_leading_edge_genes", "n_displayed_leading_edge_genes", "analysis_id", "collection", "tier", "positive_contrast", "gsea_padj_cutoff", "de_padj_cutoff", "n_mapped_sets", "n_evaluable_sets", "n_significant_sets", "n_plotted_sets", "n_available_plot_genes", "n_plotted_genes", "n_significant_set_connections", "n_significant_sets_denominator", "mean_significant_set_NES", "set_sort_label", "gene_sort_label", "selection_mode", "n_matching_sets")
  base <- as.data.frame(stats::setNames(rep(list(""), length(columns)), columns), stringsAsFactors = FALSE)
  base$figure_id <- paste0("category_evidence__", category_id)
  base$figure_type <- "category_evidence"
  base$category_id <- category_id; base$category_display_name <- cat$category_display_name
  for (name in intersect(names(evidence$metadata), columns)) base[[name]] <- evidence$metadata[[name]]
  for (name in intersect(names(cat), columns)) base[[name]] <- cat[[name]]
  base$set_sort_label <- "GSEA FDR, then |NES|, then exact ID"
  base$gene_sort_label <- "LE recurrence across all significant category sets, then exact symbol"
  base$selection_mode <- "significant sets"; base$n_matching_sets <- cat$n_significant_sets
  rows <- list(transform(base, row_type = "metadata"))
  if (nrow(ss)) for (i in seq_len(nrow(ss))) {
    row <- base; row$row_type <- "set"; row$pathway <- ss$pathway[[i]]
    row$set_order <- ss$display_order[[i]]; row$NES <- ss$NES[[i]]; row$gsea_fdr <- ss$gsea_fdr[[i]]
    row$leading_edge_state <- ss$leading_edge_state[[i]]
    row$n_leading_edge_genes <- ss$n_leading_edge_genes[[i]]
    row$n_displayed_leading_edge_genes <- sum(category_edges$pathway == ss$pathway[[i]] & category_edges$symbol %in% gg$symbol)
    rows[[length(rows) + 1L]] <- row
  }
  if (nrow(gg)) for (i in seq_len(nrow(gg))) {
    row <- base; row$row_type <- "gene"; row$symbol <- gg$symbol[[i]]; row$gene_order <- pick$display_order[[i]]
    row$n_significant_set_connections <- pick$n_significant_set_connections[[i]]
    row$n_significant_sets_denominator <- pick$n_significant_sets_denominator[[i]]
    row$log2FC <- gg$log2FC[[i]]; row$de_fdr <- gg$de_fdr[[i]]; row$de_state <- gg$de_state[[i]]
    rows[[length(rows) + 1L]] <- row
  }
  ee <- category_edges[category_edges$pathway %in% ss$pathway & category_edges$symbol %in% gg$symbol, , drop = FALSE]
  if (nrow(ee)) {
    cells <- base[rep(1L, nrow(ee)), , drop = FALSE]; cells$row_type <- "membership"
    cells$pathway <- ee$pathway; cells$symbol <- ee$symbol
    cells$set_order <- match(ee$pathway, ss$pathway); cells$gene_order <- match(ee$symbol, gg$symbol)
    rows[[length(rows) + 1L]] <- cells
  }
  rows <- lapply(rows, function(row) { row[] <- lapply(row, lisa_evidence_text); row })
  out <- do.call(rbind, rows); rownames(out) <- NULL; out$source_row_order <- seq_len(nrow(out)); out
}

# Intentionally self-contained: the identical function is copied into an offline
# base-R recipe. Source numbers and colours therefore reproduce the stored view.
lisa_category_evidence_draw <- function(source) {
  num <- function(x) suppressWarnings(as.numeric(x))
  fmt <- function(x) ifelse(is.finite(num(x)), formatC(num(x), format = "g", digits = 2), "NA")
  device_text <- function(x) gsub("[\u2013\u2014\u2212]", "-", as.character(x), perl = TRUE)
  label <- function(x, cap = 55L) ifelse(nchar(x) > cap, paste0(substr(x, 1L, cap - 3L), "..."), x)
  meta <- source[source$row_type == "metadata", , drop = FALSE][1L, , drop = FALSE]
  ss <- source[source$row_type == "set", , drop = FALSE]
  gg <- source[source$row_type == "gene", , drop = FALSE]
  ee <- source[source$row_type == "membership", , drop = FALSE]
  ss <- ss[order(num(ss$set_order)), , drop = FALSE]; gg <- gg[order(num(gg$gene_order)), , drop = FALSE]
  ns <- nrow(ss); ng <- nrow(gg)
  old <- graphics::par(mar = c(9, 23, 8, 3), xpd = NA)
  on.exit(graphics::par(old), add = TRUE)
  annotation_width <- max(10, ng * .5)
  annotation_x <- ng + annotation_width * c(.17, .5, .84)
  graphics::plot.new(); graphics::plot.window(xlim = c(.5, max(ng, 1) + annotation_width), ylim = c(.5, max(ns, 1) + 5))
  graphics::title(main = device_text(paste0(meta$category_display_name, " [", meta$category_id, "]")), line = 5.3, cex.main = 1)
  graphics::mtext(paste(meta$analysis_id, meta$collection, paste0("tier: ", meta$tier), sep = " | "), side = 3, line = 4.2, cex = .75)
  graphics::mtext(paste0("Positive contrast: ", meta$positive_contrast, "; GSEA FDR <= ", meta$gsea_padj_cutoff), side = 3, line = 3.1, cex = .72)
  graphics::mtext(paste0("Mapped/evaluable/significant sets: ", meta$n_mapped_sets, "/", meta$n_evaluable_sets, "/", meta$n_significant_sets,
    "; displayed sets ", meta$n_plotted_sets, "/", if ("n_matching_sets" %in% names(meta)) meta$n_matching_sets else meta$n_significant_sets, ", genes ", meta$n_plotted_genes, "/", meta$n_available_plot_genes), side = 3, line = 2, cex = .72)
  if (!ns || !ng) {
    graphics::text(1, 1, if (!ns) "No sets are displayed in this view.\nAll mapped set states remain in the tables." else "Displayed sets have no recorded leading-edge genes.\nSee the full set table.", adj = 0, cex = .9)
    return(invisible(NULL))
  }
  yy <- ns + 1L - seq_len(ns)
  for (i in seq_len(ns)) {
    graphics::rect(.5, yy[[i]] - .48, ng + .5, yy[[i]] + .48, col = if (i %% 2L) "#f1f5f9" else "white", border = NA)
    graphics::text(.35, yy[[i]], label(ss$pathway[[i]]), adj = 1, cex = .63)
    graphics::text(annotation_x[[1L]], yy[[i]], fmt(ss$NES[[i]]), col = if (!is.finite(num(ss$NES[[i]]))) "#9ca3af" else if (num(ss$NES[[i]]) > 0) "#b2182b" else if (num(ss$NES[[i]]) < 0) "#2166ac" else "#64748b", cex = .72)
    graphics::text(annotation_x[[2L]], yy[[i]], fmt(ss$gsea_fdr[[i]]), cex = .63)
    graphics::text(annotation_x[[3L]], yy[[i]], paste0(ss$n_displayed_leading_edge_genes[[i]], "/", ss$n_leading_edge_genes[[i]]), cex = .63)
  }
  if (nrow(ee)) {
    ii <- match(ee$pathway, ss$pathway)
    graphics::points(num(ee$gene_order), yy[ii], pch = 15, cex = .95,
      col = ifelse(num(ss$NES[ii]) > 0, "#b2182b", ifelse(num(ss$NES[ii]) < 0, "#2166ac", "#64748b")))
  }
  graphics::text(annotation_x, ns + 1.1, c("Set\nNES", "GSEA\nFDR", "LE genes\nshown / full"), cex = .68, font = 2)
  graphics::text(seq_len(ng), .35, gg$symbol, srt = 65, adj = 1, cex = .63)
  fc <- num(gg$log2FC); q <- num(gg$de_fdr)
  maxfc <- max(abs(fc[is.finite(fc)]), 1)
  colors <- vapply(fc, function(x) {
    if (!is.finite(x)) return("#d1d5db")
    rgb <- grDevices::colorRamp(c("white", if (x >= 0) "#b2182b" else "#2166ac"))(min(abs(x) / maxfc, 1))
    grDevices::rgb(rgb[[1L]], rgb[[2L]], rgb[[3L]], maxColorValue = 255)
  }, character(1L))
  graphics::rect(seq_len(ng) - .46, ns + 1.6, seq_len(ng) + .46, ns + 2.6, col = colors, border = "white")
  graphics::text(seq_len(ng), ns + 2.1, fmt(fc), srt = 90, cex = .52, col = ifelse(is.finite(fc) & abs(fc) > maxfc * .65, "white", "#111827"))
  graphics::text(seq_len(ng), ns + 3.2, fmt(q), srt = 90, cex = .5)
  graphics::text(.35, ns + c(2.1, 3.2), c("Gene DE log2FC", "Gene DE FDR"), adj = 1, cex = .7)
  if ("n_significant_set_connections" %in% names(gg)) {
    graphics::text(seq_len(ng), ns + 4.2, gg$n_significant_set_connections, cex = .53)
    graphics::text(.35, ns + 4.2, paste0("LE recurrence / ", meta$n_significant_sets, " significant sets"), adj = 1, cex = .65)
  }
  graphics::mtext("Red = positive; blue = negative; grey = missing. Dots: set NES direction. DE track: gene log2FC; gene FDR is separate.", side = 1, line = 5.5, cex = .66)
  set_label <- if ("set_sort_label" %in% names(meta)) meta$set_sort_label else "FDR, then |NES|, then ID"
  gene_label <- if ("gene_sort_label" %in% names(meta)) meta$gene_sort_label else "recurrence in ALL significant category sets, then symbol"
  graphics::mtext(paste0("Sets: ", set_label, ". Genes: ", gene_label, "."), side = 1, line = 6.7, cex = .62)
  graphics::mtext("A blank cell is not a recorded leading-edge connection; full membership may differ. Full TSVs retain all rows.", side = 1, line = 7.8, cex = .64)
  invisible(NULL)
}

lisa_category_evidence_save_figure <- function(source, stem, formats) {
  ns <- sum(source$row_type == "set"); ng <- sum(source$row_type == "gene")
  width <- max(14, min(25, 11 + ng * .3)); height <- max(8, 5 + ns * .25)
  for (format in formats) {
    lisa_guarded_write(paste0(stem, ".", format), function(target) {
      if (format == "png") grDevices::png(target, width = width, height = height, units = "in", res = 160, bg = "white")
      if (format == "pdf") grDevices::pdf(target, width = width, height = height, useDingbats = FALSE)
      if (format == "svg") grDevices::svg(target, width = width, height = height, bg = "white")
      tryCatch(lisa_category_evidence_draw(source), finally = grDevices::dev.off())
    })
  }
}

# Validate the small, presentation-only attachment contract shared with the
# offline evidence viewer. The normal FULL assembler may bind canonical
# `outputs/` or `artifacts/` references after this renderer has completed; this
# function never derives a category from an asset name.
lisa_category_evidence_category_products <- function(category_products, category_ids) {
  if (is.null(category_products)) return(list())
  if (!is.list(category_products)) stop("category_products must be a list.", call. = FALSE)
  out <- lapply(category_products, function(product) {
    if (!is.list(product) || !all(c("category_id", "product", "label", "assets") %in% names(product)))
      stop("Each category product needs category_id, product, label and assets.", call. = FALSE)
    id <- as.character(product$category_id)
    if (length(id) != 1L || is.na(id) || !id %in% category_ids)
      stop("A category product has no exact evidence category ID.", call. = FALSE)
    assets <- product$assets
    if (!is.list(assets) || !length(assets)) stop("Each category product needs at least one asset.", call. = FALSE)
    assets <- lapply(assets, function(asset) {
      if (!is.list(asset) || !all(c("format", "href") %in% names(asset)))
        stop("Each category product asset needs format and href.", call. = FALSE)
      href <- as.character(asset$href); format <- as.character(asset$format)
      if (!lisa_category_product_href_safe(href) || length(format) != 1L || is.na(format) || !nzchar(format))
        stop("A category product asset is not a safe relative file.", call. = FALSE)
      list(format = format, href = href)
    })
    product_name <- as.character(product$product); label <- as.character(product$label)
    if (length(product_name) != 1L || is.na(product_name) || !nzchar(product_name) ||
        length(label) != 1L || is.na(label) || !nzchar(label))
      stop("Each category product needs nonempty product and label.", call. = FALSE)
    list(category_id = id, product = product_name, label = label, assets = assets)
  })
  unname(out)
}

#' Render a portable, offline category-evidence view
#' @param evidence Object returned by build_lisa_category_evidence().
#' @param output_dir Destination directory.
#' @param formats Optional subset of png, pdf, svg; character() disables figures.
#' @param source_data Write source TSV sidecars for static figures.
#' @param recipes Write a shared standalone base-R reproduction script.
#' @param categories Optional exact category IDs to display; NULL displays all.
#' @param category_products Optional exact-ID product references supplied by
#'   the report assembler. Asset links must use an approved local route.
#' @param overlap_export Export displayed-set pairs (default), all pairs (explicit
#'   potentially large export), or none. Complete LE membership is always exported.
#' @return Paths to the portable HTML and its shared tables.
render_lisa_category_evidence <- function(evidence, output_dir, formats = "png",
    source_data = TRUE, recipes = TRUE, categories = NULL,
    overlap_export = c("displayed", "all", "none"), category_products = list()) {
  overlap_export <- match.arg(overlap_export)
  if (!inherits(evidence, "lisa_category_evidence")) stop("Expected lisa_category_evidence object.", call. = FALSE)
  formats <- unique(tolower(formats[nzchar(formats)]))
  if (any(!formats %in% c("png", "pdf", "svg"))) stop("Evidence formats must be png, pdf and/or svg.", call. = FALSE)
  if (!is.null(categories)) {
    unknown <- setdiff(categories, evidence$categories$category_id)
    if (length(unknown)) stop("Unknown evidence category IDs: ", paste(unknown, collapse = ", "), call. = FALSE)
    evidence$categories <- evidence$categories[evidence$categories$category_id %in% categories, , drop = FALSE]
    evidence$sets <- evidence$sets[evidence$sets$category_id %in% categories, , drop = FALSE]
    evidence$leading_edges <- evidence$leading_edges[evidence$leading_edges$category_id %in% categories, , drop = FALSE]
    evidence$gene_selection <- evidence$gene_selection[evidence$gene_selection$category_id %in% categories, , drop = FALSE]
  }
  lisa_guarded_dir_create(output_dir)
  table_dir <- file.path(output_dir, "tables"); lisa_guarded_dir_create(table_dir)
  tables <- list(categories = evidence$categories, sets = evidence$sets, genes = evidence$genes,
    leading_edges = evidence$leading_edges, gene_selection = evidence$gene_selection)
  if (evidence$metadata$complete_membership_supplied) tables$complete_memberships <- evidence$complete_memberships
  evidence$metadata$overlap_export <- overlap_export
  evidence$metadata$overlap_scope <- if (overlap_export == "displayed") "Pairs among default significant/FDR-selected displayed sets only, not all mapped sets; complete memberships are retained separately." else if (overlap_export == "all") "All within-category mapped-set pairs." else "No pair export; complete memberships and on-demand browser overlap remain available."
  evidence$metadata$n_possible_all_pairs <- sum(vapply(split(evidence$sets$pathway, evidence$sets$category_id), function(x) length(x) * (length(x) - 1) / 2, numeric(1L)))
  tables$metadata <- data.frame(key = names(evidence$metadata), value = vapply(evidence$metadata, lisa_evidence_text, character(1L)), stringsAsFactors = FALSE)
  table_paths <- character()
  for (name in names(tables)) {
    path <- file.path(table_dir, paste0(name, ".tsv"))
    if (name == "categories" && "support_schema_version" %in% names(tables[[name]]))
      lisa_write_inference_tsv(tables[[name]], path) else write_lisa_tsv(tables[[name]], path)
    table_paths[[name]] <- file.path("tables", paste0(name, ".tsv"))
  }
  # Displayed pairs are the bounded default. All-pair export is explicit and
  # streams indexed blocks, never an all-category quadratic data frame.
  if (overlap_export != "none") for (kind in c("leading_edge", if (evidence$metadata$complete_membership_supplied) "complete")) {
    name <- paste0(if (overlap_export == "displayed") "displayed_" else "all_", kind, "_overlaps")
    path <- file.path(table_dir, paste0(name, ".tsv"))
    lisa_guarded_write(path, function(target) {
      empty <- lisa_category_evidence_overlaps(evidence, "", kind)
      utils::write.table(empty, target, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
      for (id in evidence$categories$category_id) {
        lisa_category_evidence_overlaps(evidence, id, kind, displayed = overlap_export == "displayed",
          emit = function(chunk) utils::write.table(chunk, target, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE, append = TRUE, na = "NA"))
      }
    })
    table_paths[[name]] <- file.path("tables", paste0(name, ".tsv"))
  }
  figure_index <- data.frame(category_id = evidence$categories$category_id,
    anchor = paste0("category-", vapply(evidence$categories$category_id, utils::URLencode, character(1L), reserved = TRUE)),
    stem = sprintf("%03d_category_evidence", seq_len(nrow(evidence$categories))), stringsAsFactors = FALSE)
  if (length(formats)) {
    fig_dir <- file.path(output_dir, "figures"); lisa_guarded_dir_create(fig_dir)
    for (i in seq_len(nrow(figure_index))) {
      source <- lisa_category_evidence_figure_source(evidence, figure_index$category_id[[i]])
      stem <- file.path(fig_dir, figure_index$stem[[i]])
      lisa_category_evidence_save_figure(source, stem, formats)
      if (isTRUE(source_data) || isTRUE(recipes)) write_lisa_tsv(source, paste0(stem, "_source.tsv"))
    }
    if (isTRUE(recipes)) {
      recipe <- c("#!/usr/bin/env Rscript", "# Offline base-R reproduction. Usage:",
        "# Rscript --vanilla reproduce_category_evidence.R figures/001_category_evidence_source.tsv output.pdf",
        "# Source values and renderer are captured together; no server or package install is required.",
        paste0("lisa_category_evidence_draw <- ", paste(deparse(lisa_category_evidence_draw), collapse = "\n")),
        "args <- commandArgs(TRUE)", "if (length(args) != 2L) stop('Expected source TSV and output PNG/PDF/SVG paths.')",
        "source <- utils::read.delim(args[[1]], sep='\\t', quote='', comment.char='', check.names=FALSE, colClasses='character', na.strings=NULL)",
        "width <- max(14, min(25, 11 + sum(source$row_type == 'gene') * .3))", "height <- max(8, 5 + sum(source$row_type == 'set') * .25)",
        "format <- tolower(tools::file_ext(args[[2]]))", "if (!format %in% c('png','pdf','svg')) stop('Use a PNG, PDF or SVG output extension.')",
        "if (format == 'png') grDevices::png(args[[2]], width=width, height=height, units='in', res=160, bg='white')",
        "if (format == 'pdf') grDevices::pdf(args[[2]], width=width, height=height, useDingbats=FALSE)",
        "if (format == 'svg') grDevices::svg(args[[2]], width=width, height=height, bg='white')",
        "tryCatch(lisa_category_evidence_draw(source), finally=grDevices::dev.off())")
      lisa_guarded_write(file.path(output_dir, "reproduce_category_evidence.R"), function(target) writeLines(recipe, target, useBytes = TRUE))
    }
  }
  write_lisa_tsv(figure_index, file.path(table_dir, "figure_index.tsv"))
  # Standalone category exports must not advertise a nonexistent report page.
  # Report integration uses this stable offline layout; custom render callers
  # may provide an explicit relative destination in evidence metadata.
  if (is.null(evidence$metadata$gene_explorer_href)) {
    normalized_output <- gsub("\\\\", "/", output_dir)
    if (grepl("/report_pages/evidence/[^/]+/[^/]+/?$", normalized_output))
      evidence$metadata$gene_explorer_href <- "../../../gene_evidence/index.html"
  }
  # Standard evidence includes the existing member-set enrichment view at the
  # end of each sheet. Its population is independent of matrix display caps.
  # character() disables all static figure artifacts, including member-chart
  # sources/recipes. The interactive evidence and complete tables remain usable.
  member_figure_index <- if (length(formats)) {
    render_lisa_category_member_evidence(evidence, output_dir, formats)
  } else {
    data.frame(category_id = character(), stem = character(), source_tsv = character(),
      recipe_r = character(), n_significant_sets = integer(), source_sha256 = character(),
      recipe_sha256 = character(), stringsAsFactors = FALSE)
  }
  payload <- unclass(evidence)
  payload$member_figure_index <- member_figure_index
  if (!isTRUE(source_data)) payload$member_figure_index$source_tsv <- rep("", nrow(member_figure_index))
  if (!isTRUE(recipes)) payload$member_figure_index$recipe_r <- rep("", nrow(member_figure_index))
  payload$figure_index <- figure_index; payload$downloads <- as.list(table_paths)
  payload$formats <- formats; payload$source_data <- isTRUE(source_data)
  payload$recipes <- isTRUE(recipes) && length(formats) > 0L
  payload$category_products <- lisa_category_evidence_category_products(
    category_products, evidence$categories$category_id
  )
  payload$complete_memberships <- if (evidence$metadata$complete_membership_supplied) evidence$complete_memberships else NULL
  json <- jsonlite::toJSON(payload, auto_unbox = TRUE, dataframe = "rows", na = "null", null = "null", digits = 17)
  # '<' escaping protects the script-data block from identifiers like </script>.
  json <- gsub("<", "\\u003c", json, fixed = TRUE)
  asset_root <- system.file("category-evidence", package = "lisaR")
  if (!nzchar(asset_root)) stop("Installed category-evidence assets are missing.", call. = FALSE)
  html <- readLines(file.path(asset_root, "viewer.html"), warn = FALSE)
  # Keep immutable executable/style assets separate from the data-bearing HTML.
  # The report privacy gate correctly scans the latter for filesystem paths;
  # JS regex/escape syntax and CSS division are not path-bearing provenance.
  # Relative local assets work offline and their exact bytes remain bound by
  # the report's shareable SHA-256 manifest and installed-code identity ledger.
  asset_dir <- file.path(output_dir, "assets")
  lisa_guarded_dir_create(asset_dir)
  for (extension in c("css", "js")) {
    from <- file.path(asset_root, paste0("viewer.", extension))
    destination <- file.path(asset_dir, paste0("category-evidence.", extension))
    lisa_guarded_write(destination, function(target) {
      if (!file.copy(from, target, overwrite = TRUE)) stop("Cannot copy category evidence static asset.", call. = FALSE)
    })
  }
  html <- sub("<!-- EVIDENCE_STYLE -->", '<link rel="stylesheet" href="assets/category-evidence.css">', html, fixed = TRUE)
  html <- sub("<!-- EVIDENCE_DATA -->", paste0('<script type="application/json" id="evidence-data">', json, '</script>'), html, fixed = TRUE)
  html <- sub("<!-- EVIDENCE_SCRIPT -->", '<script src="assets/category-evidence.js"></script>', html, fixed = TRUE)
  html <- lisa_present_evidence_html(html, output_dir, evidence$metadata)
  path <- file.path(output_dir, "index.html")
  lisa_guarded_write(path, function(target) writeLines(html, target, useBytes = TRUE))
  invisible(list(html = path, tables = table_paths, figure_index = figure_index,
    member_figure_index = member_figure_index, status = "completed"))
}
