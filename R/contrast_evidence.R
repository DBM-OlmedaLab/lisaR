# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Descriptive comparison of existing category evidence. A and B retain their
# own enrichment, DE and leading-edge states; no differential test is created.
lisa_contrast_evidence_read <- function(x) {
  if (inherits(x, "lisa_category_evidence")) return(x)
  if (!is.character(x) || length(x) != 1L || !dir.exists(x))
    stop("Supply a category evidence object or its rendered directory.", call. = FALSE)
  names <- c("categories", "sets", "genes", "leading_edges", "metadata")
  paths <- file.path(x, "tables", paste0(names, ".tsv"))
  if (any(!file.exists(paths))) stop("Incomplete category evidence directory.", call. = FALSE)
  out <- stats::setNames(lapply(paths, lisa_evidence_read), names)
  out$metadata <- as.list(stats::setNames(out$metadata$value, out$metadata$key))
  out$input_provenance <- data.frame(table = basename(paths),
    size_bytes = as.numeric(file.info(paths)$size),
    sha256 = vapply(paths, function(path) digest::digest(file = path, algo = "sha256"), character(1L)),
    stringsAsFactors = FALSE, row.names = NULL)
  structure(out, class = "lisa_category_evidence")
}

lisa_contrast_evidence_bool <- function(x) {
  x <- tolower(as.character(x))
  if (any(!is.na(x) & !x %in% c("true", "false", "na", "")))
    stop("Malformed logical evidence state.", call. = FALSE)
  x %in% "true"
}

lisa_contrast_evidence_side <- function(x) {
  x <- lisa_contrast_evidence_read(x)
  needed <- c("category_id", "pathway", "NES", "gsea_fdr", "evaluable", "significant", "leading_edge_state")
  if (!all(needed %in% names(x$sets))) stop("Incomplete normalized set evidence.", call. = FALSE)
  for (name in c("category_id", "pathway")) x$sets[[name]] <- as.character(x$sets[[name]])
  if (anyDuplicated(x$sets[c("category_id", "pathway")]))
    stop("Duplicate exact category/pathway assignments in evidence.", call. = FALSE)
  if (anyNA(x$sets$category_id) || any(!nzchar(x$sets$category_id)) || any(x$sets$category_id == "OTHER_UNCLASSIFIED"))
    stop("Contrast evidence accepts classified categories only.", call. = FALSE)
  if (anyDuplicated(x$genes$symbol)) stop("Duplicate exact gene symbols in evidence.", call. = FALSE)
  for (name in c("NES", "gsea_fdr")) x$sets[[name]] <- lisa_evidence_num(x$sets[[name]])
  for (name in c("evaluable", "significant")) x$sets[[name]] <- lisa_contrast_evidence_bool(x$sets[[name]])
  cutoff <- lisa_evidence_num(x$metadata$gsea_padj_cutoff)
  if (length(cutoff) != 1L || !is.finite(cutoff) || cutoff < 0 || cutoff > 1)
    stop("Missing or invalid evidence GSEA cutoff.", call. = FALSE)
  evaluable <- is.finite(x$sets$NES) & is.finite(x$sets$gsea_fdr) & x$sets$gsea_fdr >= 0 & x$sets$gsea_fdr <= 1
  if (any(x$sets$evaluable & !evaluable) || any(x$sets$significant != (x$sets$evaluable & x$sets$gsea_fdr <= cutoff)))
    stop("Evidence evaluability/significance conflicts with the recorded cutoff.", call. = FALSE)
  for (name in c("log2FC", "de_fdr")) x$genes[[name]] <- lisa_evidence_num(x$genes[[name]])
  x$genes$symbol <- as.character(x$genes$symbol)
  x$leading_edges <- unique(x$leading_edges[c("pathway", "symbol")])
  x$leading_edges[] <- lapply(x$leading_edges, as.character)
  x$leading_edges <- x$leading_edges[lisa_evidence_order(x$leading_edges$pathway, x$leading_edges$symbol), , drop = FALSE]
  x
}

# Shared NES statistics, deliberately using the same helper as the three
# category plot variants. Non-evaluable rows must never enter its population.
lisa_contrast_evidence_statistics <- function(x, category_ids) {
  g <- x$sets
  g$padj <- ifelse(g$evaluable, g$gsea_fdr, NA_real_)
  lisa_category_nes_statistics(g, gsea_padj_cutoff = lisa_evidence_num(x$metadata$gsea_padj_cutoff), category_ids = category_ids)
}

build_lisa_contrast_evidence <- function(evidence_a, evidence_b, contrast_id,
    label = contrast_id, max_sets = 25L, max_genes = 40L, contrast_summary = NULL) {
  a <- lisa_contrast_evidence_side(evidence_a); b <- lisa_contrast_evidence_side(evidence_b)
  for (key in c("analysis_id", "collection", "tier")) for (x in list(a, b)) {
    value <- x$metadata[[key]]
    if (length(value) != 1L || is.na(value) || !nzchar(value)) stop("Evidence scope is incomplete: ", key, call. = FALSE)
  }
  for (key in c("collection", "tier")) if (!identical(as.character(a$metadata[[key]]), as.character(b$metadata[[key]])))
    stop("Contrast evidence requires matching ", key, "; core/expanded comparisons are not supported.", call. = FALSE)
  if (identical(a$metadata$analysis_id, b$metadata$analysis_id)) stop("A and B must identify different analyses.", call. = FALSE)
  for (key in c("gsea_padj_cutoff", "de_padj_cutoff")) if (!isTRUE(all.equal(lisa_evidence_num(a$metadata[[key]]), lisa_evidence_num(b$metadata[[key]]))))
    stop("Contrast evidence requires matching ", key, ".", call. = FALSE)
  for (x in list(contrast_id, label)) if (length(x) != 1L || is.na(x) || !nzchar(x)) stop("Contrast ID and label must be nonempty strings.", call. = FALSE)
  for (x in list(max_sets, max_genes)) if (length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x)) stop("Display caps must be positive integers.", call. = FALSE)
  category_ids <- sort(unique(c(as.character(a$categories$category_id), as.character(b$categories$category_id))), method = "radix")
  category_ids <- setdiff(category_ids[!is.na(category_ids) & nzchar(category_ids)], "OTHER_UNCLASSIFIED")
  categories <- data.frame(category_id = category_ids, stringsAsFactors = FALSE)
  for (key in c("category_display_name", "macrogroup_id", "macrogroup_name", "category_description")) {
    av <- as.character(lisa_evidence_column(a$categories, key))[match(category_ids, a$categories$category_id)]
    bv <- as.character(lisa_evidence_column(b$categories, key))[match(category_ids, b$categories$category_id)]
    categories[[key]] <- ifelse(!is.na(av) & nzchar(av), av, ifelse(!is.na(bv) & nzchar(bv), bv, if (key == "category_display_name") category_ids else ""))
  }
  paired <- merge(a$sets, b$sets, by = c("category_id", "pathway"), all = TRUE, suffixes = c("_a", "_b"), sort = FALSE)
  paired <- paired[lisa_evidence_order(paired$category_id, paired$pathway), , drop = FALSE]
  for (side in c("a", "b")) {
    paired[[paste0("present_", side)]] <- !is.na(paired[[paste0("evaluable_", side)]])
    # Missing significance flags can be false for set selection; NES/FDR and
    # missing/present states remain NA/explicit, never imputed to zero.
    for (key in c("evaluable", "significant")) paired[[paste0(key, "_", side)]] <- paired[[paste0(key, "_", side)]] %in% TRUE
  }
  paired$significance_membership <- ifelse(paired$significant_a & paired$significant_b, "both", ifelse(paired$significant_a, "A_only", ifelse(paired$significant_b, "B_only", "neither")))
  statistics <- list(a = lisa_contrast_evidence_statistics(a, category_ids), b = lisa_contrast_evidence_statistics(b, category_ids))
  # Store shared helper columns with side suffixes, without combining them.
  for (side in c("a", "b")) {
    st <- statistics[[side]]; ii <- match(category_ids, st$category_id)
    x <- if (side == "a") a else b
    if ("mean_significant_set_NES" %in% names(x$categories)) {
      observed <- lisa_evidence_num(x$categories$mean_significant_set_NES)
      si <- match(x$categories$category_id, st$category_id)
      derived <- st$mean_NES[si]
      mismatch <- xor(is.finite(observed), is.finite(derived)) |
        (is.finite(observed) & is.finite(derived) & abs(observed - derived) > 1e-10 * pmax(1, abs(derived)))
      if (any(mismatch)) stop("Recorded category mean conflicts with full significant support in side ", side, ".", call. = FALSE)
      st$mean_NES[si] <- observed
    }
    for (key in setdiff(names(st), "category_id")) categories[[paste0(key, "_", side)]] <- st[[key]][ii]
    mapped <- table(factor(x$sets$category_id, levels = category_ids))
    evaluable <- table(factor(x$sets$category_id[x$sets$evaluable], levels = category_ids))
    categories[[paste0("n_mapped_", side)]] <- as.integer(mapped)
    categories[[paste0("n_evaluable_", side)]] <- as.integer(evaluable)
  }
  for (state in c("both", "A_only", "B_only")) categories[[paste0("n_significant_", state)]] <- vapply(category_ids, function(id) sum(paired$category_id == id & paired$significance_membership == state), integer(1L))
  categories$display_mean_NES_a <- categories$mean_NES_a
  categories$display_mean_NES_b <- categories$mean_NES_b
  categories$endpoint_source_a <- ifelse(is.finite(categories$mean_NES_a), "significant_mean", "not_available")
  categories$endpoint_source_b <- ifelse(is.finite(categories$mean_NES_b), "significant_mean", "not_available")
  categories$delta_mean_NES <- categories$mean_NES_a - categories$mean_NES_b
  categories$delta_state <- ifelse(is.finite(categories$delta_mean_NES), "derived_significant_means", "not_available")
  canonical <- lisa_evidence_read(contrast_summary)
  if (!is.null(contrast_summary)) {
    required <- c("category_id", "mean_NES_A", "mean_NES_B", "display_mean_NES_A", "display_mean_NES_B", "delta_mean_NES")
    if (!all(required %in% names(canonical)) || anyDuplicated(canonical$category_id)) stop("Canonical contrast summary is incomplete or duplicated.", call. = FALSE)
    if (any(!canonical$category_id %in% category_ids)) stop("Canonical contrast contains categories absent from side evidence.", call. = FALSE)
    if ("universe" %in% names(canonical) && any(canonical$universe != a$metadata$collection)) stop("Canonical contrast universe conflicts with side evidence.", call. = FALSE)
    if ("gsea_padj_cutoff" %in% names(canonical) && any(lisa_evidence_num(canonical$gsea_padj_cutoff) != lisa_evidence_num(a$metadata$gsea_padj_cutoff))) stop("Canonical contrast cutoff conflicts with side evidence.", call. = FALSE)
    ci <- match(canonical$category_id, category_ids)
    # Category evidence may contain more categories than the existing contrast
    # product. Do not invent a delta for those absent from its canonical table.
    categories$delta_mean_NES[] <- NA_real_
    categories$delta_state[] <- "not_in_canonical_contrast"
    close <- function(u, v) !xor(is.finite(u), is.finite(v)) & (!is.finite(u) | abs(u-v) <= 1e-10*pmax(1,abs(u),abs(v)))
    for (side in c("a", "b")) {
      upper <- toupper(side); x <- if (side == "a") a else b
      supplied_mean <- lisa_evidence_num(canonical[[paste0("mean_NES_", upper)]])
      if (any(!close(supplied_mean, categories[[paste0("mean_NES_", side)]][ci]))) stop("Canonical contrast significant mean conflicts with side ", side, ".", call. = FALSE)
      value <- lisa_evidence_num(canonical[[paste0("display_mean_NES_", upper)]])
      supported <- categories[[paste0("n_genesets_significant_", side)]][ci] > 0L
      context <- vapply(canonical$category_id, function(id) {
        v <- x$sets$NES[x$sets$category_id == id & x$sets$evaluable]
        if (length(v)) mean(v) else NA_real_
      }, numeric(1L))
      if (any(supported & !close(value, supplied_mean)) || any(!supported & is.finite(value) & !close(value, context))) stop("Canonical display mean conflicts with significant/contextual support in side ", side, ".", call. = FALSE)
      categories[[paste0("display_mean_NES_", side)]][ci] <- value
      source <- ifelse(supported, "significant_mean", ifelse(is.finite(value), "contextual_mean_not_significant", "not_available"))
      categories[[paste0("endpoint_source_", side)]][ci] <- source
    }
    delta <- lisa_evidence_num(canonical$delta_mean_NES)
    calculated <- categories$display_mean_NES_a[ci] - categories$display_mean_NES_b[ci]
    if (any(!close(delta, calculated))) stop("Canonical contrast delta is not oriented A minus B.", call. = FALSE)
    categories$delta_mean_NES[ci] <- delta
    categories$delta_state[ci] <- ifelse(is.finite(delta), "canonical", "not_available")
  }
  # LE comparison is one row per unique set, not repeated per category.
  set_ids <- sort(unique(paired$pathway), method = "radix")
  ia <- match(set_ids, a$sets$pathway); ib <- match(set_ids, b$sets$pathway)
  ga <- split(a$leading_edges$symbol, a$leading_edges$pathway); gb <- split(b$leading_edges$symbol, b$leading_edges$pathway)
  conservation <- data.frame(pathway = set_ids, n_a = rep(NA_integer_, length(set_ids)), n_b = rep(NA_integer_, length(set_ids)),
    intersection_n = rep(NA_integer_, length(set_ids)), union_n = rep(NA_integer_, length(set_ids)), jaccard = rep(NA_real_, length(set_ids)),
    overlap_state = rep("not_available", length(set_ids)), stringsAsFactors = FALSE)
  for (i in seq_along(set_ids)) {
    aa <- unique(ga[[set_ids[[i]]]]); bb <- unique(gb[[set_ids[[i]]]])
    available_a <- !is.na(ia[[i]]) && a$sets$evaluable[[ia[[i]]]] && a$sets$leading_edge_state[[ia[[i]]]] %in% c("recorded", "empty")
    available_b <- !is.na(ib[[i]]) && b$sets$evaluable[[ib[[i]]]] && b$sets$leading_edge_state[[ib[[i]]]] %in% c("recorded", "empty")
    if (available_a) conservation$n_a[[i]] <- length(aa)
    if (available_b) conservation$n_b[[i]] <- length(bb)
    available <- available_a && available_b
    if (available) {
      intersection <- length(intersect(aa, bb)); union <- length(union(aa, bb))
      conservation$intersection_n[[i]] <- intersection; conservation$union_n[[i]] <- union
      if (union > 0L) { conservation$jaccard[[i]] <- intersection / union; conservation$overlap_state[[i]] <- "defined" } else conservation$overlap_state[[i]] <- "empty_union"
    }
  }
  # Remove category-repeated edge rows before serialization. DE is stored once
  # per side and linked by exact symbol; no per-category expression copies.
  genes <- merge(a$genes, b$genes, by = "symbol", all = TRUE, suffixes = c("_a", "_b"), sort = FALSE)
  genes <- genes[lisa_evidence_order(genes$symbol), , drop = FALSE]
  ea <- a$leading_edges; eb <- b$leading_edges
  ea$side <- rep("a", nrow(ea)); eb$side <- rep("b", nrow(eb))
  edges <- rbind(ea, eb)
  provenance <- data.frame(side = character(), table = character(), size_bytes = numeric(), sha256 = character(), stringsAsFactors = FALSE)
  for (side in c("a", "b")) {
    x <- if (side == "a") a else b
    if (!is.null(x$input_provenance)) provenance <- rbind(provenance, cbind(side = side, x$input_provenance, stringsAsFactors = FALSE))
  }
  if (is.character(contrast_summary) && length(contrast_summary) == 1L && file.exists(contrast_summary)) {
    provenance <- rbind(provenance, data.frame(side = "contrast", table = basename(contrast_summary),
      size_bytes = as.numeric(file.info(contrast_summary)$size), sha256 = digest::digest(file = contrast_summary, algo = "sha256"), stringsAsFactors = FALSE))
  }
  metadata <- list(schema_version = "1.0", contrast_id = as.character(contrast_id), contrast_label = as.character(label),
    analysis_a = as.character(a$metadata$analysis_id), analysis_b = as.character(b$metadata$analysis_id),
    analysis_label_a = if (nrow(canonical) && "contrast_a_label" %in% names(canonical)) as.character(canonical$contrast_a_label[[1L]]) else as.character(a$metadata$analysis_id),
    analysis_label_b = if (nrow(canonical) && "contrast_b_label" %in% names(canonical)) as.character(canonical$contrast_b_label[[1L]]) else as.character(b$metadata$analysis_id),
    collection = as.character(a$metadata$collection), tier = as.character(a$metadata$tier),
    positive_contrast_a = a$metadata$positive_contrast, positive_contrast_b = b$metadata$positive_contrast,
    gsea_padj_cutoff = lisa_evidence_num(a$metadata$gsea_padj_cutoff), de_padj_cutoff = lisa_evidence_num(a$metadata$de_padj_cutoff),
    max_sets = as.integer(max_sets), max_genes = as.integer(max_genes),
    interpretation = "A minus B is descriptive. Significant only in A does not establish a difference between A and B. Missing results are not zero.",
    support_population = "Unique significant and evaluable gene sets independently within each side; support composition may differ.",
    delta_source = if (is.null(contrast_summary)) "Significant-side means only; not a replacement for an existing contrast summary." else "Exact canonical contrast display endpoints and A-B delta; contextual non-significant endpoints retain explicit labels.",
    same_set_overlap = "Leading-edge Jaccard only when both sides have evaluable GSEA and recorded leading edges; an empty union is undefined.")
  rownames(categories) <- rownames(paired) <- rownames(genes) <- rownames(edges) <- NULL
  structure(list(metadata = metadata, categories = categories, sets = paired, genes = genes,
    leading_edges = edges, conservation = conservation, provenance = provenance), class = "lisa_contrast_evidence")
}

lisa_contrast_evidence_view <- function(evidence, category_id, selection = c("union", "intersection", "evaluable")) {
  selection <- match.arg(selection)
  if (!category_id %in% evidence$categories$category_id) stop("Unknown exact contrast category.", call. = FALSE)
  sets <- evidence$sets[evidence$sets$category_id == category_id, , drop = FALSE]
  keep <- switch(selection, union = sets$significant_a | sets$significant_b,
    intersection = sets$significant_a & sets$significant_b, evaluable = sets$evaluable_a | sets$evaluable_b)
  sets <- sets[keep, , drop = FALSE]
  fdr <- pmin(ifelse(sets$significant_a, sets$gsea_fdr_a, Inf), ifelse(sets$significant_b, sets$gsea_fdr_b, Inf))
  if (selection == "evaluable") fdr <- pmin(ifelse(sets$evaluable_a, sets$gsea_fdr_a, Inf), ifelse(sets$evaluable_b, sets$gsea_fdr_b, Inf))
  sets <- sets[lisa_evidence_order(fdr, sets$pathway), , drop = FALSE]
  n_matching <- nrow(sets)
  full <- evidence$sets[evidence$sets$category_id == category_id, , drop = FALSE]
  supported <- evidence$leading_edges[(evidence$leading_edges$side == "a" & evidence$leading_edges$pathway %in% full$pathway[full$significant_a]) |
    (evidence$leading_edges$side == "b" & evidence$leading_edges$pathway %in% full$pathway[full$significant_b]), , drop = FALSE]
  recurrence <- table(supported$symbol)
  sets <- head(sets, evidence$metadata$max_sets)
  edges <- evidence$leading_edges[evidence$leading_edges$pathway %in% sets$pathway, , drop = FALSE]
  symbols <- sort(unique(edges$symbol), method = "radix")
  counts <- as.integer(recurrence[symbols]); counts[is.na(counts)] <- 0L
  symbols <- symbols[lisa_evidence_order(-counts, symbols)]
  n_genes <- length(symbols); symbols <- head(symbols, evidence$metadata$max_genes)
  list(sets = sets, genes = evidence$genes[match(symbols, evidence$genes$symbol), , drop = FALSE],
    edges = edges[edges$symbol %in% symbols, , drop = FALSE], n_matching_sets = n_matching,
    n_available_genes = n_genes, selection = selection)
}

lisa_contrast_evidence_figure_source <- function(evidence, category_id, selection = "union") {
  view <- lisa_contrast_evidence_view(evidence, category_id, selection)
  columns <- c("row_type", "contrast_id", "contrast_label", "analysis_a", "analysis_b", "collection", "tier", "category_id", "category_display_name",
    "selection", "n_matching_sets", "n_available_genes", "max_sets", "max_genes", "side", "pathway", "symbol", "set_order", "gene_order", "NES", "gsea_fdr", "evaluable", "significant", "leading_edge_state", "log2FC", "de_fdr", "de_state")
  base <- as.data.frame(stats::setNames(rep(list(""), length(columns)), columns), stringsAsFactors = FALSE)
  for (key in intersect(columns, names(evidence$metadata))) base[[key]] <- lisa_evidence_text(evidence$metadata[[key]])
  base$category_id <- category_id
  base$category_display_name <- evidence$categories$category_display_name[match(category_id, evidence$categories$category_id)]
  base$selection <- selection; base$n_matching_sets <- view$n_matching_sets; base$n_available_genes <- view$n_available_genes
  rows <- list(transform(base, row_type = "metadata"))
  for (side in c("a", "b")) {
    for (i in seq_len(nrow(view$sets))) {
      row <- base; row$row_type <- "set"; row$side <- side; row$pathway <- view$sets$pathway[[i]]; row$set_order <- i
      for (key in c("NES", "gsea_fdr", "evaluable", "significant", "leading_edge_state")) row[[key]] <- lisa_evidence_text(view$sets[[paste0(key, "_", side)]][[i]])
      rows[[length(rows) + 1L]] <- row
    }
    for (i in seq_len(nrow(view$genes))) {
      row <- base; row$row_type <- "gene"; row$side <- side; row$symbol <- view$genes$symbol[[i]]; row$gene_order <- i
      for (key in c("log2FC", "de_fdr", "de_state")) row[[key]] <- lisa_evidence_text(view$genes[[paste0(key, "_", side)]][[i]])
      rows[[length(rows) + 1L]] <- row
    }
  }
  if (nrow(view$edges)) {
    row <- base[rep(1L, nrow(view$edges)), , drop = FALSE]; row$row_type <- "membership"
    row$side <- view$edges$side; row$pathway <- view$edges$pathway; row$symbol <- view$edges$symbol
    row$set_order <- match(row$pathway, view$sets$pathway); row$gene_order <- match(row$symbol, view$genes$symbol)
    rows[[length(rows) + 1L]] <- row
  }
  out <- do.call(rbind, lapply(rows, function(x) { x[] <- lapply(x, lisa_evidence_text); x }))
  rownames(out) <- NULL; out
}

# Self-contained base-R renderer, copied verbatim into each offline recipe.
lisa_contrast_evidence_draw <- function(source) {
  num <- function(x) suppressWarnings(as.numeric(x))
  fmt <- function(x) ifelse(is.finite(num(x)), formatC(num(x), format = "g", digits = 2), "NA")
  color <- function(x) ifelse(!is.finite(num(x)), "#9ca3af", ifelse(num(x) > 0, "#b2182b", ifelse(num(x) < 0, "#2166ac", "#64748b")))
  label <- function(x, n = 40L) ifelse(nchar(x) > n, paste0(substr(x, 1L, n - 3L), "..."), x)
  meta <- source[source$row_type == "metadata", , drop = FALSE][1L, , drop = FALSE]
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(1, 2), oma = c(4, 0, 6, 0))
  for (side in c("a", "b")) {
    ss <- source[source$row_type == "set" & source$side == side, , drop = FALSE]; ss <- ss[order(num(ss$set_order)), , drop = FALSE]
    gg <- source[source$row_type == "gene" & source$side == side, , drop = FALSE]; gg <- gg[order(num(gg$gene_order)), , drop = FALSE]
    ee <- source[source$row_type == "membership" & source$side == side, , drop = FALSE]
    ns <- nrow(ss); ng <- nrow(gg); yy <- ns + 1L - seq_len(ns)
    graphics::par(mar = c(8, if (side == "a") 20 else 4, 5, 3), xpd = NA)
    graphics::plot.new(); graphics::plot.window(xlim = c(.5, max(ng, 1) + 10), ylim = c(.5, max(ns, 1) + 5))
    graphics::title(main = paste(toupper(side), meta[[paste0("analysis_", side)]]), cex.main = .85, line = 3)
    if (!ns) { graphics::text(1, 1, "No matching gene sets", adj = 0); next }
    for (i in seq_len(ns)) {
      graphics::rect(.5, yy[[i]] - .45, max(ng, 1) + .5, yy[[i]] + .45, border = NA, col = if (i %% 2L) "#f1f5f9" else "white")
      if (side == "a") graphics::text(.35, yy[[i]], label(ss$pathway[[i]]), adj = 1, cex = .53)
      graphics::text(max(ng, 1) + 3, yy[[i]], fmt(ss$NES[[i]]), col = color(ss$NES[[i]]), cex = .63)
      graphics::text(max(ng, 1) + 7, yy[[i]], fmt(ss$gsea_fdr[[i]]), cex = .58)
    }
    graphics::text(max(ng, 1) + c(3, 7), ns + 1, c("NES", "FDR GSEA"), cex = .65)
    if (nrow(ee)) {
      ii <- match(ee$pathway, ss$pathway)
      # Match the interactive paired matrix: A is a circle and B a square.
      graphics::points(num(ee$gene_order), yy[ii], pch = if (side == "a") 16 else 15, cex = .8, col = color(ss$NES[ii]))
    }
    if (ng) {
      graphics::text(seq_len(ng), rep(.2, ng), gg$symbol, srt = 65, adj = 1, cex = .53)
      graphics::points(seq_len(ng), rep(ns + 2, ng), pch = 15, col = color(gg$log2FC), cex = 1)
      graphics::text(seq_len(ng), rep(ns + 2.2, ng), fmt(gg$log2FC), cex = .48, srt = 90, adj = c(0, .5))
    }
  }
  graphics::mtext(paste(meta$category_display_name, "|", meta$contrast_label), outer = TRUE, side = 3, line = 4, cex = 1)
  graphics::mtext(paste(meta$collection, meta$tier, "|", meta$selection, "| displayed rows / matching:", nrow(source[source$row_type == "set" & source$side == "a", ]), "/", meta$n_matching_sets), outer = TRUE, side = 3, line = 2.5, cex = .8)
  graphics::mtext("Top band: gene log2FC. Cells: leading edge. A circle; B square. Red positive, blue negative. Separate FDRs; absent is not zero.", outer = TRUE, side = 1, line = 2, cex = .7)
}

## Validate only the small, report-assembler attachment contract.  Category
## membership remains an exact supplied ID; filenames are never inspected.
lisa_contrast_evidence_category_products <- function(category_products, category_ids) {
  if (is.null(category_products)) return(list())
  if (!is.list(category_products)) stop("category_products must be a list.", call. = FALSE)
  out <- lapply(category_products, function(product) {
    if (!is.list(product) || !all(c("category_id", "product", "label", "assets") %in% names(product)))
      stop("Each category product needs category_id, product, label and assets.", call. = FALSE)
    id <- as.character(product$category_id)
    if (length(id) != 1L || is.na(id) || !id %in% category_ids)
      stop("A category product has no exact contrast category ID.", call. = FALSE)
    assets <- product$assets
    if (!is.list(assets) || !length(assets)) stop("Each category product needs at least one asset.", call. = FALSE)
    assets <- lapply(assets, function(asset) {
      if (!is.list(asset) || !all(c("format", "href") %in% names(asset))) stop("Each category product asset needs format and href.", call. = FALSE)
      href <- as.character(asset$href); format <- as.character(asset$format)
      if (!lisa_category_product_href_safe(href) || length(format) != 1L || is.na(format) || !nzchar(format)) stop("A category product asset is not a safe relative file.", call. = FALSE)
      list(format = format, href = href)
    })
    product_name <- as.character(product$product); label <- as.character(product$label)
    if (length(product_name) != 1L || is.na(product_name) || !nzchar(product_name) || length(label) != 1L || is.na(label) || !nzchar(label))
      stop("Each category product needs nonempty product and label.", call. = FALSE)
    list(category_id = id, product = product_name, label = label, assets = assets)
  })
  unname(out)
}

render_lisa_contrast_evidence <- function(evidence, output_dir, formats = "png",
    variants = c("clean", "percentages", "direction", "dispersion"), report_href = "", category_products = list()) {
  if (!inherits(evidence, "lisa_contrast_evidence")) stop("Expected lisa_contrast_evidence object.", call. = FALSE)
  formats <- unique(tolower(formats[nzchar(formats)])); variants <- unique(variants)
  if (any(!formats %in% c("png", "pdf", "svg"))) stop("Unsupported contrast figure format.", call. = FALSE)
  if (!length(variants) || any(!variants %in% c("clean", "percentages", "direction", "dispersion"))) stop("Unknown contrast navigator variant.", call. = FALSE)
  lisa_guarded_dir_create(output_dir)
  for (dir in c("tables", "assets", "figures")) lisa_guarded_dir_create(file.path(output_dir, dir))
  tables <- evidence[c("categories", "sets", "genes", "leading_edges", "conservation", "provenance")]
  tables$metadata <- data.frame(key = names(evidence$metadata), value = vapply(evidence$metadata, lisa_evidence_text, character(1L)), stringsAsFactors = FALSE)
  for (key in names(tables)) write_lisa_tsv(tables[[key]], file.path(output_dir, "tables", paste0(key, ".tsv")))
  index <- data.frame(category_id = evidence$categories$category_id, stem = sprintf("%03d_contrast_evidence", seq_len(nrow(evidence$categories))), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(index))) {
    source <- lisa_contrast_evidence_figure_source(evidence, index$category_id[[i]])
    stem <- file.path(output_dir, "figures", index$stem[[i]])
    write_lisa_tsv(source, paste0(stem, "_source.tsv"))
    width <- max(18, min(32, 14 + sum(source$row_type == "gene" & source$side == "a") * .4))
    height <- max(9, 5 + sum(source$row_type == "set" & source$side == "a") * .28)
    for (format in formats) lisa_guarded_write(paste0(stem, ".", format), function(target) {
      if (format == "png") grDevices::png(target, width = width, height = height, units = "in", res = 140, bg = "white")
      if (format == "pdf") grDevices::pdf(target, width = width, height = height, useDingbats = FALSE)
      if (format == "svg") grDevices::svg(target, width = width, height = height, bg = "white")
      tryCatch(lisa_contrast_evidence_draw(source), finally = grDevices::dev.off())
    })
  }
  write_lisa_tsv(index, file.path(output_dir, "tables", "figure_index.tsv"))
  recipe <- c("#!/usr/bin/env Rscript", "# Reproduce a paired descriptive matrix from its exact stored view; no GSEA runs.",
    paste0("lisa_contrast_evidence_draw <- ", paste(deparse(lisa_contrast_evidence_draw), collapse = "\n")),
    "args <- commandArgs(TRUE)", "if (length(args) != 2L) stop('Expected source TSV and output PNG/PDF/SVG paths.')",
    "source <- utils::read.delim(args[[1]], sep='\\t', quote='', comment.char='', check.names=FALSE, colClasses='character', na.strings=NULL)",
    "width <- max(18, min(32, 14 + sum(source$row_type == 'gene' & source$side == 'a') * .4))",
    "height <- max(9, 5 + sum(source$row_type == 'set' & source$side == 'a') * .28)",
    "format <- tolower(tools::file_ext(args[[2]]))", "if (!format %in% c('png','pdf','svg')) stop('Use PNG, PDF or SVG.')",
    "if (format == 'png') grDevices::png(args[[2]], width=width, height=height, units='in',res=140,bg='white')",
    "if (format == 'pdf') grDevices::pdf(args[[2]], width=width,height=height,useDingbats=FALSE)",
    "if (format == 'svg') grDevices::svg(args[[2]], width=width,height=height,bg='white')",
    "tryCatch(lisa_contrast_evidence_draw(source),finally=grDevices::dev.off())")
  lisa_guarded_write(file.path(output_dir, "reproduce_contrast_evidence.R"), function(target) writeLines(recipe, target, useBytes = TRUE))
  payload <- unclass(evidence); payload$variants <- variants; payload$formats <- formats; payload$figure_index <- index
  payload$category_products <- lisa_contrast_evidence_category_products(category_products, evidence$categories$category_id)
  payload$metadata$report_href <- report_href
  json <- jsonlite::toJSON(payload, auto_unbox = TRUE, dataframe = "rows", rownames = FALSE,
    na = "null", null = "null", digits = 17)
  json <- gsub("<", "\\u003c", json, fixed = TRUE)
  assets <- system.file("contrast-evidence", package = "lisaR")
  if (!nzchar(assets)) stop("Installed contrast-evidence assets are missing.", call. = FALSE)
  html <- readLines(file.path(assets, "viewer.html"), warn = FALSE)
  for (ext in c("css", "js")) lisa_guarded_write(file.path(output_dir, "assets", paste0("contrast-evidence.", ext)), function(target) {
    if (!file.copy(file.path(assets, paste0("viewer.", ext)), target, overwrite = TRUE)) stop("Cannot copy contrast evidence asset.", call. = FALSE)
  })
  html <- sub("<!-- CONTRAST_DATA -->", paste0('<script type="application/json" id="contrast-evidence-data">', json, '</script>'), html, fixed = TRUE)
  shell_metadata <- evidence$metadata
  shell_metadata$positive_contrast <- "A minus B (descriptive)"
  html <- lisa_present_evidence_html(html, output_dir, shell_metadata,
    active = "contrasts", main_id = "contrast-evidence-main")
  path <- file.path(output_dir, "index.html")
  lisa_guarded_write(path, function(target) writeLines(html, target, useBytes = TRUE))
  invisible(list(html = path, figure_index = index, status = "completed"))
}
