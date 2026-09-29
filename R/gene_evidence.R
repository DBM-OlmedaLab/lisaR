# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Global, descriptive gene evidence. Full TERM2GENE membership is deliberately
# separate from GSEA leading edges, and DE measurements are stored per analysis.
lisa_gene_evidence_read <- function(x) {
  if (is.null(x)) return(data.frame(stringsAsFactors = FALSE))
  if (is.character(x) && length(x) == 1L) return(utils::read.delim(x,
    sep = "\t", quote = "", comment.char = "", check.names = FALSE,
    colClasses = "character", na.strings = "NA", stringsAsFactors = FALSE))
  as.data.frame(x, stringsAsFactors = FALSE)
}

lisa_gene_evidence_col <- function(x, choices, default = NA_character_) {
  column <- choices[choices %in% names(x)]
  if (!length(column)) return(rep(default, nrow(x)))
  as.character(x[[column[[1L]]]])
}

lisa_gene_evidence_num <- function(x) suppressWarnings(as.numeric(x))

lisa_gene_evidence_tokens <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) return(character())
  value <- trimws(unlist(strsplit(x, "[;,|/[:space:]]+")))
  sort(unique(value[nzchar(value)]), method = "radix")
}

lisa_gene_evidence_sort <- function(x, keys) {
  if (!nrow(x)) return(x)
  x <- x[do.call(order, c(x[keys], list(na.last = TRUE, method = "radix"))), , drop = FALSE]
  rownames(x) <- NULL
  x
}

lisa_gene_evidence_de <- function(x, analysis_id, column_mapping = NULL) {
  d <- lisa_gene_evidence_read(x)
  if (nrow(d) && !"symbol" %in% names(d)) stop("Gene evidence DE table lacks symbol.", call. = FALSE)
  out <- data.frame(analysis_id = rep(analysis_id, nrow(d)),
    symbol = lisa_gene_evidence_col(d, "symbol"),
    gene_id = lisa_gene_evidence_col(d, c("gene_id", "ensembl_id"), ""),
    log2FC = lisa_gene_evidence_num(lisa_gene_evidence_col(d, c("log2FoldChange", "log2FC", "logFC"))),
    rank_value = lisa_gene_evidence_num(lisa_gene_evidence_col(d, c("rank_value", "rank_stat"))),
    statistic = lisa_gene_evidence_num(lisa_gene_evidence_col(d, c("statistic", "stat", "wald_stat", "t"))),
    pvalue = lisa_gene_evidence_num(lisa_gene_evidence_col(d, c("pvalue", "PValue", "P.Value"))),
    de_fdr = lisa_gene_evidence_num(lisa_gene_evidence_col(d, c("padj", "FDR", "de_fdr"))),
    stringsAsFactors = FALSE)
  mapping <- lisa_gene_evidence_read(column_mapping)
  source_column <- function(name) {
    if (!all(c("standard_column", "source_column") %in% names(mapping))) return("not recorded")
    value <- unique(as.character(mapping$source_column[mapping$standard_column == name]))
    value <- value[!is.na(value) & nzchar(value)]
    if (length(value) > 1L) stop("Conflicting standardised DE column provenance.", call. = FALSE)
    if (length(value)) value[[1L]] else "not recorded"
  }
  out$gene_id_source <- rep(source_column("gene_id"), nrow(out))
  out$rank_source <- rep(source_column("rank_value"), nrow(out))
  # A statistic can be recovered only from explicit column provenance, never
  # by assuming that every ranking (e.g. logFC/signed-logP) is a test statistic.
  rank_is_statistic <- tolower(out$rank_source) %in% c("stat", "statistic", "wald_stat", "wald", "t", "z")
  recover <- is.na(out$statistic) & rank_is_statistic
  out$statistic[recover] <- out$rank_value[recover]
  if (any(is.na(out$symbol) | !nzchar(out$symbol))) stop("Gene evidence DE symbols must be nonempty.", call. = FALSE)
  out$gene_id[is.na(out$gene_id)] <- ""
  out <- unique(out)
  # An ambiguous symbol/ID is retained, never silently collapsed to one row.
  lisa_gene_evidence_sort(out, c("analysis_id", "symbol", "gene_id", "log2FC", "rank_value", "pvalue", "de_fdr"))
}

build_lisa_gene_evidence <- function(scopes, memberships, tier,
    gsea_padj_cutoff = .25, de_padj_cutoff = .05) {
  if (!is.list(scopes) || !length(scopes)) stop("Gene evidence needs at least one individual-analysis scope.", call. = FALSE)
  cuts <- c(gsea_padj_cutoff, de_padj_cutoff)
  if (length(cuts) != 2L || any(!is.finite(cuts) | cuts < 0 | cuts > 1)) stop("Gene evidence FDR cutoffs must be finite in [0,1].", call. = FALSE)
  if (length(tier) != 1L || is.na(tier) || !nzchar(tier)) stop("Gene evidence tier must be explicit.", call. = FALSE)
  member <- lisa_gene_evidence_read(memberships)
  pc <- intersect(c("gs_name", "pathway", "gene_set_id"), names(member))
  gc <- intersect(c("gene_symbol", "symbol"), names(member))
  if (!length(pc) || !length(gc)) stop("Gene evidence requires full TERM2GENE pathway and gene-symbol columns.", call. = FALSE)
  member <- unique(data.frame(pathway = as.character(member[[pc[[1L]]]]), symbol = as.character(member[[gc[[1L]]]]), stringsAsFactors = FALSE))
  if (any(is.na(member$pathway) | !nzchar(member$pathway) | is.na(member$symbol) | !nzchar(member$symbol))) stop("Full gene-set memberships contain missing identifiers.", call. = FALSE)
  # Match the existing GSEA input contract, which uppercases TERM2GENE and
  # standardised ranking symbols. Preserve source spelling without inferring
  # aliases or stripping Ensembl versions/leading zeroes from identifiers.
  source_symbols <- unique(data.frame(symbol = toupper(member$symbol), source_symbol = member$symbol, stringsAsFactors = FALSE))
  member$symbol <- toupper(member$symbol)
  member <- unique(member)
  member <- lisa_gene_evidence_sort(member, c("pathway", "symbol"))
  member_index <- split(member$symbol, member$pathway)
  analyses <- list(); de_by_analysis <- list(); results <- list(); edges <- list(); scope_rows <- list()
  input_identity <- list()
  if (is.character(memberships) && length(memberships) == 1L) input_identity[[1L]] <- data.frame(kind = "TERM2GENE", analysis_id = "", collection = "", filename = basename(memberships), sha256 = digest::digest(file = memberships, algo = "sha256"), stringsAsFactors = FALSE)
  seen <- character()
  for (i in seq_along(scopes)) {
    scope <- scopes[[i]]
    for (key in c("analysis_id", "collection")) if (length(scope[[key]]) != 1L || is.na(scope[[key]]) || !nzchar(scope[[key]])) stop("Gene evidence scopes need analysis_id and collection.", call. = FALSE)
    analysis_id <- as.character(scope$analysis_id); collection <- as.character(scope$collection)
    if (!is.null(scope$analysis_type) && !identical(scope$analysis_type, "single_de")) stop("Gene evidence supports individual analyses only, not contrasts.", call. = FALSE)
    key <- paste(analysis_id, collection, sep = "\r")
    if (key %in% seen) stop("Duplicate gene evidence analysis/collection scope.", call. = FALSE)
    seen <- c(seen, key)
    g <- lisa_gene_evidence_read(scope$gsea_table)
    ledger <- lisa_gene_evidence_read(scope$universe_ledger)
    if (!"pathway" %in% names(ledger) || anyNA(ledger$pathway) || anyDuplicated(ledger$pathway)) stop("Gene evidence needs a unique complete GSEA universe ledger.", call. = FALSE)
    if (nrow(g) && !"pathway" %in% names(g)) stop("Annotated GSEA table lacks pathway.", call. = FALSE)
    for (table in list(g, ledger)) {
      if ("analysis_id" %in% names(table)) {
        value <- unique(as.character(table$analysis_id)); value <- value[!is.na(value) & nzchar(value)]
        if (length(setdiff(value, analysis_id))) stop("Gene evidence analysis scope contradicts input rows.", call. = FALSE)
      }
      for (name in intersect(c("analysis_collection", "universe"), names(table))) {
        value <- unique(as.character(table[[name]])); value <- value[!is.na(value) & nzchar(value)]
        if (length(setdiff(value, collection))) stop("Gene evidence collection scope contradicts input rows.", call. = FALSE)
      }
      if ("tier" %in% names(table)) {
        value <- unique(as.character(table$tier)); value <- value[!is.na(value) & nzchar(value)]
        allowed <- if (identical(tier, "expanded")) c("core", "expanded") else tier
        if (length(setdiff(value, allowed))) stop("Gene evidence tier scope contradicts input rows.", call. = FALSE)
      }
    }
    if (length(setdiff(as.character(g$pathway), as.character(ledger$pathway)))) stop("Annotated sets lie outside the complete gene evidence universe.", call. = FALSE)
    if (length(setdiff(as.character(ledger$pathway), names(member_index)))) stop("Full TERM2GENE membership is missing universe gene sets.", call. = FALSE)
    ledger <- lisa_gene_evidence_sort(ledger, "pathway")
    d <- lisa_gene_evidence_de(scope$de_table, analysis_id, scope$column_mapping)
    contrast <- if (length(scope$positive_contrast) == 1L && !is.na(scope$positive_contrast) && nzchar(scope$positive_contrast)) as.character(scope$positive_contrast) else "not recorded"
    if (is.null(de_by_analysis[[analysis_id]])) {
      de_by_analysis[[analysis_id]] <- d
      analyses[[analysis_id]] <- data.frame(analysis_id = analysis_id, positive_contrast = contrast, stringsAsFactors = FALSE)
    } else {
      if (!identical(de_by_analysis[[analysis_id]], d)) stop("DE evidence differs between collections of the same analysis.", call. = FALSE)
      if (!identical(analyses[[analysis_id]]$positive_contrast, contrast)) stop("Positive-contrast meaning differs within one analysis.", call. = FALSE)
    }
    by_path <- split(seq_len(nrow(g)), as.character(g$pathway))
    columns <- c("NES", "padj", "pval", "leadingEdge")
    for (name in columns) {
      if (name %in% names(g)) {
        conflict <- vapply(by_path, function(rows) length(unique(as.character(g[[name]][rows]))) > 1L, logical(1L))
        if (any(conflict)) stop("Conflicting GSEA measurements across category annotation copies.", call. = FALSE)
      }
    }
    gi <- match(as.character(ledger$pathway), as.character(g$pathway))
    value <- function(name) {
      a <- lisa_gene_evidence_col(ledger, name)
      if (name %in% names(g)) {
        b <- as.character(g[[name]][gi])
        both <- !is.na(a) & !is.na(b)
        if (any(a[both] != b[both])) stop("GSEA ledger and annotation disagree: ", name, call. = FALSE)
        a[is.na(a)] <- b[is.na(a)]
      }
      a
    }
    nes <- lisa_gene_evidence_num(value("NES")); fdr <- lisa_gene_evidence_num(value("padj"))
    eligible_text <- lisa_gene_evidence_col(ledger, "eligible_for_gsea", "true")
    eligible <- tolower(eligible_text) %in% c("true", "1")
    evaluable <- eligible & is.finite(nes) & is.finite(fdr) & fdr >= 0 & fdr <= 1
    le <- value("leadingEdge")
    le_state <- ifelse(is.na(le) | !evaluable, "unavailable", ifelse(nzchar(trimws(le)), "recorded", "empty"))
    le_list <- lapply(seq_along(le), function(j) if (le_state[[j]] == "recorded") lisa_gene_evidence_tokens(le[[j]]) else character())
    invalid <- vapply(seq_along(le_list), function(j) length(setdiff(le_list[[j]], member_index[[as.character(ledger$pathway[[j]])]])) > 0L, logical(1L))
    if (any(invalid)) {
      j <- which(invalid)[[1L]]
      absent <- setdiff(le_list[[j]], member_index[[as.character(ledger$pathway[[j]])]])
      stop("Leading-edge gene is not a member of its full TERM2GENE set: ", ledger$pathway[[j]], " / ", paste(absent, collapse = ", "), call. = FALSE)
    }
    categories <- vapply(as.character(ledger$pathway), function(pathway) {
      rows <- by_path[[pathway]]
      if (is.null(rows) || !"category_id" %in% names(g)) return("")
      ids <- unique(as.character(g$category_id[rows])); ids <- ids[!is.na(ids) & nzchar(ids) & ids != "OTHER_UNCLASSIFIED"]
      paste(sort(ids, method = "radix"), collapse = ";")
    }, character(1L))
    classification <- ifelse(nzchar(categories), "classified", "unclassified")
    if ("classification_status" %in% names(ledger)) {
      known <- as.character(ledger$classification_status)
      if (any(!is.na(known) & nzchar(known) & known != classification)) stop("Gene evidence annotation lacks the ledger's declared categories.", call. = FALSE)
    }
    results[[i]] <- data.frame(analysis_id = analysis_id, collection = collection,
      pathway = as.character(ledger$pathway), NES = nes, gsea_fdr = fdr,
      gsea_pvalue = lisa_gene_evidence_num(value("pval")), evaluable = evaluable,
      significant = evaluable & fdr <= gsea_padj_cutoff,
      classification_status = classification, category_ids = categories,
      leading_edge_state = le_state,
      eligibility_status = lisa_gene_evidence_col(ledger, c("gsea_result_status", "gsea_eligibility_status"), "not recorded"),
      stringsAsFactors = FALSE)
    lens <- lengths(le_list)
    edges[[i]] <- data.frame(analysis_id = rep(analysis_id, sum(lens)), collection = rep(collection, sum(lens)),
      pathway = rep(as.character(ledger$pathway), lens), symbol = unlist(le_list, use.names = FALSE), stringsAsFactors = FALSE)
    scope_rows[[i]] <- data.frame(analysis_id = analysis_id, collection = collection, tier = tier, n_universe_sets = nrow(ledger), stringsAsFactors = FALSE)
    for (name in c("de_table", "universe_ledger", "gsea_table", "column_mapping")) if (is.character(scope[[name]]) && length(scope[[name]]) == 1L) {
      path <- scope[[name]]
      input_identity[[length(input_identity) + 1L]] <- data.frame(kind = name, analysis_id = analysis_id, collection = collection, filename = basename(path), sha256 = digest::digest(file = path, algo = "sha256"), stringsAsFactors = FALSE)
    }
  }
  sets <- lisa_gene_evidence_sort(do.call(rbind, results), c("analysis_id", "collection", "pathway"))
  de <- lisa_gene_evidence_sort(do.call(rbind, de_by_analysis), c("analysis_id", "symbol", "gene_id", "log2FC", "rank_value", "pvalue", "de_fdr"))
  member <- member[member$pathway %in% sets$pathway, , drop = FALSE]; rownames(member) <- NULL
  identifiers <- unique(de[c("symbol", "gene_id", "gene_id_source")])
  absent <- setdiff(unique(member$symbol), identifiers$symbol)
  identifiers <- rbind(identifiers, data.frame(symbol = absent, gene_id = rep("", length(absent)), gene_id_source = rep("not recorded", length(absent)), stringsAsFactors = FALSE))
  identifiers <- lisa_gene_evidence_sort(identifiers, c("symbol", "gene_id"))
  provenance <- if (length(input_identity)) do.call(rbind, input_identity) else data.frame(kind = character(), analysis_id = character(), collection = character(), filename = character(), sha256 = character())
  structure(list(metadata = list(schema_version = "1.0", tier = tier,
    gsea_padj_cutoff = gsea_padj_cutoff, de_padj_cutoff = de_padj_cutoff,
    scope = "Individual analyses only. No contrast evidence or cross-tier comparison.",
    membership = "Complete TERM2GENE membership with the existing GSEA uppercase-symbol convention; source spelling is retained separately. No inferred aliases. Leading edge is a separate observed property.",
    missing = "Missing measurements remain unavailable; no inferred scores or expression values."),
    analyses = lisa_gene_evidence_sort(do.call(rbind, analyses), "analysis_id"),
    scopes = lisa_gene_evidence_sort(do.call(rbind, scope_rows), c("analysis_id", "collection")),
    identifiers = identifiers, de = de, sets = sets,
    leading_edges = lisa_gene_evidence_sort(do.call(rbind, edges), c("analysis_id", "collection", "pathway", "symbol")),
    memberships = member, source_symbols = lisa_gene_evidence_sort(source_symbols, c("symbol", "source_symbol")),
    provenance = lisa_gene_evidence_sort(provenance, c("kind", "analysis_id", "collection", "filename"))),
    class = "lisa_gene_evidence")
}

query_lisa_gene_evidence <- function(evidence, query, include_all = FALSE) {
  if (!inherits(evidence, "lisa_gene_evidence")) stop("Expected lisa_gene_evidence.", call. = FALSE)
  ids <- evidence$identifiers
  biological_id <- !tolower(ids$gene_id_source) %in% c("rowname", "rownames", ".rownames")
  candidates <- unique(ids$symbol[ids$symbol == query | (biological_id & nzchar(ids$gene_id) & ids$gene_id == query)])
  if (length(candidates) != 1L) return(list(status = if (length(candidates)) "ambiguous" else "not_found", candidates = candidates))
  symbol <- candidates[[1L]]
  sets <- evidence$sets[evidence$sets$pathway %in% evidence$memberships$pathway[evidence$memberships$symbol == symbol], , drop = FALSE]
  if (!isTRUE(include_all)) sets <- sets[sets$significant, , drop = FALSE]
  keys <- paste(sets$analysis_id, sets$collection, sets$pathway, sep = "\r")
  le <- evidence$leading_edges[evidence$leading_edges$symbol == symbol, , drop = FALSE]
  present <- keys %in% paste(le$analysis_id, le$collection, le$pathway, sep = "\r")
  sets$leading_edge <- ifelse(sets$leading_edge_state == "unavailable", "unavailable", ifelse(present, "yes", "no"))
  list(status = "matched", symbol = symbol, candidates = candidates,
    de = evidence$de[evidence$de$symbol == symbol, , drop = FALSE], sets = sets)
}

lisa_gene_evidence_payload <- function(evidence) {
  symbols <- sort(unique(evidence$identifiers$symbol), method = "radix")
  pathways <- sort(unique(evidence$sets$pathway), method = "radix")
  sets <- evidence$sets
  keys <- paste(sets$analysis_id, sets$collection, sets$pathway, sep = "\r")
  member <- split(match(evidence$memberships$pathway, pathways) - 1L, factor(evidence$memberships$symbol, levels = symbols), drop = FALSE)
  leading <- split(match(evidence$leading_edges$symbol, symbols) - 1L,
    factor(paste(evidence$leading_edges$analysis_id, evidence$leading_edges$collection, evidence$leading_edges$pathway, sep = "\r"), levels = keys), drop = FALSE)
  sets$pathway <- match(sets$pathway, pathways) - 1L
  de <- evidence$de; de$symbol <- match(de$symbol, symbols) - 1L
  columns <- function(x) unname(lapply(x, I))
  list(metadata = evidence$metadata, analyses = evidence$analyses, scopes = evidence$scopes,
    symbols = I(symbols), pathways = I(pathways), member_sets = unname(lapply(member, I)),
    de_columns = I(names(de)), de = columns(de), set_columns = I(names(sets)),
    sets = columns(sets), leading_edges = unname(lapply(leading, I)), provenance = evidence$provenance)
}

render_lisa_gene_evidence <- function(evidence, output_dir) {
  if (!inherits(evidence, "lisa_gene_evidence")) stop("Expected lisa_gene_evidence.", call. = FALSE)
  lisa_guarded_dir_create(output_dir)
  table_dir <- file.path(output_dir, "tables"); lisa_guarded_dir_create(table_dir)
  for (name in c("analyses", "scopes", "identifiers", "de", "sets", "leading_edges", "memberships", "source_symbols", "provenance")) {
    write_lisa_tsv(evidence[[name]], file.path(table_dir, paste0(name, ".tsv")))
  }
  json <- jsonlite::toJSON(lisa_gene_evidence_payload(evidence), auto_unbox = TRUE,
    dataframe = "rows", na = "null", null = "null", digits = 17)
  json <- gsub("<", "\\u003c", json, fixed = TRUE)
  root <- system.file("gene-evidence", package = "lisaR")
  if (!nzchar(root)) stop("Installed gene-evidence assets are missing.", call. = FALSE)
  assets <- file.path(output_dir, "assets"); lisa_guarded_dir_create(assets)
  for (extension in c("js", "css")) lisa_guarded_write(file.path(assets, paste0("gene-evidence.", extension)), function(target) {
    if (!file.copy(file.path(root, paste0("viewer.", extension)), target, overwrite = TRUE)) stop("Cannot copy gene evidence asset.", call. = FALSE)
  })
  html <- readLines(file.path(root, "viewer.html"), warn = FALSE)
  placeholder <- which(html == "<!-- GENE_EVIDENCE_DATA -->")
  if (length(placeholder) != 1L) stop("Gene evidence HTML data placeholder is invalid.", call. = FALSE)
  # Assign the complete line, rather than regex replacement, so JSON escape
  # sequences in identifiers cannot be consumed by sub() replacement rules.
  html[[placeholder]] <- paste0('<script type="application/json" id="gene-evidence-data">', json, '</script>')
  html <- lisa_present_evidence_html(html, output_dir, evidence$metadata, active = "genes")
  path <- file.path(output_dir, "index.html")
  lisa_guarded_write(path, function(target) writeLines(html, target, useBytes = TRUE))
  invisible(list(html = path, status = "completed", n_genes = length(unique(evidence$identifiers$symbol)), n_memberships = nrow(evidence$memberships)))
}
