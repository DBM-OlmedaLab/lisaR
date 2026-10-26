# A missing FULL gallery is completed-empty only when the exact saved scope
# proves zero eligible members. Do not infer this from directories or defaults.
lisa_report_selection_state <- function(project_dir, analysis, collection) {
  fail <- function(reason) stop("LISA-REPORT-EMPTY-001 ", analysis, "/", collection,
    ": ", reason, call. = FALSE)
  root <- file.path(project_dir, "outputs", "single_de", analysis,
    paste0("collection_", collection))
  one <- function(dir, prefix, suffix) {
    paths <- list.files(file.path(root, dir), full.names = TRUE)
    paths <- paths[startsWith(basename(paths), prefix) & endsWith(paths, suffix)]
    if (length(paths) != 1L) fail(paste("missing or ambiguous", dir, "evidence"))
    tryCatch(read_lisa_tsv(paths[[1L]]), error = function(e) fail(conditionMessage(e)))
  }
  summary <- one("lisa_tables", paste0(analysis, "_"), "_GSEA_category_summary.tsv")
  gsea <- one("enrichment", paste0(analysis, "_GSEA_"), "_annotated.tsv")
  if (!all(c("category_id", "n_genesets", "n_genesets_mapped",
      "n_genesets_evaluable", "n_genesets_significant", "mean_NES",
      "gsea_padj_cutoff") %in% names(summary)))
    fail("invalid category-summary schema")
  # No historical 0.25 fallback is evidence that this selection was evaluated.
  cutoffs <- numeric()
  contract_path <- file.path(project_dir, "contract_manifest.tsv")
  if (file.exists(contract_path)) {
    contract <- read_lisa_tsv(contract_path)
    if (!all(c("key", "value") %in% names(contract))) fail("invalid recorded contract")
    hit <- contract$value[contract$key == "gsea_padj_cutoff"]
    if (length(hit) != 1L) fail("missing or ambiguous recorded cutoff")
    cutoffs <- suppressWarnings(as.numeric(hit))
  }
  cutoffs <- c(cutoffs, suppressWarnings(as.numeric(summary$gsea_padj_cutoff)))
  if (!length(cutoffs) || any(!is.finite(cutoffs)) ||
      any(cutoffs < 0 | cutoffs > 1) || length(unique(cutoffs)) != 1L)
    fail("missing, invalid or inconsistent recorded cutoff")
  # Reject explicit foreign owners instead of treating their zeroes as ours.
  for (x in list(gsea, summary)) for (key in intersect(
      c("analysis_id", "analysis_collection", "collection"), names(x))) {
    expected <- if (key == "analysis_id") analysis else collection
    if (anyNA(x[[key]]) || any(as.character(x[[key]]) != expected)) fail("foreign scope evidence")
  }
  stats <- tryCatch(lisa_category_nes_statistics(gsea, cutoffs[[1L]]),
    error = function(e) fail(conditionMessage(e)))
  # Include mapped-but-unsupported categories from the complete summary too.
  ids <- as.character(summary$category_id)
  if (!all(stats$category_id %in% ids)) fail("category summary omits mapped members")
  checked <- tryCatch(lisa_category_nes_statistics(gsea, cutoffs[[1L]], ids),
    error = function(e) fail(conditionMessage(e)))
  tryCatch(lisa_nes_variant_summary(summary, checked),
    error = function(e) fail(conditionMessage(e)))
  for (key in intersect(c("n_genesets_mapped", "n_genesets_evaluable",
      "n_genesets_significant"), names(summary))) {
    count <- suppressWarnings(as.numeric(summary[[key]]))
    if (any(!is.finite(count) | count != checked[[key]]))
      fail(paste("inconsistent", key))
  }
  positive <- stats$category_id[stats$n_genesets_significant > 0L]
  if (!setequal(positive, ids[checked$n_genesets_significant > 0L]))
    fail("category summary omits significant members")
  if (length(positive)) return("missing")
  if (any(stats$n_genesets_unavailable > 0L)) fail("unavailable statistics cannot prove an empty selection")
  "completed_empty"
}
