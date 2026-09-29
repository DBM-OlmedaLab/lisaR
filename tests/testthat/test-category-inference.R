# Independent exhaustive closed-testing oracle, deliberately not the hommel
# shortcut: enumerate all local intersections and every supersets' maximum.
inference_oracle <- function(p, alpha = .05) {
  m <- length(p)
  sets <- lapply(0:(2^m - 1L), function(mask) which(as.logical(intToBits(mask))[seq_len(m)]))
  local <- vapply(sets, function(ix) {
    n <- length(ix)
    if (!n) return(1)
    min(1, n * sum(1 / seq_len(n)) * min(sort(p[ix]) / seq_len(n)))
  }, numeric(1))
  closed <- vapply(sets, function(ix) {
    if (!length(ix)) return(1)
    max(local[vapply(sets, function(j) all(ix %in% j), logical(1))])
  }, numeric(1))
  d <- vapply(sets, function(ix) {
    possible_nulls <- lengths(sets)[vapply(sets, function(j) all(j %in% ix), logical(1)) & closed > alpha]
    length(ix) - max(c(0L, possible_nulls))
  }, integer(1))
  list(sets = sets, p = closed, d = d)
}

test_that("robust closed bounds match independent exhaustive closure", {
  cases <- list(c(.001, .03, .2, 1), c(0, .01, .5), rep(0, 4), rep(1, 4),
    c(.002, .002, .04, .5, .9), c(.049, .051), c(.05))
  set.seed(311)
  cases <- c(cases, replicate(12, 10^runif(5, -5, 0), simplify = FALSE))
  for (p in cases) {
    fit <- hommel::hommel(p, simes = FALSE)
    oracle <- inference_oracle(p)
    for (i in seq_along(oracle$sets)[-1L]) {
      expect_no_warning(b <- lisaR:::lisa_inference_closed_bounds(fit, p, oracle$sets[[i]], .05))
      expect_equal(b$p, oracle$p[[i]], tolerance = 1e-12)
      expect_equal(b$discoveries, oracle$d[[i]])
      expect_identical(b$p <= .05, b$discoveries >= 1L)
    }
  }
})

test_that("inclusive category threshold and lower bound agree at machine boundaries", {
  for (m in 1:12) for (k in seq_len(m)) for (multiplier in c(1 - 1e-14, 1, 1 + 1e-14)) {
    p <- c(rep(.05 * k / (m * sum(1 / seq_len(m))) * multiplier, k), rep(1, m - k))
    b <- lisaR:::lisa_inference_closed_bounds(hommel::hommel(p, simes = FALSE), p, seq_len(k), .05)
    expect_identical(b$p <= .05, b$discoveries >= 1L)
  }
})

inference_fixture <- function() {
  ledger <- data.frame(collection = c("GOBP-C2", "GOBP-C2", "GOBP-C2", "GOBP-C2", "GOBP-C2", "PATHWAYS"),
    pathway = c("a", "b", "unclassified", "missing", "ineligible", "a"),
    hypothesis_id = c("a", "b", "u", "missing", "ineligible", "a"),
    eligible_for_gsea = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE),
    pval = c(.00001, .9, .2, NA, NA, .00002), NES = c(NA, -1, 1, 2, NA, 2))
  assignments <- data.frame(collection = c(rep("GOBP-C2", 7), "PATHWAYS"),
    pathway = c("a", "b", "missing", "ineligible", "absent", "a", "missing", "a"),
    category_id = c(rep("all", 5), "strong", "unavailable", "path"))
  categories <- data.frame(collection = c(rep("GOBP-C2", 3), "PATHWAYS"),
    category_id = c("all", "strong", "unavailable", "path"),
    category_display_name = c("All", "Strong", "Unavailable", "Path"))
  list(ledger = ledger, assignments = assignments, categories = categories)
}

test_that("family includes unclassified/missing nulls and conservative unique reruns", {
  z <- inference_fixture()
  result <- do.call(lisaR:::lisa_category_inference, z)
  expect_equal(nrow(result$family), 4L)
  expect_equal(result$family$p_for_inference[match(c("a", "missing"), result$family$hypothesis_id)], c(.00002, 1))
  all <- result$categories[result$categories$category_id == "all", ]
  expect_equal(all$n_sets_total, 5L)
  expect_equal(all$n_sets_evaluable, 2L)
  expect_equal(all$n_sets_missing_p, 1L)
  expect_equal(all$n_missing_NES, 1L)
  expect_true(all$significant)
  expect_equal(result$categories$status[result$categories$category_id == "unavailable"], "not_evaluable")
  expect_true(is.na(result$categories$category_p_adjusted[result$categories$category_id == "unavailable"]))
  z$ledger$pval[6] <- NA_real_
  result2 <- do.call(lisaR:::lisa_category_inference, z)
  expect_equal(result2$family$p_for_inference[result2$family$hypothesis_id == "a"], 1)
  expect_false(any(result2$categories$significant))
})

test_that("invalid inputs and HALLMARKS fail closed", {
  z <- inference_fixture()
  z$ledger$pval[1] <- -.1
  expect_error(do.call(lisaR:::lisa_category_inference, z), "Raw GSEA p-values")
  z <- inference_fixture()
  z$ledger$collection[1] <- "HALLMARKS"
  expect_error(do.call(lisaR:::lisa_category_inference, z), "HALLMARKS")
  z <- inference_fixture()
  z$ledger <- rbind(z$ledger, z$ledger[1, ])
  expect_error(do.call(lisaR:::lisa_category_inference, z), "Duplicate gene-set")
})

test_that("canonical TSV retains threshold bits and all categories without fake NES", {
  result <- do.call(lisaR:::lisa_category_inference, inference_fixture())
  table <- result$categories[result$categories$collection == "GOBP-C2", ]
  table$analysis_id <- "A"
  source <- lisaR:::lisa_category_inference_source(table,
    result$members[result$members$collection == "GOBP-C2", ])
  expect_setequal(unique(source$category_id), c("all", "strong", "unavailable"))
  expect_true(is.na(source$NES[source$category_id == "strong"]))
  expect_no_error(lisaR:::lisa_plot_category_inference(source[source$category_id == "strong", ]))
  source$category_p_adjusted[1] <- .05 * (1 + 1e-15)
  f <- tempfile(fileext = ".tsv")
  on.exit(unlink(f), add = TRUE)
  lisaR:::lisa_write_inference_tsv(source, f)
  stored <- lisaR:::lisa_read_figure_source_tsv(f)
  expect_identical(stored$category_p_adjusted, source$category_p_adjusted)
  p <- lisaR:::lisa_plot_category_inference(source)
  expect_true(inherits(p, "ggplot"))
  labels <- p$scales$get_scales("y")$labels
  expect_identical(labels, table$category_display_name)
  expect_false(any(grepl("*", labels, fixed = TRUE)))
  expect_false(any(grepl("ns|NE", labels)))
  png <- tempfile(fileext = ".png")
  on.exit(unlink(png), add = TRUE)
  lisaR:::lisa_category_inference_save_plot(source, png)
  bytes <- readBin(png, "raw", n = 24)
  width <- sum(as.integer(bytes[17:20]) * 256^(3:0))
  height <- sum(as.integer(bytes[21:24]) * 256^(3:0))
  expect_equal(width, source$figure_width[[1]] * source$figure_dpi[[1]])
  expect_equal(height, source$figure_height[[1]] * source$figure_dpi[[1]])
})

test_that("saved inference is copy-only, validates source bytes and reseals honestly", {
  root <- tempfile("inference-archive-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  origin <- file.path(root, "source")
  out <- file.path(root, "derived")
  collection <- file.path(origin, "outputs", "single_de", "A", "collection_GOBP-C2")
  for (d in c("qc", "inputs", "enrichment", "lisa_tables", "plots")) dir.create(file.path(collection, d), recursive = TRUE)
  dir.create(file.path(origin, "config")); dir.create(file.path(origin, "audit"))
  write <- lisaR:::write_lisa_tsv
  t2g <- file.path(root, "term2gene.tsv")
  write(data.frame(gs_name = rep("set1", 2), gene_symbol = c("G1", "G2")), t2g)
  write(data.frame(analysis_id = "A", label = "Analysis A"), file.path(origin, "config", "de_index.tsv"))
  write(data.frame(analysis_collection = "GOBP-C2"), file.path(origin, "lisa_collection_registry.tsv"))
  write(data.frame(analysis_id = "A", collection = "GOBP-C2", status = "completed"), file.path(origin, "single_de_status.tsv"))
  write(data.frame(key = "term2gene", value = t2g), file.path(origin, "lisa_pipeline_metadata.tsv"))
  write(data.frame(role = "term2gene_resource", sha256 = lisaR:::lisa_sha256_file(t2g)), file.path(origin, "resource_inventory.tsv"))
  write(data.frame(pathway = "set1", eligible_for_gsea = TRUE, ranked_gene_count = 2, pval = .001, NES = 2), file.path(collection, "qc", "A_GSEA_universe_ledger.tsv"))
  write(data.frame(symbol = c("G1", "G2"), rank_value = c(2, -1)), file.path(collection, "inputs", "A_ranked_genes.tsv"))
  write(data.frame(pathway = "set1", category_id = "category1"), file.path(collection, "enrichment", "A_GSEA_semantic_annotated.tsv"))
  write(data.frame(lisa_dictionary_path = "not_used"), file.path(collection, "qc", "A_ranking_qc.tsv"))
  write(data.frame(category_id = "category1", display_name = "Category one"), file.path(collection, "qc", "lisa_category_palette.tsv"))
  write(data.frame(gene_set_id = "set1", category_id = "category1"), file.path(collection, "qc", "lisa_category_membership.tsv"))
  manifest <- lisaR:::lisa_run_manifest(origin)
  write(manifest, file.path(origin, "run_manifest.tsv"))
  source_hash <- lisaR:::lisa_sha256_file(file.path(origin, "run_manifest.tsv"))
  add_lisa_category_inference(origin, out, plot_formats = character())
  expect_identical(lisaR:::lisa_sha256_file(file.path(origin, "run_manifest.tsv")), source_hash)
  expect_identical(lisaR:::lisa_sha256_file(file.path(out, "source_provenance", "original_run_manifest.tsv")), source_hash)
  expect_false(file.exists(file.path(out, "run_manifest.tsv")))
  table <- lisaR:::read_lisa_tsv(file.path(out, "outputs", "category_inference", "A", "category_results.tsv"))
  expect_true(table$significant)
  expect_equal(table$minimum_enriched_sets, 1)
  lisaR:::lisa_finalize_category_inference_derivation(out)
  check <- lisaR:::lisa_validate_built_run_manifest(out, lisaR:::read_lisa_tsv(file.path(out, "run_manifest.tsv")))
  expect_identical(check$gate, "PASS")
  expect_error(add_lisa_category_inference(origin, out), "must be new")
  # Legacy saved reports have no membership receipt and encode unclassified
  # annotations as blank TSV fields. Reconstruct only actual assignments,
  # while retaining the unclassified hypothesis in the multiplicity family.
  unlink(file.path(collection, "qc", "lisa_category_membership.tsv"))
  dictionary <- file.path(root, "dictionary.tsv")
  write(data.frame(universe = "GOBP-C2", gene_set_id = "set1", category_id = "category1"), dictionary)
  write(data.frame(lisa_dictionary_path = dictionary, n_lisa_gene_sets = 1), file.path(collection, "qc", "A_ranking_qc.tsv"))
  write(data.frame(gs_name = rep(c("set1", "unclassified"), each = 2), gene_symbol = rep(c("G1", "G2"), 2)), t2g)
  write(data.frame(role = c("term2gene_resource", "dictionary_resource"),
    sha256 = vapply(c(t2g, dictionary), lisaR:::lisa_sha256_file, character(1))), file.path(origin, "resource_inventory.tsv"))
  write(data.frame(pathway = c("set1", "unclassified"), eligible_for_gsea = TRUE,
    ranked_gene_count = 2, pval = c(.001, .8), NES = c(2, -.5)), file.path(collection, "qc", "A_GSEA_universe_ledger.tsv"))
  annotation_path <- file.path(collection, "enrichment", "A_GSEA_semantic_annotated.tsv")
  write(data.frame(pathway = c("set1", "unclassified"), category_id = c("category1", "")), annotation_path)
  manifest <- lisaR:::lisa_run_manifest(origin)
  write(manifest, file.path(origin, "run_manifest.tsv"))
  fallback <- file.path(root, "derived-legacy")
  add_lisa_category_inference(origin, fallback, plot_formats = character())
  legacy_table <- lisaR:::read_lisa_tsv(file.path(fallback, "outputs", "category_inference", "A", "category_results.tsv"))
  expect_equal(legacy_table$n_sets_total, 1)
  expect_equal(legacy_table$n_sets_evaluable, 1)
  expect_equal(legacy_table$family_size, 2)
  expect_true(legacy_table$significant)
  # A nonempty incompatible assignment must still fail, not be ignored.
  write(data.frame(pathway = c("set1", "unclassified"), category_id = c("category1", "unexpected_category")), annotation_path)
  manifest <- lisaR:::lisa_run_manifest(origin)
  write(manifest, file.path(origin, "run_manifest.tsv"))
  expect_error(add_lisa_category_inference(origin, file.path(root, "bad-legacy"), plot_formats = character()),
    "Saved and reconstructed category memberships disagree")
  writeLines("changed", file.path(collection, "inputs", "A_ranked_genes.tsv"))
  expect_error(add_lisa_category_inference(origin, file.path(root, "rejected")), "Source run failed integrity")
  expect_false(dir.exists(file.path(root, "rejected")))
})
