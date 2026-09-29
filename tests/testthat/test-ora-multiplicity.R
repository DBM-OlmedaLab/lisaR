ora_family_fixture <- function() {
  root <- testthat::test_path("fixtures", "ora-bh-family")
  list(
    root = root,
    de_path = file.path(root, "de.tsv"),
    term2gene_path = file.path(root, "term2gene.tsv"),
    dictionary_path = file.path(root, "dictionary.tsv"),
    category_map_path = file.path(root, "category_map.tsv"),
    de = utils::read.delim(file.path(root, "de.tsv"), sep = "\t", stringsAsFactors = FALSE),
    term2gene = utils::read.delim(file.path(root, "term2gene.tsv"), sep = "\t", stringsAsFactors = FALSE),
    dictionary = utils::read.delim(file.path(root, "dictionary.tsv"), sep = "\t", stringsAsFactors = FALSE)
  )
}

test_that("ORA BH uses every eligible unique gene set before category annotation", {
  fixture <- ora_family_fixture()
  background <- fixture$de$symbol
  query <- fixture$de$symbol[fixture$de$padj <= 0.05 & fixture$de$log2FoldChange >= 0.58]

  observed <- lisaR:::run_ora_one(
    genes = query,
    bg = background,
    term2gene = fixture$term2gene,
    lisa_dict = fixture$dictionary,
    direction = "UP",
    padj_cutoff = 0.05,
    analysis_id = "analysis_A",
    collection = "GOBP-C2"
  )

  eligible <- split(toupper(fixture$term2gene$gene_symbol), fixture$term2gene$gs_name)
  eligible <- lapply(eligible, function(set) intersect(unique(set), toupper(background)))
  eligible <- eligible[lengths(eligible) > 0L]
  overlap <- lengths(lapply(eligible, intersect, y = toupper(query)))
  raw_p <- stats::phyper(
    overlap - 1L, lengths(eligible), length(background) - lengths(eligible),
    length(query), lower.tail = FALSE
  )
  raw_p[overlap == 0L] <- 1
  expected_padj <- stats::p.adjust(raw_p, method = "BH")

  one_per_set <- observed[!duplicated(observed$pathway), , drop = FALSE]
  one_per_set <- one_per_set[match(names(eligible), one_per_set$pathway), , drop = FALSE]
  expect_identical(unname(one_per_set$pvalue), unname(raw_p))
  expect_identical(unname(one_per_set$padj), unname(expected_padj))
  expect_identical(unname(one_per_set$selected), unname(expected_padj <= 0.05))
  expect_true(all(one_per_set$pvalue[one_per_set$overlap == 0L] == 1))

  expect_identical(unique(observed$analysis_id), "analysis_A")
  expect_identical(unique(observed$collection), "GOBP-C2")
  expect_identical(unique(observed$family_id), "ORA:analysis_A:GOBP-C2:UP")
  expect_identical(unique(observed$direction), "UP")
  expect_identical(unique(observed$n_eligible), 4L)
  expect_identical(unique(observed$n_tested), 4L)
  expect_identical(unique(observed$method), "BH")

  # One tested gene set maps to two LISA categories, but the BH denominator is
  # the four unique eligible gene sets rather than the five annotated rows.
  expect_identical(sum(observed$pathway == "GOBP_ORA_HIT_STRONG"), 2L)
  expect_identical(nrow(observed), 5L)
  expect_identical(unique(observed$n_tested), length(unique(observed$pathway)))
})

test_that("ORA writes all-tested evidence and preserves the significant legacy view", {
  fixture <- ora_family_fixture()
  output_dir <- tempfile("lisaR-ora-family-")

  result <- lisaR:::run_LISA_DE(
    input = fixture$de_path,
    input_type = "de_table",
    output_dir = output_dir,
    comparison_name = "ora_fixture",
    species = "Homo sapiens",
    lisa_dictionary = "fixture",
    lisa_dictionary_path = fixture$dictionary_path,
    term2gene_path = fixture$term2gene_path,
    category_map_path = fixture$category_map_path,
    universes = "GOBP-C2",
    run_gsea = FALSE,
    run_ora = TRUE,
    outputs = c("tables", "qc"),
    plots = character(),
    export_formats = "tsv",
    file_label_prefix = "semantic",
    lisa_project_root = fixture$root,
    verbose = FALSE
  )

  all_path <- file.path(output_dir, "enrichment", "ora_fixture_ORA_all_tested_annotated.tsv")
  legacy_path <- file.path(output_dir, "enrichment", "ora_fixture_ORA_semantic_annotated.tsv")
  expect_true(file.exists(all_path))
  expect_true(file.exists(legacy_path))

  all_tested <- lisaR:::read_lisa_tsv(all_path)
  legacy <- lisaR:::read_lisa_tsv(legacy_path)
  contract <- c(
    "analysis_id", "collection", "family_id", "direction", "n_eligible",
    "n_tested", "method", "selected"
  )
  expect_true(all(contract %in% names(all_tested)))
  expect_identical(unique(all_tested$n_tested), 4L)
  expect_setequal(
    unique(all_tested$family_id),
    c("ORA:ora_fixture:GOBP-C2:UP", "ORA:ora_fixture:GOBP-C2:DOWN")
  )
  down_family <- all_tested[all_tested$direction == "DOWN", , drop = FALSE]
  expect_true(all(down_family$pvalue == 1 & down_family$padj == 1 & !down_family$selected))
  expect_true(any(all_tested$overlap == 0L & all_tested$pvalue == 1 & !all_tested$selected))
  expect_true(all(legacy$selected))
  expect_setequal(unique(legacy$pathway), unique(all_tested$pathway[all_tested$selected]))
  expect_identical(result$ora, result$ora_all_tested[result$ora_all_tested$selected, , drop = FALSE])
})
