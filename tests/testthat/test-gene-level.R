test_that("gene-level tables can be built from DE and GSEA leading edges", {
  base <- system.file("extdata/minimal", package = "lisaR")
  out <- file.path(tempdir(), "lisaR-gene-level")
  result <- build_lisa_gene_level_tables(
    de_table = file.path(base, "de_A.tsv"),
    gsea_table = file.path(base, "gsea_A.tsv"),
    output_dir = out
  )

  expect_true(file.exists(result$files[["leading_edge_gene_pathways"]]))
  expect_true(file.exists(result$files[["gene_category_contributions"]]))
  expect_true(nrow(result$leading_edge) >= 4)
  expect_true("GENE1" %in% result$gene_category$symbol)
})

test_that("gene-level tables integrate ORA overlap genes", {
  base <- system.file("extdata/minimal", package = "lisaR")
  out <- file.path(tempdir(), "lisaR-gene-level-ora")
  ora <- data.frame(
    category_id = "CAT_SIGNAL",
    pathway = "ORA_SIGNALING",
    direction = "UP",
    padj = 0.01,
    genes = "GENE1/GENE2",
    stringsAsFactors = FALSE
  )
  result <- build_lisa_gene_level_tables(
    de_table = file.path(base, "de_A.tsv"),
    gsea_table = file.path(base, "gsea_A.tsv"),
    ora_table = ora,
    output_dir = out
  )

  both <- result$gene_category[result$gene_category$symbol == "GENE1" & result$gene_category$category_id == "CAT_SIGNAL", , drop = FALSE]
  ora_only <- result$gene_category[result$gene_category$symbol == "GENE2" & result$gene_category$category_id == "CAT_SIGNAL", , drop = FALSE]

  expect_equal(nrow(both), 1)
  expect_true(isTRUE(both$in_gsea_leading_edge[[1]]))
  expect_true(isTRUE(both$in_ora_overlap[[1]]))
  expect_match(both$source_evidence[[1]], "GSEA_LE")
  expect_match(both$source_evidence[[1]], "ORA_UP")

  expect_equal(nrow(ora_only), 1)
  expect_false(isTRUE(ora_only$in_gsea_leading_edge[[1]]))
  expect_true(isTRUE(ora_only$in_ora_overlap[[1]]))
  expect_true(is.finite(ora_only$gene_contribution_score[[1]]))
})
