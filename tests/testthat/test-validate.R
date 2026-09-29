test_that("DE and contrast indexes validate", {
  base <- system.file("extdata/minimal", package = "lisaR")
  de_index <- validate_lisa_de_index(file.path(base, "de_index.tsv"), base_dir = base)
  expect_equal(nrow(de_index), 2)

  contrast_index <- validate_lisa_contrast_index(file.path(base, "contrast_index.tsv"), de_index)
  expect_equal(nrow(contrast_index), 1)
})

test_that("pipeline plan is deterministic", {
  base <- system.file("extdata/minimal", package = "lisaR")
  de_index <- validate_lisa_de_index(file.path(base, "de_index.tsv"), base_dir = base)
  contrast_index <- validate_lisa_contrast_index(file.path(base, "contrast_index.tsv"), de_index)
  plan <- lisa_pipeline_plan(de_index, contrast_index, collections = c("GOBP-C2"))
  expect_true(all(c("single_de_lisa", "category_contrasts", "root_html_report") %in% plan$stage))
  expect_false("single_de_leading_edge_gene_heatmaps" %in% plan$stage)
  expect_false("single_de_enrichmentmap" %in% plan$stage)
  expect_false(any(grepl("macrogroup_heatmap", plan$stage, fixed = TRUE)))
})

test_that("collection registry normalizes analysis and dictionary names", {
  expect_equal(
    normalize_lisa_collections(c("gobp", "lisa_pathways", "h")),
    c("GOBP-C2", "PATHWAYS", "HALLMARKS")
  )

  registry <- lisa_collection_registry(c("GOBP-C2", "PATHWAYS", "HALLMARKS"))
  expect_equal(registry$analysis_collection, c("GOBP-C2", "PATHWAYS", "HALLMARKS"))
  expect_equal(registry$dictionary_id, c("LISA_GOBP_C2", "LISA_PATHWAYS", "MSIGDB_HALLMARKS"))
  expect_false(registry$allow_missing_dictionary[registry$analysis_collection == "PATHWAYS"])
  # Input contract: HALLMARKS contrasts are enabled by the provisional testing default.
  expect_true(registry$run_contrasts[registry$analysis_collection == "HALLMARKS"])
  expect_false(any(registry$run_enrichmentmap))
})

test_that("ranked GSEA is the default and ORA requires an explicit opt-in", {
  cfg <- list(
    pipeline = list(
      schema_version = "1.0.0", profile = "targeted", evidence_mode = "full_de",
      duplicate_policies = list(
        de_table_duplicate_policy = "error", matrix_duplicate_policy = "error",
        mapped_id_collision_policy = "error"
      )
    )
  )
  contract <- lisaR:::lisa_validate_pipeline_config(cfg)
  expect_false(contract$ora$value)
  expect_identical(contract$ora$source, "default")
  cfg$pipeline$run_ora <- TRUE
  contract <- lisaR:::lisa_validate_pipeline_config(cfg)
  expect_true(contract$ora$value)
  expect_identical(contract$ora$source, "explicit")
  expect_false(formals(lisaR:::run_LISA_DE)$run_ora)
})

test_that("default pipeline includes PATHWAYS and HALLMARKS without macrogroup heatmaps", {
  base <- system.file("extdata/minimal", package = "lisaR")
  de_index <- validate_lisa_de_index(file.path(base, "de_index.tsv"), base_dir = base)
  contrast_index <- validate_lisa_contrast_index(file.path(base, "contrast_index.tsv"), de_index)
  plan <- lisa_pipeline_plan(de_index, contrast_index)

  expect_true(any(grepl("PATHWAYS", plan$collections, fixed = TRUE)))
  expect_true(any(grepl("HALLMARKS", plan$collections, fixed = TRUE)))
  expect_false("single_de_leading_edge_gene_heatmaps" %in% plan$stage)
  expect_false("single_de_enrichmentmap" %in% plan$stage)
  expect_false("single_de_hallmarks_lollipop" %in% plan$stage)
  expect_false("single_de_hallmarks_leading_edge_gene_heatmaps" %in% plan$stage)
  expect_false(any(grepl("macrogroup_heatmap", plan$stage, fixed = TRUE)))
})
