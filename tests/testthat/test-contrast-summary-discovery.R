test_that("contrast summary discovery accepts a custom file label prefix", {
  root <- tempfile("lisa-contrast-discovery-")
  collection_dir <- file.path(root, "collection_PATHWAYS", "lisa_tables")
  dir.create(collection_dir, recursive = TRUE)
  expected <- file.path(
    collection_dir,
    "responders_on_vs_pre_RIAZ_GSE91061_GSEA_category_summary.tsv"
  )
  writeLines("category_id\tmean_NES", expected)

  observed <- lisaR:::find_lisa_category_summary(
    root,
    universe = "PATHWAYS",
    analysis = "GSEA"
  )

  expect_identical(observed, normalizePath(expected, mustWork = TRUE))
})

test_that("contrast summary discovery selects the requested collection", {
  root <- tempfile("lisa-contrast-collections-")
  for (collection in c("PATHWAYS", "HALLMARKS")) {
    collection_dir <- file.path(root, paste0("collection_", collection), "lisa_tables")
    dir.create(collection_dir, recursive = TRUE)
    writeLines(
      "category_id\tmean_NES",
      file.path(
        collection_dir,
        paste0("analysis_CUSTOM_", collection, "_GSEA_category_summary.tsv")
      )
    )
  }

  observed <- lisaR:::find_lisa_category_summary(
    root,
    universe = "HALLMARKS",
    analysis = "GSEA"
  )

  expect_match(observed, "collection_HALLMARKS", fixed = TRUE)
})
