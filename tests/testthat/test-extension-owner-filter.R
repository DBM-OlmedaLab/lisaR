test_that("selected analysis and contrast owners form one explicit union", {
  catalog <- data.frame(
    unit_type = c("single_de", "single_de", "single_de", "contrast", "contrast"),
    analysis_id = c("A", "B", "C", "", ""),
    contrast_id = c("", "", "", "A_vs_B", "B_vs_C"),
    collection = "GOBP-C2", category_id = "CAT_A",
    supercategory_id = "SUPER_1",
    product = c("heatmap", "heatmap", "heatmap", "contrast_heatmap", "contrast_heatmap"),
    stringsAsFactors = FALSE
  )
  selection <- list(
    mode = "selected",
    filters = list(
      analyses = c("A", "B"), contrasts = "A_vs_B",
      collections = "GOBP-C2",
      products = c("heatmap", "contrast_heatmap")
    ),
    report = list()
  )
  selected <- lisaR:::lisa_extension_filter(catalog, selection)
  expect_identical(selected$unit_type, c("single_de", "single_de", "contrast"))
  expect_identical(selected$analysis_id, c("A", "B", ""))
  expect_identical(selected$contrast_id, c("", "", "A_vs_B"))

  selection$filters$contrasts <- "missing"
  expect_error(
    lisaR:::lisa_extension_filter(catalog, selection),
    "unknown contrasts selection.*Nearest valid ID.*A_vs_B"
  )
})
