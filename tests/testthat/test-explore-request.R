test_that("a figure request names exactly one figure", {
  request <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOBP-C2",
    category_id = "SYN_SIGNAL", product = "volcano")
  expect_s3_class(request, "lisa_figure_request")
  expect_identical(request$product, "volcano")
  expect_true(is.na(request$contrast_id))
  expect_match(lisa_request_id(request), "^req-[0-9a-f]{16}$")
})

test_that("a request refuses a vector in any field", {
  # This is the structural guarantee behind "no cross-product": a request that
  # could name two categories or two products cannot be constructed at all.
  expect_error(
    lisa_figure_request(unit_type = "single_de", analysis_id = "response_a",
                        collection = "GOBP-C2",
                        category_id = c("SYN_SIGNAL", "SYN_STRESS"),
                        product = "volcano"),
    "LISA-EXPLORE-002")
  expect_error(
    lisa_figure_request(unit_type = "single_de", analysis_id = "response_a",
                        collection = "GOBP-C2", category_id = "SYN_SIGNAL",
                        product = c("volcano", "heatmap")),
    "LISA-EXPLORE-002")
})

test_that("a request requires the owner of its own unit type and no other", {
  expect_error(
    lisa_figure_request(unit_type = "single_de", collection = "GOBP-C2",
                        category_id = "SYN_SIGNAL", product = "volcano"),
    "LISA-EXPLORE-001")
  expect_error(
    lisa_figure_request(unit_type = "single_de", analysis_id = "response_a",
                        contrast_id = "c1", collection = "GOBP-C2",
                        category_id = "SYN_SIGNAL", product = "volcano"),
    "LISA-EXPLORE-006")
  expect_error(
    lisa_figure_request(unit_type = "elsewhere", analysis_id = "a",
                        collection = "c", category_id = "k", product = "volcano"),
    "LISA-EXPLORE-004")
})

test_that("request identity distinguishes every scientific field", {
  base <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOBP-C2",
    category_id = "SYN_SIGNAL", product = "volcano")
  other_category <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOBP-C2",
    category_id = "SYN_STRESS", product = "volcano")
  other_collection <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOCC",
    category_id = "SYN_SIGNAL", product = "volcano")
  ids <- c(lisa_request_id(base), lisa_request_id(other_category),
           lisa_request_id(other_collection))
  expect_identical(length(unique(ids)), 3L)
})

test_that("one request translates into a scalar-only selection", {
  request <- lisa_figure_request(
    unit_type = "single_de", analysis_id = "response_a", collection = "GOBP-C2",
    category_id = "SYN_SIGNAL", product = "volcano")
  selection <- lisaR:::lisa_explore_selection(request)
  filters <- selection$selection
  expect_identical(filters$analyses, "response_a")
  expect_null(filters$contrasts)
  expect_true(all(lengths(filters[c("analyses", "collections", "categories",
                                    "products")]) == 1L))
  # PNG and PDF are two formats of one figure, delivered with its source table
  # and its recipe.
  expect_true(selection$report$formats$png)
  expect_true(selection$report$formats$pdf)
  expect_true(selection$report$source_data)
  expect_true(selection$report$recipes)
})
