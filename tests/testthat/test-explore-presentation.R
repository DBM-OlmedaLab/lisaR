test_that("presentation attaches an available volcano at its existing category route", {
  route <- "report_pages/category_navigation/response_a/GOBP-C2/index.html"
  presentation <- data.frame(
    key = "fig-safe", category_id = "SYN_SIGNAL", product = "volcano",
    attach_route = route, attach_anchor = "category-SYN_SIGNAL",
    attach_query = "category=SYN_SIGNAL",
    png = "explore_extensions/fig-safe/artifacts/volcano.png",
    pdf = "explore_extensions/fig-safe/artifacts/volcano.pdf",
    source_data = "explore_extensions/fig-safe/artifacts/volcano_source.tsv",
    recipe = "explore_extensions/fig-safe/artifacts/volcano_recipe.R",
    stringsAsFactors = FALSE)
  html <- "<html><body><main id=\"lisa-main\"><h1>Native category navigation</h1></main></body></html>"
  attached <- lisaR:::lisa_explore_presentation_attach_html(html, presentation, route)

  expect_match(attached, "Native category navigation", fixed = TRUE)
  expect_match(attached, "data-lisa-explore-key=\"fig-safe\"", fixed = TRUE)
  expect_match(attached, "Volcano · SYN_SIGNAL", fixed = TRUE)
  expect_match(attached, "Download PNG", fixed = TRUE)
  expect_match(attached, "Download PDF", fixed = TRUE)
  expect_match(attached, "Source data", fixed = TRUE)
  expect_match(attached, "R script", fixed = TRUE)
  expect_match(attached,
    "../../../../explore_extensions/fig-safe/artifacts/volcano.png", fixed = TRUE)
})

test_that("presentation is route-specific, idempotent, and rejects unsafe paths", {
  route <- "report_pages/category_navigation/A/H/index.html"
  presentation <- data.frame(
    key = "fig-one", category_id = "C", product = "volcano", attach_route = route,
    attach_anchor = "category-C", attach_query = "category=C", png = "explore_extensions/fig-one/a.png", pdf = "",
    source_data = "", recipe = "", stringsAsFactors = FALSE)
  html <- "<main></main>"
  once <- lisaR:::lisa_explore_presentation_attach_html(html, presentation, route)
  twice <- lisaR:::lisa_explore_presentation_attach_html(once, presentation, route)
  expect_identical(twice, once)
  expect_identical(lisaR:::lisa_explore_presentation_attach_html(html, presentation,
    "report_pages/category_navigation/A/OTHER/index.html"), html)

  unsafe <- presentation
  unsafe$png <- "../outside.png"
  expect_error(lisaR:::lisa_explore_presentation_attach_html(html, unsafe, route),
    "LISA-EXPLORE-041")
})

test_that("attachment links remain contained and escaped", {
  route <- "report_pages/category_navigation/A/H/index.html"
  presentation <- data.frame(
    key = "fig<&\"", category_id = "C<script>", product = "volcano",
    attach_route = route, attach_anchor = "C\" onmouseover=\"bad",
    attach_query = "",
    png = "explore_extensions/fig/a.png", pdf = "", source_data = "", recipe = "",
    stringsAsFactors = FALSE)
  attached <- lisaR:::lisa_explore_presentation_attach_html("<main></main>", presentation, route)
  expect_false(grepl("<script>", attached, fixed = TRUE))
  expect_false(grepl("onmouseover=\"bad", attached, fixed = TRUE))
  expect_match(attached, "C&lt;script&gt;", fixed = TRUE)
})
