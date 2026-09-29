presentation_routes_fixture <- function(report = TRUE, genes = TRUE,
    contrasts = TRUE, contract = c("expected", "plan", "none")) {
  contract <- match.arg(contract)
  root <- tempfile("presentation route lifecycle ")
  page <- file.path(root, "report_pages", "evidence", "001", "GOCC")
  dir.create(page, recursive = TRUE)
  dir.create(file.path(root, "config"))
  lisaR:::write_lisa_tsv(data.frame(analysis_id = "001", label = "Treatment versus baseline"),
    file.path(root, "config", "de_index.tsv"))
  if (contrasts) lisaR:::write_lisa_tsv(data.frame(contrast_id = "AB", output_id = "profile",
    contrast_label = "Responders versus nonresponders"), file.path(root, "config", "contrast_index.tsv"))
  stages <- c("single_de_lisa", if (report) "root_html_report",
    if (genes) "single_de_gene_evidence", if (contrasts) "category_contrasts")
  if (contract == "plan") lisaR:::write_lisa_tsv(data.frame(stage = stages),
    file.path(root, "lisa_pipeline_plan.tsv"))
  if (contract == "expected") {
    expected <- data.frame(stage = stages, expectation = "required")
    if (!report) expected <- rbind(expected,
      data.frame(stage = "root_html_report", expectation = "optional"))
    lisaR:::write_lisa_tsv(expected, file.path(root, "expected_artifacts.tsv"))
  }
  list(root = root, page = page)
}

test_that("planned gene and report routes do not depend on builder order", {
  for (contract in c("expected", "plan")) {
    f <- presentation_routes_fixture(contract = contract)
    before <- lisaR:::lisa_presentation_routes(f$root, f$page)
    expect_identical(names(before), c("overview", "analyses", "contrasts", "genes", "methods"))
    expect_identical(before$genes, "../../../gene_evidence/index.html")
    expect_identical(before$overview, "../../../../report_index.html")
    expect_identical(before$contrasts, "../../../contrasts.html")
    dir.create(file.path(f$root, "report_pages", "gene_evidence"))
    after_dir <- lisaR:::lisa_presentation_routes(f$root, f$page)
    writeLines("<html>Gene search</html>", file.path(f$root, "report_pages", "gene_evidence", "index.html"))
    after_html <- lisaR:::lisa_presentation_routes(f$root, f$page)
    expect_identical(before, after_dir)
    expect_identical(before, after_html)
    unlink(f$root, recursive = TRUE)
  }
})

test_that("evidence without a requested report never promises report pages", {
  for (contract in c("expected", "plan")) {
    f <- presentation_routes_fixture(report = FALSE, contract = contract)
    expect_identical(names(lisaR:::lisa_presentation_routes(f$root, f$page)), "genes")
    # Configuration and even unrelated stale report HTML do not override the
    # explicit current execution contract.
    writeLines("<html>Old report</html>", file.path(f$root, "report_index.html"))
    expect_identical(names(lisaR:::lisa_presentation_routes(f$root, f$page)), "genes")
    unlink(f$root, recursive = TRUE)
  }
  f <- presentation_routes_fixture(report = FALSE, genes = FALSE)
  on.exit(unlink(f$root, recursive = TRUE), add = TRUE)
  expect_identical(lisaR:::lisa_presentation_routes(f$root, f$page), list())
})

test_that("required artifact policy takes precedence over an older stage plan", {
  f <- presentation_routes_fixture(report = FALSE)
  on.exit(unlink(f$root, recursive = TRUE), add = TRUE)
  lisaR:::write_lisa_tsv(data.frame(stage = c("root_html_report", "single_de_gene_evidence")),
    file.path(f$root, "lisa_pipeline_plan.tsv"))
  expect_identical(names(lisaR:::lisa_presentation_routes(f$root, f$page)), "genes")
})

test_that("standalone viewers use actual files, not guessed destination directories", {
  f <- presentation_routes_fixture(contract = "none")
  on.exit(unlink(f$root, recursive = TRUE), add = TRUE)
  expect_identical(lisaR:::lisa_presentation_routes(NULL, f$page), list())
  expect_identical(lisaR:::lisa_presentation_routes(f$root, f$page), list())
  dir.create(file.path(f$root, "report_pages", "gene_evidence"))
  expect_identical(lisaR:::lisa_presentation_routes(f$root, f$page), list())
  writeLines("<html>Gene search</html>", file.path(f$root, "report_pages", "gene_evidence", "index.html"))
  expect_identical(lisaR:::lisa_presentation_routes(f$root, f$page),
    list(genes = "../../../gene_evidence/index.html"))
  writeLines("<html>Report</html>", file.path(f$root, "report_index.html"))
  expect_identical(names(lisaR:::lisa_presentation_routes(f$root, f$page)), c("overview", "genes"))
})

test_that("manual report assembly explicitly advertises the pages it will build", {
  f <- presentation_routes_fixture(contract = "none", contrasts = FALSE)
  on.exit(unlink(f$root, recursive = TRUE), add = TRUE)
  # The canonical empty contrast index is not a configured contrast.
  lisaR:::write_lisa_tsv(data.frame(), file.path(f$root, "config", "contrast_index.tsv"))
  routes <- lisaR:::lisa_presentation_routes(f$root, f$page, report_build = TRUE)
  expect_identical(names(routes), c("overview", "analyses", "methods"))
  lisaR:::write_lisa_tsv(data.frame(contrast_id = "AB", output_id = "profile"),
    file.path(f$root, "config", "contrast_index.tsv"))
  expect_identical(names(lisaR:::lisa_presentation_routes(f$root, f$page, report_build = TRUE)),
    c("overview", "analyses", "contrasts", "methods"))
})

test_that("human contrast labels use exact contrast IDs rather than output names", {
  f <- presentation_routes_fixture()
  on.exit(unlink(f$root, recursive = TRUE), add = TRUE)
  index <- data.frame(contrast_id = c("AB", "001"), output_id = c("profile", "AB"),
    contrast_label = c("Responders versus nonresponders", "Second contrast"))
  lisaR:::write_lisa_tsv(index, file.path(f$root, "config", "contrast_index.tsv"))
  expect_identical(lisaR:::lisa_presentation_label(f$root, "AB", contrast = TRUE),
    "Responders versus nonresponders")
  expect_identical(lisaR:::lisa_presentation_label(f$root, "001", contrast = TRUE), "Second contrast")
  expect_identical(lisaR:::lisa_presentation_label(f$root, "profile", contrast = TRUE), "profile")
  expect_identical(lisaR:::lisa_presentation_label(f$root, "001"), "Treatment versus baseline")
  html <- lisaR:::lisa_present_evidence_html("<html><head></head><body><main>Evidence</main></body></html>",
    f$page, metadata = list(contrast_id = "AB", analysis_a = "001", analysis_b = "002",
      collection = "GOCC", tier = "core", positive_contrast = "A minus B (descriptive)"),
    active = "contrasts")
  expect_match(html, "Responders versus nonresponders", fixed = TRUE)
  expect_match(html, "contrast_id", fixed = TRUE)
  expect_match(html, "analysis_a", fixed = TRUE)
  expect_match(html, "analysis_b", fixed = TRUE)
  expect_match(html, '<code>001</code>', fixed = TRUE)
  expect_match(html, 'data-lisa-route="contrasts" href="../../../contrasts.html" aria-current="page"', fixed = TRUE)
})
