test_that("standalone navigation uses exact persisted contexts and portable links", {
  root <- tempfile("report nav "); dir.create(root)
  dir.create(file.path(root, "report_pages"))
  dir.create(file.path(root, "report_pages", "evidence", "001", "GOMF"), recursive = TRUE)
  page <- file.path(root, "report_pages", "evidence", "001", "GOMF")
  inventory <- list(version = 1L, contexts = list(list(id = "analysis:001", scientific_id = "001",
    label = "A < B", kind = "analysis", collections = list(list(id = "GOMF", label = "GOMF",
      sections = list(list(id = "actual-id", kind = "category-evidence", label = "Evidence",
        href = "report_pages/single_de.html#actual-id")))))))
  jsonlite::write_json(inventory, file.path(root, "report_pages", "navigation_inventory.json"), auto_unbox = TRUE)
  html <- '<html><body><main id="main">kept</main></body></html>'
  got <- lisaR:::lisa_attach_evidence_navigation(html, root, page,
    list(analysis_id = "001", collection = "GOMF"))
  expect_match(got, '"context":"analysis:001"', fixed = TRUE)
  expect_match(got, '../../../single_de.html#actual-id', fixed = TRUE)
  expect_match(got, 'A \\u003c B', fixed = TRUE)
  again <- lisaR:::lisa_attach_evidence_navigation(got, root, page,
    list(analysis_id = "001", collection = "GOMF"))
  expect_identical(again, got)
  expect_match(got, '<main id="main">kept</main>', fixed = TRUE)
  inventory$contexts[[1]]$collections[[1]]$sections[[1]]$href <- '../outside.html'
  jsonlite::write_json(inventory, file.path(root, "report_pages", "navigation_inventory.json"), auto_unbox = TRUE)
  expect_error(lisaR:::lisa_attach_evidence_navigation(html, root, page), "contained relative")
})
