test_that("contrast evidence executes its inline selected-category behavior", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node is optional for offline viewer verification")
  assets <- system.file("contrast-evidence", package = "lisaR")
  skip_if(!nzchar(assets), "Installed contrast assets are required")
  probe <- test_path("fixtures", "contrast-navigator-presentation.js")
  result <- suppressWarnings(system2(node, c(shQuote(probe),
    shQuote(file.path(assets, "viewer.js")), shQuote(file.path(assets, "viewer.html"))),
    stdout = TRUE, stderr = TRUE))
  status <- attr(result, "status")
  expect_true(is.null(status) || identical(status, 0L), info = paste(result, collapse = "\n"))
  expect_match(paste(result, collapse = "\n"), "PASS contrast viewer behavior", fixed = TRUE)
})
