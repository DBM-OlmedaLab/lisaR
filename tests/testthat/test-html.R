test_that("HTML index is generated", {
  output_file <- file.path(tempdir(), "lisaR-index.html")
  rows <- data.frame(
    label = c("Example A", "Example contrast"),
    href = c("example_A.html", "C01_example_A_vs_B.html"),
    stringsAsFactors = FALSE
  )
  path <- build_lisa_html_index(rows, output_file)
  expect_true(file.exists(path))
  html <- readLines(path, warn = FALSE)
  expect_true(any(grepl("Example A", html, fixed = TRUE)))
})
