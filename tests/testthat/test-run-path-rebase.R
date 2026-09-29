test_that("transaction text paths are rebased from staging to final root", {
  root <- tempfile("lisa-rebase-")
  dir.create(file.path(root, "nested"), recursive = TRUE)
  from <- file.path(root, "run.staging-example")
  to <- file.path(root, "run")
  dir.create(from)
  path <- file.path(from, "nested.tsv")
  writeLines(
    c("path", file.path(from, "source_data", "example.tsv")),
    path
  )

  changed <- lisaR:::lisa_rebase_run_text_paths(from, from, to)

  to_portable <- normalizePath(to, winslash = "/", mustWork = FALSE)
  expect_identical(
    changed, normalizePath(path, winslash = "/", mustWork = TRUE)
  )
  expect_match(readLines(path)[[2]], to_portable, fixed = TRUE)
  expect_false(grepl(from, readLines(path)[[2]], fixed = TRUE))
})

test_that("transaction path rebasing accepts a non-canonical staging spelling", {
  root <- tempfile("lisa-rebase-dotdot-")
  dir.create(file.path(root, "config"), recursive = TRUE)
  dir.create(file.path(root, "results"), recursive = TRUE)
  staging <- file.path(root, "results", "run.staging-example")
  staging_spelling <- file.path(root, "config", "..", "results", "run.staging-example")
  final <- file.path(root, "results", "run")
  dir.create(staging)
  path <- file.path(staging, "index.tsv")
  writeLines(c("path", file.path(staging_spelling, "outputs", "table.tsv")), path)

  lisaR:::lisa_rebase_run_text_paths(staging, staging_spelling, final)

  observed <- readLines(path)[[2]]
  expect_match(
    observed, normalizePath(final, winslash = "/", mustWork = FALSE),
    fixed = TRUE
  )
  expect_false(grepl(".staging-", observed, fixed = TRUE))
})
