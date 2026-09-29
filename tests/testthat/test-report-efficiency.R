report_efficiency_env <- function() {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  env
}

test_that("standalone report sanitization needs no inherited package helpers", {
  # Test environments usually inherit the loaded lisaR namespace through the
  # search path. Isolate the standalone helpers from that package lookup.
  e <- new.env(parent = baseenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = e)
  expect_identical(parent.env(e), baseenv())
  # Newer R also supplies this operator in base. Trap its use so this check
  # still protects supported older R versions where it is unavailable.
  e[["%||%"]] <- function(...) stop("Standalone helpers must support older R without %||%.")
  root <- withr::local_tempdir(pattern = "portable-standalone-")
  internal <- file.path(root, "inputs", "de.tsv")
  # Exercise both the scalar fallback and vectorized pre-normalized-root path.
  expect_identical(e$report_portable_text_one(internal, root), "inputs/de.tsv")
  expect_identical(e$report_portable_text(
    c(internal, "/private/source.tsv", "file:///private/source.tsv",
      "1.2345678901234567", "https://example.org/public.tsv", NA_character_), root),
    c("inputs/de.tsv", "source.tsv", "source.tsv",
      "1.2345678901234567", "https://example.org/public.tsv", NA_character_))
  source <- file.path(root, "source.tsv")
  destination <- file.path(root, "portable.tsv.gz")
  writeLines(c("symbol\tvalue\tinput_path", paste("G1", "1.2345678901234567", internal, sep = "\t")), source)
  e$report_portable_table_copy(source, destination, root)
  portable <- e$report_read_delimited_text(destination)
  expect_identical(portable$value, "1.2345678901234567")
  expect_identical(portable$input_path, "inputs/de.tsv")
})

test_that("vectorized portable tables retain old path semantics and exact numbers", {
  e <- report_efficiency_env()
  root <- tempfile("portable-vector-"); dir.create(root)
  values <- c(NA_character_, "", "1.2345678901234567", "ENSG000001", "A/B/C",
    "https://example.org/a/b", "/private/source.tsv", "C:\\private\\source.tsv",
    paste0(root, "/inputs/de.tsv"), "see /private/source.tsv here",
    "file:///private/source.tsv", "a duplicated name", "a duplicated name")
  reference <- vapply(values, e$report_portable_text_one, character(1),
    project_dir = root, USE.NAMES = FALSE)
  expect_identical(e$report_portable_text(values, root), reference)
  expect_identical(e$report_portable_text(values, root, force_path = TRUE),
    vapply(values, e$report_portable_text_one, character(1),
      project_dir = root, force_path = TRUE, USE.NAMES = FALSE))
  numeric <- data.frame(x = c(pi, NA_real_, NaN, Inf, -Inf),
    int = c(1L, NA_integer_, 3L, 4L, 5L), flag = c(TRUE, FALSE, NA, TRUE, FALSE))
  expect_identical(e$report_portable_data_frame(numeric, root), numeric)
})

test_that("download reuse is content-scoped and never trusts a stale output", {
  e <- report_efficiency_env()
  root <- tempfile("portable-reuse-")
  dir.create(file.path(root, "report_pages"), recursive = TRUE)
  page <- file.path(root, "report_pages", "downloads.html"); file.create(page)
  a <- file.path(root, "one.tsv"); b <- file.path(root, "two.tsv")
  writeLines(c("gene\tvalue\tinput_path", "G1\t1.2345678901234567\t/private/de.tsv"), a)
  file.copy(a, b)
  calls <- 0L; original <- e$report_portable_table_copy
  e$report_portable_table_copy <- function(...) { calls <<- calls + 1L; original(...) }
  first <- e$short_file_copy(a, page)
  expect_identical(e$short_file_copy(b, page), first)
  expect_equal(calls, 1L)
  target <- file.path(dirname(page), first)
  con <- gzfile(target, "rt"); lines <- readLines(con); close(con)
  expect_true(any(grepl("1.2345678901234567", lines, fixed = TRUE)))
  expect_false(any(grepl("/private/", lines, fixed = TRUE)))
  writeLines("tampered", target)
  expect_identical(e$short_file_copy(b, page), first)
  expect_equal(calls, 2L)
  # A fresh invocation must sanitize once, even with an existing same-name file.
  fresh <- report_efficiency_env()
  writeLines("unsanitized /private/leak.tsv", target)
  fresh$short_file_copy(a, page)
  con <- gzfile(target, "rt"); restored <- readLines(con); close(con)
  expect_identical(restored, lines)
  writeLines(c("gene\tvalue", "G2\t2"), b)
  expect_false(identical(e$short_file_copy(b, page), first))
})
