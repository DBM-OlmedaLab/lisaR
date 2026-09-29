resolve_simple_settings <- function(...) {
  args <- list(collections = NULL, title = NULL, workers = NULL, dry_run = NULL,
    flags = list(), resources = NULL, pipeline = NULL, report = NULL)
  args[names(list(...))] <- list(...)
  do.call(lisaR:::lisa_run_resolve_settings, args)
}

test_that("direct settings cannot be silently overridden, even by equal values", {
  cases <- list(
    list(workers = 1L, pipeline = list(workers = 8L)),
    list(workers = 1L, pipeline = list(workers = 1L)),
    list(dry_run = TRUE, pipeline = list(dry_run = FALSE)),
    list(dry_run = TRUE, pipeline = list(dry_run = TRUE)),
    list(title = "Study", pipeline = list(report_title = "Other")),
    list(resources = "core", pipeline = list(dictionary_resource = "custom@1")),
    list(resources = list(dictionary = "custom@1"), pipeline = list(lisa_dictionary = "core")),
    list(resources = list(term2gene = "custom@1"), pipeline = list(term2gene_resource = "custom@1")),
    list(resources = list(category_map = "custom@1"), pipeline = list(category_map_resource = "custom@2")),
    list(flags = list(svg = FALSE), report = list(formats = list(svg = TRUE))))
  for (args in cases) expect_error(do.call(resolve_simple_settings, args), "LISA-RUN-ARGS-004")
  expect_error(resolve_simple_settings(pipeline = list(workers = 1, workers = 2)),
    "LISA-RUN-ARGS-017")
  expect_identical(resolve_simple_settings(dry_run = TRUE)$pipeline$dry_run, TRUE)
  expect_identical(resolve_simple_settings(pipeline = list(workers = 2L))$pipeline$workers, 2L)
  settings <- resolve_simple_settings(resources = "core", pipeline = list(workers = 2L))
  expect_identical(settings$pipeline$lisa_dictionary, "core")
})

test_that("a conflicting dry run is rejected before any project is created", {
  target <- tempfile("conflicting-run-")
  expect_error(run_lisa(de = "nonexistent.tsv", output_dir = target,
    species = "Homo sapiens", dry_run = TRUE, pipeline = list(dry_run = FALSE)),
    "LISA-RUN-ARGS-004")
  expect_false(dir.exists(target))
})
