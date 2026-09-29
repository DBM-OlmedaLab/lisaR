# Static contract checks over the Shiny shell's own source.
#
# testthat runs with tests/testthat as the working directory, so these files are
# addressed relative to the package root rather than to the cwd. When the suite
# runs against an installed package without its sources beside it, the checks
# skip instead of failing on a path that legitimately is not there. The JS asset
# IS installed, so it is read from the installed package when present.

explore_pkg_source <- function(...) {
  candidates <- c(testthat::test_path("..", "..", ...), file.path(...))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) {
    testthat::skip(paste("package source not available:",
                         paste(..., sep = "/")))
  }
  paste(readLines(hit[[1L]], warn = FALSE), collapse = "\n")
}

test_that("the Shiny entry point is optional and keeps generation explicit", {
  app <- explore_pkg_source("R", "explore_app.R")
  expect_match(app, "requireNamespace\\(\"shiny\"", perl = TRUE)
  expect_match(app, "lisa_explore_shiny_request(payload)", fixed = TRUE)
  expect_match(app, "lisa_explore_submit(ws, request, background = TRUE)", fixed = TRUE)
  expect_match(app, "background = TRUE", fixed = TRUE)
  expect_match(app, "lisa_explore_submit_export(ws, destination", fixed = TRUE)
  expect_match(app, "explore_export_cancel", fixed = TRUE)
  expect_match(app, "lisa-explore-report", fixed = TRUE)
  # The parent shell names its installed stylesheet so the iframe injector can
  # reuse that exact local asset URL; an iframe is a separate document and does
  # not inherit parent CSS.
  expect_match(app, "id = \"lisa-explore-stylesheet\"", fixed = TRUE)
  # The shell must not invent its own route. It delegates to the engine's
  # attachment rule, which resolves an existing page of the native report.
  # (The delivered test looked for the literal "category_navigation" in this
  # file; it is not there and never was, because the route is built in
  # explore_export.R. Asserting the delegation is what the check meant.)
  expect_match(app, "lisa_explore_attachment(", fixed = TRUE)
  # that rule resolves the category evidence sheet -- the page the reader
  # opens for one category -- and carries the category as a query route, so a
  # block can be scoped to the open category instead of stacked on a collection.
  resolved <- lisaR:::lisa_explore_attachment(lisa_figure_request(
    unit_type = "single_de", analysis_id = "a", collection = "H",
    category_id = "C", product = "volcano"))
  expect_identical(resolved$route, "report_pages/evidence/a/H/index.html")
  expect_identical(resolved$query, "category=C")
  expect_identical(resolved$anchor, "category-C")
  # Both observers that can raise are wrapped, so an engine error reports itself
  # instead of closing the Shiny session.
  expect_match(app, "lisa_explore_submit_export", fixed = TRUE)
  expect_match(app, "lisa_explore_submit(ws, request, background = TRUE)", fixed = TRUE)
  expect_identical(
    length(gregexpr("tryCatch", app, fixed = TRUE)[[1L]]) >= 3L, TRUE)
  expect_false(grepl("run_lisa_de|run_lisa_contrast|complete-missing", app))
})

test_that("UI assets bind controls to an explicit Shiny event", {
  installed <- system.file("shiny_assets", "explore.js", package = "lisaR")
  js <- if (nzchar(installed) && file.exists(installed)) {
    paste(readLines(installed, warn = FALSE), collapse = "\n")
  } else {
    explore_pkg_source("inst", "shiny_assets", "explore.js")
  }
  expect_match(js, "Generate figure", fixed = TRUE)
  expect_match(js, "kegg_pathway_map: 'KEGG pathway map'", fixed = TRUE)
  expect_match(js, "de_recurrent_genes: 'Recurrent genes'", fixed = TRUE)
  expect_match(js, "contrast_paired_heatmap: 'Paired gene heatmaps'", fixed = TRUE)
  expect_match(js, "contrast_gene_category_network: 'Gene-category network'",
               fixed = TRUE)
  expect_match(js, "contrast_kegg_map: 'KEGG maps / painted pathways'",
               fixed = TRUE)
  expect_match(js, "row.label", fixed = TRUE)
  expect_match(js, "Native KEGG pathway (title and ID)", fixed = TRUE)
  expect_match(js, "row.entity_title", fixed = TRUE)
  expect_match(js, "groupedRows(shown)", fixed = TRUE)
  expect_match(js, "setInputValue\\('explore_generate'", perl = TRUE)
  # It reads the report's own routes rather than inventing any.
  expect_match(js, "report_pages", fixed = TRUE)
  expect_match(js, "container.querySelector('section.panel')", fixed = TRUE)
  expect_match(js, "function containerFor(doc, main, row)", fixed = TRUE)
  # A collection-wide block targets the exact owner and collection DOM context,
  # never the first panel on the page.
  expect_match(js, "data-lisa-nav-context=", fixed = TRUE)
  expect_match(js, "data-lisa-nav-collection=", fixed = TRUE)
  expect_match(js, "if (signature === applied) { return; }", fixed = TRUE)
  # a refused submission is reported, and never reported as success.
  expect_match(js, "lisa-explore-error", fixed = TRUE)
  # The native report remains untouched. The child document gets the same
  # installed stylesheet exactly once after each iframe navigation.
  expect_match(js, "function installStyles(doc)", fixed = TRUE)
  expect_match(js, "lisa-explore-iframe-stylesheet", fixed = TRUE)
  expect_match(js, "parentStyle.href", fixed = TRUE)
  expect_match(js, "installStyles(doc);", fixed = TRUE)
  expect_false(grepl("Generation requested", js, fixed = TRUE))
  expect_false(grepl("fetch\\(|XMLHttpRequest", js))
})

test_that("the static bundle script reveals only the open category", {
  installed <- system.file("shiny_assets", "explore-static.js", package = "lisaR")
  js <- if (nzchar(installed) && file.exists(installed)) {
    paste(readLines(installed, warn = FALSE), collapse = "\n")
  } else {
    explore_pkg_source("inst", "shiny_assets", "explore-static.js")
  }
  # It reads the page and toggles visibility; it must never reach the network or
  # generate anything, because it runs in an exported bundle with R closed.
  expect_match(js, "data-lisa-explore-category", fixed = TRUE)
  expect_match(js, "MutationObserver", fixed = TRUE)
  expect_false(grepl("fetch\\(|XMLHttpRequest|document.write", js))
})

test_that("the Shiny entry point is declared and optional in the package", {
  # The integration the UI delivery asked the parent to make: shiny stays a
  # Suggests, and the entry point is exported.
  expect_true("explore_lisa_run" %in% getNamespaceExports("lisaR"))
  description <- read.dcf(system.file("DESCRIPTION", package = "lisaR"))
  suggests <- if ("Suggests" %in% colnames(description)) description[1, "Suggests"] else ""
  imports <- if ("Imports" %in% colnames(description)) description[1, "Imports"] else ""
  expect_match(suggests, "shiny")
  expect_false(grepl("\\bshiny\\b", imports))
  # The assets the shell serves are actually installed.
  expect_true(nzchar(system.file("shiny_assets", "explore.css", package = "lisaR")))
  expect_true(nzchar(system.file("shiny_assets", "explore.js", package = "lisaR")))
})
