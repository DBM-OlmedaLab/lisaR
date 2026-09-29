test_that("production single-DE caller forwards category evidence context", {
  e <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), e)
  calls <- list()
  visit <- function(x) {
    if (missing(x)) return(invisible(NULL))
    if (!is.call(x) && !is.expression(x)) return(invisible(NULL))
    if (is.call(x) && identical(x[[1L]], as.name("report_lisa_categories_section")))
      calls[[length(calls) + 1L]] <<- x
    for (child in as.list(x)) visit(child)
  }
  visit(body(e$main))
  expect_length(calls, 2L)
  for (call in calls)
    expect_identical(tail(as.list(call), 1L)[[1L]], as.name("evidence_index"))
})

test_that("PRE and ON LISA downloads cannot select the member-set NES figure", {
  e <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), e)
  root <- tempfile("download-identity-")
  on.exit(unlink(root, recursive = TRUE))
  page <- file.path(root, "report_pages", "single_de.html")
  dir.create(dirname(page), recursive = TRUE)
  e$report_image_formats <- function() c("png", "svg", "pdf")
  e$short_file_copy <- e$short_media_copy <- function(path, ...) basename(path)
  e$figure_source_contract <- function(path, ...) list(
    source_tsv = paste0(tools::file_path_sans_ext(basename(path)), "_source.tsv"),
    recipe_r = paste0(tools::file_path_sans_ext(basename(path)), "_recipe.R"))
  e$report_inline_category_navigation <- function(...) ""
  for (owner in c("PRE", "ON")) {
    collection <- file.path(root, "outputs", "single_de", owner, "collection_GOCC")
    dir.create(file.path(collection, "plots"), recursive = TRUE)
    dir.create(file.path(collection, "lisa_tables"))
    utils::write.table(data.frame(n_sets_evaluable = 2L),
      file.path(collection, "lisa_tables", paste0(owner, "_LISA_category_inference.tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE)
    for (kind in c("GSEA_lollipop_hommel_support", "LISA_category_inference"))
      for (format in c("png", "svg", "pdf"))
        file.create(file.path(collection, "plots", paste0(owner, "_", kind, ".", format)))
    evidence <- file.path(root, "report_pages", "evidence", owner, "index.html")
    dir.create(dirname(evidence), recursive = TRUE); file.create(evidence)
    section <- e$report_lisa_categories_section(owner, collection, character(), page,
      evidence_index = evidence)
    expect_false(grepl("LISA_category_inference|category-results-downloads", section))
    for (format in c("png", "svg", "pdf"))
      expect_match(section, paste0('download href="', owner, '_GSEA_lollipop_hommel_support.', format, '"'), fixed = TRUE)
    for (suffix in c("_source.tsv", "_recipe.R"))
      expect_match(section, paste0(owner, "_GSEA_lollipop_hommel_support", suffix), fixed = TRUE)
    cover <- e$report_category_evidence_section(owner, file.path(root, "nav"), evidence,
      page, collection_dir = collection)
    expect_match(cover, "Member-set NES distribution (not the LISA category plot)", fixed = TRUE)
    expect_match(cover, 'data-figure-kind="member-set-nes-distribution"', fixed = TRUE)
    expect_match(cover, 'data-category-results-tsv', fixed = TRUE)
    expect_match(cover, paste0(owner, "_LISA_category_inference.svg"), fixed = TRUE)
    expect_false(grepl("Saved category figure", paste0(section, cover), fixed = TRUE))
  }
})
