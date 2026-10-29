full_anchor_report_env <- function() {
  e <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), e)
  e$report_output_policy <- lisaR:::lisa_output_policy_table("svg", TRUE, TRUE)
  e$report_package_dir <- system.file(package = "lisaR")
  e$report_gallery_assets <- list()
  e
}

full_anchor_html_values <- function(html, attribute) {
  pattern <- paste0('\\b', attribute, '="[^"]*"')
  hits <- regmatches(html, gregexpr(pattern, html, perl = TRUE))[[1L]]
  sub('"$', "", sub(paste0('^', attribute, '="'), "", hits))
}

test_that("FULL overview links resolve to real rendered positive sections, never empty sections", {
  e <- full_anchor_report_env()
  root <- tempfile("full-anchor-render-")
  dir.create(file.path(root, "config"), recursive = TRUE)
  dir.create(file.path(root, "report_pages"))
  e$report_project_dir <- root
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  lisaR:::write_lisa_tsv(data.frame(analysis_id = "A"),
    file.path(root, "config", "de_index.tsv"))
  lisaR:::write_lisa_tsv(data.frame(contrast_id = "C", output_id = "paired"),
    file.path(root, "config", "contrast_index.tsv"))
  navigation <- list()
  input_files <- character()
  for (kind in c("category", "contrast")) {
    owner <- if (kind == "category") "A" else "C_paired"
    page <- if (kind == "category") "single_de.html" else "contrasts.html"
    page_file <- file.path(root, "report_pages", page)
    evidence_dir <- file.path(root, "report_pages",
      if (kind == "category") "evidence" else "contrast_evidence", owner, "GOCC")
    dir.create(file.path(evidence_dir, "tables"), recursive = TRUE)
    evidence_index <- file.path(evidence_dir, "index.html")
    writeLines('<!doctype html><html><body>Saved synthetic evidence</body></html>', evidence_index)
    metadata <- data.frame(key = c(if (kind == "category") "analysis_id" else "contrast_id",
      "collection", "tier"), value = c(owner, "GOCC", "synthetic"))
    lisaR:::write_lisa_tsv(metadata, file.path(evidence_dir, "tables", "metadata.tsv"))
    art <- file.path("artifacts", owner, "GOCC", "volcano", "CAT_POSITIVE")
    dir.create(dirname(file.path(root, art)), recursive = TRUE)
    writeLines('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"><circle cx="5" cy="5" r="2"/></svg>',
      paste0(file.path(root, art), ".svg"))
    writeLines(c("category_id\tvalue", "CAT_POSITIVE\t1"),
      paste0(file.path(root, art), "_source.tsv"))
    writeLines("plot(1)", paste0(file.path(root, art), "_recipe.R"))
    paths <- paste0(art, c(".svg", "_source.tsv", "_recipe.R"))
    input_files <- c(input_files, file.path(root, paths))
    attachment <- list(category_products = list(list(product = "volcano",
      category_id = "CAT_POSITIVE", label = "Volcano overlays",
      assets = Map(function(path, format) list(source_path = path, format = format),
        paths, c("svg", "tsv", "r")))), unclassified = character())
    title <- paste(owner, "GOCC", "FULL figures", sep = " - ")
    before <- tools::md5sum(file.path(root, paths))
    positive <- e$report_full_products_section(title, attachment, evidence_index,
      page_file, kind, owner, "GOCC", project_dir = root)
    # Render both absence paths. Neither is a positive FULL destination.
    empty_title <- paste(owner, "GOMF", "FULL figures", sep = " - ")
    empty <- e$report_full_products_section(empty_title,
      list(category_products = list(), unclassified = character()), evidence_index,
      page_file, kind, owner, "GOMF", project_dir = root)
    no_sheet_empty <- e$report_full_gallery_section(
      paste(owner, "HALLMARKS", "FULL figures", sep = " - "), list(), page_file, root)
    expect_identical(empty, "")
    expect_identical(no_sheet_empty, "")
    expect_match(positive, "Volcano overlays", fixed = TRUE)
    expect_match(positive, "CAT_POSITIVE.svg", fixed = TRUE)
    expect_identical(tools::md5sum(file.path(root, paths)), before)
    category_title <- paste(owner, "GOMF", "LISA categories", sep = " - ")
    empty_notice <- e$empty_completed_section(category_title,
      "No eligible member gene sets met the recorded GSEA FDR cutoff.", "")
    contexts <- list(list(id = owner, label = paste("Synthetic", owner),
      collections = list(
        list(id = "GOCC", label = "GOCC", sections = list(
          e$report_navigation_section(title, "FULL figures"))),
        list(id = "GOMF", label = "GOMF", sections = list(
          e$report_navigation_section(category_title, "LISA categories"))),
        list(id = "HALLMARKS", label = "HALLMARKS", sections = list()))))
    navigation[[kind]] <- contexts
    e$write_page(page_file, e$page_shell("Public synthetic mixed report",
      if (kind == "category") "single_de" else "contrasts",
      paste0(positive, empty_notice), root, page_file, navigation = contexts))
  }
  before <- tools::md5sum(input_files)
  overview <- e$report_full_overview_section(navigation$category, navigation$contrast)
  hrefs <- full_anchor_html_values(overview, "href")
  expect_length(hrefs, 2L)
  expect_true(any(startsWith(hrefs, "report_pages/single_de.html#")))
  expect_true(any(startsWith(hrefs, "report_pages/contrasts.html#")))
  expect_false(any(grepl("gomf|hallmarks", hrefs)))
  for (href in hrefs) {
    pieces <- strsplit(href, "#", fixed = TRUE)[[1L]]
    expect_length(pieces, 2L)
    target <- file.path(root, pieces[[1L]])
    expect_true(file.exists(target))
    html <- paste(readLines(target, warn = FALSE), collapse = "\n")
    expect_true(pieces[[2L]] %in% full_anchor_html_values(html, "id"))
    expect_match(html, "No eligible member gene sets met", fixed = TRUE)
  }
  expect_identical(tools::md5sum(input_files), before)
  expect_identical(e$report_full_overview_section(list(), list()), "")
})

test_that("FULL overview retains real no-evidence galleries while omitting empty galleries", {
  e <- full_anchor_report_env()
  root <- tempfile("full-gallery-anchor-")
  e$report_project_dir <- root
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  art <- file.path(root, "artifacts", "A", "HALLMARKS", "gene_cards")
  dir.create(art, recursive = TRUE)
  writeLines('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10"/>',
    file.path(art, "H_SIGNAL.svg"))
  writeLines(c("category_id", "H_SIGNAL"), file.path(art, "H_SIGNAL_source.tsv"))
  writeLines("plot(1)", file.path(art, "H_SIGNAL_recipe.R"))
  page <- file.path(root, "report_pages", "single_de.html")
  dir.create(dirname(page), recursive = TRUE)
  title <- "A - HALLMARKS - FULL figures"
  html <- e$report_full_gallery_section(title, list(GeneCards = list(path = art)), page, root)
  expect_true(nzchar(html))
  expect_identical(e$report_full_gallery_section("A - GOMF - FULL figures", list(), page, root), "")
  e$write_page(page, html)
  context <- list(list(id = "A", label = "A", collections = list(
    list(id = "HALLMARKS", label = "HALLMARKS", sections = list(
      e$report_navigation_section(title, "FULL figures"))))))
  entry <- e$report_full_overview_section(context, list())
  expect_identical(full_anchor_html_values(entry, "href"),
    "report_pages/single_de.html#a-hallmarks-full-figures")
  expect_true("a-hallmarks-full-figures" %in% full_anchor_html_values(html, "id"))
  expect_match(html, "H_SIGNAL.svg", fixed = TRUE)
})
