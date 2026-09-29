navigation_fixture <- function() {
  data.frame(category_id = c("001", "CAT&/positive", "CAT_NEG", "CAT_NONE", "OTHER_UNCLASSIFIED"),
    category_display_name = c("Zero mean with mixed support", "Positive <label>", "Negative mean", "No support", "Residual"),
    mean_NES = c(0, 2.25, -1.8, NA, 3),
    n_genesets_significant = c(2L, 3L, 2L, 0L, 1L),
    n_pos_genesets = c(1L, 3L, 0L, 0L, 1L), n_neg_genesets = c(1L, 0L, 2L, 0L, 0L),
    macrogroup_name = "A group", gsea_padj_cutoff = .25, analysis = "GSEA", stringsAsFactors = FALSE)
}

test_that("category navigation retains source metric, order, scope, and absence", {
  x <- navigation_fixture()
  nav <- lisaR:::build_lisa_category_navigation(x, "A_1", "GOBP-C2", "core")
  expect_identical(nav$categories$category_id, x$category_id[1:4])
  expect_identical(nav$categories$mean_NES, x$mean_NES[1:4])
  expect_equal(nav$categories$n_genesets_significant, c(2, 3, 2, 0))
  expect_identical(nav$categories$support_state[[4L]], "no_significant_member_sets")
  expect_true(is.na(nav$categories$mean_NES[[4L]]))
  expect_identical(nav$categories$mean_NES[[1L]], 0)
  expect_match(nav$metadata$inference, "no category-level p-value or FDR", fixed = TRUE)
  links <- lisaR:::lisa_category_navigation_links(nav, "../../../evidence/A_1/GOBP-C2/index.html")
  expect_identical(links[[2L]], "../../../evidence/A_1/GOBP-C2/index.html?category=CAT%26%2Fpositive&analysis_id=A_1&collection=GOBP-C2&tier=core")
  expect_false(any(grepl("OTHER_UNCLASSIFIED", links)))
})

test_that("navigator renders ordered supercategory bands without changing category evidence", {
  x <- navigation_fixture()[1:4, ]
  x$macrogroup_id <- c("IMMUNE", "IMMUNE", "METABOLISM", "METABOLISM")
  x$macrogroup_name <- c("Immune response and inflammatory signalling", "Immune response and inflammatory signalling", "Metabolism", "Metabolism")
  x$macrogroup_order <- c(1L, 1L, 2L, 2L)
  x$category_order_within_macrogroup <- c(1L, 2L, 1L, 2L)
  x$color <- c("#1b9e77", "#1b9e77", "#7570b3", "#7570b3")
  nav <- lisaR:::build_lisa_category_navigation(x, "A_1", "GOBP-C2", "core")
  out <- tempfile("grouped-category-navigation-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  paths <- lisaR:::render_lisa_category_navigation(nav, out)
  svg <- paste(readLines(paths$svg, warn = FALSE), collapse = "\n")
  html <- paste(readLines(file.path(out, "index.html"), warn = FALSE), collapse = "\n")
  expect_match(html, "this historical navigator does not show category-level adjusted P", fixed = TRUE)
  expect_false(grepl("Stars show adjusted category P", html, fixed = TRUE))
  source <- lisaR:::lisa_navigation_read(paths$source)
  expect_equal(nrow(source), nrow(x))
  expect_identical(source$category_id, x$category_id)
  expect_identical(source$n_genesets_significant, as.character(x$n_genesets_significant))
  expect_identical(source$evidence_url, lisaR:::lisa_category_navigation_links(nav, "../../../evidence/A_1/GOBP-C2/index.html"))
  expect_true(all(c("macrogroup_id", "macrogroup_name", "macrogroup_order", "category_order_within_macrogroup", "category_color") %in% names(source)))
  expect_equal(length(gregexpr('class="supercategory-band"', svg, fixed = TRUE)[[1L]]), 2L)
  expect_match(svg, "Immune response and", fixed = TRUE)
  expect_match(svg, "inflammatory signalling", fixed = TRUE)
  expect_match(svg, 'stroke="#1b9e77"', fixed = TRUE)
  expect_match(svg, 'data-category-id="CAT&amp;/positive"', fixed = TRUE)
})

test_that("navigation rejects wrong scope, contrasts and inconsistent source summaries", {
  make <- function(x, ...) lisaR:::build_lisa_category_navigation(x, "A", "GOBP-C2", "core", ...)
  x <- navigation_fixture()
  expect_error(make(transform(x, analysis_id = "B")), "scope mismatch")
  expect_error(make(transform(x, universe = "GOCC")), "scope mismatch")
  expect_error(make(transform(x, tier = "expanded")), "scope mismatch")
  expect_error(make(transform(x, analysis = "ORA")), "GSEA summaries only")
  expect_error(make(transform(x, mean_NES_A = mean_NES)), "individual analyses")
  expect_error(make(x, gsea_padj_cutoff = .05), "cutoff differs")
  expect_error(make(rbind(x, x[1L, ])), "Duplicate category IDs")
  y <- x; y$mean_NES[[4L]] <- 0
  expect_error(make(y), "not zero")
  y <- x; y$mean_NES[[1L]] <- NA_real_
  expect_error(make(y), "finite existing mean_NES")
  y <- x; y$n_genesets_significant[[1L]] <- 1.5
  expect_error(make(y), "nonnegative integers")
  y <- x; y$n_pos_genesets[[1L]] <- 3
  expect_error(make(y), "Directional set counts")
  expect_error(lisaR:::build_lisa_category_navigation(x, "../contrast", "GOBP-C2", "core"), "path component")
})

test_that("portable navigator provides accessible exact links and a reproducible SVG", {
  nav <- lisaR:::build_lisa_category_navigation(navigation_fixture(), "A_1", "GOBP-C2", "core", positive_contrast = "Treatment versus control")
  out <- tempfile("category-navigation-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  paths <- lisaR:::render_lisa_category_navigation(nav, out)
  expect_true(all(file.exists(unlist(paths))))
  html <- paste(readLines(paths$html, warn = FALSE), collapse = "\n")
  svg <- paste(readLines(paths$svg, warn = FALSE), collapse = "\n")
  expect_false(grepl("Accessible category table", html, fixed = TRUE))
  expect_match(svg, 'target="_top" tabindex="0"', fixed = TRUE)
  expect_match(svg, 'role="group"', fixed = TRUE)
  expect_match(svg, 'Positive &lt;label&gt;', fixed = TRUE)
  expect_match(svg, "Mean NES of significant member sets", fixed = TRUE)
  expect_match(svg, "Red: positive mean NES", fixed = TRUE)
  expect_match(svg, "2.00 (1.00 / 1.00)", fixed = TRUE)
  expect_match(svg, 'category=001&amp;analysis_id=A_1', fixed = TRUE)
  expect_match(svg, '#b2182b', fixed = TRUE); expect_match(svg, '#2166ac', fixed = TRUE)
  expect_match(svg, 'not available', fixed = TRUE)
  expect_false(grepl("OTHER_UNCLASSIFIED", html, fixed = TRUE))
  expect_false(grepl("fetch(", html, fixed = TRUE))
  table <- lisaR:::lisa_navigation_read(paths$source)
  expect_identical(table$category_id[[1L]], "001")
  metadata <- jsonlite::read_json(paths$metadata, simplifyVector = TRUE)
  expect_identical(metadata$source_sha256, digest::digest(file = paths$source, algo = "sha256"))
  expect_identical(lisaR:::lisa_category_navigation_svg(table, metadata$title, metadata$subtitle), svg)
  recipe <- readLines(file.path(out, "reproduce_category_navigation.R"), warn = FALSE)
  expect_true(any(grepl('colClasses="character"', recipe, fixed = TRUE)))
  expect_false(any(grepl("library(", recipe, fixed = TRUE)))
  expect_false(any(grepl(out, recipe, fixed = TRUE)))
  reproduced <- file.path(out, "reproduced.svg")
  log <- file.path(out, "recipe.log")
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla",
    shQuote(file.path(out, "reproduce_category_navigation.R")), shQuote(paths$source), shQuote(reproduced)),
    stdout = log, stderr = log)
  expect_equal(status, 0L, info = paste(readLines(log, warn = FALSE), collapse = "\n"))
  expect_identical(readBin(reproduced, "raw", file.info(reproduced)$size), readBin(paths$svg, "raw", file.info(paths$svg)$size))
  expect_error(lisaR:::render_lisa_category_navigation(nav, tempfile(), evidence_base = "javascript:alert(1)"), "relative HTML")
})

test_that("the clickable SVG passes the actual report privacy gate", {
  skip_if_not(nzchar(system.file("scripts", "build_LISA_report.R", package = "lisaR")))
  nav <- lisaR:::build_lisa_category_navigation(navigation_fixture(), "A", "GOCC", "core")
  out <- tempfile("navigation-privacy-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  paths <- lisaR:::render_lisa_category_navigation(nav, out)
  report <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = report)
  expect_no_error(report$report_assert_no_private_paths(out, paths$html))
  nav$metadata$positive_contrast <- "/srv/private-study/design.tsv"
  unsafe <- lisaR:::render_lisa_category_navigation(nav, out)
  expect_error(report$report_assert_no_private_paths(out, unsafe$html), "LISA-REPORT-PRIVACY-004")
})

test_that("empty and residual-only summaries stay navigable without invented means", {
  for (x in list(navigation_fixture()[FALSE, ], navigation_fixture()[5L, ])) {
    nav <- lisaR:::build_lisa_category_navigation(x, "A", "GOCC", "core")
    expect_equal(nrow(nav$categories), 0L)
    out <- tempfile("empty-navigation-")
    paths <- lisaR:::render_lisa_category_navigation(nav, out)
    html <- paste(readLines(paths$html, warn = FALSE), collapse = "\n")
    svg <- paste(readLines(paths$svg, warn = FALSE), collapse = "\n")
    expect_match(html, "No classified categories", fixed = TRUE)
    expect_false(grepl('supercategory-band', svg, fixed = TRUE))
    expect_equal(nrow(lisaR:::lisa_navigation_read(paths$source)), 0L)
    unlink(out, recursive = TRUE)
  }
})

test_that("navigation remains usable when its optional downloads are disabled", {
  nav <- lisaR:::build_lisa_category_navigation(navigation_fixture(), "A_1", "GOBP-C2", "core")
  out <- tempfile("navigation-policy-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  paths <- lisaR:::render_lisa_category_navigation(nav, out,
    download_policy = list(svg = FALSE, source_data = FALSE, recipes = FALSE))
  html <- paste(readLines(paths$html, warn = FALSE), collapse = "\n")
  expect_true(all(file.exists(unlist(paths))))
  expect_match(html, '<div class="diagram">', fixed = TRUE)
  expect_match(html, 'Download provenance', fixed = TRUE)
  expect_false(grepl('Download SVG|Download source table|>R script<', html))
  source <- lisaR:::lisa_navigation_read(paths$source)
  expect_identical(source$category_id, navigation_fixture()$category_id[1:4])
  expect_equal(as.numeric(source$mean_NES), navigation_fixture()$mean_NES[1:4])
  expect_match(html, 'category=CAT%26%2Fpositive', fixed = TRUE)
  paths <- lisaR:::render_lisa_category_navigation(nav, out,
    download_policy = list(svg = TRUE, source_data = FALSE, recipes = TRUE))
  html <- paste(readLines(paths$html, warn = FALSE), collapse = "\n")
  expect_match(html, 'Download SVG', fixed = TRUE)
  expect_match(html, '>R script<', fixed = TRUE)
  expect_false(grepl('Download source table', html, fixed = TRUE))
})
