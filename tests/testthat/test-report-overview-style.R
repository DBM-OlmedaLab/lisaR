overview_style_env <- function() {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), env)
  env
}

test_that("overview exposes results and Analysis map before collapsed detail", {
  env <- overview_style_env()
  root <- tempfile("overview-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  de <- data.frame(analysis_id = c("A", "B", "C"),
    label = c("Before: treatment & control", "After <24h>", "Third"))
  co <- data.frame(contrast_id = c("c1", "c2"), output_id = c("one", "two"),
    label = c("First difference", "Second difference"))
  map <- '<section id="analysis-map"><h2>Analysis map</h2></section>'
  details <- '<section id="top-category-signals"><table><tr><td>0.05</td></tr></table></section>'
  for (mode in c("standard", "full")) {
    html <- env$report_overview_entry("Trial <A&B>", mode, de, co, root, map, details)
    expect_match(html, "Trial &lt;A&amp;B&gt;", fixed = TRUE)
    expect_match(html, "<strong>3</strong> DE analyses", fixed = TRUE)
    expect_match(html, "<strong>2</strong> contrasts", fixed = TRUE)
    expect_equal(length(regmatches(html, gregexpr('class="style-result-card ', html, fixed = TRUE))[[1]]), 5L)
    expect_match(html, "treatment &amp; control", fixed = TRUE)
    expect_match(html, "After &lt;24h&gt;", fixed = TRUE)
    expect_match(html, 'href="report_pages/single_de.html#', fixed = TRUE)
    expect_match(html, 'href="report_pages/contrasts.html#', fixed = TRUE)
    expect_match(html, '<details class="style-report-details" id="style-report-details">', fixed = TRUE)
    expect_lt(regexpr(map, html, fixed = TRUE)[[1]], regexpr('<details ', html, fixed = TRUE)[[1]])
    expect_gt(regexpr(details, html, fixed = TRUE)[[1]], regexpr('<summary>Report details</summary>', html, fixed = TRUE)[[1]])
    expect_match(html, paste0('style-mode-label">', toupper(mode)), fixed = TRUE)
  }
  html <- env$report_overview_entry("Single study", "standard", de[1, ], co[FALSE, ], root, map, details)
  expect_match(html, "<strong>1</strong> DE analysis", fixed = TRUE)
  expect_match(html, "<strong>0</strong> contrasts", fixed = TRUE)
  expect_false(grepl('style-result-card--contrast', html, fixed = TRUE))
})

test_that("theme is offline, idempotent and does not rewrite plot or script content", {
  env <- overview_style_env()
  root <- tempfile("theme root "); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  source <- '<!doctype html><html><head><style>svg{color:red}</style></head><body class="viewer"><svg id="plot"><path fill="#123456"/></svg><script>const value = 0.0123;</script><a download href="figure.svg">SVG</a></body></html>'
  cases <- c("report_index.html" = "overview", "report_pages/single_de.html" = "analysis",
    "report_pages/contrasts.html" = "contrast", "report_pages/evidence/A/GOCC/index.html" = "analysis",
    "report_pages/contrast_evidence/C/GOCC/index.html" = "contrast", "report_pages/gene_evidence/index.html" = "neutral")
  for (path in names(cases)) {
    page <- file.path(root, path)
    html <- env$report_apply_theme(source, page, root)
    expect_match(html, paste0('data-style-scope="', cases[[path]], '"'), fixed = TRUE)
    expect_match(html, '<svg id="plot"><path fill="#123456"/></svg>', fixed = TRUE)
    expect_match(html, '<script>const value = 0.0123;</script>', fixed = TRUE)
    expect_match(html, '<a download href="figure.svg">SVG</a>', fixed = TRUE)
    expect_identical(env$report_apply_theme(html, page, root), html)
    expect_false(grepl(root, html, fixed = TRUE))
    if (grepl('/GOCC/', path, fixed = TRUE))
      expect_match(html, '../../../../report_assets/lisa_pastel.css', fixed = TRUE)
  }
  env$copy_assets(system.file(package = "lisaR"), root)
  expect_true(file.exists(file.path(root, "report_assets", "lisa_pastel.css")))
  expect_identical(tools::md5sum(file.path(root, "report_assets", "lisa_pastel.css"))[[1]],
    tools::md5sum(system.file("report_assets", "lisa_pastel.css", package = "lisaR"))[[1]])
})
