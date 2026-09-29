contrast_navigation_fixture <- function() {
  data.frame(category_id = c("DEV_MORPH", "APOPTOSIS", "NO_SUPPORT", "OTHER_UNCLASSIFIED"),
    category_display_name = c("Development and morphogenesis", "Apoptosis and regulated cell death",
      "No support", "Residual"),
    macrogroup_id = c("DEV", "DEV", "DEV", "DEV"), macrogroup_name = "Developmental programs",
    color = "#E88BA1",
    mean_NES_A = c(1.2, -0.4, NA, 3), mean_NES_B = c(0.3, 0.9, NA, 1),
    n_genesets_significant_A = c(12L, 4L, 0L, 1L), n_genesets_significant_B = c(6L, 9L, 0L, 1L),
    delta_mean_NES = c(0.9, -1.3, NA, 2), universe = "GOBP-C2", stringsAsFactors = FALSE)
}

test_that("contrast navigation retains canonical A/B values without recomputing enrichment", {
  x <- contrast_navigation_fixture()
  nav <- lisaR:::build_lisa_contrast_navigation(x, contrast_id = "response_dynamics",
    scope = "response_dynamics_response_dynamics", analysis_a = "responders_on_vs_pre",
    analysis_b = "pd_on_vs_pre", collection = "GOBP-C2", tier = "core",
    label = "Treatment dynamics: responders minus PD")
  expect_identical(nav$categories$category_id, x$category_id[1:3])
  expect_identical(nav$categories$mean_NES, x$delta_mean_NES[1:3])
  expect_true(is.na(nav$categories$mean_NES[[3L]]))
  expect_equal(nav$categories$n_pos_genesets, x$n_genesets_significant_A[1:3])
  expect_equal(nav$categories$n_neg_genesets, x$n_genesets_significant_B[1:3])
  expect_true(all(is.na(nav$categories$n_genesets_significant)))
  expect_match(nav$metadata$metric, "Existing mean_NES_A and mean_NES_B", fixed = TRUE)
  expect_identical(nav$categories$mean_NES_A,x$mean_NES_A[1:3])
  expect_identical(nav$categories$mean_NES_B,x$mean_NES_B[1:3])
  expect_match(nav$metadata$inference, "no category-level p-value, FDR or new test", fixed = TRUE)
  expect_false(any(grepl("OTHER_UNCLASSIFIED", nav$categories$category_id)))
})

test_that("contrast navigation rejects a single-analysis table and duplicate/unsafe scope", {
  x <- contrast_navigation_fixture()
  make <- function(y, ...) lisaR:::build_lisa_contrast_navigation(y, contrast_id = "response_dynamics",
    scope = "response_dynamics_response_dynamics", analysis_a = "responders_on_vs_pre",
    analysis_b = "pd_on_vs_pre", collection = "GOBP-C2", tier = "core", ...)
  single_shaped <- x[, setdiff(names(x), c("mean_NES_A", "mean_NES_B", "delta_mean_NES"))]
  single_shaped$mean_NES <- c(1, 2, NA, 3)
  single_shaped$n_genesets_significant <- c(1L, 2L, 0L, 1L)
  expect_error(make(single_shaped), "category-contrast table")
  expect_error(make(rbind(x, x[1L, ])), "Duplicate category IDs")
  expect_error(make(transform(x, universe = "GOCC")), "scope mismatch")
  expect_error(lisaR:::build_lisa_contrast_navigation(x, contrast_id = "../x",
    scope = "response_dynamics_response_dynamics", analysis_a = "a", analysis_b = "b",
    collection = "GOBP-C2", tier = "core"), "path component")
})

test_that("contrast navigation links reuse the exact scope so the contrast-evidence viewer selects the intended category", {
  nav <- lisaR:::build_lisa_contrast_navigation(contrast_navigation_fixture(), contrast_id = "response_dynamics",
    scope = "response_dynamics_response_dynamics", analysis_a = "responders_on_vs_pre",
    analysis_b = "pd_on_vs_pre", collection = "GOBP-C2", tier = "core")
  links <- lisaR:::lisa_contrast_navigation_links(nav,
    "../../../contrast_evidence/response_dynamics_response_dynamics/GOBP-C2/index.html")
  expect_identical(links[[1L]], paste0(
    "../../../contrast_evidence/response_dynamics_response_dynamics/GOBP-C2/index.html",
    "?category_id=DEV_MORPH&contrast_id=response_dynamics&analysis_a=responders_on_vs_pre",
    "&analysis_b=pd_on_vs_pre&collection=GOBP-C2&tier=core"))
})

test_that("the rendered navigator shows each side, not the difference, and exact A/B counts", {
  nav <- lisaR:::build_lisa_contrast_navigation(contrast_navigation_fixture(), contrast_id = "response_dynamics",
    scope = "response_dynamics_response_dynamics", analysis_a = "responders_on_vs_pre",
    analysis_b = "pd_on_vs_pre", collection = "GOBP-C2", tier = "core",
    label = "Treatment dynamics: responders minus PD")
  out <- tempfile("contrast-navigation-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  paths <- lisaR:::render_lisa_contrast_navigation(nav, out)
  expect_true(all(file.exists(unlist(paths))))
  svg <- paste(readLines(paths$svg, warn = FALSE), collapse = "\n")
  expect_match(svg, "Significant sets (A / B)", fixed = TRUE)
  expect_match(svg, "Mean NES: A and B</text>", fixed = TRUE)
  expect_false(grepl("A minus B delta mean NES of significant member sets", svg, fixed = TRUE))
  expect_match(svg, "Blue circles: A", fixed = TRUE)
  expect_match(svg, 'data-profile="A" data-mean-nes="1.2"', fixed = TRUE)
  expect_match(svg, "A 12 / B 6", fixed = TRUE)
  expect_false(grepl("NA (12 / 6)", svg, fixed = TRUE))
  expect_match(svg, "A: unavailable", fixed = TRUE)
  expect_match(svg, 'target="_top" tabindex="0"', fixed = TRUE)
  expect_match(svg, "contrast_evidence/response_dynamics_response_dynamics/GOBP-C2/index.html", fixed = TRUE)
  expect_match(svg, "category_id=DEV_MORPH", fixed = TRUE)
  expect_false(grepl("category_id=OTHER_UNCLASSIFIED", svg, fixed = TRUE))
  html <- paste(readLines(paths$html, warn = FALSE), collapse = "\n")
  expect_match(html, "<svg", fixed = TRUE)
  expect_match(html, "Both analysis profiles are shown separately", fixed = TRUE)
  recipe_path <- file.path(out, "reproduce_contrast_navigation.R")
  expect_true(file.exists(recipe_path))
  reproduced <- tempfile("contrast-navigation-reproduced-", fileext = ".svg")
  log <- tempfile("contrast-navigation-recipe-", fileext = ".log")
  status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla",
    shQuote(recipe_path), shQuote(paths$source), shQuote(reproduced)), stdout = log, stderr = log)
  expect_identical(status, 0L, info = paste(readLines(log, warn = FALSE), collapse = "\n"))
  expect_true(file.exists(reproduced))
  if (file.exists(reproduced))
    expect_identical(readLines(reproduced, warn = FALSE), readLines(paths$svg, warn = FALSE))
})

test_that("contrast navigator preserves missing significant support rather than plotting a contextual value as zero", {
  x <- contrast_navigation_fixture()[1L, ]
  x$mean_NES_A <- NA_real_
  x$mean_NES_B <- 0.2
  x$display_mean_NES_A <- 0.7
  x$display_mean_NES_B <- 0.2
  x$delta_mean_NES <- 0.5
  x$n_genesets_significant_A <- 0L
  x$n_genesets_significant_B <- 2L
  nav <- lisaR:::build_lisa_contrast_navigation(x, contrast_id = "response_dynamics",
    scope = "response_dynamics_response_dynamics", analysis_a = "responders_on_vs_pre",
    analysis_b = "pd_on_vs_pre", collection = "GOBP-C2", tier = "core")
  out <- tempfile("contrast-contextual-navigation-"); on.exit(unlink(out, recursive = TRUE), add = TRUE)
  svg <- paste(readLines(lisaR:::render_lisa_contrast_navigation(nav, out)$svg, warn = FALSE), collapse = "\n")
  expect_match(svg, "A: unavailable", fixed = TRUE)
  expect_match(svg, "A 0 / B 2", fixed = TRUE)
  expect_false(grepl("A minus B delta mean NES unavailable", svg, fixed = TRUE))
  expect_false(grepl("delta mean NES of significant member sets", svg, fixed = TRUE))
})

test_that("the contrast navigator is registered in the code-identity and CLI-argument contracts", {
  closures <- lisaR:::lisa_post_script_dependencies()
  expect_true("build_contrast_navigation.R" %in% names(closures))
  interfaces <- lisaR:::lisa_post_script_interfaces()
  expect_true(all(c("contrast-id", "scope", "analysis-a", "analysis-b", "universe", "tier",
    "gsea-padj-cutoff", "output-dir") %in% interfaces[["build_contrast_navigation.R"]]))
})

test_that("contrast navigation preserves both profiles when optional downloads are disabled", {
  nav <- lisaR:::build_lisa_contrast_navigation(contrast_navigation_fixture(),
    contrast_id = "response", scope = "response_pair", analysis_a = "on",
    analysis_b = "pre", collection = "GOBP-C2", tier = "core")
  out <- tempfile("contrast-navigation-policy-")
  on.exit(unlink(out, recursive = TRUE), add = TRUE)
  off <- lisaR:::render_lisa_contrast_navigation(nav, out,
    download_policy = list(svg = FALSE, source_data = FALSE, recipes = FALSE))
  original <- readBin(off$svg, "raw", n = file.info(off$svg)$size)
  html <- paste(readLines(off$html, warn = FALSE), collapse = "\n")
  expect_match(html, '<div class="diagram">', fixed = TRUE)
  expect_match(html, 'data-profile="A"', fixed = TRUE)
  expect_match(html, 'data-profile="B"', fixed = TRUE)
  expect_false(grepl('Download SVG|Download source table|>R script<', html))
  on <- lisaR:::render_lisa_contrast_navigation(nav, out,
    download_policy = list(svg = TRUE, source_data = FALSE, recipes = TRUE))
  html <- paste(readLines(on$html, warn = FALSE), collapse = "\n")
  expect_match(html, 'Download SVG', fixed = TRUE)
  expect_match(html, '>R script<', fixed = TRUE)
  expect_false(grepl('Download source table', html, fixed = TRUE))
  expect_identical(readBin(on$svg, "raw", n = file.info(on$svg)$size), original)
})

test_that("report_inline_category_navigation embeds the contrast navigator with links relative to contrasts.html, mirroring the analyses category-evidence-cover", {
  report <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = report)

  project_dir <- tempfile("contrast-navigator-integration-")
  dir.create(project_dir, recursive = TRUE)
  contrasts_file <- file.path(project_dir, "report_pages", "contrasts.html")
  dir.create(dirname(contrasts_file), recursive = TRUE)

  scope <- "response_dynamics_response_dynamics"; collection <- "GOBP-C2"
  evidence_dir <- file.path(project_dir, "report_pages", "contrast_evidence", scope, collection)
  dir.create(evidence_dir, recursive = TRUE)
  evidence_index <- file.path(evidence_dir, "index.html")
  writeLines("<!doctype html><html><body>stub contrast evidence page</body></html>", evidence_index)

  navigation_dir <- file.path(project_dir, "report_pages", "category_navigation", scope, collection)
  nav <- lisaR:::build_lisa_contrast_navigation(contrast_navigation_fixture(), contrast_id = "response_dynamics",
    scope = scope, analysis_a = "responders_on_vs_pre", analysis_b = "pd_on_vs_pre",
    collection = collection, tier = "core", label = "Treatment dynamics: responders minus PD")
  lisaR:::render_lisa_contrast_navigation(nav, navigation_dir)

  embedded <- report$report_inline_category_navigation(navigation_dir, contrasts_file)
  expect_match(embedded, '<div class="category-navigator"', fixed = TRUE)
  expect_match(embedded, "<svg", fixed = TRUE)
  # The raw navigator link is 3 levels up from category_navigation/<scope>/<collection>;
  # rebased against contrasts.html (one level up from report_pages), it becomes
  # a same-directory relative path, not the original "../../../" prefix.
  expect_match(embedded, paste0("contrast_evidence/", scope, "/", collection, "/index.html?category_id=DEV_MORPH"), fixed = TRUE)
  expect_false(grepl('href="\\.\\./\\.\\./\\.\\./contrast_evidence', embedded))

  # Absence must stay silent (older/standard reports without a navigator yet).
  expect_identical(report$report_inline_category_navigation(file.path(project_dir, "no-such-dir"), contrasts_file), "")
})
