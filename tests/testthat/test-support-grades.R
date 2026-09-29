support_fixture <- function() {
  n <- c(579L, 10L, 254L, 8L, 1L, 10L, 0L)
  d <- c(1L, 1L, 104L, 5L, 1L, 0L, NA_integer_)
  data.frame(analysis_id = "A", collection = "GOCC", category_id = paste0("C", seq_along(n)),
    category_display_name = c("Duplicate label", "Duplicate label", "Broad support", "Half or more", "Single set", "Zero bound", "Not evaluable"),
    n_sets_evaluable = n, minimum_enriched_sets = d, minimum_enriched_pct = 100 * d/n,
    n_sets_total = pmax(n, 1L), n_sets_missing_p = as.integer(n == 0), n_sets_not_eligible_or_absent = 0L,
    confidence_level = .95, category_p_adjusted = c(rep(.001, 5), 1, NA_real_),
    significant = !is.na(d) & d > 0, status = c(rep("significant", 5), "not_significant", "not_evaluable"),
    n_positive_NES = pmax(n - 1L, 0L), n_negative_NES = as.integer(n > 0),
    n_zero_NES = 0L, n_missing_NES = 0L, stringsAsFactors = FALSE)
}

support_lollipop_fixture <- function(variant = "plain", grouped = TRUE) {
  x <- support_fixture()
  base <- data.frame(category_id = x$category_id, category_display_name = x$category_display_name,
    macrogroup_id = "M", macrogroup_name = "Group", macrogroup_order = 1L,
    category_order_within_macrogroup = seq_len(nrow(x)),
    mean_NES = c(NA, 1, -1, 2, -2, NA, NA), n_genesets = c(0, 1, 30, 5, 1, 0, 0),
    same_direction_pct = c(NA, 100, 60, 80, 100, NA, NA),
    consistency = c(NA, 1, .6, .8, 1, NA, NA), color = "#445566")
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  lisaR:::lisa_write_lollipop_figure_source(base, path, grouped, "supracategory",
    "Support fixture", "Descriptive NES", "Positive favours group A", annotation_variant = variant,
    figure_dpi = 72)
  lisaR:::lisa_read_figure_source_tsv(path)
}

test_that("editorial boundaries use exact counts and never upward rounding", {
  d <- c(0, 1, 9, 10, 24, 25, 49, 50, 74, 75, 100, 1, 104, 5, 1, 1,
    214748364, 214748365, 536870911, 536870912, 1073741823, 1073741824, 1610612735, 1610612736)
  n <- c(rep(100, 11), 579, 254, 8, 1, 2147483647, rep(2147483647, 8))
  z <- lisaR:::lisa_support_grade(d, n)
  expect_identical(z$support_grade, c("", "S1", "S1", "S2", "S2", "S3", "S3", "S4", "S4", "S5", "S5", "S1", "S3", "S4", "S5", "S1", "S1", "S2", "S2", "S3", "S3", "S4", "S4", "S5"))
  expect_identical(z$minimum_support_pct_text[c(12, 13, 14, 15, 16)], c("0.17%", "40.94%", "62.5%", "100%", "0.00000004%"))
  shown <- as.numeric(sub("%", "", z$minimum_support_pct_text))
  expect_true(all(shown <= 100 * d/n))
  expect_true(all(shown[d > 0] > 0))
  expect_identical(lisaR:::lisa_support_grade(1, 3)$minimum_support_pct_text, "33.33%")
  expect_identical(lisaR:::lisa_support_grade(1, 101)$minimum_support_pct_text, "0.99%")
  expect_identical(lisaR:::lisa_support_grade(999, 10000)$support_grade, "S1")
  expect_identical(lisaR:::lisa_support_grade(999, 10000)$minimum_support_pct_text, "9.99%")
  expect_error(lisaR:::lisa_support_grade(1, 0), "integer counts")
  expect_error(lisaR:::lisa_support_grade(NA_real_, 1), "integer counts")
  expect_error(lisaR:::lisa_support_grade(.5, 1), "integer counts")
  expect_error(lisaR:::lisa_support_grade(1, .Machine$integer.max + 1), "integer counts")
})

test_that("zero and no evaluable members stay distinct; scientific values remain unchanged", {
  x <- support_fixture()
  z <- lisaR:::lisa_category_support_presentation(x)
  expect_identical(z[names(x)], x)
  expect_identical(z$support_grade[6:7], c("", ""))
  expect_identical(z$support_annotation[6:7], c("", ""))
  expect_identical(z$minimum_support_pct_text[6:7], c("0%", NA_character_))
  expect_identical(z$support_count[6:7], c("0/10", NA_character_))
  expect_match(z$support_label[6], "No minimum greater than zero")
  expect_match(z$support_label[7], "Not evaluable")
  expect_match(z$observed_NES_pattern[1], "both signs")
  expect_identical(z$observed_NES_pattern[5], "Observed negative NES")
  display <- lisaR:::lisa_category_inference_display(x)
  expect_identical(display[["Minimum enriched / evaluable sets (d/N)"]][6:7], c("0/10", "Not evaluable"))
  bad <- x; bad$collection <- "HALLMARKS"
  expect_error(lisaR:::lisa_category_support_presentation(bad), "HALLMARKS")
  bad <- x; bad$significant[1] <- FALSE
  expect_error(lisaR:::lisa_category_support_presentation(bad), "significance disagree")
  expect_error(lisaR:::lisa_category_support_presentation(rbind(x, x[1, ])), "Duplicate")
})

test_that("support lollipops join identity not labels and do not fabricate means", {
  x <- support_fixture()
  original <- support_lollipop_fixture()
  z <- lisaR:::lisa_support_lollipop_source(original, x, "A", "GOCC")
  expect_identical(z$mean_NES, original$mean_NES)
  expect_identical(z$n_genesets, original$n_genesets)
  expect_identical(z$color, original$color)
  expect_identical(z$support_grade[1:2], c("S1", "S2"))
  expect_false(z$has_lollipop_geometry[1])
  expect_match(z$support_annotation[1], "1/579", fixed = TRUE)
  missing <- lisaR:::lisa_support_lollipop_source(original[-1, ], x, "A", "GOCC")
  expect_true(is.na(missing$mean_NES[1]))
  expect_identical(missing$support_annotation[1], z$support_annotation[1])
  x$analysis_id[1] <- "B"
  expect_error(lisaR:::lisa_support_lollipop_source(original, x, "A", "GOCC"), "different analysis")
})

test_that("both support lollipop recipes retain exact labels and missing geometry", {
  skip_if_not_installed("png")
  for (variant in c("plain", "direction")) for (grouped in c(TRUE, FALSE)) {
    z <- lisaR:::lisa_support_lollipop_source(support_lollipop_fixture(variant, grouped), support_fixture(), "A", "GOCC")
    z$figure_dpi <- 30
    path <- tempfile(fileext = ".tsv")
    on.exit(unlink(path), add = TRUE)
    lisaR:::lisa_write_inference_tsv(z, path)
    saved <- lisaR:::lisa_read_figure_source_tsv(path)
    rendered <- lisaR:::lisa_lollipop_from_source(saved)
    build <- ggplot2::ggplot_build(rendered$plot)
    labels <- build$data[[4]]$label
    expect_setequal(labels, z$support_annotation)
    expect_false(any(grepl("^S[1-5]", labels[!nzchar(labels)])))
    points <- build$data[[7]]
    expect_equal(nrow(points), 4L)
    expect_setequal(points$x, c(1, -1, 2, -2))
    expected <- tempfile(fileext = ".png"); replay <- tempfile(fileext = ".png")
    on.exit(unlink(c(expected, replay)), add = TRUE)
    lisaR:::lisa_support_save_lollipop(saved, expected)
    e <- new.env(parent = globalenv()); e$output_png <- replay
    for (expr in parse(system.file("scripts", "reproduce_lisa_figure.R", package = "lisaR")))
      if (is.call(expr) && identical(expr[[1]], as.name("<-")) && identical(expr[[2]], as.name("render_lisa_lollipop"))) eval(expr, e)
    e$render_lisa_lollipop(saved)
    expect_equal(png::readPNG(replay), png::readPNG(expected), tolerance = 1/255)
  }
})

test_that("an entirely absent mean-NES summary still displays support without fake points", {
  for (variant in c("plain", "direction")) {
    original <- support_lollipop_fixture(variant)
    original$mean_NES <- NA_real_
    original$n_genesets <- 0L
    original$same_direction_pct <- original$consistency <- NA_real_
    z <- lisaR:::lisa_support_lollipop_source(original, support_fixture(), "A", "GOCC")
    expect_no_error(rendered <- lisaR:::lisa_lollipop_from_source(z))
    expect_no_error(build <- ggplot2::ggplot_build(rendered$plot))
    expect_equal(nrow(build$data[[7]]), 0L)
    expect_setequal(build$data[[4]]$label, z$support_annotation)
  }
})

test_that("presentation refresh is copy-only and does not recompute inference", {
  root <- tempfile("support-copy-"); origin <- file.path(root, "source"); out <- file.path(root, "new")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  family <- file.path(origin, "outputs", "category_inference", "A")
  collection <- file.path(origin, "outputs", "single_de", "A", "collection_GOCC")
  for (path in c(family, file.path(collection, c("plots", "lisa_tables")), file.path(origin, "audit"))) dir.create(path, recursive = TRUE)
  x <- support_fixture(); write <- lisaR:::lisa_write_inference_tsv
  write(x, file.path(family, "category_results.tsv"))
  write(x, file.path(collection, "lisa_tables", "A_LISA_category_inference.tsv"))
  writeLines("untouched family", file.path(family, "hypothesis_family.tsv"))
  writeLines("untouched raw results", file.path(family, "raw_results.tsv"))
  writeLines('{"original":"method"}', file.path(family, "method.json"))
  members <- data.frame(analysis_id = "A", collection = "GOCC", category_id = x$category_id[1:6], pathway = paste0("set", 1:6),
    NES = 1, pval = .001, hypothesis_id = paste0("set", 1:6), included_in_family = TRUE, raw_p_available = TRUE)
  fig <- lisaR:::lisa_category_inference_source(x, members)
  fig <- fig[, !names(fig) %in% c(names(lisaR:::lisa_support_grade(1, 1)), "observed_NES_pattern", "observed_NES_counts")]
  write(fig, file.path(collection, "plots", "A_LISA_category_inference_source.tsv"))
  for (variant in c("plain", "direction")) write(support_lollipop_fixture(variant), file.path(collection, "plots",
    paste0("A_GSEA_lollipop", if (variant == "direction") "_direction_stats" else "", "_source.tsv")))
  manifest <- lisaR:::lisa_run_manifest(origin)
  lisaR:::write_lisa_tsv(manifest, file.path(origin, "run_manifest.tsv"))
  hashes <- vapply(file.path(origin, manifest$path), lisaR:::lisa_sha256_file, character(1))
  testthat::local_mocked_bindings(lisa_category_inference = function(...) stop("Must not recalculate"), .package = "lisaR")
  refresh_lisa_support_grades(origin, out, plot_formats = character())
  expect_identical(vapply(file.path(origin, manifest$path), lisaR:::lisa_sha256_file, character(1)), hashes)
  expect_false(file.exists(file.path(out, "run_manifest.tsv")))
  y <- lisaR:::read_lisa_tsv(file.path(out, "outputs", "category_inference", "A", "category_results.tsv"))
  expect_equal(y[names(x)], x, tolerance = 0)
  for (name in c("hypothesis_family.tsv", "raw_results.tsv", "method.json")) expect_identical(
    lisaR:::lisa_sha256_file(file.path(family, name)), lisaR:::lisa_sha256_file(file.path(out, "outputs", "category_inference", "A", name)))
  expect_true(file.exists(file.path(out, "outputs", "single_de", "A", "collection_GOCC", "plots", "A_GSEA_lollipop_support_source.tsv")))
  expect_true(file.exists(file.path(out, "outputs", "single_de", "A", "collection_GOCC", "plots", "A_GSEA_lollipop_hommel_support_source.tsv")))
  expect_true(lisaR:::lisa_finalize_support_grades_derivation(out))
  expect_identical(lisaR:::lisa_validate_built_run_manifest(out, lisaR:::read_lisa_tsv(file.path(out, "run_manifest.tsv")))$gate, "PASS")
  # The supported input to a second UI revision is already a derived run.
  # Keep the original support sources and its archive byte-for-byte.
  old_files <- list.files(out, recursive = TRUE, full.names = TRUE)
  old_files <- old_files[grepl("support_source[.]tsv$|historical_v1|source_provenance/", old_files) & !grepl("hommel_support_source", old_files)]
  old_hashes <- vapply(old_files, lisaR:::lisa_sha256_file, character(1))
  second <- file.path(root, "second")
  prior_inference <- file.path(out, "outputs", "single_de", "A", "collection_GOCC", "plots", "A_LISA_category_inference_source.tsv")
  prior_inference_sha <- lisaR:::lisa_sha256_file(prior_inference)
  refresh_lisa_support_grades(out, second, plot_formats = character())
  copied <- file.path(second, substring(old_files, nchar(out) + 2L))
  expect_identical(unname(vapply(copied, lisaR:::lisa_sha256_file, character(1))), unname(old_hashes))
  expect_true(file.exists(file.path(second, "source_provenance", "hommel_support_ui_v1_2", "parent_run_manifest.tsv")))
  expect_identical(lisaR:::lisa_sha256_file(file.path(second, "outputs", "single_de", "A", "collection_GOCC", "plots",
    "A_LISA_category_inference_historical_support_v1_source.tsv")), prior_inference_sha)
  expect_true(lisaR:::lisa_finalize_support_grades_derivation(second))
  expect_error(refresh_lisa_support_grades(origin, out), "must be new")
  expect_error(refresh_lisa_support_grades(origin, file.path(origin, "inside")), "outside")
  writeLines("corrupt", file.path(family, "raw_results.tsv"))
  expect_error(refresh_lisa_support_grades(origin, file.path(root, "reject")), "failed integrity")
  expect_false(dir.exists(file.path(root, "reject")))
})

test_that("compact result downloads preserve historical figures without a second results section", {
  e <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = e)
  root <- tempfile("support-html-"); collection <- file.path(root, "outputs", "single_de", "A", "collection_GOCC")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(collection, "plots"), recursive = TRUE)
  dir.create(file.path(collection, "lisa_tables"))
  table <- lisaR:::lisa_category_support_presentation(support_fixture())
  lisaR:::lisa_write_inference_tsv(table, file.path(collection, "lisa_tables", "A_LISA_category_inference.tsv"))
  for (name in c("A_LISA_category_inference.png", "A_GSEA_lollipop.png", "A_GSEA_lollipop_direction_stats.png",
    "A_GSEA_lollipop_support.png", "A_GSEA_lollipop_direction_stats_support.png", "A_GSEA_lollipop_support_source.tsv"))
    file.create(file.path(collection, "plots", name))
  e$report_image_formats <- function() "png"
  e$figure_card <- function(path, ...) paste0('<figure>', basename(path), '</figure>')
  e$table_html <- function(x, ...) paste(unlist(x), collapse = "|")
  e$short_file_copy <- function(path, ...) basename(path)
  e$short_media_copy <- function(path, ...) basename(path)
  e$figure_source_contract <- function(...) list(source_tsv="saved-figure.tsv.gz", recipe_r="saved-figure.R")
  html <- paste0(e$report_category_inference_downloads(collection, file.path(root, "report_pages", "single_de.html")),
    e$report_support_lollipops(collection, file.path(root, "report_pages", "single_de.html")))
  for (text in c('class="support-lollipop"', 'class="support-direction-lollipop"',
    'class="support-historical-lollipops"', 'not all evaluated sets',
    'data-category-results-tsv', 'Complete category results (TSV)', 'saved-figure.tsv.gz', 'saved-figure.R')) expect_match(html, text, fixed = TRUE)
  expect_false(grepl('\\* Significant category', html))
  expect_false(grepl('hommel-support-technical|support-results-table|Category enrichment: minimum supported extent', html))
  expect_identical(e$report_category_inference_downloads(file.path(root, "collection_HALLMARKS"), "unused"), "")
})

test_that("portable support lollipop downloads preserve every scientific numeric cell", {
  e <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = e)
  root <- tempfile("support-download-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "outputs"), recursive = TRUE)
  dir.create(file.path(root, "report_pages"))
  figure <- file.path(root, "outputs", "A_GSEA_lollipop_support.png")
  file.create(figure)
  x <- lisaR:::lisa_support_lollipop_source(support_lollipop_fixture(), support_fixture(), "A", "GOCC")
  x$category_p_adjusted[1] <- 1.0143905644500296e-07
  x$minimum_enriched_pct[1] <- 1.1173184357541899
  source <- sub("[.]png$", "_source.tsv", figure)
  lisaR:::lisa_write_inference_tsv(x, source)
  e$report_requested <- function(key, default = TRUE) key != "recipes"
  result <- e$figure_source_contract(figure, file.path(root, "report_pages", "single_de.html"))
  saved <- lisaR:::lisa_read_figure_source_tsv(file.path(root, "report_pages", result$source_tsv))
  for (field in names(x)[vapply(x, is.numeric, logical(1))]) expect_equal(saved[[field]], x[[field]], tolerance = 0, info = field)
  expect_identical(saved$category_p_adjusted, x$category_p_adjusted)
  expect_identical(saved$minimum_enriched_pct, x$minimum_enriched_pct)
  expect_identical(saved$support_annotation, x$support_annotation)
})
