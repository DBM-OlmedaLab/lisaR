test_that("external verification checks declared sizes as well as content hashes", {
  final <- tempfile("audit-byte-metadata-")
  tx <- lisaR:::lisa_transaction_begin(final, lisa_test_run_contract(), run_id = "auditbytes")
  writeLines("unchanged content", file.path(tx$staging_dir, "artifact.txt"))
  lisaR:::lisa_transaction_promote(tx)
  expect_identical(lisaR:::verify_run(final)$gate, "PASS")
  manifest_path <- file.path(final, "run_manifest.tsv")
  original <- lisaR:::read_lisa_tsv(manifest_path)
  for (invalid in c("99999", "not-a-size", NA_character_)) {
    manifest <- original
    manifest$bytes <- as.character(manifest$bytes)
    manifest$bytes[manifest$path == "artifact.txt"] <- invalid
    lisaR:::write_lisa_tsv(manifest, manifest_path)
    observed <- lisaR:::verify_run(final)
    expect_identical(observed$gate, "FAIL")
    expect_true("resized:artifact.txt" %in% observed$findings)
    expect_false("modified:artifact.txt" %in% observed$findings)
  }
})

test_that("analysis path parser explicitly requests Windows-compatible slashes", {
  script <- system.file("scripts", "build_single_de_category_count_heatmaps.R", package = "lisaR")
  code <- parse(script)
  e <- new.env(parent = baseenv())
  # Emulate the documented Windows normalizePath return convention, not an
  # OS certification. The real parser must request '/' to recover the ID.
  e$normalizePath <- function(path, winslash = "\\", mustWork = FALSE) {
    gsub("/", winslash, "C:/study/outputs/single_de/response/collection_GOCC/lisa_tables/summary.tsv", fixed = TRUE)
  }
  for (expr in code) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        identical(expr[[2L]], as.name("parse_analysis_id"))) eval(expr, e)
  }
  expect_identical(e$parse_analysis_id("ignored"), "response")
})

lollipop_audit_fixture <- function() data.frame(
  category_id = c("CAT_B", "CAT_EMPTY", "CAT_A"),
  category_display_name = c("B programme", "No support", "A programme"),
  macrogroup_id = c("M2", "M1", "M1"),
  macrogroup_name = c("Second group", "First group", "First group"),
  macrogroup_order = c(2, 1, 1), category_order_within_macrogroup = c(1, 2, 1),
  mean_NES = c(-1.2, NA_real_, 1.8), n_genesets = c(2L, 0L, 5L),
  same_direction_pct = c(50, NA_real_, 80), consistency = c(.5, NA_real_, .8),
  color = c("#115599", "#777777", "#BB5500"), stringsAsFactors = FALSE)

test_that("lollipop recipes retain geometry support groups order and context", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("png")
  fixture <- lollipop_audit_fixture()
  script <- system.file("scripts", "reproduce_lisa_figure.R", package = "lisaR")
  exprs <- parse(script)
  for (variant in c("plain", "direction")) for (grouped in c(TRUE, FALSE)) {
    source <- tempfile(fileext = ".tsv")
    original <- tempfile(fileext = ".png")
    replay <- tempfile(fileext = ".png")
    title <- "Response %00"
    subtitle <- "On treatment\nminus pretreatment"
    lisaR:::lisa_write_lollipop_figure_source(fixture, source, grouped,
      "mean_NES", title, subtitle, "Positive favours on treatment",
      annotation_variant = variant, figure_width = 6, figure_height = 4,
      figure_dpi = 72)
    e <- new.env(parent = globalenv())
    e$output_png <- replay
    for (expr in exprs) if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        identical(expr[[2L]], as.name("render_lisa_lollipop"))) eval(expr, e)
    written <- utils::read.delim(source, check.names = FALSE, quote = "")
    e$render_lisa_lollipop(written)
    renderer <- if (variant == "plain") lisaR:::plot_lisa_gsea_lollipop else lisaR:::plot_lisa_gsea_direction_lollipop
    plot <- renderer(fixture, grouped, "mean_NES")
    plot <- lisaR:::lisa_apply_plot_context(plot, title, subtitle, "Positive favours on treatment")
    ggplot2::ggsave(original, plot, width = 6, height = 4, dpi = 72,
      bg = "white", device = "png", limitsize = FALSE)
    expect_equal(png::readPNG(replay), png::readPNG(original), tolerance = 1 / 255)
    expect_identical(written$selected_for_plot, c(TRUE, FALSE, TRUE))
    expect_identical(written$color, fixture$color)
  }
})

test_that("lollipop source handles absent titles and rejects legacy incomplete contracts", {
  fixture <- lollipop_audit_fixture()
  source <- tempfile(fileext = ".tsv")
  lisaR:::lisa_write_lollipop_figure_source(fixture, source, TRUE, "supracategory")
  written <- utils::read.delim(source, check.names = FALSE, quote = "")
  expect_s3_class(lisaR:::lisa_lollipop_from_source(written)$plot, "ggplot")
  expect_error(lisaR:::lisa_lollipop_from_source(fixture), "metadata")
  written$group_by_supracategory[[1L]] <- FALSE
  expect_error(lisaR:::lisa_lollipop_from_source(written), "consistent metadata")
})
