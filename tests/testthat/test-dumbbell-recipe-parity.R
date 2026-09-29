dumbbell_parity_fixture <- function() {
  data.frame(
    category_id = c("CAT_SIGNAL", "CAT_BLANK", "CAT_CONTEXT"),
    category_display_name = c("Signal category", "Blank category", "Context category"),
    macrogroup_name = c("Immune programmes", "Immune programmes", "Matrix programmes"),
    color = c("#336699", "#777777", "#AA5500"),
    mean_NES_A = c(1.4, NA, -0.8),
    mean_NES_B = c(-1.1, NA, NA),
    mean_NES_contextual_A = c(1.4, NA, -0.8),
    mean_NES_contextual_B = c(-1.1, NA, 0.2),
    display_mean_NES_A = c(1.4, NA, -0.8),
    display_mean_NES_B = c(-1.1, NA, 0.2),
    n_genesets_A = c(3L, 0L, 2L),
    n_genesets_B = c(4L, 0L, 0L),
    n_genesets_evaluable_A = c(8L, 6L, 5L),
    n_genesets_evaluable_B = c(9L, 7L, 7L),
    same_direction_pct_A = c(75, NA, 100),
    same_direction_pct_B = c(50, NA, NA),
    has_significant_support_A = c(TRUE, FALSE, TRUE),
    has_significant_support_B = c(TRUE, FALSE, FALSE),
    plot_has_any_significant_support = c(TRUE, FALSE, TRUE),
    gsea_padj_cutoff = 0.05,
    stringsAsFactors = FALSE
  )
}

run_dumbbell_renderer <- function(renderer, source, output) {
  log <- system2(
    lisaR:::lisa_rscript_executable(),
    c("--vanilla", shQuote(renderer), shQuote(source), shQuote(output)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(log, "status")
  if (is.null(status)) status <- 0L
  list(status = status, log = log)
}

parse_dumbbell_record <- function(log) {
  record <- log[grepl("^LISA_DUMBBELL_REPRODUCTION\\t", log)]
  expect_equal(length(record), 1L, info = paste(log, collapse = "\n"))
  fields <- strsplit(record, "\t", fixed = TRUE)[[1L]][-1L]
  stats::setNames(sub("^[^=]*=", "", fields), sub("=.*$", "", fields))
}

test_that("formula-derived dumbbell dimensions survive TSV round trips", {
  fixture <- dumbbell_parity_fixture()
  observed <- fixture[rep(seq_len(nrow(fixture)), length.out = 40L), , drop = FALSE]
  observed$category_id <- sprintf("ROUNDTRIP_%02d", seq_len(nrow(observed)))

  for (variant in c("plain", "annotated")) {
    source <- tempfile(paste0("dumbbell-roundtrip-", variant, "-"), fileext = ".tsv")
    base_height <- lisaR:::lisa_contrast_plot_height(nrow(observed))
    expected_height <- if (identical(variant, "annotated")) {
      base_height * 1.08
    } else {
      base_height
    }
    expected_width <- if (identical(variant, "annotated")) 15.5 else 12

    expect_no_error(lisaR:::lisa_write_contrast_figure_source(
      observed, source, "roundtrip", "A minus B",
      "Formula-derived height round trip", "all", variant,
      "A", "B", 0.05, group_by_supracategory = TRUE
    ))
    written <- lisaR:::read_lisa_tsv(source)
    expect_equal(unique(written$figure_width), expected_width)
    expect_equal(
      unique(written$figure_height), expected_height,
      tolerance = .Machine$double.eps^0.5
    )
    expect_equal(unique(written$figure_dpi), 300)
  }
})

test_that("dumbbell recipe matches canonical grouped and ungrouped plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("png")

  renderer <- system.file("scripts", "reproduce_lisa_figure.R", package = "lisaR")
  expect_true(file.exists(renderer))
  fixture <- dumbbell_parity_fixture()

  run_case <- function(variant, grouped) {
    annotated <- identical(variant, "annotated")
    width <- if (annotated) 8 else 6
    height <- 4
    dpi <- 72
    source <- tempfile(paste0("dumbbell-", variant, "-"), fileext = ".tsv")
    output <- tempfile(paste0("dumbbell-", variant, "-"), fileext = ".png")
    canonical <- tempfile(paste0("dumbbell-canonical-", variant, "-"), fileext = ".png")

    lisaR:::lisa_write_contrast_figure_source(
      fixture, source, "fixture", "A minus B fixture",
      "Responder minus progressor (A - B)\nPositive values favour A.", "all", variant,
      "Responder", "Progressor", 0.05,
      group_by_supracategory = grouped,
      figure_width = width, figure_height = height, figure_dpi = dpi
    )
    plot <- lisaR:::plot_lisa_contrast_dumbbell(
      fixture, "Responder", "Progressor", "Responder", "Progressor",
      "A minus B fixture", "Responder minus progressor (A - B)\nPositive values favour A.",
      group_by_supracategory = grouped,
      annotate_gene_sets = annotated
    )
    ggplot2::ggsave(
      canonical, plot, width = width, height = height, dpi = dpi,
      bg = "white", limitsize = FALSE
    )
    result <- run_dumbbell_renderer(renderer, source, output)
    expect_equal(result$status, 0L, info = paste(result$log, collapse = "\n"))
    expect_true(file.exists(output))

    observed <- png::readPNG(output)
    expected <- png::readPNG(canonical)
    expect_identical(dim(observed), dim(expected))
    expect_equal(observed, expected, tolerance = 1 / 255)
    list(metadata = parse_dumbbell_record(result$log), dimensions = dim(observed))
  }

  annotated <- run_case("annotated", TRUE)
  expect_identical(annotated$metadata[c(
    "plot_set", "variant", "grouped", "facets", "panel_rows", "categories",
    "blank", "segments", "supported", "contextual", "labels"
  )], c(
    plot_set = "all", variant = "annotated", grouped = "true", facets = "2",
    panel_rows = "2,1", categories = "3", blank = "1", segments = "2",
    supported = "3", contextual = "1", labels = "4"
  ))
  expect_identical(annotated$dimensions[2:1], c(8L * 72L, 4L * 72L))

  plain <- run_case("plain", FALSE)
  expect_identical(plain$metadata[["variant"]], "plain")
  expect_identical(plain$metadata[["grouped"]], "false")
  expect_identical(plain$metadata[["facets"]], "1")
  expect_identical(plain$metadata[["panel_rows"]], "3")
  expect_identical(plain$metadata[["labels"]], "0")
  expect_identical(plain$dimensions[2:1], c(6L * 72L, 4L * 72L))
})

test_that("empty dumbbell recipes preserve canonical variant dimensions", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("png")
  renderer <- system.file("scripts", "reproduce_lisa_figure.R", package = "lisaR")
  empty <- dumbbell_parity_fixture()[FALSE, , drop = FALSE]

  for (variant in c("plain", "annotated")) {
    source <- tempfile(paste0("empty-", variant, "-"), fileext = ".tsv")
    output <- tempfile(paste0("empty-", variant, "-"), fileext = ".png")
    width <- if (identical(variant, "annotated")) 15.5 else 12
    height <- if (identical(variant, "annotated")) 4.8 * 1.08 else 4.8
    lisaR:::lisa_write_contrast_figure_source(
      empty, source, "fixture", "Empty fixture", "No matching categories",
      "opposite_direction", variant, "A", "B", 0.05,
      group_by_supracategory = TRUE, figure_dpi = 72
    )
    result <- run_dumbbell_renderer(renderer, source, output)
    expect_equal(result$status, 0L, info = paste(result$log, collapse = "\n"))
    metadata <- parse_dumbbell_record(result$log)
    expect_identical(metadata[["categories"]], "0")
    expect_equal(as.numeric(metadata[["width"]]), width)
    expect_equal(as.numeric(metadata[["height"]]), height)
    image <- png::readPNG(output)
    expect_equal(dim(image)[[2L]], round(width * 72), tolerance = 1)
    expect_equal(dim(image)[[1L]], round(height * 72), tolerance = 1)
  }
})

test_that("dumbbell recipe rejects aliases and inconsistent metadata before writing", {
  skip_if_not_installed("ggplot2")
  renderer <- system.file("scripts", "reproduce_lisa_figure.R", package = "lisaR")
  fixture <- dumbbell_parity_fixture()
  source <- tempfile("dumbbell-safe-", fileext = ".tsv")
  lisaR:::lisa_write_contrast_figure_source(
    fixture, source, "fixture", "A minus B", "Fixture", "all", "plain",
    "A", "B", 0.05, figure_dpi = 72
  )
  source_sha <- lisaR:::lisa_sha256_file(source)

  if (.Platform$OS.type != "windows") {
    alias <- tempfile("dumbbell-alias-", fileext = ".png")
    expect_true(file.symlink(source, alias))
    result <- suppressWarnings(run_dumbbell_renderer(renderer, source, alias))
    expect_false(identical(result$status, 0L))
    expect_match(paste(result$log, collapse = "\n"), "symbolic link", ignore.case = TRUE)
    expect_identical(lisaR:::lisa_sha256_file(source), source_sha)
  }

  original <- lisaR:::read_lisa_tsv(source)
  for (field in c("plot_set", "annotation_variant", "group_by_supracategory")) {
    mixed <- original
    mixed[[field]][[2L]] <- switch(
      field,
      plot_set = "same_direction",
      annotation_variant = "annotated",
      group_by_supracategory = FALSE
    )
    mixed_source <- tempfile(paste0("dumbbell-mixed-", field, "-"), fileext = ".tsv")
    output <- tempfile(paste0("dumbbell-mixed-", field, "-"), fileext = ".png")
    utils::write.table(
      mixed, mixed_source, sep = "\t", quote = FALSE, row.names = FALSE, na = ""
    )
    result <- suppressWarnings(run_dumbbell_renderer(renderer, mixed_source, output))
    expect_false(identical(result$status, 0L))
    expect_match(paste(result$log, collapse = "\n"), "inconsistent metadata")
    expect_false(file.exists(output))
  }
})
