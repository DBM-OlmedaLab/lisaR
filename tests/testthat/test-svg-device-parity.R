# Public synthetic summaries, real devices and complete-file comparisons.
# No enrichment, external fixture paths, mocked renderer or SVG normalization.
# Bounded: 12 history/plot SVG children + one path check + one failure child.

svg_parity_require_cairo <- function() {
  if (!isTRUE(capabilities("cairo"))) {
    if (identical(Sys.info()[["sysname"]], "Linux"))
      stop("Linux SVG parity requires native Cairo; this is not a skip.")
    skip("Native Cairo is unavailable on this non-Linux platform.")
  }
}

svg_parity_plots <- function() {
  # Explicit public controls: one significant mixed-direction category and
  # one evaluated category with no significant members. No study data.
  summary <- data.frame(
    category_id = c("CAT_BP", "CAT_MF"),
    category_display_name = c("Positive control", "Mixed empty eligibility"),
    macrogroup_id = c("SUPER_BP", "SUPER_MF"),
    macrogroup_name = c("Positive control", "Eligibility control"),
    macrogroup_order = c(1L, 2L), category_order_within_macrogroup = 1L,
    color = "#E16A86", n_genesets = c(2L, 0L),
    n_genesets_mapped = c(2L, 4L), n_genesets_evaluable = c(2L, 1L),
    n_genesets_significant = c(2L, 0L), n_pos_genesets = c(1L, 0L),
    n_neg_genesets = c(1L, 0L), mean_NES = c(0.0364306891523489, NA_real_),
    mean_NES_contextual = c(0.0364306891523489, 0.589926131970361),
    mean_pos_NES = c(3.96542567802728, NA_real_),
    mean_neg_NES = c(-3.89256429972258, NA_real_),
    median_NES = c(0.0364306891523489, NA_real_),
    consistency = c(0.5, NA_real_), min_padj = c(3.71106344204542e-54, NA_real_),
    min_padj_mapped = c(3.71106344204542e-54, 0.975378787878788),
    has_significant_support = c(TRUE, FALSE), gsea_padj_cutoff = 0.05,
    same_direction_pct = c(50, NA_real_),
    mean_NES_direction = c("Positive NES", NA_character_),
    analysis = "GSEA", stringsAsFactors = FALSE
  )
  contextualize <- function(plot) lisaR:::lisa_apply_plot_context(plot,
    "A", "Synthetic A minus control", "Higher synthetic values in A")
  empty <- summary[2L, , drop = FALSE]
  filled <- summary[1L, , drop = FALSE]
  list(
    empty = contextualize(lisaR:::plot_lisa_gsea_lollipop(empty, TRUE, "supracategory")),
    direction = contextualize(lisaR:::plot_lisa_gsea_direction_lollipop(empty, TRUE, "supracategory")),
    filled = contextualize(lisaR:::plot_lisa_gsea_lollipop(filled, TRUE, "supracategory"))
  )
}

svg_parity_bytes <- function(path) {
  stopifnot(file.exists(path), !is.na(file.info(path)$size), file.info(path)$size > 0)
  readBin(path, "raw", n = file.info(path)$size)
}

test_that("real SVG bytes ignore filled SVG, 180 dpi PNG and active PDF histories", {
  svg_parity_require_cairo()
  output <- tempfile("lisa svg device parity ")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  old_options <- options(lisaR.run_root = NULL)
  on.exit(options(old_options), add = TRUE)
  plots <- svg_parity_plots()
  references <- list()
  for (history in c("cold", "filled-svg", "png180-after-svg", "caller-pdf")) {
    dest <- file.path(output, history)
    dir.create(dest)
    # Render in THIS process to reproduce the mixed-device history. Calling
    # save_plot_multi for this primer would isolate it and invalidate the test.
    if (history != "cold") ggplot2::ggsave(file.path(dest, "prior.svg"),
      plots$filled, width = 11, height = 12, device = grDevices::svg,
      limitsize = FALSE, bg = "white")
    if (history == "png180-after-svg") ggplot2::ggsave(
      file.path(dest, "prior180.png"), plots$filled, width = 4, height = 3,
      dpi = 180, device = "png", limitsize = FALSE, bg = "white")
    caller <- NULL
    if (history == "caller-pdf") {
      grDevices::pdf(file.path(dest, "caller.pdf"), width = 5, height = 4)
      caller <- grDevices::dev.cur()
      on.exit({
        if (caller %in% grDevices::dev.list()) grDevices::dev.off(caller)
      }, add = TRUE, after = FALSE)
      graphics::plot(1, 1)
    }
    for (kind in names(plots)) {
      stem <- file.path(dest, kind)
      before_devices <- grDevices::dev.list()
      before_current <- grDevices::dev.cur()
      lisaR:::save_plot_multi(plots[[kind]], stem, "svg",
        width = if (kind == "direction") 11.5 else 11, height = 12)
      expect_identical(grDevices::dev.list(), before_devices)
      expect_identical(grDevices::dev.cur(), before_current)
      svg <- paste0(stem, ".svg")
      expect_true(file.exists(svg))
      header <- paste(readLines(svg, n = 3L, warn = FALSE), collapse = "\n")
      expect_match(header, "<svg")
      expect_match(header, if (kind == "direction") 'width="828pt"' else 'width="792pt"', fixed = TRUE)
      expect_match(header, 'height="864pt"', fixed = TRUE)
      bytes <- svg_parity_bytes(svg)
      if (history == "cold") references[[kind]] <- bytes else
        expect_identical(bytes, references[[kind]])
    }
    if (!is.null(caller)) {
      expect_identical(grDevices::dev.cur(), caller)
      grDevices::dev.off(caller)
    }
  }
})

test_that("PNG retains direct renderer bytes and SVG handles space-containing paths", {
  svg_parity_require_cairo()
  output <- tempfile("lisa svg path with spaces ")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  old_options <- options(lisaR.run_root = NULL)
  on.exit(options(old_options), add = TRUE)
  plot <- svg_parity_plots()$empty
  direct <- file.path(output, "direct.png")
  ggplot2::ggsave(direct, plot, width = 11, height = 12, dpi = 300,
    device = "png", limitsize = FALSE, bg = "white")
  lisaR:::save_plot_multi(plot, file.path(output, "package"), "png", 11, 12)
  expect_identical(svg_parity_bytes(file.path(output, "package.png")),
    svg_parity_bytes(direct))
  stem <- file.path(output, "name with spaces")
  before_devices <- grDevices::dev.list()
  before_current <- grDevices::dev.cur()
  lisaR:::save_plot_multi(plot, stem, "svg", 11, 12)
  expect_identical(grDevices::dev.list(), before_devices)
  expect_identical(grDevices::dev.cur(), before_current)
  expect_match(paste(readLines(paste0(stem, ".svg"), n = 3L,
    warn = FALSE), collapse = "\n"), "<svg")
  expect_gt(length(svg_parity_bytes(paste0(stem, ".svg"))), 0L)
})

test_that("failed SVG child preserves the guarded destination and caller device", {
  svg_parity_require_cairo()
  output <- tempfile("lisa-svg-render-failure-")
  dir.create(output)
  on.exit(unlink(output, recursive = TRUE), add = TRUE)
  old_options <- options(lisaR.run_root = NULL)
  on.exit(options(old_options), add = TRUE)
  target <- file.path(output, "preserved.svg")
  writeLines("existing-sentinel", target)
  before <- svg_parity_bytes(target)
  grDevices::pdf(file.path(output, "caller.pdf"), width = 5, height = 4)
  caller <- grDevices::dev.cur()
  on.exit({
    if (caller %in% grDevices::dev.list()) grDevices::dev.off(caller)
  }, add = TRUE, after = FALSE)
  graphics::plot(1, 1)
  before_devices <- grDevices::dev.list()
  before_files <- sort(list.files(output, all.files = TRUE, no.. = TRUE))
  before_bundles <- sort(list.files(tempdir(), pattern = "^lisa-svg-device-",
    full.names = TRUE, all.files = TRUE))
  bad <- structure(list(), class = "lisa_svg_parity_invalid_plot")
  expect_error(lisaR:::save_plot_multi(bad,
    sub("[.]svg$", "", target), "svg", 11, 12),
    "Isolated SVG figure rendering failed")
  expect_identical(svg_parity_bytes(target), before)
  expect_identical(grDevices::dev.list(), before_devices)
  expect_identical(grDevices::dev.cur(), caller)
  expect_identical(sort(list.files(output, all.files = TRUE, no.. = TRUE)), before_files)
  expect_identical(sort(list.files(tempdir(), pattern = "^lisa-svg-device-",
    full.names = TRUE, all.files = TRUE)), before_bundles)
})
