lisa_nes_variant_fixture <- function() {
  x <- data.frame(category_id = c(rep("MIX", 5), "ONE", "NONE", "NONE", "ZERO", "ZERO", "OTHER_UNCLASSIFIED"),
    pathway = c("s1", "s2", "s3", "s4", "s5", "001", "n1", "n2", "z1", "z2", "u1"),
    NES = c(-2, -1, 0, 1, 4, 1.2345678901234567, 2, 8, -1, 1, 99),
    padj = c(.01, .02, .03, .04, .05, .001, .3, NA, .01, .01, .001),
    stringsAsFactors = FALSE)
  # One repeated dictionary assignment is not another enriched set.
  rbind(x, x[1L, , drop = FALSE])
}

lisa_nes_variant_summary_fixture <- function() {
  data.frame(category_id = c("MIX", "ONE", "NONE", "ZERO", "EMPTY"),
    category_display_name = c("Mixed directions", "Single member", "No significant support", "Exactly zero mean", "No members"),
    mean_NES = c(.4, 1.2345678901234567, NA, 0, NA),
    n_genesets = c(5L, 1L, 0L, 2L, 0L), stringsAsFactors = FALSE)
}

test_that("numeric NES inputs retain exact doubles and inclusive cutoff membership", {
  cutoff <- 0.12345678901234567
  effect <- 1.2345678901234567
  numeric_values <- c(effect, -effect, cutoff, 1e-300, NA_real_, Inf)
  expect_identical(lisaR:::lisa_nes_number(numeric_values, "NES"), numeric_values)
  # TSV fields and factors are still parsed by their labels, not factor codes.
  exact_text <- sprintf("%.17g", c(effect, cutoff))
  expect_identical(lisaR:::lisa_nes_number(exact_text, "NES"), c(effect, cutoff))
  expect_identical(lisaR:::lisa_nes_number(factor(exact_text), "NES"), c(effect, cutoff))
  expect_error(lisaR:::lisa_nes_number("1,2", "NES"), "Invalid numeric")
  members <- data.frame(category_id = "PRECISE", pathway = "001",
    NES = effect, padj = cutoff, stringsAsFactors = FALSE)
  s <- lisaR:::lisa_category_nes_statistics(members, cutoff)
  expect_identical(s$n_genesets_significant, 1L)
  expect_identical(s$mean_NES, effect)
  expect_identical(s$median_NES, effect)
  expect_identical(s$gsea_padj_cutoff, cutoff)
})

test_that("NES distributions count unique finite significant members with explicit denominators", {
  s <- lisaR:::lisa_category_nes_statistics(lisa_nes_variant_fixture(), .05,
    lisa_nes_variant_summary_fixture()$category_id)
  expect_identical(s$category_id, lisa_nes_variant_summary_fixture()$category_id)
  mix <- s[s$category_id == "MIX", ]
  expect_equal(mix$n_genesets_mapped, 5)
  expect_equal(mix$n_genesets_significant, 5)
  expect_equal(mix$mean_NES, .4)
  expect_equal(mix$median_NES, 0)
  expect_equal(c(mix$p25_NES, mix$p75_NES), c(-1, 1))
  expect_equal(c(mix$positive_pct, mix$negative_pct, mix$zero_pct), c(40, 40, 20))
  expect_equal(mix$same_direction_pct, 40)
  expect_equal(mix$quantile_type, 7)
  expect_false("OTHER_UNCLASSIFIED" %in% s$category_id)
  expect_true(all(s$positive_pct + s$negative_pct + s$zero_pct == 100, na.rm = TRUE))
  # Cutoff is inclusive and does not use precomputed significant flags.
  lower <- lisaR:::lisa_category_nes_statistics(lisa_nes_variant_fixture(), .04, "MIX")
  expect_equal(lower$n_genesets_significant, 4)
  expect_equal(lower$mean_NES, -.5)
})

test_that("NES distributions preserve n=0/1, unavailability, and exactly zero means", {
  s <- lisaR:::lisa_category_nes_statistics(lisa_nes_variant_fixture(), .05,
    lisa_nes_variant_summary_fixture()$category_id)
  none <- s[s$category_id == "NONE", ]
  expect_equal(none$n_genesets_mapped, 2)
  expect_equal(none$n_genesets_evaluable, 1)
  expect_equal(none$n_genesets_unavailable, 1)
  expect_equal(none$n_genesets_significant, 0)
  expect_true(all(is.na(none[c("mean_NES", "median_NES", "p25_NES", "p75_NES", "positive_pct", "same_direction_pct")])))
  one <- s[s$category_id == "ONE", ]
  expect_equal(unname(unlist(one[c("mean_NES", "median_NES", "p25_NES", "p75_NES")])),
    rep(1.2345678901234567, 4))
  zero <- s[s$category_id == "ZERO", ]
  expect_equal(zero$mean_NES, 0)
  expect_true(is.na(zero$same_direction_pct))
  expect_equal(zero$mean_NES_direction, "Zero mean NES")
  empty <- s[s$category_id == "EMPTY", ]
  expect_equal(empty$n_genesets_mapped, 0)
  expect_true(is.na(empty$mean_NES))
  z <- data.frame(category_id = "ALL_ZERO", pathway = "zero", NES = 0, padj = 0)
  only_zero <- lisaR:::lisa_category_nes_statistics(z, 0)
  expect_equal(only_zero$zero_pct, 100)
  expect_true(is.na(only_zero$same_direction_pct))
})

test_that("NES statistics reject conflicting duplicate results and malformed identity", {
  x <- lisa_nes_variant_fixture()
  x$NES[nrow(x)] <- 100
  expect_error(lisaR:::lisa_category_nes_statistics(x, .05), "Conflicting duplicate")
  x <- lisa_nes_variant_fixture(); x$padj[[1L]] <- -1
  expect_error(lisaR:::lisa_category_nes_statistics(x, .05), "\\[0,1\\]")
  x <- lisa_nes_variant_fixture(); x$NES[[1L]] <- "not-a-number"
  expect_error(lisaR:::lisa_category_nes_statistics(x, .05), "Invalid numeric")
  expect_error(lisaR:::lisa_category_nes_statistics(lisa_nes_variant_fixture(), .05,
    c("MIX", "MIX")), "unique classified")
  expect_error(lisaR:::lisa_category_nes_statistics(lisa_nes_variant_fixture(), .05,
    "OTHER_UNCLASSIFIED"), "unique classified")
  x <- lisa_nes_variant_fixture(); names(x)[names(x) == "padj"] <- "gsea_fdr"
  expect_equal(lisaR:::lisa_category_nes_statistics(x, .05, "MIX")$mean_NES, .4)
})

test_that("preparing plot variants validates but never substitutes the canonical mean", {
  summary <- lisa_nes_variant_summary_fixture()
  before <- summary
  v <- lisaR:::lisa_prepare_category_nes_variants(summary, lisa_nes_variant_fixture(), .05)
  expect_identical(summary, before)
  expect_identical(v$source$mean_NES, summary$mean_NES)
  expect_identical(v$source$category_id, summary$category_id)
  expect_identical(v$source$source_row_order, seq_len(nrow(summary)))
  bad <- summary; bad$mean_NES[[1L]] <- .5
  expect_error(lisaR:::lisa_prepare_category_nes_variants(bad, lisa_nes_variant_fixture(), .05), "mean_NES disagrees")
  bad <- summary; bad$n_genesets[[1L]] <- 6
  expect_error(lisaR:::lisa_prepare_category_nes_variants(bad, lisa_nes_variant_fixture(), .05), "support count")
  bad <- summary; bad$mean_NES[[3L]] <- 0
  expect_error(lisaR:::lisa_prepare_category_nes_variants(bad, lisa_nes_variant_fixture(), .05), "mean_NES disagrees")
})

test_that("four variants share scale, order and mean geometry", {
  skip_if_not_installed("ggplot2")
  v <- lisaR:::lisa_prepare_category_nes_variants(lisa_nes_variant_summary_fixture(), lisa_nes_variant_fixture(), .05)
  variants <- c("clean", "percentages", "direction", "dispersion")
  plots <- lapply(variants, function(variant) lisaR:::plot_lisa_category_nes_variant(v, variant))
  sources <- lapply(plots, attr, which = "lisa_nes_source")
  expect_identical(sources[[1]], sources[[2]])
  expect_identical(sources[[2]], sources[[3]])
  expect_identical(sources[[3]], sources[[4]])
  built <- lapply(plots, ggplot2::ggplot_build)
  xrange <- lapply(built, function(p) p$layout$panel_params[[1]]$x.range)
  expect_equal(xrange[[1]], xrange[[2]])
  expect_equal(xrange[[2]], xrange[[3]])
  yrange <- lapply(built, function(p) p$layout$panel_params[[1]]$y.range)
  expect_equal(yrange[[1]], yrange[[3]])
  expect_identical(v$source$display_color[v$source$category_id == "MIX"], "#b2182b")
  expect_identical(v$source$display_color[v$source$category_id == "ZERO"], "#687787")
  expect_equal(length(plots[[2]]$layers), length(plots[[1]]$layers) + 1)
  expect_equal(length(plots[[3]]$layers), length(plots[[2]]$layers) + 1)
  expect_equal(length(plots[[4]]$layers), length(plots[[3]]$layers) + 1)
  expect_match(plots[[4]]$labels$caption, "not[[:space:]]+confidence")
})

test_that("contrast variants keep independent supports and the unchanged existing delta", {
  a <- lisa_nes_variant_summary_fixture()
  b <- a
  b$mean_NES <- -b$mean_NES
  members_b <- lisa_nes_variant_fixture(); members_b$NES <- -members_b$NES
  contrast <- data.frame(category_id = a$category_id,
    display_mean_NES_A = c(a$mean_NES[1:2], .3, 0, NA),
    display_mean_NES_B = c(b$mean_NES[1:2], -.2, 0, NA),
    delta_mean_NES = c(.8, 2 * a$mean_NES[[2]], .5, 0, NA))
  before <- contrast
  v <- lisaR:::lisa_prepare_contrast_nes_variants(contrast, a, b,
    lisa_nes_variant_fixture(), members_b, .05)
  expect_identical(contrast, before)
  expect_identical(v$source$delta_mean_NES[1:5], contrast$delta_mean_NES)
  expect_false(any(grepl("delta.*(p25|p75|quantile)", names(v$source))))
  sa <- v$source[v$source$side == "A", ]; sb <- v$source[v$source$side == "B", ]
  expect_equal(sa$p25_NES, -sb$p75_NES)
  expect_equal(sa$positive_pct, sb$negative_pct)
  expect_identical(sa$endpoint_state[[3L]], "contextual_not_significant")
  expect_true(is.na(sa$median_NES[[3L]]))
  expect_true(is.na(sb$p25_NES[[3L]]))
  expect_no_warning(ggplot2::ggplot_build(lisaR:::plot_lisa_category_nes_variant(v, "dispersion")))
  # JSON serialization changes the representation of named side labels. The
  # legend order must not depend on ggplot's implicit guide hash ordering.
  skip_if_not_installed("png")
  folder <- tempfile("paired-nes-recipe-")
  on.exit(unlink(folder, recursive = TRUE), add = TRUE)
  files <- lisaR:::render_lisa_category_nes_variants(v, folder, formats = "png",
    width = 10, height = 5, dpi = 90)
  for (i in seq_len(nrow(files))) {
    target <- file.path(folder, paste0("replayed-", files$variant[[i]], ".png"))
    lisaR:::lisa_reproduce_category_nes_variant(files$source[[i]], files$settings[[i]], target)
    expect_identical(png::readPNG(files$path[[i]]), png::readPNG(target))
  }
})

test_that("default selection produces four variants, explicit subset produces only its files", {
  skip_if_not_installed("ggplot2")
  expect_identical(lisaR:::lisa_nes_variants_validate(), c("clean", "percentages", "direction", "dispersion"))
  for (bad in list(character(), c("clean", "clean"), "all", NA_character_))
    expect_error(lisaR:::lisa_nes_variants_validate(bad), "nonempty unique")
  v <- lisaR:::lisa_prepare_category_nes_variants(lisa_nes_variant_summary_fixture(), lisa_nes_variant_fixture(), .05)
  dir <- tempfile("nes-selected-")
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  out <- lisaR:::render_lisa_category_nes_variants(v, dir, variants = "direction",
    formats = "png", width = 10, height = 5, dpi = 90)
  expect_equal(nrow(out), 1)
  expect_identical(out$variant, "direction")
  expect_true(all(file.exists(c(out$path, out$source, out$settings))))
  expect_length(list.files(dir, pattern = "clean|dispersion"), 0)
  settings <- jsonlite::read_json(out$settings[[1]], simplifyVector = TRUE)
  expect_identical(settings$variant, "direction")
  expect_identical(settings$presentation_version, "2.0")
  expect_identical(settings$render_context$profile, "fresh_r_process_v1")
  expect_equal(settings$x_limits, v$metadata$x_limits)
  expect_match(settings$quantiles, "type 7")
})

test_that("exact-source recipe reproduces pixels and rejects content changes", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("png")
  summary <- lisa_nes_variant_summary_fixture()
  summary$category_id[[2L]] <- "001"
  x <- lisa_nes_variant_fixture(); x$category_id[x$category_id == "ONE"] <- "001"
  v <- lisaR:::lisa_prepare_category_nes_variants(summary, x, .05)
  dir <- tempfile("nes-recipe-"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  out <- lisaR:::render_lisa_category_nes_variants(v, dir, variants = "dispersion",
    formats = "png", width = 10, height = 5, dpi = 90)
  target <- file.path(dir, "reproduced.png")
  expect_no_warning(lisaR:::lisa_reproduce_category_nes_variant(out$source[[1]], out$settings[[1]], target))
  expect_identical(png::readPNG(out$path[[1]]), png::readPNG(target))
  stored <- utils::read.delim(out$source[[1]], colClasses = "character", quote = '"', check.names = FALSE)
  expect_true("001" %in% stored$category_id)
  expect_identical(as.numeric(stored$mean_NES[[2]]), summary$mean_NES[[2]])
  cat("\n", file = out$source[[1]], append = TRUE)
  expect_error(lisaR:::lisa_reproduce_category_nes_variant(out$source[[1]], out$settings[[1]], target), "hash mismatch")
})

test_that("empty classified universes remain drawable without fake support", {
  skip_if_not_installed("ggplot2")
  x <- lisa_nes_variant_fixture()[FALSE, , drop = FALSE]
  s <- lisa_nes_variant_summary_fixture()[FALSE, , drop = FALSE]
  v <- lisaR:::lisa_prepare_category_nes_variants(s, x, .05)
  expect_equal(nrow(v$source), 0)
  expect_identical(v$source$macrogroup_name, character())
  expect_identical(lisaR:::lisa_nes_direction_label(v$source), character())
  for (variant in c("clean", "percentages", "direction", "dispersion")) {
    expect_no_warning(ggplot2::ggplot_build(lisaR:::plot_lisa_category_nes_variant(v, variant)))
  }
  # Empty contrast summaries need a stable paired schema as well.
  contrast <- data.frame(category_id = character(), display_mean_NES_A = numeric(),
    display_mean_NES_B = numeric(), delta_mean_NES = numeric())
  paired <- lisaR:::lisa_prepare_contrast_nes_variants(contrast, s, s, x, x, .05)
  expect_equal(nrow(paired$source), 0)
  expect_identical(paired$metadata$scope, "contrast")
  for (variant in c("clean", "percentages", "direction", "dispersion")) {
    p <- lisaR:::plot_lisa_category_nes_variant(paired, variant)
    expect_no_warning(ggplot2::ggplot_build(p))
    expect_null(p$scales$get_scales("shape"))
  }
})

test_that("categories without significant support have no recycled direction annotation", {
  skip_if_not_installed("ggplot2")
  summary <- lisa_nes_variant_summary_fixture()
  summary <- summary[summary$category_id %in% c("NONE", "EMPTY"), , drop = FALSE]
  v <- lisaR:::lisa_prepare_category_nes_variants(summary, lisa_nes_variant_fixture(), .05)
  expect_identical(v$source$n_genesets_significant, c(0L, 0L))
  expect_true(all(is.na(v$source$median_NES)))
  expect_true(all(is.na(v$source$positive_pct)))
  p <- lisaR:::plot_lisa_category_nes_variant(v, "direction")
  expect_no_warning(ggplot2::ggplot_build(p))
  # A saved empty/supported-zero figure remains an executable ordinary recipe.
  out_dir <- tempfile("nes-no-support-")
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  files <- lisaR:::render_lisa_category_nes_variants(v, out_dir,
    variants = "direction", formats = "png", width = 10, height = 5, dpi = 72)
  reproduced <- file.path(out_dir, "reproduced.png")
  expect_no_warning(lisaR:::lisa_reproduce_category_nes_variant(
    files$source[[1L]], files$settings[[1L]], reproduced))
  expect_true(file.exists(reproduced))
})


test_that("direction labels name counts and percentages without false negative signs", {
  x <- data.frame(n_genesets_significant = c(3L, 0L), positive_pct = c(100, NA),
    negative_pct = c(0, NA), n_zero_genesets = c(0L, 0L), zero_pct = c(0, NA),
    same_direction_pct = c(100, NA))
  before <- x
  labels <- lisaR:::lisa_nes_direction_label(x)
  expect_identical(labels[[1L]], "n=3 | Pos 100% / Neg 0% | Agree 100%")
  expect_match(labels[[2L]], "Agree NA", fixed = TRUE)
  expect_identical(x, before)
  expect_identical(lisaR:::lisa_nes_direction_label(x[FALSE, ]), character())
})

test_that("count size legends have integer breaks without changing size mapping", {
  v <- lisaR:::lisa_prepare_category_nes_variants(lisa_nes_variant_summary_fixture(), lisa_nes_variant_fixture(), .05)
  p <- lisaR:::plot_lisa_category_nes_variant(v, "dispersion")
  scale <- p$scales$get_scales("size")
  expect_identical(scale$breaks(c(2, 3)), c(2, 3))
  expect_identical(scale$breaks(c(3, 3)), 3)
  expect_length(scale$breaks(c(NA_real_, NA_real_)), 0L)
  b <- ggplot2::ggplot_build(p)
  trained <- b$plot$scales$get_scales("size")
  breaks <- trained$get_breaks()
  expect_true(all(breaks == floor(breaks)))
  expect_identical(trained$get_labels(breaks), format(breaks, scientific = FALSE, trim = TRUE))
  expect_identical(attr(p, "lisa_nes_source")$n_genesets_significant, v$source$n_genesets_significant)
})

test_that("supercategory bands and canonical colours survive every NES view", {
  skip_if_not_installed("ggplot2")
  summary <- lisa_nes_variant_summary_fixture()
  summary$macrogroup_name <- c("Immune biology", "Immune biology", "Other programmes", "Other programmes", "Other programmes")
  summary$color <- c("#00AABB", "#226677", "#AB3322", "#882211", "#662211")
  prepared <- lisaR:::lisa_prepare_category_nes_variants(summary, lisa_nes_variant_fixture(), .05)
  for (variant in c("clean", "percentages", "direction", "dispersion")) {
    p <- lisaR:::plot_lisa_category_nes_variant(prepared, variant)
    expect_s3_class(p$facet, "FacetGrid")
    b <- ggplot2::ggplot_build(p)
    expect_identical(as.character(b$layout$layout$macrogroup_name), unique(summary$macrogroup_name))
    expect_identical(attr(p, "lisa_nes_source")$category_color, summary$color)
    expect_identical(attr(p, "lisa_nes_source")$category_id, summary$category_id)
    # Band labels are horizontal on the left, with no new right-hand column.
    expect_identical(p$theme$strip.text.y.left$angle, 0)
    expect_lt(as.numeric(p$theme$plot.margin)[[2L]], 30)
    if (variant != "clean") {
      text_layers <- Filter(function(layer) inherits(layer$geom, "GeomText"), p$layers)
      expect_length(text_layers, 1L)
      rows <- text_layers[[1L]]$data
      expect_true(all(abs(rows$label_x - rows$display_mean_NES) < diff(prepared$metadata$x_limits) * .05))
      expect_false(any(grepl("Pos |Neg |Agree ", rows$point_label)))
      expect_true(any(grepl("Direction undefined", rows$point_label)))
    }
  }
})

test_that("contrast subsets preserve canonical membership, endpoints and delta", {
  a <- lisa_nes_variant_summary_fixture(); b <- a
  b$mean_NES[[2L]] <- -b$mean_NES[[2L]]
  members_b <- lisa_nes_variant_fixture()
  members_b$NES[members_b$category_id == "ONE"] <- -members_b$NES[members_b$category_id == "ONE"]
  contrast <- data.frame(category_id = a$category_id,
    display_mean_NES_A = c(a$mean_NES[1:2], .3, 0, NA),
    display_mean_NES_B = c(b$mean_NES[1:2], -.2, 0, NA),
    delta_mean_NES = c(0, 2 * a$mean_NES[[2L]], .5, 0, NA),
    is_same_direction = c(TRUE, FALSE, FALSE, FALSE, FALSE),
    is_opposite_direction = c(FALSE, TRUE, FALSE, FALSE, FALSE))
  prepared <- lisaR:::lisa_prepare_contrast_nes_variants(contrast, a, b,
    lisa_nes_variant_fixture(), members_b, .05)
  before <- prepared
  same <- lisaR:::lisa_nes_select_plot_set(prepared, "same_direction")
  opposite <- lisaR:::lisa_nes_select_plot_set(prepared, "opposite_direction")
  expect_identical(same$source$category_id, c("MIX", "MIX"))
  expect_identical(opposite$source$category_id, c("ONE", "ONE"))
  expect_identical(same$source$delta_mean_NES, c(0, 0))
  expect_identical(opposite$source$delta_mean_NES, rep(contrast$delta_mean_NES[[2L]], 2))
  expect_identical(same$metadata$x_limits, prepared$metadata$x_limits)
  expect_identical(opposite$metadata$x_limits, prepared$metadata$x_limits)
  expect_identical(prepared, before)
  missing <- prepared; missing$source$is_same_direction <- NULL
  expect_error(lisaR:::lisa_nes_select_plot_set(missing, "same_direction"), "canonical contrast membership")
  expect_error(lisaR:::lisa_nes_plot_sets_validate(character(), "contrast"), "Invalid nonempty")
  expect_error(lisaR:::lisa_nes_plot_sets_validate("opposite_direction", "single"), "Invalid nonempty")
  skip_if_not_installed("ggplot2")
  directory <- tempfile("nes-subsets-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  files <- lisaR:::render_lisa_category_nes_variants(prepared, directory, variants = "percentages",
    formats = "png", width = 10, height = 5, dpi = 72,
    plot_sets = c("all", "same_direction", "opposite_direction"))
  expect_identical(files$plot_set, c("all", "same_direction", "opposite_direction"))
  expect_length(unique(files$path), 3L)
  expect_identical(basename(files$path[[1L]]), "category_nes_percentages.png")
  expect_identical(basename(files$path[[2L]]), "category_nes_same_direction_percentages.png")
  for (i in seq_len(nrow(files))) {
    settings <- jsonlite::read_json(files$settings[[i]], simplifyVector = TRUE)
    expect_identical(settings$plot_set, files$plot_set[[i]])
    expect_identical(settings$source_sha256, digest::digest(file = files$source[[i]], algo = "sha256"))
  }
})

test_that("archived settings keep the legacy rendering meaning", {
  skip_if_not_installed("ggplot2")
  old <- lisaR:::lisa_prepare_category_nes_variants(lisa_nes_variant_summary_fixture(), lisa_nes_variant_fixture(), .05)
  old$metadata$presentation_version <- NULL; old$metadata$render_context <- NULL
  old$metadata$x_limits <- c(-4.48, 4.48)
  for (variant in c("clean", "direction", "dispersion")) {
    p <- lisaR:::plot_lisa_category_nes_variant(old, variant)
    expected <- lisaR:::plot_lisa_category_nes_variant_legacy(old, variant)
    expect_identical(ggplot2::ggplot_build(p)$data, ggplot2::ggplot_build(expected)$data)
  }
  expect_error(lisaR:::plot_lisa_category_nes_variant(old, "percentages"), "Legacy NES settings")
})

test_that("fresh-device origin and replay ignore preceding 180 and 300 dpi plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("png")
  skip_if_not(capabilities("cairo"))
  s <- lisa_nes_variant_summary_fixture()
  s$category_display_name <- paste(s$category_display_name, "with a deliberately long scientific label")
  s$macrogroup_name <- c("First long supercategory", "First long supercategory", rep("Second supercategory", 3))
  a <- s; b <- s; b$mean_NES <- -b$mean_NES
  members_b <- lisa_nes_variant_fixture(); members_b$NES <- -members_b$NES
  contrast <- data.frame(category_id = a$category_id,
    display_mean_NES_A = c(a$mean_NES[1:2], .3, 0, NA),
    display_mean_NES_B = c(b$mean_NES[1:2], -.2, 0, NA),
    delta_mean_NES = c(.8, 2 * a$mean_NES[[2L]], .5, 0, NA))
  prepared <- lisaR:::lisa_prepare_contrast_nes_variants(contrast, a, b,
    lisa_nes_variant_fixture(), members_b, .05,
    label_a = "Responders: ON-treatment minus PRE-treatment",
    label_b = "Progressors: ON-treatment minus PRE-treatment",
    subtitle = "Two independently adjusted analyses\nPositive delta indicates a larger change in A.")
  directory <- tempfile("nes-device-history-"); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  previous <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_text(label = "Text-rich previous drawing") + ggplot2::theme_minimal()
  for (history in c(180, 300)) {
    ggplot2::ggsave(file.path(directory, paste0("previous-", history, ".png")), previous,
      width = 4, height = 3, dpi = history)
    files <- lisaR:::render_lisa_category_nes_variants(prepared, file.path(directory, paste0("origin-", history)),
      formats = "png", width = 13, height = 7, dpi = 180)
    for (i in seq_len(nrow(files))) {
      target <- file.path(directory, paste0("replay-", history, "-", files$variant[[i]], ".png"))
      lisaR:::lisa_reproduce_category_nes_variant(files$source[[i]], files$settings[[i]], target)
      expect_identical(png::readPNG(files$path[[i]]), png::readPNG(target))
      if (history == 300) expect_identical(png::readPNG(target),
        png::readPNG(file.path(directory, paste0("replay-180-", files$variant[[i]], ".png"))))
    }
  }
})

test_that("median and percentile lanes do not share the mean or percentage rows", {
  skip_if_not_installed("ggplot2")
  summary <- lisa_nes_variant_summary_fixture()
  summary$macrogroup_name <- c("First", "First", rep("Second", 3L))
  single <- lisaR:::lisa_prepare_category_nes_variants(summary, lisa_nes_variant_fixture(), .05)
  opposite <- summary; opposite$mean_NES <- -opposite$mean_NES
  members_b <- lisa_nes_variant_fixture(); members_b$NES <- -members_b$NES
  canonical <- data.frame(category_id = summary$category_id,
    display_mean_NES_A = c(summary$mean_NES[1:2], .3, 0, NA),
    display_mean_NES_B = c(opposite$mean_NES[1:2], -.2, 0, NA),
    delta_mean_NES = c(.8, 2 * summary$mean_NES[[2L]], .5, 0, NA))
  paired <- lisaR:::lisa_prepare_contrast_nes_variants(canonical, summary, opposite,
    lisa_nes_variant_fixture(), members_b, .05)
  key <- function(rows) paste(rows$category_id, rows$side, sep = "\r")
  mapping_is <- function(layer, field) identical(rlang::get_expr(layer$mapping$x), as.name(field))
  for (prepared in list(single, paired)) {
    before <- prepared$source
    expect_identical(prepared$metadata$distribution_layout, "separate_lanes_v1")
    mean_reference <- NULL
    for (variant in c("clean", "percentages", "direction", "dispersion")) {
      plot <- lisaR:::plot_lisa_category_nes_variant(prepared, variant)
      built <- ggplot2::ggplot_build(plot)
      mean_indices <- which(vapply(plot$layers, function(layer)
        inherits(layer$geom, "GeomPoint") && mapping_is(layer, "display_mean_NES"), logical(1L)))
      means <- do.call(rbind, lapply(mean_indices, function(i) data.frame(
        key = key(plot$layers[[i]]$data), x = built$data[[i]]$x,
        y = as.numeric(built$data[[i]]$y), stringsAsFactors = FALSE)))
      rownames(means) <- NULL
      if (is.null(mean_reference)) mean_reference <- means else expect_identical(means, mean_reference)
      if (variant %in% c("direction", "dispersion")) {
        text_indices <- which(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1L)))
        labels <- do.call(rbind, lapply(text_indices, function(i) data.frame(
          key = key(plot$layers[[i]]$data), y = as.numeric(built$data[[i]]$y), stringsAsFactors = FALSE)))
        median_indices <- which(vapply(plot$layers, function(layer)
          inherits(layer$geom, "GeomPoint") && mapping_is(layer, "median_NES"), logical(1L)))
        for (i in median_indices) {
          rows <- plot$layers[[i]]$data
          yy <- as.numeric(built$data[[i]]$y)
          base_y <- means$y[match(key(rows), means$key)]
          label_y <- labels$y[match(key(rows), labels$key)]
          expected <- ifelse(rows$side == "B", -.4, .4)
          expect_equal(yy - base_y, expected)
          expect_true(all(abs(yy - label_y) >= .19))
          expect_equal(built$data[[i]]$x, rows$median_NES)
        }
        if (variant == "dispersion") {
          interval_indices <- which(vapply(plot$layers, function(layer)
            inherits(layer$geom, "GeomSegment") && mapping_is(layer, "p25_NES"), logical(1L)))
          expect_length(interval_indices, length(median_indices))
          for (i in interval_indices) {
            rows <- plot$layers[[i]]$data
            base_y <- means$y[match(key(rows), means$key)]
            expect_equal(as.numeric(built$data[[i]]$y) - base_y, ifelse(rows$side == "B", -.4, .4))
            expect_equal(built$data[[i]]$x, rows$p25_NES)
            expect_equal(built$data[[i]]$xend, rows$p75_NES)
          }
        }
      }
      expect_identical(prepared$source, before)
      expect_identical(plot$scales$get_scales("x")$limits, prepared$metadata$x_limits)
    }
    # Already emitted pre-correction 2.0 recipes retain their recorded layout.
    archived <- prepared; archived$metadata$distribution_layout <- NULL
    expect_equal(lisaR:::lisa_nes_distribution_lane_offset(archived$metadata, "single"), 0)
    expect_equal(lisaR:::lisa_nes_distribution_lane_offset(archived$metadata, "A"), .14)
    expect_equal(lisaR:::lisa_nes_distribution_lane_offset(archived$metadata, "B"), -.14)
    invalid <- prepared; invalid$metadata$distribution_layout <- "unrecognized"
    expect_error(lisaR:::plot_lisa_category_nes_variant(invalid, "clean"), "Unsupported NES distribution layout")
  }
})
