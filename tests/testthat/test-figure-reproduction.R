test_that("figure-specific source contracts regenerate category, volcano and GeneCard PNGs", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  renderer <- system.file("scripts", "reproduce_lisa_figure.R", package = "lisaR")
  expect_true(file.exists(renderer))

  run_recipe <- function(df, stem) {
    source <- file.path(tempdir(), paste0(stem, "_source.tsv"))
    output <- file.path(tempdir(), paste0(stem, ".png"))
    utils::write.table(df, source, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    status <- system2(lisaR:::lisa_rscript_executable(),
      c("--vanilla", shQuote(renderer), shQuote(source), shQuote(output)),
      stdout = TRUE, stderr = TRUE)
    exit_status <- attr(status, "status")
    if (is.null(exit_status)) exit_status <- 0L
    expect_equal(exit_status, 0L, info = paste(status, collapse = "\n"))
    expect_true(file.exists(output))
    expect_gt(file.info(output)$size, 1000)
    if (identical(stem, "category_contract")) {
      skip_if_not_installed("png")
      expect_equal(dim(png::readPNG(output))[1:2], c(1350L, 3750L))
    }
  }

  category <- data.frame(
    figure_id = "cat_1", figure_type = "lisa_category_gene_sets", source_row_order = 1:2,
    selected_for_plot = TRUE, highlighted = TRUE, labelled = TRUE,
    category_id = "IMMUNE", category_display_name = "Immune signalling",
    macrogroup_name = "Immunity", category_color = "#336699", plot_subtitle = "fixture",
    pathway = c("SET_A", "SET_B"), pathway_label = c("Set A", "Set B"),
    NES = c(-1.4, 1.8), padj = c(0.01, 0.02), plot_x_nes = c(-1.4, 1.8),
    plot_point_size_neg_log10_fdr = -log10(c(0.01, 0.02)), plot_x_min = -2, plot_x_max = 2,
    plotted_order = 1:2, figure_width = 12.5, figure_height = 4.5,
    figure_dpi = 300, stringsAsFactors = FALSE)
  run_recipe(category, "category_contract")

  volcano <- data.frame(
    figure_id = "volcano_1", figure_type = "volcano_overlay", source_row_order = 1:4,
    selected_for_plot = TRUE, highlighted = c(TRUE, TRUE, FALSE, FALSE), labelled = c(TRUE, TRUE, FALSE, FALSE),
    category_id = "IMMUNE", category_display_name = "Immune signalling", macrogroup_name = "Immunity",
    category_color = "#336699", plot_subtitle = "fixture", symbol = c("A", "B", "C", "D"),
    log2FC = c(-1.2, 1.4, 0.1, -0.2), padj = c(0.01, 0.02, 0.8, 0.6),
    plot_x_log2FC = c(-1.2, 1.4, 0.1, -0.2), plot_y_neg_log10_fdr = -log10(c(0.01, 0.02, 0.8, 0.6)),
    de_class = c("DE down", "DE up", "not DE", "not DE"), threshold_abs_log2FC = 0.58,
    threshold_de_fdr = 0.05, stringsAsFactors = FALSE)
  run_recipe(volcano, "volcano_contract")

  card <- data.frame(
    figure_id = "card_1", figure_type = "gene_card", source_row_order = 1:3,
    selected_for_plot = TRUE, highlighted = TRUE, labelled = TRUE, plotted_order = 1:3,
    category_id = "IMMUNE", category_display_name = "Immune signalling", macrogroup_name = "Immunity",
    category_color = "#336699", color = "#336699", symbol = c("A", "B", "C"),
    log2FC = c(-1.2, 1.4, 0.3), rank_value = c(-3.4, 2.9, 0.5),
    padj = c(0.01, 0.02, 0.2), gene_contribution_score = c(2.1, 1.8, 0.8),
    de_is_significant = c(TRUE, TRUE, FALSE), in_gsea_leading_edge = c(TRUE, TRUE, FALSE),
    in_ora_overlap = c(FALSE, TRUE, TRUE), category_recurrence = c(3, 2, 1),
    source_geneset_count = c(4, 3, 1), stringsAsFactors = FALSE)
  run_recipe(card, "card_contract")
})

test_that("KEGG source contract regenerates without network access", {
  skip_if_not_installed("png")
  renderer <- system.file("scripts", "reproduce_lisa_figure.R", package = "lisaR")
  base <- file.path(tempdir(), "kegg_fixture_base.png")
  png::writePNG(array(1, dim = c(120, 180, 3)), base)
  source <- file.path(tempdir(), "kegg_fixture_source.tsv")
  output <- file.path(tempdir(), "kegg_fixture.png")
  x <- data.frame(
    figure_id = "kegg_1", figure_type = "kegg_single", source_row_order = 1,
    selected_for_plot = TRUE, highlighted = TRUE, labelled = FALSE,
    base_image_file = basename(base), x = 90, y = 60, width = 50, height = 24,
    log2FC = 1.1, max_abs_log2fc = 1.5, color_power = 1,
    figure_title = "mmu00000 fixture", figure_subtitle = "offline fixture", stringsAsFactors = FALSE)
  utils::write.table(x, source, sep = "\t", quote = FALSE, row.names = FALSE)
  status <- system2(lisaR:::lisa_rscript_executable(),
    c("--vanilla", shQuote(renderer), shQuote(source), shQuote(output)), stdout = TRUE, stderr = TRUE)
  exit_status <- attr(status, "status")
  if (is.null(exit_status)) exit_status <- 0L
  expect_equal(exit_status, 0L, info = paste(status, collapse = "\n"))
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 1000)
})

test_that("selected contrast heatmap recipe retains the native size and labelled axis", {
  skip_if_not_installed("png")
  source <- tempfile(fileext = ".tsv"); output <- tempfile(fileext = ".png")
  canonical <- tempfile(fileext = ".png")
  x <- data.frame(figure_type = "lisa_contrast_heatmap", figure_title = "Signal category",
    category_display_name = "Signal category", contrast_a_label = "ON", contrast_b_label = "PRE",
    display_mean_NES_A = 1.6, display_mean_NES_B = -.5,
    mean_NES_A = 1.6, mean_NES_B = -.5,
    has_significant_support_A = TRUE, has_significant_support_B = FALSE,
    plot_has_any_significant_support = TRUE, gsea_padj_cutoff = .25,
    figure_width = 5.5, figure_height = 3.8, figure_dpi = 300)
  lisaR:::write_lisa_tsv(x, source)
  # Both the real renderer and the recipe run in fresh R processes. Use that
  # same boundary for the reference: the test session's graphics device and
  # font state must not become part of the expected figure contract.
  reference <- quote({
  values <- data.frame(side = factor(c("ON", "PRE"), levels = c("ON", "PRE")),
    mean_NES = c(1.6, -.5), has_significant_support = c(TRUE, FALSE))
  p <- ggplot2::ggplot(values, ggplot2::aes(x = side, y = "Signal category", fill = mean_NES,
    alpha = has_significant_support)) + ggplot2::geom_tile(color = "white", linewidth = .5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f%s", mean_NES,
      ifelse(has_significant_support, "", " (NS)"))), size = 4) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = .35), guide = "none") +
    ggplot2::labs(title = "Signal category", x = NULL, y = NULL, fill = "Mean NES") +
    ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(panel.grid = ggplot2::element_blank())
  ggplot2::ggsave(canonical, p, width = 5.5, height = 3.8, dpi = 300, bg = "white")
  })
  reference_script <- tempfile(fileext = ".R")
  writeLines(c(paste0("canonical <- ", encodeString(canonical, quote = '"')),
    deparse(reference)), reference_script)
  reference_status <- system2(lisaR:::lisa_rscript_executable(),
    c("--vanilla", shQuote(reference_script)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(reference_status, "status"),
    info = paste(reference_status, collapse = "\n"))
  status <- system2(lisaR:::lisa_rscript_executable(), c("--vanilla",
    shQuote(system.file("scripts", "reproduce_lisa_figure.R", package = "lisaR")),
    shQuote(source), shQuote(output)), stdout = TRUE, stderr = TRUE)
  expect_null(attr(status, "status"), info = paste(status, collapse = "\n"))
  a <- png::readPNG(canonical); b <- png::readPNG(output)
  expect_identical(dim(b), c(1140L, 1650L, 3L))
  expect_true(identical(a, b), info = paste("Mean absolute pixel difference:", mean(abs(a - b))))
})
