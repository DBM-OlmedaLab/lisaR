lisa_ab_category_map <- function() {
  data.frame(
    category_id = c("BOTH_SIG", "A_ONLY", "B_ONLY", "NEITHER", "SUPPORTED_ZERO"),
    display_name = c("Both significant", "A only", "B only", "Neither", "Supported zero"),
    macrogroup_id = "TEST",
    macrogroup_name = "Test categories",
    macrogroup_order = 1,
    category_order_within_macrogroup = 1:5,
    color = c("#0072B2", "#009E73", "#D55E00", "#737373", "#CC79A7"),
    stringsAsFactors = FALSE
  )
}

lisa_ab_gsea_fixture <- function(side) {
  if (identical(side, "A")) {
    data.frame(
      category_id = c("BOTH_SIG", "A_ONLY", "B_ONLY", "NEITHER", "SUPPORTED_ZERO", "SUPPORTED_ZERO"),
      pathway = c("A_BOTH", "A_ONLY_SIG", "A_BONLY_NS", "A_NEITHER_NS", "A_ZERO_POS", "A_ZERO_NEG"),
      NES = c(2.0, 1.5, -0.4, 0.6, 1.0, -1.0),
      padj = c(0.05, 0.01, 0.20, 0.20, 0.02, 0.03),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      category_id = c("BOTH_SIG", "A_ONLY", "B_ONLY", "NEITHER", "SUPPORTED_ZERO"),
      pathway = c("B_BOTH", "B_AONLY_NS", "B_ONLY_SIG", "B_NEITHER_NS", "B_ZERO_OTHER"),
      NES = c(1.0, 0.4, -0.8, -0.6, 0.5),
      padj = c(0.04, 0.20, 0.03, 0.40, 0.01),
      stringsAsFactors = FALSE
    )
  }
}

lisa_ab_summary_fixture <- function(side, cutoff = 0.05) {
  lisaR:::build_lisa_summaries(
    gsea_all = lisa_ab_gsea_fixture(side),
    ora_all = data.frame(),
    category_map = lisa_ab_category_map(),
    include_empty_categories = TRUE,
    gsea_padj_cutoff = cutoff,
    ora_padj_cutoff = 0.05
  )$gsea
}

test_that("single-DE summaries separate significant and contextual category support", {
  a <- lisa_ab_summary_fixture("A")
  b <- lisa_ab_summary_fixture("B")

  expect_identical(as.character(a$category_id), lisa_ab_category_map()$category_id)
  expect_true(all(a$n_genesets == a$n_genesets_significant))
  expect_true(all(a$gsea_padj_cutoff == 0.05))

  both_a <- a[a$category_id == "BOTH_SIG", ]
  expect_identical(both_a$n_genesets_significant, 1L)
  expect_true(both_a$has_significant_support)
  expect_equal(both_a$mean_NES, 2)

  a_only_b <- b[b$category_id == "A_ONLY", ]
  expect_identical(a_only_b$n_genesets_mapped, 1L)
  expect_identical(a_only_b$n_genesets_evaluable, 1L)
  expect_identical(a_only_b$n_genesets_significant, 0L)
  expect_false(a_only_b$has_significant_support)
  expect_true(is.na(a_only_b$mean_NES))
  expect_equal(a_only_b$mean_NES_contextual, 0.4)

  neither <- a[a$category_id == "NEITHER", ]
  expect_false(neither$has_significant_support)
  expect_true(is.na(neither$mean_NES))
  expect_equal(neither$mean_NES_contextual, 0.6)

  zero <- a[a$category_id == "SUPPORTED_ZERO", ]
  expect_true(zero$has_significant_support)
  expect_identical(zero$n_genesets_significant, 2L)
  expect_equal(zero$mean_NES, 0)
})

test_that("A/B contrast retains independent member sets and approved display states", {
  a_raw <- lisa_ab_gsea_fixture("A")
  b_raw <- lisa_ab_gsea_fixture("B")
  expect_length(intersect(
    a_raw$pathway[a_raw$category_id == "BOTH_SIG"],
    b_raw$pathway[b_raw$category_id == "BOTH_SIG"]
  ), 0L)

  a <- lisa_ab_summary_fixture("A")
  b <- lisa_ab_summary_fixture("B")
  observed <- lisaR:::build_lisa_contrast_summary(
    a, b, "A", "B", "PATHWAYS", TRUE, 0,
    gsea_padj_cutoff = 0.05
  )

  expect_identical(as.character(observed$category_id), lisa_ab_category_map()$category_id)
  expect_true(all(observed$gsea_padj_cutoff == 0.05))

  a_only <- observed[observed$category_id == "A_ONLY", ]
  expect_true(a_only$has_significant_support_A)
  expect_false(a_only$has_significant_support_B)
  expect_equal(a_only$mean_NES_A, 1.5)
  expect_true(is.na(a_only$mean_NES_B))
  expect_equal(a_only$display_mean_NES_A, 1.5)
  expect_equal(a_only$display_mean_NES_B, 0.4)
  expect_identical(a_only$endpoint_source_B, "contextual_mean_not_significant")
  expect_equal(a_only$delta_mean_NES, 1.1)

  b_only <- observed[observed$category_id == "B_ONLY", ]
  expect_false(b_only$has_significant_support_A)
  expect_true(b_only$has_significant_support_B)
  expect_equal(b_only$display_mean_NES_A, -0.4)
  expect_equal(b_only$display_mean_NES_B, -0.8)

  neither <- observed[observed$category_id == "NEITHER", ]
  expect_false(neither$plot_has_any_significant_support)
  expect_true(is.na(neither$mean_NES_A))
  expect_true(is.na(neither$mean_NES_B))
  expect_true(is.na(neither$display_mean_NES_A))
  expect_true(is.na(neither$display_mean_NES_B))
  expect_true(is.na(neither$delta_mean_NES))
  expect_identical(neither$direction_class, "no_significant_support")

  zero <- observed[observed$category_id == "SUPPORTED_ZERO", ]
  expect_true(zero$has_significant_support_A)
  expect_equal(zero$display_mean_NES_A, 0)
  expect_identical(zero$direction_A, "zero")

  swapped <- lisaR:::build_lisa_contrast_summary(
    b, a, "B", "A", "PATHWAYS", TRUE, 0,
    gsea_padj_cutoff = 0.05
  )
  expect_equal(swapped$display_mean_NES_A, observed$display_mean_NES_B)
  expect_equal(swapped$display_mean_NES_B, observed$display_mean_NES_A)
  expect_equal(swapped$delta_mean_NES, -observed$delta_mean_NES)
  expect_identical(swapped$has_significant_support_A, observed$has_significant_support_B)
  expect_identical(swapped$has_significant_support_B, observed$has_significant_support_A)
})

test_that("all-category dumbbell keeps blank rows and distinguishes contextual endpoints", {
  skip_if_not_installed("ggplot2")
  observed <- lisaR:::build_lisa_contrast_summary(
    lisa_ab_summary_fixture("A"), lisa_ab_summary_fixture("B"),
    "A", "B", "PATHWAYS", TRUE, 0, gsea_padj_cutoff = 0.05
  )
  all_rows <- lisaR:::filter_lisa_contrast_plot_set(observed, "all")
  expect_equal(nrow(all_rows), 5L)

  plot <- lisaR:::plot_lisa_contrast_dumbbell(
    all_rows, "A", "B", "A", "B", "A minus B", "Fixture",
    group_by_supracategory = FALSE, annotate_gene_sets = TRUE
  )
  built <- ggplot2::ggplot_build(plot)
  expect_equal(nrow(built$data[[1]]), 5L) # blank geometry trains every row
  expect_equal(nrow(built$data[[3]]), 4L) # no segment for NEITHER
  expect_equal(nrow(built$data[[4]]), 6L) # supported endpoints, including zero
  expect_equal(nrow(built$data[[5]]), 2L) # contextual non-significant endpoints
  expect_true(any(built$data[[4]]$x == 0))
  expect_true(all(c("Both significant", "A only", "B only", "Neither", "Supported zero") %in%
    built$layout$panel_params[[1]]$y$get_labels()))
  expect_true(grepl("\n", plot$labels$caption, fixed = TRUE))
  expect_identical(plot$theme$plot.caption$hjust, 0)
})

test_that("faceted dumbbells keep each category in exactly one macrogroup panel", {
  skip_if_not_installed("ggplot2")
  observed <- lisaR:::build_lisa_contrast_summary(
    lisa_ab_summary_fixture("A"), lisa_ab_summary_fixture("B"),
    "A", "B", "PATHWAYS", TRUE, 0, gsea_padj_cutoff = 0.05
  )
  observed$macrogroup_name <- c("Group one", "Group one", "Group two", "Group two", "Group two")

  plot <- lisaR:::plot_lisa_contrast_dumbbell(
    observed, "A", "B", "A", "B", "A minus B", "Fixture",
    group_by_supracategory = TRUE, annotate_gene_sets = TRUE
  )
  built <- ggplot2::ggplot_build(plot)
  panel_breaks <- lapply(built$layout$panel_params, function(panel) {
    as.character(panel$y$get_breaks())
  })

  expect_equal(length(panel_breaks), 2L)
  expect_equal(sort(lengths(panel_breaks)), c(2L, 3L))
  expect_equal(sum(lengths(panel_breaks)), nrow(observed))
  expect_equal(length(unique(unlist(panel_breaks, use.names = FALSE))), nrow(observed))
  expect_true(any(vapply(panel_breaks, function(x) {
    any(grepl("^NEITHER__", x))
  }, logical(1))))
})

test_that("gsea_padj_cutoff is one typed inclusive config value", {
  config <- function(value = NULL) {
    pipeline <- list(
      schema_version = "1.0.0", profile = "targeted", evidence_mode = "full_de",
      duplicate_policies = list(
        de_table_duplicate_policy = "error", matrix_duplicate_policy = "error",
        mapped_id_collision_policy = "error"
      )
    )
    if (!is.null(value)) pipeline$gsea_padj_cutoff <- value
    list(pipeline = pipeline)
  }
  expect_equal(lisaR:::lisa_validate_pipeline_config(config())$gsea_padj_cutoff, 0.25)
  expect_equal(lisaR:::lisa_validate_pipeline_config(config(0))$gsea_padj_cutoff, 0)
  expect_equal(lisaR:::lisa_validate_pipeline_config(config(1))$gsea_padj_cutoff, 1)
  for (value in list("0.05", NA_real_, Inf, -0.01, 1.01, c(0.05, 0.1), TRUE)) {
    expect_error(
      lisaR:::lisa_validate_pipeline_config(config(value)),
      "LISA-CONFIG-GSEA-001"
    )
  }

  schema <- jsonlite::fromJSON(
    system.file("schema", "lisa-config.schema.json", package = "lisaR"),
    simplifyVector = FALSE
  )
  cutoff <- schema$properties$pipeline$properties$gsea_padj_cutoff
  expect_identical(cutoff$type, "number")
  expect_equal(cutoff$minimum, 0)
  expect_equal(cutoff$maximum, 1)
})

test_that("historical empty aliases are read as missing, never invented zero", {
  legacy <- data.frame(
    category_id = c("EMPTY", "SUPPORTED"),
    n_genesets = c(0, 1), mean_NES = c(0, 0),
    stringsAsFactors = FALSE
  )
  normalized <- lisaR:::normalize_lisa_gsea_summary_support(
    legacy, 0.05, "legacy fixture"
  )
  expect_true(is.na(normalized$mean_NES[normalized$category_id == "EMPTY"]))
  expect_false(normalized$has_significant_support[normalized$category_id == "EMPTY"])
  expect_equal(normalized$mean_NES[normalized$category_id == "SUPPORTED"], 0)
  expect_true(normalized$has_significant_support[normalized$category_id == "SUPPORTED"])
})

test_that("A/B comparison fails closed when contextual geometry or category maps disagree", {
  a <- lisa_ab_summary_fixture("A")
  b <- lisa_ab_summary_fixture("B")

  missing_context <- b
  row <- missing_context$category_id == "A_ONLY"
  missing_context$n_genesets_mapped[row] <- 0L
  missing_context$n_genesets_evaluable[row] <- 0L
  missing_context$mean_NES_contextual[row] <- NA_real_
  expect_error(
    lisaR:::build_lisa_contrast_summary(
      a, missing_context, "A", "B", "PATHWAYS", TRUE, 0,
      gsea_padj_cutoff = 0.05
    ),
    "LISA-CONTRAST-SUPPORT-002"
  )

  incompatible_map <- b
  incompatible_map$category_order_within_macrogroup[
    incompatible_map$category_id == "BOTH_SIG"
  ] <- 99L
  expect_error(
    lisaR:::build_lisa_contrast_summary(
      a, incompatible_map, "A", "B", "PATHWAYS", TRUE, 0,
      gsea_padj_cutoff = 0.05
    ),
    "LISA-CONTRAST-SUPPORT-001"
  )
})

test_that("saved dumbbell source reproduces the approved support geometry", {
  skip_if_not_installed("ggplot2")
  observed <- lisaR:::build_lisa_contrast_summary(
    lisa_ab_summary_fixture("A"), lisa_ab_summary_fixture("B"),
    "A", "B", "PATHWAYS", TRUE, 0, gsea_padj_cutoff = 0.05
  )
  source <- tempfile("lisa-ab-source-", fileext = ".tsv")
  output <- tempfile("lisa-ab-reproduced-", fileext = ".png")
  lisaR:::lisa_write_contrast_figure_source(
    observed, source, "fixture", "A minus B", "Fixture", "all", "plain",
    "A", "B", 0.05
  )
  written <- lisaR:::read_lisa_tsv(source)
  expect_true(all(written$plot_set == "all"))
  expect_true(all(written$annotation_variant == "plain"))
  expect_true(all(written$group_by_supracategory))
  expect_true(all(written$figure_width == 12))
  expect_true(all(written$figure_height == lisaR:::lisa_contrast_plot_height(nrow(observed))))
  expect_true(all(written$figure_dpi == 300))
  expect_error(
    lisaR:::lisa_write_contrast_figure_source(
      observed, tempfile(fileext = ".tsv"), "fixture", "A minus B", "Fixture",
      NA_character_, "plain", "A", "B", 0.05
    ),
    "LISA-FIGURE-SOURCE-001"
  )
  recipe <- system.file("scripts", "reproduce_lisa_figure.R", package = "lisaR")
  log <- system2(lisaR:::lisa_rscript_executable(),
    c(shQuote(recipe), shQuote(source), shQuote(output)), stdout = TRUE, stderr = TRUE)
  status <- attr(log, "status")
  if (is.null(status)) status <- 0L
  expect_equal(status, 0L, info = paste(log, collapse = "\n"))
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
})

test_that("empty directional dumbbell subsets preserve reproducible metadata only", {
  skip_if_not_installed("ggplot2")
  observed <- lisaR:::build_lisa_contrast_summary(
    lisa_ab_summary_fixture("A"), lisa_ab_summary_fixture("B"),
    "A", "B", "PATHWAYS", TRUE, 0, gsea_padj_cutoff = 0.05
  )
  empty <- observed[FALSE, , drop = FALSE]
  source <- tempfile("lisa-ab-empty-source-", fileext = ".tsv")
  expect_no_error(lisaR:::lisa_write_contrast_figure_source(
    empty, source, "fixture", "A minus B", "Fixture",
    "opposite_direction", "plain", "A", "B", 0.05
  ))
  written <- lisaR:::read_lisa_tsv(source)
  expect_identical(nrow(written), 1L)
  expect_false(written$selected_for_plot[[1L]])
  expect_identical(written$figure_record_kind[[1L]], "empty_plot_metadata")
  expect_true(is.na(written$category_id[[1L]]))
  expect_true(all(c(
    "figure_id", "figure_type", "source_row_order", "selected_for_plot",
    "figure_record_kind", "figure_title", "figure_subtitle", "plot_set", "annotation_variant",
    "group_by_supracategory", "figure_width", "figure_height", "figure_dpi",
    "contrast_a_label", "contrast_b_label", "gsea_padj_cutoff"
  ) %in% colnames(written)))

  neither <- observed[observed$category_id == "NEITHER", , drop = FALSE]
  plot <- lisaR:::plot_lisa_contrast_dumbbell(
    neither, "A", "B", "A", "B", "A minus B", "Fixture",
    group_by_supracategory = FALSE, annotate_gene_sets = TRUE
  )
  expect_no_warning(ggplot2::ggplot_build(plot))

  recipe <- system.file("scripts", "reproduce_lisa_figure.R", package = "lisaR")
  output <- tempfile("lisa-ab-empty-reproduced-", fileext = ".png")
  log <- system2(
    lisaR:::lisa_rscript_executable(),
    c(shQuote(recipe), shQuote(source), shQuote(output)),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(log, "status")
  if (is.null(status)) status <- 0L
  expect_equal(status, 0L, info = paste(log, collapse = "\n"))
  expect_false(any(grepl("Warning", log, fixed = TRUE)))
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
})
