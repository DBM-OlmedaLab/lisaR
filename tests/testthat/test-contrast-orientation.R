test_that("LISA category deltas are oriented as contrast A minus contrast B", {
  make_summary <- function(mean_nes) {
    data.frame(
      category_id = "SYN_CATEGORY_A",
      category_display_name = "T-cell receptor signaling",
      n_genesets = 4,
      mean_NES = mean_nes,
      stringsAsFactors = FALSE
    )
  }

  observed <- lisaR:::build_lisa_contrast_summary(
    summary_a = make_summary(1.75),
    summary_b = make_summary(0.50),
    contrast_a_label = "Responders",
    contrast_b_label = "PD",
    universe = "PATHWAYS",
    include_missing_categories = TRUE,
    direction_epsilon = 0
  )

  expect_equal(observed$delta_mean_NES, 1.25)
  expect_equal(observed$abs_delta_mean_NES, 1.25)
})
