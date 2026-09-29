test_that("calibration sensitivity grid preserves individual genes and weights", {
  evidence <- data.frame(
    priority_id = c("A", "A", "B", "B", "C"),
    gene_set_id = c("S1", "S2", "S1", "S3", "S4"),
    category_id = c("C1", "C2", "C1", "C3", "C4"),
    effect = c(2, 2, 1, 1, 3),
    padj = c(.01, .01, .02, .02, .03)
  )
  grid <- lisa_gps_sensitivity_grid(evidence)
  expect_equal(nrow(grid$weights), 4L)
  expect_equal(nrow(grid$scores), 12L)
  expect_setequal(grid$scores$priority_id, c("A", "B", "C"))
  expect_true(all(grid$scores$rank[grid$scores$applicability == "applicable"] >= 1))
  expect_equal(grid$weights$scenario_id[[1]], "equal")
  expect_equal(grid$weights[grid$weights$scenario_id == "equal", c("effect_weight", "confidence_weight", "breadth_weight")], data.frame(effect_weight = 1 / 3, confidence_weight = 1 / 3, breadth_weight = 1 / 3))
})

test_that("calibration stability and concordance do not pool shRNA contrasts", {
  scores <- data.frame(
    scenario_id = rep(c("equal", "effect_emphasis"), each = 3),
    priority_id = rep(c("A", "B", "C"), 2),
    rank = c(1L, 2L, 3L, 2L, 1L, 3L),
    lisa_gps = c(.9, .6, .2, .7, .8, .1),
    stringsAsFactors = FALSE
  )
  stability <- lisa_gps_rank_stability(scores)
  expect_equal(nrow(stability), 3L)
  expect_true(all(c("rank_min", "rank_max", "rank_sd") %in% names(stability)))
  sh2 <- data.frame(priority_id = c("A", "B", "C"), rank = c(1L, 2L, 3L), lisa_gps = c(.8, .6, .2))
  sh3 <- data.frame(priority_id = c("A", "B", "C"), rank = c(2L, 1L, 3L), lisa_gps = c(.7, .9, .1))
  out <- lisa_gps_concordance(sh2, sh3, top_n = 2)
  expect_equal(out$n_shared, 3L)
  expect_equal(out$top_n_overlap, 2L)
  expect_equal(out$comparison, "separate_inputs_no_pooling")
})

test_that("calibration evidence spaces retain collection-qualified identities", {
  evidence <- data.frame(
    priority_id = c("A", "A"),
    gene_set_id = c("GOBP::SHARED_SET", "GOMF::SHARED_SET"),
    category_id = c("GOBP::CAT_PROCESS", "GOMF::CAT_FUNCTION"),
    effect = c(2, 2),
    padj = c(.01, .01),
    stringsAsFactors = FALSE
  )
  out <- lisa_gps_sensitivity_grid(evidence)
  equal <- out$scores[out$scores$scenario_id == "equal", , drop = FALSE]
  expect_equal(equal$priority_id, "A")
  expect_equal(equal$independent_category_breadth, 2L)
  expect_equal(equal$raw_gene_set_count, 2L)
})
