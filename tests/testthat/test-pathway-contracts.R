test_that("lisa-gps produces one auditable row per priority ID", {
  evidence <- data.frame(
    priority_id = c("P1", "P1", "P1", "P2", "P2", "P3"),
    gene_set_id = c("GS1", "GS1", "GS2", "GS3", "GS4", "GS5"),
    category_id = c("CAT_A", "CAT_A", "CAT_B", "CAT_A", "CAT_C", "CAT_D"),
    effect = c(1, 1, 1, 2, 2, 3),
    padj = c(0.1, 0.1, 0.1, 0.01, 0.01, 0),
    stringsAsFactors = FALSE
  )

  out <- lisa_gps_prioritization(evidence)
  expect_identical(out$priority_id, c("P1", "P2", "P3"))
  expect_identical(out$raw_gene_set_count, c(2L, 2L, 1L))
  expect_identical(out$independent_category_breadth, c(2L, 2L, 1L))
  expect_identical(out$gene_set_ids, c("GS1;GS2", "GS3;GS4", "GS5"))
  expect_identical(out$category_ids, c("CAT_A;CAT_B", "CAT_A;CAT_C", "CAT_D"))
  expect_true(out$zero_padj_clamped[[3]])
  expect_true(is.finite(out$de_confidence_raw[[3]]))
  expect_equal(out$lisa_gps, (out$effect_strength_normalized * out$de_confidence_normalized * out$independent_category_breadth_normalized)^(1 / 3))
  expect_match(lisa_gps_development_contract()$status, "provisional")
})

test_that("lisa-gps requires explicit IDs and never resolves inconsistent DE statistics", {
  complete <- data.frame(priority_id = c("P1", "P1"), gene_set_id = c("GS1", "GS2"), category_id = c("CAT_A", "CAT_B"), effect = c(1, 1), padj = c(0.1, 0.1))
  expect_error(lisa_gps_prioritization(complete[-2]), "gene_set_id")
  expect_error(lisa_gps_prioritization(transform(complete, category_id = "")), "must be non-empty")
  expect_error(lisa_gps_prioritization(transform(complete, effect = c(1, 2))), "inconsistent effect or padj.*will not select or aggregate DE statistics")
  expect_error(lisa_gps_prioritization(transform(complete, padj = c(0.1, 0.2))), "inconsistent effect or padj")
  expect_error(lisa_gps_prioritization(transform(complete, padj = -0.1)), "within \\[0, 1\\]")
})

test_that("lisa-gps preserves legitimate many-to-many dictionary assignments", {
  evidence <- data.frame(
    priority_id = c("P1", "P1"),
    gene_set_id = c("GS1", "GS1"),
    category_id = c("CAT_A", "CAT_B"),
    effect = c(1, 1),
    padj = c(0.1, 0.1),
    stringsAsFactors = FALSE
  )
  out <- lisa_gps_prioritization(evidence)
  expect_equal(out$raw_gene_set_count, 1L)
  expect_equal(out$independent_category_breadth, 2L)
  expect_equal(out$gene_set_ids, "GS1")
  expect_equal(out$category_ids, "CAT_A;CAT_B")
})

test_that("lisa-gps handles missing values, ties, and non-full evidence deterministically", {
  x <- data.frame(
    priority_id = c("z", "a", "n", "zero"),
    gene_set_id = c("G1", "G2", "G3", "G4"),
    category_id = c("C1", "C2", "C3", "C4"),
    effect = c(1, 1, NA, 1),
    padj = c(0.1, 0.1, 0.2, 0),
    stringsAsFactors = FALSE
  )
  out <- lisa_gps_prioritization(x)
  expect_identical(out$priority_id, c("a", "n", "z", "zero"))
  expect_equal(out$lisa_gps[out$priority_id %in% c("a", "z")], rep(out$lisa_gps[out$priority_id == "a"], 2))
  expect_true(is.na(out$lisa_gps[out$priority_id == "n"]))
  expect_true(all(is.na(out$pareto_front[out$priority_id %in% c("a", "z")])))
  expect_equal(out$pareto_front[out$priority_id == "zero"], 1L)
  rank_only <- lisa_gps_prioritization(x, "rank_only")
  expect_true(all(rank_only$applicability == "not_applicable"))
  expect_true(all(is.na(rank_only$lisa_gps)))
  expect_true(all(is.na(rank_only$pareto_front)))
  expect_true(all(rank_only$pareto_role == "not_applicable"))
})

test_that("Pareto first-front membership matches brute force with ties", {
  set.seed(7)
  values <- matrix(sample(0:9, 300, replace = TRUE), ncol = 3)
  x <- data.frame(
    effect_strength_normalized = values[, 1],
    de_confidence_normalized = values[, 2],
    independent_category_breadth_normalized = values[, 3],
    applicability = "applicable"
  )
  observed <- !is.na(lisa_gps_pareto_front(x))
  expected <- vapply(seq_len(nrow(values)), function(i) {
    !any(vapply(seq_len(nrow(values)), function(j) {
      j != i && all(values[j, ] >= values[i, ]) && any(values[j, ] > values[i, ])
    }, logical(1)))
  }, logical(1))
  expect_identical(observed, expected)
})
