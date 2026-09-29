make_g3_deseq2_dds <- function() {
  set.seed(1729)
  condition <- factor(
    rep(c("control", "treated"), each = 3L),
    levels = c("control", "treated")
  )
  baseline <- seq(40, 240, length.out = 120L)
  means <- vapply(
    seq_along(condition),
    function(index) {
      treatment_effect <- if (condition[[index]] == "treated") {
        rep(c(2, 0.5, 1), length.out = length(baseline))
      } else {
        rep(1, length(baseline))
      }
      baseline * treatment_effect
    },
    numeric(length(baseline))
  )
  counts <- matrix(
    stats::rnbinom(length(means), mu = as.vector(means), size = 8),
    nrow = length(baseline),
    dimnames = list(
      paste0("GENE", seq_along(baseline)),
      paste0("sample", seq_along(condition))
    )
  )
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = data.frame(condition = condition, row.names = colnames(counts)),
    design = ~ condition
  )
  DESeq2::DESeq(dds, fitType = "mean", quiet = TRUE)
}

test_that("DESeqDataSet accepts one named coefficient selector", {
  skip_if_not_installed("DESeq2")
  dds <- make_g3_deseq2_dds()
  coefficient <- "condition_treated_vs_control"

  extracted <- lisaR:::extract_de_table(
    dds,
    "deseq2_dds",
    deseq2_contrast = NULL,
    deseq2_name = coefficient
  )
  expected <- as.data.frame(DESeq2::results(dds, name = coefficient))

  expect_equal(extracted, expected)
  expect_identical(rownames(extracted), rownames(dds))
})

test_that("DESeqDataSet accepts one three-part contrast selector", {
  skip_if_not_installed("DESeq2")
  dds <- make_g3_deseq2_dds()
  contrast <- c("condition", "treated", "control")

  extracted <- lisaR:::extract_de_table(
    dds,
    "deseq2_dds",
    deseq2_contrast = contrast,
    deseq2_name = NULL
  )
  expected <- as.data.frame(DESeq2::results(dds, contrast = contrast))

  expect_equal(extracted, expected)
  expect_identical(rownames(extracted), rownames(dds))
})

test_that("DESeqDataSet requires exactly one result selector", {
  skip_if_not_installed("DESeq2")
  dds <- make_g3_deseq2_dds()

  expect_error(
    lisaR:::extract_de_table(dds, "deseq2_dds", NULL, NULL),
    "LISA-INPUT-DESEQ2-001"
  )
  expect_error(
    lisaR:::extract_de_table(
      dds,
      "deseq2_dds",
      c("condition", "treated", "control"),
      "condition_treated_vs_control"
    ),
    "LISA-INPUT-DESEQ2-002"
  )
})

test_that("DESeqDataSet selectors fail closed when malformed or unknown", {
  skip_if_not_installed("DESeq2")
  dds <- make_g3_deseq2_dds()

  expect_error(
    lisaR:::extract_de_table(dds, "deseq2_dds", NULL, ""),
    "LISA-INPUT-DESEQ2-003"
  )
  expect_error(
    lisaR:::extract_de_table(dds, "deseq2_dds", NULL, "missing_coefficient"),
    "LISA-INPUT-DESEQ2-004"
  )
  expect_error(
    lisaR:::extract_de_table(dds, "deseq2_dds", c("condition", "treated"), NULL),
    "LISA-INPUT-DESEQ2-005"
  )
})
