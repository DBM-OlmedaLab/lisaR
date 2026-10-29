test_that("real path translation reconciles only preverified FULL target hashes", {
  f <- receipt_fixture()
  before <- receipt_tree_hashes(f$root)
  changed <- lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$final)
  ledger <- lisaR:::read_lisa_tsv(f$marker)
  observed <- vapply(file.path(f$root, ledger$path), lisaR:::lisa_sha256_file, character(1))
  expect_identical(ledger$sha256, unname(observed))
  expect_false(identical(ledger$sha256[[1]], f$ledger$sha256[[1]]))
  expect_identical(ledger$sha256[-1], f$ledger$sha256[-1])
  expect_identical(ledger$input_fingerprint, f$ledger$input_fingerprint)
  expect_identical(lisaR:::lisa_sha256_file(file.path(f$root, "products/data.tsv")), unname(before[["products/data.tsv"]]))
  expect_identical(lisaR:::lisa_sha256_file(file.path(f$root, "products/figure.pdf")), unname(before[["products/figure.pdf"]]))
  expect_match(readLines(file.path(f$root, "products/index.tsv"))[[2]], f$final, fixed = TRUE)
  expect_setequal(changed, c(file.path(f$root, "products/index.tsv"), f$marker))
})

test_that("all current FULL families use the same checked rebase contract", {
  f <- receipt_fixture()
  families <- c("prioritized_genes", "volcano", "gene_heatmap", "recurrent_genes",
    "kegg_maps", "paired_gene_products", "gene_category_network")
  for (family in families) lisaR:::write_lisa_tsv(f$ledger,
    file.path(f$root, "audit/full_products", paste0("A__TEST__", family, ".tsv")))
  lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$final)
  for (family in families) {
    ledger <- lisaR:::read_lisa_tsv(file.path(f$root, "audit/full_products", paste0("A__TEST__", family, ".tsv")))
    expect_identical(ledger$sha256, unname(vapply(file.path(f$root, ledger$path), lisaR:::lisa_sha256_file, character(1))))
    expect_identical(ledger$input_fingerprint, f$ledger$input_fingerprint)
  }
})

test_that("invalid legacy receipts fail before any fixture byte is changed", {
  cases <- c("tampered_sha", "missing", "escape", "absolute", "drive", "backslash",
    "duplicate", "dot_alias", "double_slash", "empty_path", "bad_sha", "bad_fingerprint",
    "mixed_fingerprint", "extra_column", "duplicate_header", "missing_column", "cycle", "pending", "late_bad_receipt")
  for (case in cases) {
    f <- receipt_fixture()
    ledger <- f$ledger
    if (case == "tampered_sha") ledger$sha256[[1]] <- strrep("0", 64)
    if (case == "missing") unlink(file.path(f$root, ledger$path[[1]]))
    if (case == "escape") ledger$path[[1]] <- "../outside.tsv"
    if (case == "absolute") ledger$path[[1]] <- file.path(f$root, ledger$path[[1]])
    if (case == "drive") ledger$path[[1]] <- "C:/products/index.tsv"
    if (case == "backslash") ledger$path[[1]] <- "products\\index.tsv"
    if (case == "duplicate") ledger <- rbind(ledger, ledger[1L, ])
    if (case == "dot_alias") ledger$path[[1]] <- "products/./index.tsv"
    if (case == "double_slash") ledger$path[[1]] <- "products//index.tsv"
    if (case == "empty_path") ledger$path[[1]] <- ""
    if (case == "bad_sha") ledger$sha256[[1]] <- "broken"
    if (case == "bad_fingerprint") ledger$input_fingerprint[] <- "broken"
    if (case == "mixed_fingerprint") ledger$input_fingerprint[[1]] <- strrep("2", 64)
    if (case == "extra_column") ledger$unknown <- "not_a_known_schema"
    if (case == "duplicate_header") names(ledger)[[3]] <- "sha256"
    if (case == "missing_column") ledger$sha256 <- NULL
    if (case == "cycle") ledger$path[[1]] <- "audit/full_products/A__TEST__recurrent_genes.tsv"
    if (case == "pending") writeLines("interrupted", paste0(f$marker, ".pending"))
    if (case == "late_bad_receipt") {
      ledger$sha256[[1]] <- strrep("0", 64)
      lisaR:::write_lisa_tsv(ledger, file.path(f$root, "audit/full_products/Z__TEST__volcano.tsv"))
      ledger <- f$ledger
    }
    lisaR:::write_lisa_tsv(ledger, f$marker)
    before <- receipt_tree_hashes(f$root)
    expect_error(lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$final), "LISA-RUN-009", info = case)
    expect_identical(receipt_tree_hashes(f$root), before, info = case)
  }
})

test_that("same-root no-op preserves bytes without waiving corrupt receipts", {
  f <- receipt_fixture()
  before <- receipt_tree_hashes(f$root)
  expect_identical(lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$root), character())
  expect_identical(receipt_tree_hashes(f$root), before)
  ledger <- f$ledger
  ledger$sha256[[1]] <- strrep("0", 64)
  lisaR:::write_lisa_tsv(ledger, f$marker)
  before <- receipt_tree_hashes(f$root)
  expect_error(lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$root), "LISA-RUN-009")
  expect_identical(receipt_tree_hashes(f$root), before)
})

test_that("copy-only saved derivations rebase their own validated targets", {
  f <- receipt_fixture()
  before <- receipt_tree_hashes(f$root)
  copy <- paste0(f$root, "-copy")
  fs::dir_copy(f$root, copy)
  lisaR:::lisa_rebase_run_text_paths(copy, f$root, copy)
  expect_identical(receipt_tree_hashes(f$root), before)
  ledger <- lisaR:::read_lisa_tsv(file.path(copy, "audit/full_products", basename(f$marker)))
  expect_identical(ledger$sha256, unname(vapply(file.path(copy, ledger$path), lisaR:::lisa_sha256_file, character(1))))
  expect_identical(ledger$input_fingerprint, f$ledger$input_fingerprint)
  expect_match(readLines(file.path(copy, "products/index.tsv"))[[2]], copy, fixed = TRUE)
})

test_that("ordinary STANDARD text trees without a FULL ledger retain existing behavior", {
  root <- tempfile("ordinary-rebase-")
  dir.create(root)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  writeLines(c("path", paste0(root, "/output.tsv")), file.path(root, "index.tsv"))
  final <- paste0(root, "-final")
  expect_identical(lisaR:::lisa_rebase_run_text_paths(root, root, final), file.path(root, "index.tsv"))
  expect_identical(readLines(file.path(root, "index.tsv")), c("path", paste0(final, "/output.tsv")))
  expect_false(dir.exists(file.path(root, "audit/full_products")))
})
