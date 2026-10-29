inventory_fixture <- function(f) {
  rows <- data.frame(owner = c("A", "A", "A", "A"), collection = "TEST",
    family = c("recurrent_genes", "volcano", "gene_heatmap", "prioritized_genes"),
    state = c("present", "completed_empty", "empty_with_evidence", "image_formats_not_requested"),
    path = c("products/index.tsv", "", "", ""),
    sha256 = c(f$ledger$sha256[[1L]], "", "", ""), stringsAsFactors = FALSE)
  path <- file.path(f$root, "audit/full_products/required_family_inventory.tsv")
  lisaR:::write_lisa_tsv(rows, path)
  list(path = path, rows = rows)
}

test_that("aggregate inventory binds real targets without acquiring marker semantics", {
  f <- receipt_fixture(); inv <- inventory_fixture(f)
  before <- receipt_tree_hashes(f$root)
  expect_identical(lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$root), character())
  expect_identical(receipt_tree_hashes(f$root), before)
  lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$final)
  observed <- utils::read.delim(inv$path, colClasses = "character")
  expect_identical(names(observed), c("owner", "collection", "family", "state", "path", "sha256"))
  expect_identical(observed$sha256[[1L]], lisaR:::lisa_sha256_file(file.path(f$root, observed$path[[1L]])))
  expect_identical(observed[-1, ], inv$rows[-1, ])
  expect_false(identical(observed$sha256[[1L]], inv$rows$sha256[[1L]]))
  copy <- paste0(f$root, "-copy"); fs::dir_copy(f$root, copy)
  lisaR:::lisa_rebase_run_text_paths(copy, f$final, copy)
  observed_copy <- utils::read.delim(file.path(copy, "audit/full_products/required_family_inventory.tsv"), colClasses = "character")
  expect_identical(observed_copy$sha256[[1L]], lisaR:::lisa_sha256_file(file.path(copy, "products/index.tsv")))
})

test_that("unknown TSV and malformed aggregate controls reject without mutation", {
  for (case in c("unknown", "schema", "state", "blank_hash", "stale", "empty_has_path", "traversal", "duplicate", "inventory_proof")) {
    f <- receipt_fixture(); inv <- inventory_fixture(f); rows <- inv$rows
    if (case == "unknown") lisaR:::write_lisa_tsv(f$ledger, file.path(f$root, "audit/full_products/unknown.tsv"))
    if (case == "schema") rows$input_fingerprint <- strrep("a", 64)
    if (case == "state") rows$state[[2]] <- "unreviewed"
    if (case == "blank_hash") rows$sha256[[1]] <- ""
    if (case == "stale") rows$sha256[[1]] <- strrep("0", 64)
    if (case == "empty_has_path") rows$path[[2]] <- "products/data.tsv"
    if (case == "traversal") rows$path[[1]] <- "../outside.tsv"
    if (case == "duplicate") rows <- rbind(rows, rows[1, ])
    if (case == "inventory_proof") saveRDS(list(schema = "not a family"), paste0(inv$path, ".inputs.rds"))
    lisaR:::write_lisa_tsv(rows, inv$path)
    before <- receipt_tree_hashes(f$root)
    expect_error(lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$final), "LISA-RUN-009", info = case)
    expect_identical(receipt_tree_hashes(f$root), before, info = case)
  }
})
