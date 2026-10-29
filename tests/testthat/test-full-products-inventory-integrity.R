inventory_tx_fixture <- function() {
  final <- tempfile("lisa-final-inventory-long-name-")
  tx <- lisaR:::lisa_transaction_begin(final, lisa_test_run_contract(), run_id = "inventoryfinalization")
  folder <- file.path(tx$staging_dir, "outputs/gene_level/A")
  dir.create(folder, recursive = TRUE)
  writeLines(c("path", file.path(tx$staging_dir, "source_data/table.tsv")), file.path(folder, "index.tsv"))
  writeLines(c("gene\tstat", "G1\t2"), file.path(tx$staging_dir, "science.tsv"))
  lisaR:::write_lisa_tsv(data.frame(key = "report_mode", value = "full"), file.path(tx$staging_dir, "lisa_pipeline_metadata.tsv"))
  old <- options(lisaR.run_root = tx$staging_dir)
  on.exit(options(old), add = TRUE)
  lisaR:::lisa_write_full_products_inventory(tx$staging_dir)
  tx
}

inventory_read <- function(root) lisaR:::read_lisa_tsv(file.path(root, "artifacts/full_products_inventory.tsv"))

inventory_seal_outer_manifest <- function(root) {
  # Eager evaluation: snapshot before write_lisa_tsv creates its own temp file.
  manifest <- lisaR:::lisa_run_manifest(root)
  lisaR:::write_lisa_tsv(manifest, file.path(root, "run_manifest.tsv"))
}

test_that("promotion finalizes embedded FULL byte counts after root rebasing", {
  tx <- inventory_tx_fixture()
  before <- lisaR:::lisa_sha256_file(file.path(tx$staging_dir, "science.tsv"))
  lisaR:::lisa_transaction_promote(tx)
  rows <- inventory_read(tx$output_dir)
  expect_equal(rows$bytes, as.numeric(file.info(file.path(tx$output_dir, rows$path))$size))
  expect_identical(lisaR:::lisa_sha256_file(file.path(tx$output_dir, "science.tsv")), before)
  expect_identical(verify_lisa_run(tx$output_dir)$gate, "PASS")
  expect_false(any(grepl(tx$staging_dir, readLines(file.path(tx$output_dir, rows$path[[1]])), fixed = TRUE)))
})

test_that("public verification rejects stale embedded sizes even with a valid outer manifest", {
  tx <- inventory_tx_fixture(); lisaR:::lisa_transaction_promote(tx)
  rows <- inventory_read(tx$output_dir); rows$bytes[[1]] <- rows$bytes[[1]] + 1
  lisaR:::write_lisa_tsv(rows, file.path(tx$output_dir, "artifacts/full_products_inventory.tsv"))
  inventory_seal_outer_manifest(tx$output_dir)
  observed <- verify_lisa_run(tx$output_dir)
  expect_identical(observed$gate, "FAIL")
  expect_true(any(grepl("full_products_inventory", observed$findings, fixed = TRUE)))
  expect_false(any(grepl("^(absent|undeclared|modified|resized):", observed$findings)))
})

test_that("embedded inventory corruption is never blessed during path rebasing", {
  for (case in c("bytes", "duplicate", "omitted", "format", "traversal", "schema", "NA")) {
    tx <- inventory_tx_fixture(); rows <- inventory_read(tx$staging_dir)
    if (case == "bytes") rows$bytes[[1]] <- rows$bytes[[1]] + 1
    if (case == "duplicate") rows <- rbind(rows, rows)
    if (case == "omitted") rows <- rows[FALSE, ]
    if (case == "format") rows$format[[1]] <- "pdf"
    if (case == "traversal") rows$path[[1]] <- "../outside.tsv"
    if (case == "schema") rows$format <- NULL
    if (case == "NA") rows$bytes[[1]] <- NA_real_
    lisaR:::write_lisa_tsv(rows, file.path(tx$staging_dir, "artifacts/full_products_inventory.tsv"))
    paths <- list.files(tx$staging_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE)
    before <- vapply(paths, lisaR:::lisa_sha256_file, character(1))
    expect_error(lisaR:::lisa_rebase_run_text_paths(tx$staging_dir, tx$staging_dir, tx$output_dir), "LISA-RUN-010", info = case)
    expect_identical(vapply(paths, lisaR:::lisa_sha256_file, character(1)), before, info = case)
  }
})

test_that("FULL inventory is required by recorded FULL mode but not standard", {
  tx <- inventory_tx_fixture(); lisaR:::lisa_transaction_promote(tx)
  unlink(file.path(tx$output_dir, "artifacts/full_products_inventory.tsv"))
  inventory_seal_outer_manifest(tx$output_dir)
  observed <- verify_lisa_run(tx$output_dir)
  expect_identical(observed$gate, "FAIL")
  expect_true("full_products_inventory_absent" %in% observed$findings)
  lisaR:::write_lisa_tsv(data.frame(key = "report_mode", value = "standard"), file.path(tx$output_dir, "lisa_pipeline_metadata.tsv"))
  inventory_seal_outer_manifest(tx$output_dir)
  expect_identical(verify_lisa_run(tx$output_dir)$gate, "PASS")
})

test_that("contract FULL mode survives missing or conflicting pipeline metadata", {
  for (case in c("missing", "conflicting", "malformed")) {
    tx <- inventory_tx_fixture(); lisaR:::lisa_transaction_promote(tx)
    lisaR:::write_lisa_tsv(data.frame(key = "report_mode", value = "full"), file.path(tx$output_dir, "contract_manifest.tsv"))
    unlink(file.path(tx$output_dir, "artifacts/full_products_inventory.tsv"))
    metadata <- file.path(tx$output_dir, "lisa_pipeline_metadata.tsv")
    if (case == "missing") unlink(metadata)
    if (case == "conflicting") lisaR:::write_lisa_tsv(data.frame(key = "report_mode", value = "standard"), metadata)
    if (case == "malformed") writeLines("bad_header\nbad_value", metadata)
    inventory_seal_outer_manifest(tx$output_dir)
    observed <- verify_lisa_run(tx$output_dir)
    expect_identical(observed$gate, "FAIL", info = case)
    expect_true("full_products_inventory_absent" %in% observed$findings, info = case)
    if (case == "conflicting") expect_true("full_products_inventory_conflicting_modes" %in% observed$findings)
    if (case == "malformed") expect_true(any(grepl("invalid_mode", observed$findings, fixed = TRUE)))
    expect_false(any(grepl("^(absent|undeclared|modified|resized):", observed$findings)), info = case)
  }
})

test_that("inventory excludes only itself and normalizes equivalent root spellings", {
  tx <- inventory_tx_fixture()
  product <- file.path(tx$staging_dir, "outputs/gene_level/A/full_products_inventory.tsv")
  writeLines("legitimate same-basename product", product)
  old <- options(lisaR.run_root = tx$staging_dir)
  on.exit(options(old), add = TRUE)
  lisaR:::lisa_write_full_products_inventory(paste0(tx$staging_dir, "/"))
  rows <- inventory_read(tx$staging_dir)
  expect_true("outputs/gene_level/A/full_products_inventory.tsv" %in% rows$path)
  expect_false("artifacts/full_products_inventory.tsv" %in% rows$path)
  expect_identical(lisaR:::lisa_validate_full_products_inventory(paste0(tx$staging_dir, "/")), character())
  lisaR:::lisa_transaction_promote(tx)
  expect_identical(verify_lisa_run(tx$output_dir)$gate, "PASS")
})
