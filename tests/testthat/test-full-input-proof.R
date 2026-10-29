proof_fixture <- function() {
  f <- receipt_fixture()
  input <- "products/index.tsv"
  payload <- list(inputs = setNames(lisaR:::lisa_sha256_file(file.path(f$root, input)), input),
    builder = strrep("b", 64), args = c("--project-dir", f$root))
  fp <- digest::digest(payload, algo = "sha256")
  f$ledger$input_fingerprint <- fp
  lisaR:::write_lisa_tsv(f$ledger, f$marker)
  lisaR:::lisa_write_full_input_proof(f$marker, payload, fp, f$root)
  f$payload <- payload
  f$fp <- fp
  f$proof <- paste0(f$marker, ".inputs.rds")
  f
}

test_that("passive proof roundtrip preserves the exact old typed digest object", {
  f <- proof_fixture()
  proof <- lisaR:::lisa_read_full_input_proof(f$marker, f$fp, f$root)
  expect_identical(proof$payload, f$payload)
  expect_identical(digest::digest(proof$payload, algo = "sha256"), f$fp)
  expect_identical(digest::digest(proof$payload, algo = "sha256",
    serializeVersion = proof$digest_serialize_version), f$fp)
  expect_identical(proof$r_version, R.version.string)
  expect_identical(proof$digest_version, as.character(packageVersion("digest")))
  expect_error(lisaR:::lisa_write_full_input_proof(f$marker, f$payload, f$fp, f$root), "overwrite")
  payload <- f$payload; payload$args[[2L]] <- f$final
  expect_false(identical(digest::digest(payload, algo = "sha256"), f$fp))
})

test_that("rebase changes current target receipts but not immutable prepromotion proof", {
  f <- proof_fixture()
  before_proof <- lisaR:::lisa_sha256_file(f$proof)
  before_input <- f$payload$inputs[[1L]]
  lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$final)
  ledger <- lisaR:::read_lisa_tsv(f$marker)
  expect_identical(lisaR:::lisa_sha256_file(f$proof), before_proof)
  expect_identical(lisaR:::lisa_read_full_input_proof(f$marker, f$fp, f$root)$payload, f$payload)
  expect_identical(unique(ledger$input_fingerprint), f$fp)
  expect_identical(ledger$sha256[[1L]], lisaR:::lisa_sha256_file(file.path(f$root, ledger$path[[1L]])))
  expect_false(identical(ledger$sha256[[1L]], before_input))
  # An inverse ROOT transform proves old bytes only when the exact old hash agrees.
  current <- readBin(file.path(f$root, "products/index.tsv"), "raw", n = 1e6)
  restored <- charToRaw(gsub(f$final, f$root, rawToChar(current), fixed = TRUE))
  expect_identical(digest::digest(restored, algo = "sha256", serialize = FALSE), before_input)
  restored[[length(restored)]] <- as.raw(0)
  expect_false(identical(digest::digest(restored, algo = "sha256", serialize = FALSE), before_input))
})

test_that("malformed, tampered, unsafe, orphan and required missing proofs fail shut", {
  cases <- c("hash", "builder", "args", "input_name", "duplicate_names", "input_type",
    "argument_type", "list_order", "fingerprint", "serialize", "schema", "unreadable",
    "traversal", "absolute", "input_attribute", "orphan")
  for (case in cases) {
    f <- proof_fixture(); proof <- readRDS(f$proof)
    if (case == "hash") proof$payload$inputs[[1]] <- strrep("0", 64)
    if (case == "builder") proof$payload$builder <- strrep("0", 64)
    if (case == "args") proof$payload$args[[2]] <- f$final
    if (case == "input_name") names(proof$payload$inputs) <- "products/data.tsv"
    if (case == "duplicate_names") proof$payload$inputs <- rep(proof$payload$inputs, 2)
    if (case == "input_type") proof$payload$inputs <- as.list(proof$payload$inputs)
    if (case == "argument_type") proof$payload$args <- as.list(proof$payload$args)
    if (case == "list_order") proof$payload <- proof$payload[c(3, 1, 2)]
    if (case == "fingerprint") proof$input_fingerprint <- strrep("0", 64)
    if (case == "serialize") proof$digest_serialize_version <- 99L
    if (case == "schema") proof$schema <- "unknown"
    if (case == "traversal") names(proof$payload$inputs) <- "../outside.tsv"
    if (case == "absolute") names(proof$payload$inputs) <- file.path(f$root, "products/index.tsv")
    if (case == "input_attribute") attr(proof$payload$inputs, "unexpected") <- TRUE
    saveRDS(proof, f$proof, version = 2L)
    if (case == "unreadable") writeLines("not rds", f$proof)
    if (case == "orphan") unlink(f$marker)
    before <- receipt_tree_hashes(f$root)
    expect_error(lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$final), "(LISA-FULL-PROOF|LISA-RUN-009)", info = case)
    expect_identical(receipt_tree_hashes(f$root), before, info = case)
  }
  f <- proof_fixture(); unlink(f$proof)
  expect_error(lisaR:::lisa_read_full_input_proof(f$marker, f$fp, f$root), "missing FULL")
  expect_null(lisaR:::lisa_read_full_input_proof(f$marker, f$fp, f$root, required = FALSE))
  # Legacy compatibility cannot prove that a missing sidecar never existed.
  expect_silent(lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$root))
  f <- proof_fixture(); target <- paste0(f$root, "-outside.rds")
  file.copy(f$proof, target); unlink(f$proof); file.symlink(target, f$proof)
  expect_error(lisaR:::lisa_read_full_input_proof(f$marker, f$fp, f$root), "(link|Symbolic)")
})

completion_fixture <- function() {
  root <- tempfile("completion-proof-"); dir.create(root)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  dir.create(file.path(root, "config"))
  lisaR:::write_lisa_tsv(data.frame(analysis_id = "A", species = "Homo sapiens"), file.path(root, "config/de_index.tsv"))
  lisaR:::write_lisa_tsv(data.frame(contrast_id = character(), output_id = character(), contrast_a = character(), contrast_b = character()), file.path(root, "config/contrast_index.tsv"))
  dir.create(file.path(root, "outputs/single_de/A/collection_TEST"), recursive = TRUE)
  gene <- file.path(root, "outputs/gene_level/single_de/A/collection_TEST")
  dir.create(gene, recursive = TRUE)
  lisaR:::write_lisa_tsv(data.frame(symbol = "PUBLIC_SYNTHETIC"), file.path(gene, "A_TEST_gene_level_gene_category_contributions.tsv"))
  list(root = root, image = file.path(gene, "category_gene_cards/public.png"),
    marker = file.path(root, "audit/full_products/A__TEST__prioritized_genes.tsv"))
}

test_that("new family lifecycle writes proof and existing reuse validates it", {
  f <- completion_fixture()
  # Isolate ONLY the external presentation renderer. Real completion, hashes,
  # guards, proof construction, marker writes and reused validation run unchanged.
  e <- new.env(parent = environment(lisaR:::lisa_complete_full_products))
  e$system2 <- function(...) {
    dir.create(dirname(f$image), recursive = TRUE, showWarnings = FALSE)
    writeBin(as.raw(c(137, 80, 78, 71, 0)), f$image)
    0L
  }
  complete <- lisaR:::lisa_complete_full_products; environment(complete) <- e
  result <- complete(f$root, formats = "png", families = "prioritized_genes")
  expect_identical(result$status, "generated_missing")
  ledger <- lisaR:::read_lisa_tsv(f$marker)
  proof <- lisaR:::lisa_read_full_input_proof(f$marker, unique(ledger$input_fingerprint), f$root)
  expect_identical(digest::digest(proof$payload, algo = "sha256"), unique(ledger$input_fingerprint))
  before <- lisaR:::lisa_sha256_file(paste0(f$marker, ".inputs.rds"))
  result <- complete(f$root, formats = "png", families = "prioritized_genes")
  expect_identical(result$status, "reused_verified")
  expect_identical(lisaR:::lisa_sha256_file(paste0(f$marker, ".inputs.rds")), before)
  proof$payload$builder <- strrep("0", 64)
  saveRDS(proof, paste0(f$marker, ".inputs.rds"), version = 2L)
  expect_error(complete(f$root, formats = "png", families = "prioritized_genes"), "LISA-FULL-PROOF")
})

test_that("inverse root translation proves an inherited ledger's OLD target SHA exactly", {
  f <- receipt_fixture()
  # The downloaded receipt is not regenerated/blessed after promotion. Its old
  # target hash remains evidence to check against reconstructed old target bytes.
  downloaded <- file.path(f$root, "downloaded_receipt.tsv")
  lisaR:::write_lisa_tsv(f$ledger, downloaded)
  original_ledger_hash <- lisaR:::lisa_sha256_file(downloaded)
  old_index_sha <- f$ledger$sha256[[1L]]
  lisaR:::lisa_rebase_run_text_paths(f$root, f$root, f$final)
  current <- file.path(f$root, "products/index.tsv")
  raw <- readBin(current, "raw", n = file.info(current)$size)
  inverse <- charToRaw(gsub(f$final, f$root, rawToChar(raw), fixed = TRUE))
  expect_identical(digest::digest(inverse, algo = "sha256", serialize = FALSE), old_index_sha)
  expect_identical(lisaR:::lisa_sha256_file(downloaded), original_ledger_hash)
  expect_false(identical(lisaR:::lisa_sha256_file(current), old_index_sha))
})
