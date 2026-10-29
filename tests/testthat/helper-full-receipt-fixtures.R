receipt_fixture <- function() {
  root <- tempfile("receipt-rebase-")
  dir.create(root)
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  dir.create(file.path(root, "products"))
  dir.create(file.path(root, "audit/full_products"), recursive = TRUE)
  writeLines(c("path", file.path(root, "products/data.tsv")), file.path(root, "products/index.tsv"))
  writeLines(c("gene\tvalue", "PUBLIC_SYNTHETIC\t0.050000000000000003"), file.path(root, "products/data.tsv"))
  writeBin(as.raw(c(37, 80, 68, 70, 45, 0, 255, 10)), file.path(root, "products/figure.pdf"))
  paths <- c("products/index.tsv", "products/data.tsv", "products/figure.pdf")
  marker <- file.path(root, "audit/full_products/A__TEST__recurrent_genes.tsv")
  ledger <- data.frame(path = paths,
    sha256 = vapply(file.path(root, paths), lisaR:::lisa_sha256_file, character(1)),
    input_fingerprint = strrep("a", 64), stringsAsFactors = FALSE)
  lisaR:::write_lisa_tsv(ledger, marker)
  list(root = root, final = paste0(root, "-final"), marker = marker, ledger = ledger)
}

receipt_tree_hashes <- function(root) {
  files <- sort(list.files(root, recursive = TRUE, full.names = TRUE, all.files = TRUE))
  setNames(vapply(files, lisaR:::lisa_sha256_file, character(1)), substring(files, nchar(root) + 2L))
}

