content_cache_args <- function(root = tempfile("lisa-content-cache-"), mode = "readwrite") {
  list(
    pathways = list(SET_A = c("G1", "G2"), SET_B = c("G2", "G3", "OUTSIDE_RANKS")),
    ranks = c(G1 = 3, G2 = -1, G3 = -2), min_gs_size = 1L,
    max_gs_size = 600L, n_threads = 1L, fgsea_nperm = 32L, fgsea_eps = 0,
    random_seed = 1729L, task_id = "gsea:fixture:GOBP-C2",
    universe = list(collection = "GOBP-C2", c2_sources = "REACTOME"),
    cache_options = lisaR:::lisa_content_cache_options(root, mode, 1024^2)
  )
}

content_cache_result <- function() {
  result <- data.frame(pathway = c("SET_A", "SET_B"), pval = c(0.12345678901234566, 1),
                       padj = c(0.24691357802469132, 1), ES = c(0.5, -0.25),
                       NES = c(1.1234567890123457, -0.7654321098765432),
                       size = c(2L, 2L), stringsAsFactors = FALSE)
  result$leadingEdge <- list(c("G1", "G2"), "G3")
  attr(result, "lisa_notes") <- c("Synthetic fixture note")
  result
}

test_that("passive cache encoding preserves exact numeric bytes and primitive types", {
  result <- data.frame(
    numeric = c(0, -0, .Machine$double.xmin, .Machine$double.xmax,
                .Machine$double.xmin * .Machine$double.eps, pi, NA_real_, NaN, Inf, -Inf),
    integer = c(NA_integer_, seq_len(9)), logical = rep(c(TRUE, FALSE, NA, TRUE, FALSE), 2),
    character = c(NA_character_, "", "Unicode: \u03b1", rep("a/b\t\n", 7)),
    stringsAsFactors = FALSE
  )
  result$leadingEdge <- rep(list(c("G1", "G2"), character()), 5)
  attr(result, "lisa_notes") <- "note"
  encoded <- lisaR:::lisa_cache_json(lisaR:::lisa_cache_encode_result(result))
  decoded <- lisaR:::lisa_cache_decode_result(jsonlite::fromJSON(encoded, simplifyVector = FALSE))
  expect_identical(decoded, result)
  expect_identical(writeBin(decoded$numeric, raw()), writeBin(result$numeric, raw()))
  expect_identical(lisaR:::lisa_cache_decode_result(jsonlite::fromJSON(
    lisaR:::lisa_cache_json(lisaR:::lisa_cache_encode_result(data.frame())),
    simplifyVector = FALSE)), data.frame())
  expect_error(lisaR:::lisa_cache_encode_result(structure(result, class = c("unsafe", "data.frame"))), "plain")
  expect_error(lisaR:::lisa_cache_encode_vector(list(new.env())), "Unsupported")
})

test_that("every enrichment-affecting input changes the content contract", {
  base <- content_cache_args()
  base$cache_options <- NULL
  base$versions <- list(R = "test R", fgsea = "test fgsea", scoreType = "std", gseaParam = 1)
  key <- function(args) do.call(lisaR:::lisa_fgsea_cache_contract, args)$key
  original <- key(base)
  variants <- list()
  changed <- base; changed$ranks[[1]] <- changed$ranks[[1]] + 1e-12; variants$ranks <- changed
  changed <- base; names(changed$ranks)[[1]] <- "OTHER"; variants$names <- changed
  changed <- base; changed$ranks <- rev(changed$ranks); variants$order <- changed
  changed <- base; changed$pathways$SET_A <- c("G1", "G3"); variants$membership <- changed
  changed <- base; changed$pathways$SET_B[[3]] <- "OTHER_OUTSIDE_RANKS"; variants$unmatched_membership <- changed
  changed <- base; changed$pathways <- changed$pathways[1]; variants$complete_family <- changed
  changed <- base; changed$universe$collection <- "GOCC"; variants$universe <- changed
  for (field in c("min_gs_size", "max_gs_size", "n_threads", "fgsea_nperm", "fgsea_eps", "random_seed")) {
    changed <- base; changed[[field]] <- changed[[field]] + 1; variants[[field]] <- changed
  }
  changed <- base; changed$task_id <- "gsea:other:GOBP-C2"; variants$task_id <- changed
  for (field in names(base$versions)) {
    changed <- base; changed$versions[[field]] <- "changed"; variants[[paste0("version_", field)]] <- changed
  }
  expect_identical(key(base), original)
  for (variant in variants) expect_false(identical(key(variant), original))
  snapshot <- lisaR:::lisa_rng_snapshot()
  on.exit(lisaR:::lisa_restore_rng(snapshot), add = TRUE)
  RNGkind(normal.kind = "Box-Muller")
  expect_false(identical(key(base), original))
})

test_that("a hit reuses exact output without engine work or caller RNG changes", {
  calls <- 0L
  testthat::local_mocked_bindings(
    lisa_fgsea_cache_versions = function() list(fgsea = "fixture"),
    run_fgsea_lisa = function(...) {
      calls <<- calls + 1L
      result <- content_cache_result()
      result$pval[[1L]] <- stats::runif(1L)
      result
    }, .package = "lisaR"
  )
  args <- content_cache_args()
  set.seed(713L)
  before <- lisaR:::lisa_rng_snapshot()
  miss <- do.call(lisaR:::lisa_cached_fgsea, args)
  expect_identical(lisaR:::lisa_rng_snapshot(), before)
  hit <- do.call(lisaR:::lisa_cached_fgsea, args)
  expect_identical(lisaR:::lisa_rng_snapshot(), before)
  expect_identical(miss$qc$status, "miss")
  expect_identical(miss$qc$write_status, "stored")
  expect_identical(hit$qc$status, "hit")
  expect_identical(hit$result, miss$result)
  expect_identical(calls, 1L)
  expect_match(hit$qc$artifact_sha256, "^[a-f0-9]{64}$")

  args$cache_options$mode <- "off"
  disabled <- do.call(lisaR:::lisa_cached_fgsea, args)
  expect_identical(disabled$result, hit$result)
  expect_identical(disabled$qc$status, "disabled")
  expect_identical(calls, 2L)
  expect_identical(lisaR:::lisa_rng_snapshot(), before)
  args$cache_options$mode <- "refresh"
  refreshed <- do.call(lisaR:::lisa_cached_fgsea, args)
  expect_identical(refreshed$result, hit$result)
  expect_identical(refreshed$qc$reason, "refresh_requested")
  expect_identical(calls, 3L)
})

test_that("corrupt records are misses and do not bypass artifact digest validation", {
  calls <- 0L
  testthat::local_mocked_bindings(
    lisa_fgsea_cache_versions = function() list(fgsea = "fixture"),
    run_fgsea_lisa = function(...) { calls <<- calls + 1L; content_cache_result() }, .package = "lisaR"
  )
  args <- content_cache_args()
  first <- do.call(lisaR:::lisa_cached_fgsea, args)
  slot <- list.files(args$cache_options$dir, pattern = "^slot-.*json$", full.names = TRUE, recursive = TRUE)
  record <- jsonlite::fromJSON(slot, simplifyVector = FALSE)
  record$payload <- sub("SET_A", "SET_Z", record$payload, fixed = TRUE)
  writeLines(lisaR:::lisa_cache_json(record), slot)
  repaired <- do.call(lisaR:::lisa_cached_fgsea, args)
  expect_identical(repaired$qc$status, "miss")
  expect_identical(repaired$qc$reason, "digest_mismatch")
  expect_identical(repaired$result, first$result)
  expect_identical(calls, 2L)
  writeLines("malformed JSON", slot)
  malformed <- do.call(lisaR:::lisa_cached_fgsea, args)
  expect_identical(malformed$qc$status, "miss")
  expect_identical(malformed$qc$reason, "malformed_or_unreadable")
  expect_identical(calls, 3L)
})

test_that("disabled and readonly cache policies never create a cache", {
  testthat::local_mocked_bindings(
    lisa_fgsea_cache_versions = function() list(fgsea = "fixture"),
    run_fgsea_lisa = function(...) content_cache_result(), .package = "lisaR"
  )
  expect_error(lisaR:::lisa_content_cache_options(NULL, "readwrite"), "explicit local")
  expect_error(lisaR:::lisa_content_cache_options(NULL, "off", 1.5), "whole")
  for (mode in c("off", "readonly")) {
    args <- content_cache_args(mode = mode)
    result <- do.call(lisaR:::lisa_cached_fgsea, args)
    expect_false(dir.exists(args$cache_options$dir))
    expect_identical(result$qc$write_status, "not_requested")
  }
})

test_that("bounded direct-mapped slots reject collisions and oversized entries", {
  options <- lisaR:::lisa_content_cache_options(tempfile("lisa-cache-bound-"), "readwrite", 131072)
  first <- list(key = paste0("a", paste(rep("1", 63), collapse = "")), parts = list(schema = "fixture"))
  second <- first; second$key <- paste0("a", paste(rep("2", 63), collapse = ""))
  slot <- lisaR:::lisa_cache_slot(options, first$key, create = TRUE)
  expect_identical(lisaR:::lisa_cache_write(slot, first, content_cache_result()), "stored")
  expect_true(lisaR:::lisa_cache_read(slot, first)$hit)
  expect_false(lisaR:::lisa_cache_read(slot, second)$hit)
  expect_identical(lisaR:::lisa_cache_write(slot, second, content_cache_result()), "stored")
  expect_false(lisaR:::lisa_cache_read(slot, first)$hit)
  expect_true(lisaR:::lisa_cache_read(slot, second)$hit)
  oversized <- content_cache_result(); oversized$leadingEdge[[1]] <- rep("G1", 5000)
  expect_identical(lisaR:::lisa_cache_write(slot, second, oversized), "entry_exceeds_slot_capacity")
  paths <- list.files(slot$root, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  expect_lte(sum(file.info(paths)$size), options$max_bytes)
  expect_length(paths, 1L)
  pending <- paste0(slot$path, ".pending")
  dir.create(pending)
  expect_identical(lisaR:::lisa_cache_write(slot, second, content_cache_result()), "slot_writer_busy_or_stale")
  expect_true(lisaR:::lisa_cache_read(slot, second)$hit)
})

test_that("cache annotations are regenerated independently from the shared enrichment", {
  testthat::local_mocked_bindings(
    lisa_fgsea_cache_versions = function() list(fgsea = "fixture"),
    run_fgsea_lisa = function(...) content_cache_result(), .package = "lisaR"
  )
  args <- content_cache_args()
  first <- do.call(lisaR:::lisa_cached_fgsea, args)
  core <- g2_dictionary("core", 4)[1L, , drop = FALSE]
  core$gene_set_id <- "SET_A"
  core$category_id <- "CORE"
  core$source_family <- "fixture"
  expanded <- core
  expanded$category_id <- "EXPANDED"
  core_annotated <- lisaR:::annotate_lisa(first$result, core, by_col = "pathway")
  hit <- do.call(lisaR:::lisa_cached_fgsea, args)
  expanded_annotated <- lisaR:::annotate_lisa(hit$result, expanded, by_col = "pathway")
  expect_identical(hit$qc$status, "hit")
  expect_false("category_id" %in% names(hit$result))
  expect_true("CORE" %in% core_annotated$category_id)
  expect_true("EXPANDED" %in% expanded_annotated$category_id)
  expect_identical(hit$result, first$result)
})

test_that("real small fgsea output is bit-exact on cache hit and disabled execution", {
  skip_if_not_installed("fgsea")
  args <- content_cache_args()
  args$ranks <- stats::setNames(seq(3, -3, length.out = 60), paste0("G", seq_len(60)))
  args$pathways <- list(UP = paste0("G", 1:12), DOWN = paste0("G", 45:60))
  first <- do.call(lisaR:::lisa_cached_fgsea, args)
  hit <- do.call(lisaR:::lisa_cached_fgsea, args)
  args$cache_options$mode <- "off"
  off <- do.call(lisaR:::lisa_cached_fgsea, args)
  expect_identical(first$qc$write_status, "stored")
  expect_identical(hit$qc$status, "hit")
  expect_identical(hit$result, first$result)
  expect_identical(hit$result, off$result)
})
