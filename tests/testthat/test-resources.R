test_that("Resource dictionary resources resolve only from a managed cache", {
  root <- tempfile("lisa-dictionary-cache-")
  artifact <- file.path(root, "fixture_dictionary", "1.0.0", "dictionary.tsv")
  dir.create(dirname(artifact), recursive = TRUE)
  writeLines("fixture dictionary", artifact)
  registry <- data.frame(
    logical_id = "fixture_dictionary", version = "1.0.0", species = "Homo sapiens",
    modality = "transcriptomic/genomic", schema = "lisa_dictionary@2",
    sha256 = lisaR:::lisa_sha256_file(artifact), compatibility = "lisaR>=0.5.0",
    approved_origin = "local test fixture", artifact = "dictionary.tsv", stringsAsFactors = FALSE
  )
  resolved <- lisa_resolve_dictionary_resource("fixture_dictionary", "1.0.0", "Homo sapiens", registry = registry, cache_root = root)
  # lisaR reports managed paths with forward slashes on every platform;
  # `normalizePath()` defaults to backslashes on Windows.
  expect_identical(resolved$path, normalizePath(artifact, winslash = "/"))
  writeLines("corrupt", artifact)
  expect_error(lisa_resolve_dictionary_resource("fixture_dictionary", "1.0.0", "Homo sapiens", registry = registry, cache_root = root), "SHA-256 mismatch")
  expect_error(lisa_resolve_dictionary_resource("absent", registry = registry, cache_root = root), "did not resolve")
})

test_that("Resource species contract prevents human/mouse mismatches", {
  human <- lisa_species_contract("Homo sapiens")
  mouse <- lisa_species_contract("Mus musculus")
  expect_identical(human$kegg_code, "hsa")
  expect_identical(human$annotation_db, "org.Hs.eg.db")
  expect_identical(mouse$kegg_code, "mmu")
  expect_identical(mouse$annotation_db, "org.Mm.eg.db")
  expect_error(lisa_species_contract("Homo sapiens", kegg_code = "mmu"), "conflicts")
  expect_error(lisa_species_contract("Mus musculus", msigdb_mode = "human"), "conflicts")
})

test_that("Resource KEGG snapshots are immutable and cache-only", {
  root <- tempfile("lisa-kegg-cache-")
  source <- tempfile("lisa-kegg-fixture-")
  writeLines("local KEGG fixture", source)
  first <- lisa_store_kegg_snapshot(source, root, "snapshot-1", "hsa", "pathway_links", "hsa:1")
  second <- lisa_resolve_kegg_resource("cache_only", root, "snapshot-1", "hsa", "pathway_links", "hsa:1")
  expect_identical(first$metadata$sha256[[1]], second$metadata$sha256[[1]])
  expect_identical(first$metadata$key[[1]], "hsa:1")
  expect_identical(first$path, second$path)
  expect_false(grepl(":", basename(first$path), fixed = TRUE))
  expect_silent(lisaR:::lisa_safe_id(basename(first$path)))
  expect_error(lisa_store_kegg_snapshot(source, root, "snapshot-1", "hsa", "pathway_links", "hsa:1"), "already exists")
  expect_error(lisa_resolve_kegg_resource("cache_only", root, "missing", "hsa", "pathway_links", "hsa:1"), "absent")
})

test_that("Resource requires species before analysis", {
  de <- data.frame(analysis_id = "x", de_path = "missing.tsv", stringsAsFactors = FALSE)
  expect_error(run_lisa_pipeline(de, dictionary_dir = "unused", term2gene = "unused", output_dir = tempfile("lisa-contract-"), dry_run = TRUE), "every analysis must declare species")
})

test_that("Resource rejects per-run resource roots", {
  cfg <- list(pipeline = list(
    schema_version = "1.0.0", profile = "transcriptomic/genomic", evidence_mode = "full_de",
    duplicate_policies = list(
      de_table_duplicate_policy = "error", matrix_duplicate_policy = "error",
      mapped_id_collision_policy = "error"
    ),
    output_dir = tempfile("lisa-contract-config-"), dictionary_cache_root = "/tmp/free-root"
  ))
  path <- tempfile(fileext = ".json")
  jsonlite::write_json(cfg, path, auto_unbox = TRUE)
  expect_error(run_lisa_pipeline_from_config(path), "per-run resource root")
})

test_that("config execution restores the run root after an error", {
  old <- options("lisaR.run_root")
  on.exit(options(old), add = TRUE)
  options(lisaR.run_root = "sentinel-run-root")
  cfg <- list(pipeline = list(
    schema_version = "1.0.0", profile = "transcriptomic/genomic", evidence_mode = "full_de",
    duplicate_policies = list(
      de_table_duplicate_policy = "error", matrix_duplicate_policy = "error",
      mapped_id_collision_policy = "error"
    ),
    output_dir = tempdir(), dictionary_cache_root = "/tmp/free-root"
  ))
  path <- tempfile(fileext = ".json")
  jsonlite::write_json(cfg, path, auto_unbox = TRUE)

  expect_error(run_lisa_pipeline_from_config(path), "per-run resource root")
  expect_identical(getOption("lisaR.run_root"), "sentinel-run-root")
})

test_that("KEGG painters use explicit external or cache-only access", {
  scripts <- system.file("scripts", c(
    "build_single_de_kegg_pathway_painter.R",
    "build_contrast_kegg_pathway_painter.R",
    "kegg_snapshot_helpers.R"
  ), package = "lisaR")
  expect_true(all(file.exists(scripts)))
  code <- paste(unlist(lapply(scripts, readLines, warn = FALSE)), collapse = "\n")
  expect_false(grepl("species = \\\"mmu\\\"", code, fixed = TRUE))
  expect_false(grepl("(?:/home|/Users)/[^/[:space:]]+/", code, perl = TRUE))
  expect_match(code, "cache_only", fixed = TRUE)
  expect_match(code, "prefer_cache", fixed = TRUE)
  expect_match(code, "lisa_fetch_kegg_resource", fixed = TRUE)
  expect_match(code, "org.Hs.eg.db", fixed = TRUE)
  expect_match(code, "org.Mm.eg.db", fixed = TRUE)
})

test_that("external KEGG retrieval rejects unsupported requests before network access", {
  expect_error(
    lisaR:::lisa_fetch_kegg_resource("eco", "pathway_list", "all"),
    "unsupported KEGG organism"
  )
  expect_error(
    lisaR:::lisa_fetch_kegg_resource("hsa", "unsupported", "all"),
    "unsupported KEGG resource type"
  )
  expect_error(
    lisaR:::lisa_fetch_kegg_resource("hsa", "kgml", "mmu00010"),
    "does not match the declared organism"
  )
})

test_that("binary KEGG image retrieval records only the known text-encoding warning", {
  exact <- lisaR:::lisa_kegg_binary_response(function() {
    warning("No encoding supplied: defaulting to UTF-8.")
    array(1, dim = c(2, 2, 3))
  })
  expect_identical(dim(exact$value), c(2L, 2L, 3L))
  expect_match(exact$encoding_note, "binary PNG response", fixed = TRUE)

  expect_warning(
    lisaR:::lisa_kegg_binary_response(function() {
      warning("unexpected transport warning")
      array(1, dim = c(2, 2, 3))
    }),
    "unexpected transport warning"
  )
})

test_that("Resource painters require an explicit project directory", {
  scripts <- system.file("scripts", c(
    "build_single_de_kegg_pathway_painter.R",
    "build_contrast_kegg_pathway_painter.R"
  ), package = "lisaR")
  expect_true(all(file.exists(scripts)))
  for (script in scripts) {
    output <- suppressWarnings(system2(lisaR:::lisa_rscript_executable(), c(script), stdout = TRUE, stderr = TRUE))
    expect_true(!is.null(attr(output, "status")) && attr(output, "status") != 0L)
    expect_match(paste(output, collapse = "\n"), "--project-dir")
  }
})

test_that("Resource KEGG painters collapse multiline KGML and reject unpainted renders", {
  scripts <- system.file("scripts", c(
    "build_single_de_kegg_pathway_painter.R",
    "build_contrast_kegg_pathway_painter.R"
  ), package = "lisaR")
  expect_true(all(file.exists(scripts)))
  code <- paste(unlist(lapply(scripts, readLines, warn = FALSE)), collapse = "\n")
  expect_match(code, 'paste(kgml, collapse = "\\n")', fixed = TRUE)
  expect_match(code, "LISA-KEGG-017", fixed = TRUE)
  expect_match(code, "LISA-KEGG-018", fixed = TRUE)
  expect_match(code, "LISA-KEGG-019", fixed = TRUE)
  expect_match(code, "stale_maps", fixed = TRUE)
  expect_match(code, "LISA-KEGG-016", fixed = TRUE)
})

test_that("enabled EnrichmentMap stages render archived preview graphs", {
  pipeline_path <- testthat::test_path("..", "..", "R", "lisa_pipeline.R")
  skip_if_not(file.exists(pipeline_path), "Source-only pipeline code is not installed in the built package")
  pipeline <- readLines(pipeline_path, warn = FALSE)
  code <- paste(pipeline, collapse = "\n")
  expect_match(code, '"--render-graphs", "true"', fixed = TRUE)
})

test_that("Resource cache-only snapshots are identical offline and fail closed", {
  root <- tempfile("lisa-kegg-offline-")
  source <- tempfile(fileext = ".rds")
  saveRDS(c("hsa:1" = "path:hsa00010"), source)
  lisa_store_kegg_snapshot(source, root, "frozen", "hsa", "pathway_links", "all", retrieved_at = "2026-07-12T00:00:00Z")
  a <- readRDS(lisa_resolve_kegg_resource("cache_only", root, "frozen", "hsa", "pathway_links", "all")$path)
  b <- readRDS(lisa_resolve_kegg_resource("cache_only", root, "frozen", "hsa", "pathway_links", "all")$path)
  expect_identical(a, b)
  expect_error(lisa_resolve_kegg_resource("cache_only", root, "frozen", "mmu", "pathway_links", "all"), "absent")
  expect_error(lisa_resolve_kegg_resource("prefer_cache", root, "missing", "hsa", "pathway_links", "all"), "explicit fetch")
})
