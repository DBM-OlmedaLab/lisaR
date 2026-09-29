# Focused, installed tests for the general example-acquisition entry
# `install_lisa_example_bundle()`. These exercise transport and integrity only:
# no DE/limma, no GSEA, no report, no KEGG, and no scientific or external API
# call. Network transport is tested against a controlled local endpoint
# (`file://` and a loopback HTTP server serving a temporary directory).
#
# Optional integration fixtures are supplied explicitly by the caller.
# Tests that require a distribution archive skip if it is not available.

example_distribution_archives <- function() {
  list(
    riaz = Sys.getenv("LISAR_RIAZ_ARCHIVE", unset = ""),
    cptac = Sys.getenv("LISAR_CPTAC_ARCHIVE", unset = "")
  )
}

skip_without_archive <- function(path) {
  skip_if_not(
    file.exists(path) && !dir.exists(path),
    paste0("reviewed distribution archive is not present in this environment: ", path)
  )
}

# Serve one directory over loopback HTTP with the interpreter that is already
# present; the server is a bounded child of this test and is always stopped.
local_http_endpoint <- function(directory, envir = parent.frame()) {
  python <- Sys.which("python3")
  skip_if_not(nzchar(python), "python3 is not available for the local HTTP endpoint test")
  port <- 24000L + sample.int(2000L, 1L)
  log <- tempfile("lisar-http-endpoint-")
  pid_file <- tempfile("lisar-http-endpoint-pid-")
  # The server is a bounded child of this test; its own pid file is the only
  # thing this suite ever stops.
  system2(
    "/bin/sh",
    c("-c", shQuote(sprintf(
      "echo $$ > %s; exec %s -m http.server %d --bind 127.0.0.1 --directory %s",
      shQuote(pid_file), shQuote(python), port, shQuote(directory)
    ))),
    stdout = log, stderr = log, wait = FALSE
  )
  withr::defer({
    if (file.exists(pid_file)) {
      pid <- suppressWarnings(as.integer(readLines(pid_file, warn = FALSE)[[1]]))
      if (!is.na(pid)) try(tools::pskill(pid), silent = TRUE)
    }
  }, envir = envir)
  base <- sprintf("http://127.0.0.1:%d/", port)
  ready <- FALSE
  for (attempt in seq_len(80L)) {
    probe <- suppressWarnings(tryCatch({
      connection <- url(base, open = "rb")
      on.exit(close(connection), add = TRUE)
      readBin(connection, "raw", 1L)
      TRUE
    }, error = function(error) FALSE))
    if (isTRUE(probe)) {
      ready <- TRUE
      break
    }
    Sys.sleep(0.1)
  }
  skip_if_not(ready, "the local HTTP endpoint did not become reachable")
  base
}

test_that("no public endpoint is invented when none is published", {
  contracts <- lisaR:::lisa_example_distribution_contracts()
  for (key in names(contracts)) contracts[[key]]$source_url <- NA_character_
  testthat::local_mocked_bindings(
    lisa_example_distribution_contracts = function() contracts, .package = "lisaR"
  )
  root <- withr::local_tempdir("lisar-example-gate-")
  destination <- file.path(root, "bundle")
  expect_error(
    install_lisa_example_bundle(
      "riaz-gse91061", destination,
      cache_root = file.path(root, "cache")
    ),
    "LISA-EXAMPLE-031.*no public distribution endpoint is published yet"
  )
  # The gate must not be reported by silently creating anything.
  expect_false(file.exists(destination))
})

test_that("unknown examples are rejected before any filesystem work", {
  root <- withr::local_tempdir("lisar-example-unknown-")
  expect_error(
    install_lisa_example_bundle(
      "not-an-example", file.path(root, "bundle"),
      cache_root = file.path(root, "cache")
    ),
    "LISA-EXAMPLE-030 unknown example"
  )
  expect_false(file.exists(file.path(root, "bundle")))
})

test_that("a reviewed archive is obtained over file:// and promoted atomically", {
  archives <- example_distribution_archives()
  skip_without_archive(archives$riaz)
  root <- withr::local_tempdir("lisar-example-file-")
  destination <- file.path(root, "riaz-bundle")
  cache <- file.path(root, "cache")

  result <- install_lisa_example_bundle(
    "riaz-gse91061", destination,
    source = paste0("file://", normalizePath(archives$riaz, mustWork = TRUE)),
    cache_root = cache
  )

  expect_true(dir.exists(destination))
  expect_identical(result$example, "riaz-gse91061")
  expect_true(result$acquired)
  expect_identical(
    result$archive_sha256,
    "983efbe6fce2bfbbd0a3df691d60603c53b452d1e2cdc3d5ab08642df808b144"
  )
  expect_identical(result$scientific_files, 17L)
  expect_identical(result$files, 22L)
  expect_identical(result$public_endpoint, lisaR:::lisa_example_distribution_contract(result$example)$source_url)
  # The promoted bundle is exactly what the existing preparer expects.
  expect_true(file.exists(file.path(destination, "prepared_bundle_manifest.tsv")))
  expect_true(file.exists(file.path(destination, "prepared_bundle_provenance.tsv")))
  expect_true(file.exists(file.path(
    destination, "outputs/lisa_inputs/responders_on_vs_pre.tsv"
  )))
  expect_match(result$next_command, "LISAR_RIAZ_PREPARED_BUNDLE=")
  expect_match(result$next_command, "11_prepare_prepared_project\\.R")
  # The command must be runnable as printed: a raw <placeholder> would be read
  # by the shell as a redirection and would leave the project variable empty.
  expect_false(grepl("[<>]", result$next_command))
  expect_match(
    result$next_command,
    paste0("LISAR_RIAZ_PROJECT_DIR=", shQuote(paste0(destination, "-project"))),
    fixed = TRUE
  )
  # This preparer never calls run_lisa(), so it carries no prepare-only switch.
  expect_false(grepl("PREPARE_ONLY", result$next_command))
  # No staging residue is left beside the destination.
  expect_length(
    list.files(root, pattern = "bundle-staging", all.files = TRUE), 0L
  )

  # A second install reuses the verified cached archive with no transport at
  # all: an unreachable source must not be consulted.
  second <- install_lisa_example_bundle(
    "riaz-gse91061", file.path(root, "riaz-bundle-2"),
    source = "https://127.0.0.1:1/never-contacted.tar.gz",
    cache_root = cache
  )
  expect_false(second$acquired)
  expect_identical(second$bundle_manifest_sha256, result$bundle_manifest_sha256)
})

test_that("a reviewed archive is obtained over a controlled local HTTP endpoint", {
  archives <- example_distribution_archives()
  skip_without_archive(archives$cptac)
  skip_if_not_installed("withr")
  served <- withr::local_tempdir("lisar-example-served-")
  file.copy(archives$cptac, file.path(served, basename(archives$cptac)))
  base <- local_http_endpoint(served)

  root <- withr::local_tempdir("lisar-example-http-")
  destination <- file.path(root, "cptac-bundle")
  result <- install_lisa_example_bundle(
    "cptac-ccrcc", destination,
    source = paste0(base, basename(archives$cptac)),
    cache_root = file.path(root, "cache")
  )

  expect_true(dir.exists(destination))
  expect_true(result$acquired)
  expect_identical(result$scientific_files, 6L)
  expect_true(file.exists(file.path(destination, "cptac_ccrcc_bundle_manifest.tsv")))
  expect_true(file.exists(file.path(
    destination, "de/de_primary_conservative_80pairs.tsv"
  )))
  expect_match(result$next_command, "LISAR_CPTAC_CCRCC_BUNDLE=")
  expect_identical(result$public_endpoint, lisaR:::lisa_example_distribution_contract(result$example)$source_url)
  # Runnable as printed, with a concrete project directory.
  expect_false(grepl("[<>]", result$next_command))
  expect_match(
    result$next_command,
    paste0("LISAR_CPTAC_CCRCC_PROJECT_DIR=",
           shQuote(paste0(destination, "-project"))),
    fixed = TRUE
  )
  # The proteomic preparer continues into run_lisa() unless pinned, so
  # obtaining a bundle must never hand the recipient a scientific run.
  expect_match(result$next_command, "LISAR_CPTAC_CCRCC_PREPARE_ONLY=1",
               fixed = TRUE)
})

test_that("a cached archive with the same version but a different identity is replaced", {
  archives <- example_distribution_archives()
  skip_without_archive(archives$cptac)
  root <- withr::local_tempdir("lisar-example-cache-identity-")
  cache <- file.path(root, "cache")
  contract <- lisaR:::lisa_example_distribution_contract("cptac-ccrcc")
  cached <- lisaR:::lisa_example_distribution_archive_path(cache, contract)
  dir.create(dirname(cached), recursive = TRUE)
  bytes <- readBin(archives$cptac, "raw", n = file.size(archives$cptac))
  middle <- length(bytes) %/% 2L
  bytes[[middle]] <- as.raw(bitwXor(as.integer(bytes[[middle]]), 1L))
  writeBin(bytes, cached)

  result <- install_lisa_example_bundle(
    "cptac-ccrcc", file.path(root, "bundle"),
    source = archives$cptac, cache_root = cache
  )
  expect_true(result$acquired)
  expect_identical(lisa_sha256_file(cached), contract$archive_sha256)
  expect_identical(result$archive_sha256, contract$archive_sha256)
  expect_true(file.exists(file.path(result$bundle, contract$bundle_manifest)))
})

test_that("a truncated transfer leaves the destination and cache untouched", {
  archives <- example_distribution_archives()
  skip_without_archive(archives$cptac)
  root <- withr::local_tempdir("lisar-example-truncated-")
  broken <- file.path(root, basename(archives$cptac))
  bytes <- readBin(archives$cptac, "raw", n = 4096L)
  writeBin(bytes, broken)
  destination <- file.path(root, "bundle")
  cache <- file.path(root, "cache")

  expect_error(
    install_lisa_example_bundle(
      "cptac-ccrcc", destination,
      source = paste0("file://", normalizePath(broken, mustWork = TRUE)),
      cache_root = cache
    ),
    "LISA-EXAMPLE-033"
  )
  expect_false(file.exists(destination))
  expect_false(file.exists(file.path(
    cache, "example-archives", "cptac-1.0.0.tar.gz"
  )))
})

test_that("a same-size altered archive is refused on digest alone", {
  archives <- example_distribution_archives()
  skip_without_archive(archives$cptac)
  root <- withr::local_tempdir("lisar-example-substituted-")
  # Correct name, correct byte count, one changed byte: only the pinned digest
  # can catch this, so the byte-size bound must not be what rejects it.
  altered <- file.path(root, basename(archives$cptac))
  bytes <- readBin(archives$cptac, "raw", n = file.size(archives$cptac))
  middle <- length(bytes) %/% 2L
  bytes[[middle]] <- as.raw(bitwXor(as.integer(bytes[[middle]]), 1L))
  writeBin(bytes, altered)
  expect_identical(file.size(altered), file.size(archives$cptac))
  destination <- file.path(root, "bundle")

  expect_error(
    install_lisa_example_bundle(
      "cptac-ccrcc", destination,
      source = paste0("file://", normalizePath(altered, mustWork = TRUE)),
      cache_root = file.path(root, "cache")
    ),
    "LISA-EXAMPLE-033.*SHA-256 does not match"
  )
  expect_false(file.exists(destination))
  expect_false(file.exists(file.path(
    root, "cache", "example-archives", "cptac-1.0.0.tar.gz"
  )))
})

test_that("an oversized archive is refused by the bounded acquisition policy", {
  archives <- example_distribution_archives()
  skip_without_archive(archives$riaz)
  root <- withr::local_tempdir("lisar-example-oversized-")
  # The 19 MiB Riaz archive served under the proteomic example's name exceeds
  # that example's bound and is refused before it is ever read as an archive.
  oversized <- file.path(root, basename(
    lisaR:::lisa_example_distribution_contract("cptac-ccrcc")$archive_name
  ))
  file.copy(archives$riaz, oversized)
  destination <- file.path(root, "bundle")

  expect_error(
    install_lisa_example_bundle(
      "cptac-ccrcc", destination,
      source = paste0("file://", normalizePath(oversized, mustWork = TRUE)),
      cache_root = file.path(root, "cache")
    ),
    "LISA-EXAMPLE-032.*bounded acquisition policy"
  )
  expect_false(file.exists(destination))
})

test_that("an existing destination is never overwritten", {
  archives <- example_distribution_archives()
  skip_without_archive(archives$riaz)
  root <- withr::local_tempdir("lisar-example-existing-")
  destination <- file.path(root, "bundle")
  dir.create(destination)
  writeLines("keep me", file.path(destination, "sentinel.txt"))

  expect_error(
    install_lisa_example_bundle(
      "riaz-gse91061", destination,
      source = paste0("file://", normalizePath(archives$riaz, mustWork = TRUE)),
      cache_root = file.path(root, "cache")
    ),
    "LISA-EXAMPLE-030 destination already exists"
  )
  expect_identical(readLines(file.path(destination, "sentinel.txt")), "keep me")
  expect_length(list.files(destination), 1L)
})

test_that("unsafe and duplicated archive members are refused before extraction", {
  root <- withr::local_tempdir("lisar-example-members-")
  contract <- lisaR:::lisa_example_distribution_contract("cptac-ccrcc")

  build <- function(name, paths) {
    staging <- file.path(root, name)
    dir.create(staging)
    for (path in paths) {
      target <- file.path(staging, path)
      dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
      writeLines("value", target)
    }
    archive <- file.path(root, paste0(name, ".tar.gz"))
    old <- setwd(staging)
    on.exit(setwd(old), add = TRUE)
    utils::tar(archive, list.files(".", recursive = TRUE), compression = "gzip")
    archive
  }

  safe <- build("safe", c(
    "cptac_ccrcc_bundle_manifest.tsv", "CHECKSUMS.sha256", "README.md",
    "de/de_primary_conservative_80pairs.tsv"
  ))
  members <- lisaR:::lisa_example_distribution_members(safe, contract)
  expect_true("cptac_ccrcc_bundle_manifest.tsv" %in% members$files)

  incomplete <- build("incomplete", c("README.md", "de/only.tsv"))
  expect_error(
    lisaR:::lisa_example_distribution_members(incomplete, contract),
    "LISA-EXAMPLE-034.*missing required bundle-root"
  )

  # Both remaining cases are packaging defects that a post-extraction tree scan
  # cannot see, so they are built with the system archiver and rejected from the
  # listing alone.
  system_tar <- Sys.which("tar")
  skip_if_not(nzchar(system_tar), "system tar is not available for the packaging-defect cases")
  base <- file.path(root, "safe")

  unsafe <- file.path(root, "unsafe.tar.gz")
  escape <- file.path(root, "escape.tsv")
  writeLines("value", escape)
  expect_identical(system2(system_tar, c(
    "-czf", shQuote(unsafe), "-P", "-C", shQuote(base), shQuote("../escape.tsv")
  )), 0L)
  expect_error(
    lisaR:::lisa_example_distribution_members(unsafe, contract),
    "LISA-EXAMPLE-034.*absolute, parent-relative"
  )

  # Two entries, one unique name: exactly the historical duplicate/hard-link
  # packaging defect.
  duplicate <- file.path(root, "duplicate.tar.gz")
  expect_identical(system2(system_tar, c(
    "-czf", shQuote(duplicate), "-C", shQuote(base),
    "README.md", "README.md"
  )), 0L)
  expect_error(
    lisaR:::lisa_example_distribution_members(duplicate, contract),
    "LISA-EXAMPLE-034.*same path more than once"
  )
})

test_that("an acquired bundle is consumed unchanged by the existing prepare-only preparer", {
  archives <- example_distribution_archives()
  skip_without_archive(archives$cptac)
  cache <- Sys.getenv("LISAR_TEST_RESOURCE_CACHE", unset = "")
  skip_if_not(
    dir.exists(cache),
    "the reviewed normal resource cache is not present in this checkout"
  )
  root <- withr::local_tempdir("lisar-example-chain-")
  bundle <- file.path(root, "cptac-bundle")
  acquired <- install_lisa_example_bundle(
    "cptac-ccrcc", bundle,
    source = paste0("file://", normalizePath(archives$cptac, mustWork = TRUE)),
    cache_root = file.path(root, "cache")
  )

  script <- system.file(
    "examples/cptac-ccrcc/prepare_and_run_cptac_ccrcc.R", package = "lisaR",
    mustWork = TRUE
  )
  project <- file.path(root, "prepared-project")
  # Export rather than `system2(env = )`: on Windows that argument puts the
  # NAME=value tokens on Rscript's command line, where the first one is taken
  # as the script to execute (?system2). This test is currently skipped on
  # Windows, so the defect is latent here, but the mechanism is identical.
  out <- lisaR:::lisa_with_child_environment(
    c(
      LISAR_CPTAC_CCRCC_BUNDLE = acquired$bundle,
      LISAR_CPTAC_CCRCC_PROJECT_DIR = project,
      LISAR_CPTAC_CCRCC_RESOURCE_CACHE = cache,
      LISAR_CPTAC_CCRCC_PREPARE_ONLY = "1",
      R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)
    ),
    suppressWarnings(system2(
      file.path(R.home("bin"), "Rscript"), c("--vanilla", script),
      stdout = TRUE, stderr = TRUE
    ))
  )
  expect_true(any(grepl("CPTAC_CCRCC_LOCAL_ROUTE=PREPARE_ONLY_PASS", out)),
              info = paste(utils::tail(out, 20L), collapse = "\n"))
  expect_true(dir.exists(project))
  # Prepare-only really stopped before any analysis.
  expect_false(dir.exists(file.path(project, "runs")))
  # The acquired bundle itself was not modified by the preparer.
  expect_identical(
    lisaR:::lisa_sha256_file(file.path(acquired$bundle, acquired$bundle_manifest)),
    acquired$bundle_manifest_sha256
  )
})

test_that("a tampered scientific file inside the bundle is refused", {
  archives <- example_distribution_archives()
  skip_without_archive(archives$cptac)
  root <- withr::local_tempdir("lisar-example-tampered-")
  staging <- file.path(root, "staging")
  dir.create(staging, mode = "0700")
  expect_identical(
    as.integer(utils::untar(archives$cptac, exdir = staging, tar = "internal")),
    0L
  )
  contract <- lisaR:::lisa_example_distribution_contract("cptac-ccrcc")
  members <- lisaR:::lisa_example_distribution_members(archives$cptac, contract)
  staging <- normalizePath(staging, winslash = "/", mustWork = TRUE)

  # The intact tree verifies.
  intact <- lisaR:::lisa_example_distribution_verify_tree(staging, members, contract)
  expect_identical(nrow(intact$manifest), 6L)

  tampered <- file.path(staging, "audit/data_quality_summary.tsv")
  lines <- readLines(tampered, warn = FALSE)
  writeLines(c(lines, "tampered\textra\trow"), tampered)
  expect_error(
    lisaR:::lisa_example_distribution_verify_tree(staging, members, contract),
    "LISA-EXAMPLE-036.*does not match its own reviewed manifest"
  )
})


test_that("published examples resolve their pinned versioned assets without network discovery", {
  for (contract in lisaR:::lisa_example_distribution_contracts()) {
    source <- lisaR:::lisa_example_distribution_source(NULL, contract)
    expect_true(source$pinned)
    expect_identical(source$url, paste0(
      "https://github.com/DBM-OlmedaLab/lisaR-example-data/releases/download/",
      "v1.0.0/", contract$archive_name))
  }
})
