test_that("install_lisa_resource copies and persistently registers local TSVs", {
  root <- tempfile("lisa-resource-install-")
  dir.create(root)
  source <- file.path(root, "source-term2gene.tsv")
  utils::write.table(
    data.frame(
      gs_collection = "C5", gs_subcollection = "GO:BP",
      gs_name = "GS_A", gs_exact_source = "TEST", gene_symbol = "GENE1",
      stringsAsFactors = FALSE
    ),
    source, sep = "\t", quote = FALSE, row.names = FALSE
  )
  cache <- file.path(root, "cache")
  withr::local_options(list(
    lisaR.dictionary_registry = NULL,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))

  installed <- install_lisa_resource(
    source = source,
    logical_id = "fixture_term2gene",
    version = "1.0.0",
    species = "Homo sapiens",
    modality = "transcriptomic/genomic",
    schema = "term2gene@1",
    approved_origin = "authorised local test fixture",
    expected_sha256 = lisaR:::lisa_sha256_file(source)
  )

  expect_true(file.exists(installed$path))
  expect_true(file.exists(installed$registry_path))
  expect_identical(installed$resource_id, "fixture_term2gene@1.0.0")
  expect_identical(
    lisaR:::lisa_sha256_file(installed$path),
    lisaR:::lisa_sha256_file(source)
  )
  expect_identical(
    lisaR:::lisa_active_resource_context()$registry,
    installed$registry_path
  )
  resolved <- lisaR:::lisa_resolve_dictionary_resource(
    "fixture_term2gene", "1.0.0", "Homo sapiens",
    modality = "transcriptomic/genomic",
    registry = installed$registry_path,
    cache_root = cache
  )
  expect_identical(resolved$path, installed$path)

  repeated <- install_lisa_resource(
    source = source,
    logical_id = "fixture_term2gene",
    version = "1.0.0",
    species = "Homo sapiens",
    modality = "transcriptomic/genomic",
    schema = "term2gene@1",
    approved_origin = "authorised local test fixture",
    expected_sha256 = lisaR:::lisa_sha256_file(source)
  )
  expect_identical(repeated$resource_id, installed$resource_id)
  expect_identical(repeated$path, installed$path)
})

test_that("install_lisa_resource enforces immutable versions and accepts reviewed source families", {
  root <- tempfile("lisa-resource-immutable-")
  dir.create(root)
  source <- file.path(root, "source-term2gene.tsv")
  utils::write.table(
    data.frame(
      gs_collection = "C5", gs_subcollection = "GO:BP",
      gs_name = "GS_A", gs_exact_source = "TEST", gene_symbol = "GENE1",
      stringsAsFactors = FALSE
    ),
    source, sep = "\t", quote = FALSE, row.names = FALSE
  )
  cache <- file.path(root, "cache")
  args <- list(
    source = source, logical_id = "fixture_term2gene", version = "1.0.0",
    species = "Homo sapiens", modality = "transcriptomic/genomic",
    schema = "term2gene@1",
    approved_origin = "authorised local test fixture",
    expected_sha256 = lisaR:::lisa_sha256_file(source), cache_root = cache
  )
  do.call(install_lisa_resource, args)
  utils::write.table(
    data.frame(
      gs_collection = c("C5", "C5"),
      gs_subcollection = c("GO:BP", "GO:BP"),
      gs_name = c("GS_A", "GS_B"),
      gs_exact_source = c("TEST", "TEST"),
      gene_symbol = c("GENE1", "GENE2"),
      stringsAsFactors = FALSE
    ),
    source, sep = "\t", quote = FALSE, row.names = FALSE
  )
  args$expected_sha256 <- lisaR:::lisa_sha256_file(source)
  expect_error(
    do.call(install_lisa_resource, args),
    "LISA-RESOURCE-024.*different metadata or SHA-256.*new version"
  )

  reviewed <- file.path(root, "reviewed-dictionary.tsv")
  utils::write.table(
    data.frame(
      universe = c("GOBP-C2", "GOBP-C2"),
      gene_set_id = c("BIOCARTA_REVIEWED", "KEGG_REVIEWED"),
      gene_set_name = c("BioCarta reviewed", "KEGG Legacy reviewed"),
      source_id = c("BIOCARTA", "KEGG_LEGACY"),
      category_id = "CAT_A", category_display_name = "Category A",
      tier = "custom", stringsAsFactors = FALSE
    ),
    reviewed, sep = "\t", quote = FALSE, row.names = FALSE
  )
  expect_no_error(
    install_lisa_resource(
      source = reviewed, logical_id = "fixture_reviewed_dictionary", version = "1.0.0",
      species = "Homo sapiens", modality = "transcriptomic/genomic",
      schema = "lisa_dictionary@2",
      approved_origin = "reviewed local source-family fixture",
      expected_sha256 = lisaR:::lisa_sha256_file(reviewed),
      cache_root = cache
    )
  )
})

test_that("prepare_lisa_msigdb_resource rejects an unpinned ZIP before staging", {
  root <- tempfile("lisa-msigdb-negative-")
  dir.create(root)
  source_zip <- file.path(root, "not-msigdb.zip")
  writeLines("synthetic unpinned ZIP fixture", source_zip, useBytes = TRUE)
  cache <- file.path(root, "cache")

  expect_error(
    prepare_lisa_msigdb_resource(source_zip, cache_root = cache),
    "LISA-RESOURCE-028.*SHA-256 does not match"
  )
  expect_false(file.exists(cache))
})

test_that("MSigDB acquisition requires recipient consent before network or cache writes", {
  root <- tempfile("lisa-msigdb-consent-")
  dir.create(root)
  cache <- file.path(root, "cache")
  testthat::local_mocked_bindings(
    lisa_msigdb_2026_1_download_zip = function(...) {
      stop("network must not be reached without consent")
    },
    .package = "lisaR"
  )

  expect_error(
    prepare_lisa_msigdb_resource(cache_root = cache),
    "LISA-RESOURCE-028.*recipient has not explicitly accepted.*MSigDB terms.*KEGG terms.*accept_terms = TRUE"
  )
  expect_false(file.exists(cache))
})

test_that("MSigDB acquisition failures do not promote corrupt archives or registry rows", {
  root <- tempfile("lisa-msigdb-download-failure-")
  dir.create(root)
  cache <- file.path(root, "cache")
  archive <- file.path(cache, "archives", "msigdb.2026.1.zip")
  testthat::local_mocked_bindings(
    lisa_msigdb_2026_1_download_file = function(url, destination, contract) {
      writeLines("corrupt synthetic download", destination, useBytes = TRUE)
      invisible(destination)
    },
    .package = "lisaR"
  )

  expect_error(
    prepare_lisa_msigdb_resource(accept_terms = TRUE, cache_root = cache),
    "LISA-RESOURCE-028.*source_zip SHA-256"
  )
  expect_false(file.exists(archive))
  expect_false(file.exists(file.path(cache, "resource_registry.tsv")))
})

test_that("MSigDB acquisition leaves an existing failed archive untouched", {
  root <- tempfile("lisa-msigdb-archive-preserve-")
  dir.create(root)
  cache <- file.path(root, "cache")
  archive <- file.path(cache, "archives", "msigdb.2026.1.zip")
  dir.create(dirname(archive), recursive = TRUE)
  writeLines("old corrupt archive", archive, useBytes = TRUE)
  original <- readBin(archive, "raw", n = file.info(archive)$size)
  testthat::local_mocked_bindings(
    lisa_msigdb_2026_1_download_file = function(...) {
      stop("synthetic transfer failure")
    },
    .package = "lisaR"
  )

  expect_error(
    prepare_lisa_msigdb_resource(accept_terms = TRUE, cache_root = cache),
    "synthetic transfer failure"
  )
  expect_identical(readBin(archive, "raw", n = file.info(archive)$size), original)
  expect_false(file.exists(file.path(cache, "resource_registry.tsv")))
})

test_that("MSigDB acquisition reuses a previously verified archive", {
  root <- tempfile("lisa-msigdb-archive-hit-")
  dir.create(root)
  cache <- file.path(root, "cache")
  archive <- file.path(cache, "archives", "msigdb.2026.1.zip")
  dir.create(dirname(archive), recursive = TRUE)
  writeLines("fixture verified by the narrowly mocked verifier", archive,
             useBytes = TRUE)
  contract <- lisaR:::lisa_msigdb_2026_1_term2gene_contract()
  testthat::local_mocked_bindings(
    lisa_assert_msigdb_2026_1_zip = function(source_zip, contract) source_zip,
    lisa_msigdb_2026_1_download_file = function(...) {
      stop("a verified archive must not download again")
    },
    .package = "lisaR"
  )

  expect_identical(
    lisaR:::lisa_msigdb_2026_1_download_zip(cache, contract),
    archive
  )
})

test_that("prepare_lisa_msigdb_resource reuses only the exact registered identity", {
  root <- tempfile("lisa-msigdb-idempotent-")
  dir.create(root)
  cache <- file.path(root, "cache")
  dir.create(cache)
  registry <- file.path(cache, "resource_registry.tsv")
  contract <- lisaR:::lisa_msigdb_2026_1_term2gene_contract()
  row <- data.frame(
    logical_id = "msigdb_term2gene", version = "2026.1",
    species = "Homo sapiens", modality = "all", schema = "term2gene@1",
    sha256 = contract$output_sha256, compatibility = "lisaR>=0.6.0",
    approved_origin = lisaR:::lisa_msigdb_2026_1_term2gene_origin(contract),
    artifact = "homo_sapiens/all/term2gene.tsv", stringsAsFactors = FALSE
  )
  utils::write.table(row, registry, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  resolve_calls <- 0L
  testthat::local_mocked_bindings(
    lisa_msigdb_2026_1_download_zip = function(...) {
      stop("existing registered resource must not download")
    },
    lisa_resolve_dictionary_resource = function(...) {
      resolve_calls <<- resolve_calls + 1L
      list(path = "/synthetic/already-verified-term2gene.tsv",
           resource_id = "msigdb_term2gene@2026.1")
    },
    lisa_msigdb_2026_1_term2gene_table = function(...) {
      stop("existing resource must not be transformed")
    },
    .package = "lisaR"
  )

  reused <- prepare_lisa_msigdb_resource(
    cache_root = cache, registry_path = registry
  )
  expect_identical(reused$resource_id, "msigdb_term2gene@2026.1")
  expect_identical(resolve_calls, 1L)
  expect_identical(reused$source_zip_sha256, contract$source_zip_sha256)

  row$approved_origin <- paste(row$approved_origin, "different provenance")
  utils::write.table(row, registry, sep = "\t", quote = FALSE,
                     row.names = FALSE)
  expect_error(
    prepare_lisa_msigdb_resource(
      cache_root = cache, registry_path = registry
    ),
    "LISA-RESOURCE-024.*different metadata or SHA-256"
  )
})

test_that("install_lisa_resource isolates species and supports modality fallback", {
  root <- tempfile("lisa-resource-dimensions-")
  dir.create(root)
  source <- file.path(root, "term2gene.tsv")
  utils::write.table(
    data.frame(
      gs_collection = "C5", gs_subcollection = "GO:BP",
      gs_name = "GS_A", gs_exact_source = "TEST", gene_symbol = "GENE1",
      stringsAsFactors = FALSE
    ),
    source, sep = "\t", quote = FALSE, row.names = FALSE
  )
  cache <- file.path(root, "cache")
  digest <- lisaR:::lisa_sha256_file(source)
  install_one <- function(species, modality) {
    install_lisa_resource(
      source = source,
      logical_id = "fixture_multidimensional",
      version = "1.0.0",
      species = species,
      modality = modality,
      schema = "term2gene@1",
      approved_origin = "authorised local test fixture",
      expected_sha256 = digest,
      cache_root = cache
    )
  }

  human_all <- install_one("Homo sapiens", "all")
  mouse_all <- install_one("Mus musculus", "all")
  human_exact <- install_one("Homo sapiens", "targeted")

  expect_false(identical(human_all$path, mouse_all$path))
  expect_false(identical(human_all$path, human_exact$path))
  expect_match(human_all$path, "homo_sapiens/all/term2gene[.]tsv$")
  expect_match(mouse_all$path, "mus_musculus/all/term2gene[.]tsv$")
  expect_match(human_exact$path, "homo_sapiens/targeted/term2gene[.]tsv$")

  resolve <- function(species, modality = NULL) {
    lisaR:::lisa_resolve_dictionary_resource(
      "fixture_multidimensional", "1.0.0", species,
      modality = modality, registry = human_all$registry_path,
      cache_root = cache
    )
  }
  expect_identical(resolve("Homo sapiens", "targeted")$modality, "targeted")
  expect_identical(
    resolve("Homo sapiens", "global proteomic")$modality, "all"
  )
  expect_identical(resolve("Homo sapiens")$modality, "all")
  expect_identical(resolve("Mus musculus", "targeted")$modality, "all")
})

test_that("install_lisa_resource supports nested registries and literal NA text", {
  root <- tempfile("lisa-resource-nested-registry-")
  dir.create(root)
  source <- file.path(root, "term2gene.tsv")
  utils::write.table(
    data.frame(
      gs_collection = "C5", gs_subcollection = "GO:BP",
      gs_name = "GS_A", gs_exact_source = "TEST", gene_symbol = "GENE1",
      stringsAsFactors = FALSE
    ),
    source, sep = "\t", quote = FALSE, row.names = FALSE
  )
  cache <- file.path(root, "cache")
  registry <- file.path(cache, "registries", "resource_registry.tsv")
  unrelated_run_root <- file.path(root, "active-run")
  dir.create(unrelated_run_root)
  withr::local_options(list(lisaR.run_root = unrelated_run_root))
  installed <- install_lisa_resource(
    source = source,
    logical_id = "fixture_nested_registry",
    version = "1.0.0",
    species = "Homo sapiens",
    modality = "all",
    schema = "term2gene@1",
    approved_origin = "NA",
    expected_sha256 = lisaR:::lisa_sha256_file(source),
    cache_root = cache,
    registry_path = registry
  )

  reread <- lisaR:::lisa_read_dictionary_registry(installed$registry_path)
  expect_identical(reread$approved_origin, "NA")
  expect_false(file.exists(paste0(installed$registry_path, ".lock")))
  expect_identical(
    lisaR:::lisa_resolve_dictionary_resource(
      "fixture_nested_registry", "1.0.0", "Homo sapiens",
      modality = "targeted", registry = installed$registry_path,
      cache_root = cache
    )$path,
    installed$path
  )
})

test_that("install_lisa_resource fails closed on digest and registry locks", {
  root <- tempfile("lisa-resource-lock-")
  dir.create(root)
  source <- file.path(root, "term2gene.tsv")
  utils::write.table(
    data.frame(
      gs_collection = "C5", gs_subcollection = "GO:BP",
      gs_name = "GS_A", gs_exact_source = "TEST", gene_symbol = "GENE1",
      stringsAsFactors = FALSE
    ),
    source, sep = "\t", quote = FALSE, row.names = FALSE
  )
  cache <- file.path(root, "cache")
  digest <- lisaR:::lisa_sha256_file(source)
  expect_error(
    install_lisa_resource(
      source = source, logical_id = "fixture_lock", version = "1.0.0",
      species = "Homo sapiens", modality = "all",
      schema = "term2gene@1", approved_origin = "local negative test",
      expected_sha256 = paste0(substr(digest, 1L, 63L),
                                if (substr(digest, 64L, 64L) == "0") "1" else "0"),
      cache_root = cache
    ),
    "LISA-RESOURCE-025.*does not match expected_sha256"
  )

  dir.create(cache, recursive = TRUE)
  lock <- file.path(cache, "resource_registry.tsv.lock")
  dir.create(lock)
  expect_error(
    install_lisa_resource(
      source = source, logical_id = "fixture_lock", version = "1.0.0",
      species = "Homo sapiens", modality = "all",
      schema = "term2gene@1", approved_origin = "local negative test",
      expected_sha256 = digest, cache_root = cache
    ),
    "LISA-RESOURCE-026.*locked"
  )
})

test_that("install_lisa_resource rejects linked cache components", {
  skip_on_os("windows")
  root <- tempfile("lisa-resource-linked-cache-")
  dir.create(root)
  source <- file.path(root, "term2gene.tsv")
  utils::write.table(
    data.frame(
      gs_collection = "C5", gs_subcollection = "GO:BP",
      gs_name = "GS_A", gs_exact_source = "TEST", gene_symbol = "GENE1",
      stringsAsFactors = FALSE
    ),
    source, sep = "\t", quote = FALSE, row.names = FALSE
  )
  cache <- file.path(root, "cache")
  outside <- file.path(root, "outside")
  dir.create(cache)
  dir.create(outside)
  expect_true(file.symlink(outside, file.path(cache, "fixture_linked")))

  expect_error(
    install_lisa_resource(
      source = source, logical_id = "fixture_linked", version = "1.0.0",
      species = "Homo sapiens", modality = "all",
      schema = "term2gene@1", approved_origin = "local negative test",
      expected_sha256 = lisaR:::lisa_sha256_file(source),
      cache_root = cache
    ),
    "LISA-RESOURCE-008.*symbolic link"
  )
  expect_length(list.files(outside, all.files = TRUE, no.. = TRUE), 0L)
  expect_false(file.exists(file.path(cache, "resource_registry.tsv.lock")))
})
