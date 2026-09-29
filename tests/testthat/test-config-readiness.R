readiness_fixture <- function() {
  root <- tempfile("lisa-readiness-")
  dir.create(root)
  input <- file.path(root, "de.tsv")
  matrix <- file.path(root, "matrix.tsv")
  writeLines("symbol\tlog2FC\tpadj\nA\t1\t0.01", input)
  writeLines("symbol\tsample_a\nA\t1", matrix)

  cache <- file.path(root, "cache")
  ids <- c(
    dictionary_resource = "fixture_dictionary",
    term2gene_resource = "fixture_term2gene",
    category_map_resource = "fixture_category_map"
  )
  artifacts <- c("dictionary.tsv", "term2gene.tsv", "category_map.tsv")
  schemas <- c("lisa_dictionary@2", "term2gene@1", "category_map@1")
  paths <- file.path(cache, ids, "1.0.0", artifacts)
  tables <- list(g2_dictionary(), g2_term2gene(), g2_category_map())
  for (i in seq_along(paths)) {
    g2_write_tsv(tables[[i]], paths[[i]])
  }
  registry <- data.frame(
    logical_id = unname(ids),
    version = "1.0.0",
    species = "Homo sapiens",
    modality = "transcriptomic/genomic",
    schema = schemas,
    sha256 = vapply(paths, lisaR:::lisa_sha256_file, character(1)),
    compatibility = "lisaR>=0.6.0",
    approved_origin = "local test fixture",
    artifact = artifacts,
    stringsAsFactors = FALSE
  )
  cfg <- list(
    pipeline = list(
      schema_version = "1.0.0",
      profile = "transcriptomic/genomic",
      evidence_mode = "full_de",
      output_dir = file.path(root, "results"),
      workers = 1L,
      dry_run = TRUE,
      dictionary_resource = paste0(ids[["dictionary_resource"]], "@1.0.0"),
      term2gene_resource = paste0(ids[["term2gene_resource"]], "@1.0.0"),
      category_map_resource = paste0(ids[["category_map_resource"]], "@1.0.0"),
      run_kegg_maps = FALSE,
      duplicate_policies = list(
        de_table_duplicate_policy = "error",
        matrix_duplicate_policy = "error",
        mapped_id_collision_policy = "error"
      )
    ),
    report = list(
      mode = "standard",
      formats = list(png = TRUE, svg = FALSE, pdf = FALSE),
      source_data = TRUE,
      recipes = FALSE
    ),
    collections = list("GOBP-C2"),
    single_de = list(list(
      analysis_id = "analysis_a",
      de_path = basename(input),
      expression_matrix_path = basename(matrix),
      matrix_feature_col = "symbol",
      species = "Homo sapiens",
      symbol_col = "symbol",
      logfc_col = "log2FC",
      padj_col = "padj"
    ))
  )
  config <- file.path(root, "study.json")
  jsonlite::write_json(cfg, config, auto_unbox = TRUE, pretty = TRUE)
  list(
    root = root, cache = cache, registry = registry, cfg = cfg,
    config = config, input = input, matrix = matrix
  )
}

test_that("path resolvers recognise portable absolute-path forms", {
  absolute <- c(
    "/srv/study/input.tsv",
    "C:/study/input.tsv",
    "C:\\study\\input.tsv",
    "//server/share/input.tsv",
    "\\\\server\\share\\input.tsv",
    "///srv/study/input.tsv"
  )
  relative <- "study/input.tsv"

  expect_true(all(vapply(
    absolute, lisaR:::lisa_is_absolute_path, logical(1)
  )))
  expect_false(any(vapply(
    relative, lisaR:::lisa_is_absolute_path, logical(1)
  )))
  expect_error(
    lisaR:::lisa_is_absolute_path("\\\\?\\C:\\study\\input.tsv"),
    "device-path"
  )
  expect_error(
    lisaR:::lisa_is_absolute_path("\\\\.\\C:\\study\\input.tsv"),
    "device-path"
  )
  expect_error(
    lisaR:::lisa_is_absolute_path("C:study/input.tsv"),
    "drive-relative"
  )
  expect_error(
    lisaR:::lisa_is_absolute_path("//server"),
    "Incomplete UNC"
  )
  expect_error(
    lisaR:::lisa_config_path("C:study/input.tsv", "D:/base"),
    "drive-relative"
  )
  expect_error(
    lisaR:::lisa_norm_path("//server", "D:/base"),
    "Incomplete UNC"
  )

  for (path in absolute) {
    expect_identical(lisaR:::lisa_config_path(path, "D:/base"), path)
    expect_identical(lisaR:::lisa_norm_path(path, "D:/base"), path)
  }
  expect_identical(
    lisaR:::lisa_config_path("study/input.tsv", "D:/base"),
    normalizePath("D:/base/study/input.tsv", winslash = "/", mustWork = FALSE)
  )
})

test_that("configuration readiness resolves all registered resources read-only", {
  fixture <- readiness_fixture()
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))

  validation <- validate_lisa_config(fixture$config, check_files = TRUE)

  expect_true(validation$structural_valid)
  expect_true(validation$inputs_ready)
  expect_true(validation$dependencies_ready)
  expect_true(validation$resources_ready)
  expect_true(validation$execution_ready)
  expect_true(validation$valid)
  expect_identical(
    validation$readiness$stage,
    c("structure", "inputs", "dependencies", "resources", "execution")
  )
  expect_true(all(validation$readiness$ready))
  expect_setequal(
    validation$resources$config_key,
    c("dictionary_resource", "term2gene_resource", "category_map_resource")
  )
  expect_identical(
    validation$resources$schema,
    c("lisa_dictionary@2", "term2gene@1", "category_map@1")
  )
  expect_true(all(validation$resources$ready))
  expect_identical(validation$resources$selection_source, rep("explicit", 3L))
  expect_true(all(file.exists(validation$resources$path)))
  expect_true(all(vapply(validation$resources$sha256, lisaR:::lisa_sha256_is_valid, logical(1))))
  expect_setequal(validation$inputs$input, c("de_path", "expression_matrix_path"))
  expect_false(dir.exists(fixture$cfg$pipeline$output_dir))
})

test_that("an unknown versioned resource cannot validate as ready", {
  fixture <- readiness_fixture()
  fixture$cfg$pipeline$term2gene_resource <- "missing_term2gene@999"
  jsonlite::write_json(
    fixture$cfg, fixture$config, auto_unbox = TRUE, pretty = TRUE
  )
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))

  validation <- validate_lisa_config(fixture$config, check_files = TRUE)

  expect_true(validation$structural_valid)
  expect_true(validation$inputs_ready)
  expect_false(validation$resources_ready)
  expect_false(validation$execution_ready)
  expect_false(validation$valid)
  missing <- validation$resources[
    validation$resources$config_key == "term2gene_resource", , drop = FALSE
  ]
  expect_false(missing$ready)
  expect_identical(missing$status, "not_ready")
  expect_match(missing$message, "did not resolve", fixed = TRUE)
  expect_error(
    run_lisa(fixture$config),
    "LISA-READINESS-001.*resources"
  )
  expect_false(dir.exists(fixture$cfg$pipeline$output_dir))
})

test_that("structural-only validation cannot claim execution readiness", {
  fixture <- readiness_fixture()
  fixture$cfg$pipeline$dictionary_resource <- "missing_dictionary@999"
  fixture$cfg$single_de[[1]]$de_path <- "missing-de.tsv"
  jsonlite::write_json(
    fixture$cfg, fixture$config, auto_unbox = TRUE, pretty = TRUE
  )

  validation <- validate_lisa_config(fixture$config, check_files = FALSE)

  expect_true(validation$structural_valid)
  expect_true(is.na(validation$inputs_ready))
  expect_true(is.na(validation$resources_ready))
  expect_false(validation$execution_ready)
  expect_false(validation$valid)
  expect_true(all(validation$inputs$status == "not_checked"))
  expect_true(all(validation$resources$status == "not_checked"))
  expect_true("selection_source" %in% names(validation$resources))
})

test_that("missing referenced inputs are returned as a non-ready stage", {
  fixture <- readiness_fixture()
  fixture$cfg$single_de[[1]]$expression_matrix_path <- "missing-matrix.tsv"
  jsonlite::write_json(
    fixture$cfg, fixture$config, auto_unbox = TRUE, pretty = TRUE
  )
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))

  validation <- validate_lisa_config(fixture$config, check_files = TRUE)

  expect_false(validation$inputs_ready)
  expect_true(validation$resources_ready)
  expect_false(validation$execution_ready)
  expect_false(validation$valid)
  missing <- validation$inputs[!validation$inputs$ready, , drop = FALSE]
  expect_identical(missing$input, "expression_matrix_path")
  expect_identical(missing$status, "not_ready")
})

test_that("source_data allowlist failures are part of input readiness", {
  fixture <- readiness_fixture()
  source <- file.path(fixture$root, "source.tsv")
  writeLines("value\nsource", source)
  fixture$cfg$source_data <- list(list(
    path = basename(source), role = "readiness fixture",
    target_subdir = "registered_sources"
  ))
  fixture$cfg$allowlisted_source_paths <- "different.tsv"
  jsonlite::write_json(
    fixture$cfg, fixture$config, auto_unbox = TRUE, pretty = TRUE
  )
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))

  blocked <- validate_lisa_config(fixture$config, check_files = TRUE)
  expect_false(blocked$inputs_ready)
  expect_true(blocked$resources_ready)
  expect_false(blocked$execution_ready)
  source_row <- blocked$inputs[blocked$inputs$input == "source_data", , drop = FALSE]
  expect_identical(source_row$status, "not_ready")
  expect_match(source_row$message, "exact allowlist match", fixed = TRUE)

  fixture$cfg$allowlisted_source_paths <- basename(source)
  jsonlite::write_json(
    fixture$cfg, fixture$config, auto_unbox = TRUE, pretty = TRUE
  )
  ready <- validate_lisa_config(fixture$config, check_files = TRUE)
  expect_true(ready$inputs_ready)
  expect_true(ready$execution_ready)
  source_row <- ready$inputs[
    ready$inputs$input == "source_data[1]", , drop = FALSE
  ]
  expect_identical(source_row$status, "ready")
  expect_identical(source_row$path, normalizePath(source, winslash = "/"))
})

test_that("duplicate-policy mapping files are part of input readiness", {
  fixture <- readiness_fixture()
  fixture$cfg$pipeline$duplicate_policies$de_table_duplicate_policy <- list(
    type = "mapping_file", mapping_file = "missing-mapping.tsv"
  )
  jsonlite::write_json(
    fixture$cfg, fixture$config, auto_unbox = TRUE, pretty = TRUE
  )
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))

  validation <- validate_lisa_config(fixture$config, check_files = TRUE)

  expect_false(validation$inputs_ready)
  expect_true(validation$resources_ready)
  expect_false(validation$execution_ready)
  mapping <- validation$inputs[
    grepl("mapping_file$", validation$inputs$input), , drop = FALSE
  ]
  expect_identical(mapping$status, "not_ready")
  expect_match(mapping$message, "mapping file", fixed = TRUE)
})

test_that("missing dependencies are returned as a non-ready stage", {
  fixture <- readiness_fixture()
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))
  testthat::local_mocked_bindings(
    lisa_config_dependency_requirements = function(...) data.frame(
      package = "fixtureMissingPackage",
      reason = "readiness fixture",
      manager = "CRAN",
      available = FALSE,
      stringsAsFactors = FALSE
    ),
    .package = "lisaR"
  )

  validation <- validate_lisa_config(fixture$config, check_files = TRUE)

  expect_true(validation$inputs_ready)
  expect_false(validation$dependencies_ready)
  expect_true(validation$resources_ready)
  expect_false(validation$execution_ready)
  expect_false(validation$valid)
  expect_identical(validation$dependencies$package, "fixtureMissingPackage")
  expect_identical(validation$dependencies$available, FALSE)
})

test_that("check_files is one explicit non-missing logical value", {
  fixture <- readiness_fixture()
  for (value in list("false", 0L, NA, c(TRUE, FALSE))) {
    expect_error(
      validate_lisa_config(fixture$config, check_files = value),
      "LISA-CONFIG-BOOL-001.*field=check_files"
    )
  }
})
