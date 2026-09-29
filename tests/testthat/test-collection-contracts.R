# Collection and resource validation contracts.
#
g1_all_collections <- function() {
  c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS", "HALLMARKS")
}

g1_example_collections <- function() {
  c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS")
}

g1_package_path <- function(...) {
  normalizePath(
    system.file(..., package = "lisaR"),
    mustWork = TRUE
  )
}

g1_pipeline_config <- function(collections = NULL, run_hallmarks = NULL) {
  pipeline <- list(
    schema_version = "1.0.0",
    profile = "targeted",
    evidence_mode = "full_de",
    duplicate_policies = list(
      de_table_duplicate_policy = "error",
      matrix_duplicate_policy = "error",
      mapped_id_collision_policy = "error"
    )
  )
  if (!is.null(run_hallmarks)) pipeline$run_hallmarks <- run_hallmarks
  out <- list(pipeline = pipeline)
  if (!is.null(collections)) out$collections <- as.list(collections)
  out
}

g1_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = ""
  )
}

# `score` is retained in the signature and ignored: these fixtures now use the
# canonical score-free `lisa_dictionary@2` contract , and tier
# alone carries what the score used to encode.
g1_dictionary <- function(tier = "custom", score = NA_real_) {
  data.frame(
    universe = c("GOBP-C2", "GOBP-C2", "PATHWAYS"),
    gene_set_id = c("GS_A", "GS_A", "GS_B"),
    gene_set_name = c("Gene set A", "Gene set A", "Gene set B"),
    source_id = c("TEST", "TEST", "TEST"),
    category_id = c("CAT_A", "CAT_B", "CAT_A"),
    category_display_name = c("Category A", "Category B", "Category A"),
    tier = rep(tier, 3L),
    stringsAsFactors = FALSE
  )
}

g1_category_map <- function() {
  data.frame(
    category_id = c("CAT_A", "CAT_B"),
    display_name = c("Category A", "Category B"),
    macrogroup_id = c("SUPER_A", "SUPER_B"),
    macrogroup_name = c("Supercategory A", "Supercategory B"),
    macrogroup_order = c(1L, 2L),
    category_order_within_macrogroup = c(1L, 1L),
    notes = c("Test fixture", "Test fixture"),
    stringsAsFactors = FALSE
  )
}

g1_term2gene <- function() {
  data.frame(
    gs_collection = c("C5", "C5", "C2", "C2"),
    gs_subcollection = c("GO:BP", "GO:BP", "CP", "CP"),
    gs_name = c("GS_A", "GS_A", "GS_B", "GS_B"),
    gs_exact_source = c("TEST", "TEST", "TEST", "TEST"),
    gene_symbol = c("GENE1", "GENE2", "GENE2", "GENE3"),
    stringsAsFactors = FALSE
  )
}

g1_resource_fixture <- function(species = "Homo sapiens") {
  root <- tempfile("lisa-g1-contract-")
  cache <- file.path(root, "cache")

  # These rows stand in for the registered scientific defaults, so they must
  # carry the identities lisa_pipeline_resource_defaults() actually selects.
  # Registering the retired identities instead would leave the real bundled
  # dictionary resolving against this fixture's tiny category map.
  definitions <- list(
    lisa_dictionary_core = list(
      version = "1.0.0", schema = "lisa_dictionary@2",
      artifact = "dictionary.tsv", data = g1_dictionary("core", 4)
    ),
    lisa_dictionary_expanded = list(
      version = "1.0.0", schema = "lisa_dictionary@2",
      artifact = "dictionary.tsv", data = g1_dictionary("expanded", 3)
    ),
    msigdb_term2gene = list(
      version = "2026.1", schema = "term2gene@1",
      artifact = "term2gene.tsv", data = g1_term2gene()
    ),
    lisa_category_map = list(
      version = "1.0.0", schema = "category_map@1",
      artifact = "category_map.tsv", data = g1_category_map()
    ),
    lab_dictionary = list(
      version = "1.0.0", schema = "lisa_dictionary@2",
      artifact = "dictionary.tsv", data = g1_dictionary()
    ),
    lab_term2gene = list(
      version = "1.0.0", schema = "term2gene@1",
      artifact = "term2gene.tsv", data = g1_term2gene()
    ),
    lab_category_map = list(
      version = "1.0.0", schema = "category_map@1",
      artifact = "category_map.tsv", data = g1_category_map()
    )
  )

  rows <- lapply(names(definitions), function(logical_id) {
    item <- definitions[[logical_id]]
    path <- file.path(cache, logical_id, item$version, item$artifact)
    g1_write_tsv(item$data, path)
    data.frame(
      logical_id = logical_id,
      version = item$version,
      species = species,
      modality = "transcriptomic/genomic",
      schema = item$schema,
      sha256 = lisaR:::lisa_sha256_file(path),
      compatibility = "lisaR>=0.6.0",
      approved_origin = "Local test fixture",
      artifact = item$artifact,
      stringsAsFactors = FALSE
    )
  })
  registry <- do.call(rbind, rows)
  rownames(registry) <- NULL

  refs <- stats::setNames(
    paste0(registry$logical_id, "@", registry$version),
    registry$logical_id
  )
  paths <- stats::setNames(
    file.path(cache, registry$logical_id, registry$version, registry$artifact),
    registry$logical_id
  )

  list(
    root = root,
    cache = cache,
    species = species,
    registry = registry,
    refs = refs,
    paths = paths
  )
}

g1_replace_resource_table <- function(fixture, logical_id, table) {
  g1_write_tsv(table, fixture$paths[[logical_id]])
  row <- fixture$registry$logical_id == logical_id
  fixture$registry$sha256[row] <- lisaR:::lisa_sha256_file(
    fixture$paths[[logical_id]]
  )
  fixture
}

g1_resolve <- function(pipeline, fixture) {
  if (is.null(pipeline$profile)) {
    pipeline$profile <- "transcriptomic/genomic"
  }
  lisaR:::lisa_resolve_pipeline_resources(
    pipeline = pipeline,
    species = fixture$species,
    registry = fixture$registry,
    cache_root = fixture$cache,
    shared_root = NULL
  )
}

g1_validate_custom <- function(fixture, species = fixture$species) {
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))
  validator <- getExportedValue("lisaR", "validate_lisa_custom_resources")
  validator(
    dictionary_resource = fixture$refs[["lab_dictionary"]],
    category_map_resource = fixture$refs[["lab_category_map"]],
    term2gene_resource = fixture$refs[["lab_term2gene"]],
    species = species
  )
}

test_that("Configuration has one canonical five-collection default", {
  expected <- g1_all_collections()

  expect_identical(
    lisaR:::lisa_builtin_collection_registry()$analysis_collection,
    expected
  )
  expect_identical(lisaR:::lisa_default_collections(), expected)
  expect_identical(lisaR:::normalize_lisa_collections("all"), expected)
  expect_identical(lisaR:::normalize_lisa_universes("all"), expected)
  expect_identical(lisaR:::lisa_scientific_report_collection_order(), expected)

  functions <- list(
    lisaR::run_lisa_de,
    lisaR:::run_LISA_DE,
    lisaR::run_lisa_contrast,
    lisaR:::run_LISA_contrast
  )
  for (fn in functions) {
    expect_identical(
      eval(formals(fn)$universes, envir = asNamespace("lisaR")),
      expected
    )
  }

  normalizer_body <- paste(
    deparse(body(lisaR:::normalize_lisa_universes)), collapse = "\n"
  )
  expect_match(normalizer_body, "normalize_lisa_collections", fixed = TRUE)
})

test_that("Configuration validates and records canonical collection selections", {
  expected <- g1_all_collections()

  schema <- jsonlite::read_json(
    g1_package_path("schema", "lisa-config.schema.json")
  )
  collection_schema <- schema$properties$collections
  expect_identical(unlist(collection_schema$items$enum), expected)
  expect_identical(unlist(collection_schema$default), expected)
  expect_identical(collection_schema$minItems, 1L)
  expect_true(collection_schema$uniqueItems)

  validated <- lisaR:::lisa_validate_pipeline_config(g1_pipeline_config())
  expect_identical(validated$collections, expected)

  without_hallmarks <- lisaR:::lisa_validate_pipeline_config(
    g1_pipeline_config(run_hallmarks = FALSE)
  )
  expect_identical(without_hallmarks$collections, g1_example_collections())

  expect_error(
    lisaR:::lisa_validate_pipeline_config(g1_pipeline_config(character())),
    "LISA-COLLECTION.*at least one.*Repair:"
  )
  expect_error(
    lisaR:::lisa_validate_pipeline_config(g1_pipeline_config("UNKNOWN")),
    "LISA-COLLECTION.*UNKNOWN.*Repair:"
  )
  expect_error(
    lisaR:::lisa_validate_pipeline_config(
      g1_pipeline_config(c("GOMF", "GOMF"))
    ),
    "LISA-COLLECTION.*duplicate.*Repair:"
  )
})

test_that("Example configurations select four bounded collections", {
  relative <- c(
    "examples/minimal-study.yaml",
    "examples/minimal-study.json",
    "examples/quick-start/study.yml",
    "examples/riaz-gse91061/config/riaz-gse91061.yml",
    "examples/riaz-gse91061/config/riaz-gse91061.json"
  )
  for (file in relative) {
    path <- do.call(g1_package_path, as.list(strsplit(file, "/", fixed = TRUE)[[1L]]))
    config <- if (grepl("[.]json$", file)) {
      jsonlite::read_json(path)
    } else {
      yaml::read_yaml(path)
    }
    expect_identical(
      unname(unlist(config$collections, use.names = FALSE)),
      g1_example_collections(),
      info = file
    )
  }
})

test_that("low-level semantic defaults resolve the registered resource bundle", {
  fixture <- g1_resource_fixture()
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))

  resolved <- lisaR:::lisa_resolve_low_level_resource_paths(
    species = fixture$species,
    lisa_dictionary = "core",
    universes = g1_example_collections()
  )
  expect_identical(
    resolved$lisa_dictionary_path,
    normalizePath(fixture$paths[["lisa_dictionary_core"]], winslash = "/")
  )
  expect_identical(
    resolved$term2gene_path,
    normalizePath(fixture$paths[["msigdb_term2gene"]], winslash = "/")
  )
  expect_identical(
    resolved$category_map_path,
    normalizePath(fixture$paths[["lisa_category_map"]], winslash = "/")
  )
  expect_true(resolved$registered_category_map)
  expect_identical(resolved$source, "registered_defaults")

  expect_error(
    lisaR:::lisa_resolve_low_level_resource_paths(
      species = fixture$species,
      lisa_dictionary = "core",
      universes = "GOBP-C2",
      term2gene_path = fixture$paths[["msigdb_term2gene"]]
    ),
    "LISA-RESOURCE-027.*all three"
  )
})
