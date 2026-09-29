test_that("G2 omitted and explicit scientific defaults resolve identically", {
  expected_ids <- c(
    dictionary_resource = "fixture_lisa_core@0.1.0",
    term2gene_resource = "fixture_msigdb_term2gene@2026.1",
    category_map_resource = "fixture_lisa_category_map@0.1.1"
  )

  for (species in c("Homo sapiens", "Mus musculus")) {
    fixture <- g2_resource_fixture(species)
    omitted <- g2_resolve_resources(fixture = fixture)
    explicit <- g2_resolve_resources(list(
      dictionary_resource = expected_ids[["dictionary_resource"]],
      term2gene_resource = expected_ids[["term2gene_resource"]],
      category_map_resource = expected_ids[["category_map_resource"]]
    ), fixture)

    identity_columns <- c(
      "config_key", "resource_id", "species", "modality", "schema",
      "sha256", "path", "compatibility"
    )
    expect_identical(
      omitted$status[, identity_columns],
      explicit$status[, identity_columns]
    )
    expect_identical(
      stats::setNames(omitted$status$resource_id, omitted$status$config_key),
      expected_ids
    )
    expect_identical(omitted$status$selection_source, rep("explicit", 3L))
    expect_identical(explicit$status$selection_source, rep("explicit", 3L))
    expect_true(all(omitted$status$ready))
  }
})

test_that("G2 validation exposes effective default IDs in the normalized config", {
  fixture <- g2_resource_fixture()
  de_path <- file.path(fixture$root, "de.tsv")
  g2_write_tsv(data.frame(
    symbol = "GENE1", log2FC = 1, padj = 0.01,
    stringsAsFactors = FALSE
  ), de_path)
  config <- list(
    pipeline = list(
      schema_version = "1.0.0",
      profile = "transcriptomic/genomic",
      evidence_mode = "full_de",
      dictionary_resource = "fixture_lisa_core@0.1.0",
      term2gene_resource = "fixture_msigdb_term2gene@2026.1",
      category_map_resource = "fixture_lisa_category_map@0.1.1",
      output_dir = file.path(fixture$root, "result"),
      duplicate_policies = list(
        de_table_duplicate_policy = "error",
        matrix_duplicate_policy = "error",
        mapped_id_collision_policy = "error"
      )
    ),
    collections = list("GOBP-C2"),
    single_de = list(list(
      analysis_id = "default_resources",
      de_path = de_path,
      species = "Homo sapiens",
      symbol_col = "symbol",
      logfc_col = "log2FC",
      padj_col = "padj"
    ))
  )
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))

  result <- validate_lisa_config(config, check_files = TRUE)
  expect_true(result$resources_ready)
  expect_identical(
    unlist(result$config$pipeline[c(
      "dictionary_resource", "term2gene_resource", "category_map_resource"
    )], use.names = TRUE),
    c(
      dictionary_resource = "fixture_lisa_core@0.1.0",
      term2gene_resource = "fixture_msigdb_term2gene@2026.1",
      category_map_resource = "fixture_lisa_category_map@0.1.1"
    )
  )
  expect_identical(result$resources$selection_source, rep("explicit", 3L))
})

test_that("G2 configured execution passes resolved resources to the engine", {
  fixture <- g2_resource_fixture()
  de_path <- file.path(fixture$root, "de.tsv")
  g2_write_tsv(data.frame(
    symbol = paste0("GENE", seq_len(120L)),
    log2FC = seq(-1, 1, length.out = 120L),
    padj = seq(0.001, 0.9, length.out = 120L),
    stringsAsFactors = FALSE
  ), de_path)
  output <- file.path(fixture$root, "configured-result")
  config <- list(
    pipeline = list(
      schema_version = "1.0.0",
      profile = "transcriptomic/genomic",
      evidence_mode = "full_de",
      dictionary_resource = "fixture_lisa_core@0.1.0",
      term2gene_resource = "fixture_msigdb_term2gene@2026.1",
      category_map_resource = "fixture_lisa_category_map@0.1.1",
      output_dir = output,
      dry_run = TRUE,
      duplicate_policies = list(
        de_table_duplicate_policy = "error",
        matrix_duplicate_policy = "error",
        mapped_id_collision_policy = "error"
      )
    ),
    collections = list("GOBP-C2"),
    single_de = list(list(
      analysis_id = "resource_propagation",
      de_path = de_path,
      species = "Homo sapiens",
      symbol_col = "symbol",
      logfc_col = "log2FC",
      padj_col = "padj"
    ))
  )
  config_path <- file.path(fixture$root, "study.json")
  jsonlite::write_json(config, config_path, auto_unbox = TRUE)
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))
  captured <- new.env(parent = emptyenv())
  testthat::local_mocked_bindings(
    lisa_assert_config_dependencies = function(...) data.frame(),
    run_lisa_pipeline = function(...) {
      captured$args <- list(...)
      list(status = "captured")
    },
    .package = "lisaR"
  )

  result <- suppressWarnings(lisaR:::run_lisa_pipeline_from_config(config_path))
  expect_identical(result$gate, "PLAN_PASS")
  expect_true(captured$args$registered_category_map)
  expect_identical(
    captured$args$dictionary_path,
    normalizePath(fixture$paths[["fixture_lisa_core"]], winslash = "/")
  )
  expect_identical(
    captured$args$term2gene,
    normalizePath(fixture$paths[["fixture_msigdb_term2gene"]], winslash = "/")
  )
  expect_identical(
    captured$args$category_map_path,
    normalizePath(fixture$paths[["fixture_lisa_category_map"]], winslash = "/")
  )
})

test_that("G2 dictionary tier precedence is deterministic", {
  fixture <- g2_resource_fixture()
  expanded <- g2_resolve_resources(list(
    dictionary_resource = "fixture_lisa_expanded@0.1.0"
  ), fixture)
  expect_identical(
    expanded$resolved$dictionary_resource$resource_id,
    "fixture_lisa_expanded@0.1.0"
  )
  expect_identical(
    expanded$status$selection_source[
      expanded$status$config_key == "dictionary_resource"
    ],
    "explicit"
  )

  expect_no_error(g2_resolve_resources(list(
      dictionary_resource = "fixture_lisa_core@0.1.0"
  ), fixture))
  expect_no_error(g2_resolve_resources(list(
    dictionary_resource = "fixture_lisa_expanded@0.1.0"
  ), fixture))
  expect_no_error(g2_resolve_resources(list(
    dictionary_resource = "lab_dictionary@1.0.0",
    term2gene_resource = "lab_term2gene@1.0.0",
    category_map_resource = "lab_category_map@1.0.0"
  ), fixture))
  expect_error(
    g2_resolve_resources(list(
      lisa_dictionary = "core",
      dictionary_resource = "lab_dictionary@1.0.0",
      term2gene_resource = "lab_term2gene@1.0.0",
      category_map_resource = "lab_category_map@1.0.0"
    ), fixture),
    "LISA-RESOURCE-018.*custom.*Repair:"
  )
  expect_error(
    g2_resolve_resources(list(
      dictionary_resource = "lab_dictionary",
      term2gene_resource = "lab_term2gene@1.0.0",
      category_map_resource = "lab_category_map@1.0.0"
    ), fixture),
    "LISA-RESOURCE-019.*version.*Repair:"
  )

  future_core <- fixture
  current <- future_core$registry$logical_id == "fixture_lisa_core"
  row <- future_core$registry[current, , drop = FALSE]
  row$version <- "0.3.0"
  source <- future_core$paths[["fixture_lisa_core"]]
  target <- file.path(
    future_core$cache, "fixture_lisa_core", row$version, row$artifact
  )
  g2_write_tsv(read.delim(
    source, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
  ), target)
  row$sha256 <- lisaR:::lisa_sha256_file(target)
  future_core$registry <- rbind(future_core$registry, row)
  future <- g2_resolve_resources(list(
    dictionary_resource = "fixture_lisa_core@0.3.0"
  ), future_core)
  expect_true(future$ready)
  expect_identical(future$dictionary_tier, "custom")
})

test_that("G2 missing default TERM2GENE fails locally with one repair action", {
  fixture <- g2_resource_fixture()
  fixture$registry <- fixture$registry[
    fixture$registry$logical_id != "fixture_msigdb_term2gene", , drop = FALSE
  ]

  expect_error(
    g2_resolve_resources(fixture = fixture),
    "LISA-RESOURCE-006.*fixture_msigdb_term2gene@2026[.]1.*Repair:"
  )
  status <- g2_resolve_resources(
    fixture = fixture, stop_on_error = FALSE
  )$status
  missing <- status[
    status$config_key == "term2gene_resource", , drop = FALSE
  ]
  expect_false(missing$ready)
  expect_identical(missing$selection_source, "explicit")
  expect_match(missing$message, "fixture_msigdb_term2gene@2026[.]1")
  expect_match(missing$message, "Repair:", fixed = TRUE)
  expect_false(grepl("download|http", missing$message, ignore.case = TRUE))
})

test_that("G2 missing dictionary default has one repair across selection forms", {
  fixture <- g2_resource_fixture()
  fixture$registry <- fixture$registry[
    fixture$registry$logical_id != "fixture_lisa_core", , drop = FALSE
  ]
  selections <- list(
    automatic = list(),
    explicit = list(dictionary_resource = "fixture_lisa_core@0.1.0")
  )

  for (source in names(selections)) {
    expect_error(
      g2_resolve_resources(selections[[source]], fixture),
      "LISA-RESOURCE-006.*fixture_lisa_core@0[.]1[.]0.*Repair:"
    )
    status <- g2_resolve_resources(
      selections[[source]], fixture, stop_on_error = FALSE
    )$status
    missing <- status[
      status$config_key == "dictionary_resource", , drop = FALSE
    ]
    expect_false(missing$ready)
    expect_identical(missing$selection_source, "explicit")
    expect_match(missing$message, "LISA-RESOURCE-006", fixed = TRUE)
  }
})

test_that("G2 validates registered many-to-many custom bundles offline", {
  for (species in c("Homo sapiens", "Mus musculus")) {
    result <- g2_validate_custom_resources(g2_resource_fixture(species))
    expect_true(result$valid)
    expect_identical(result$species, species)
    expect_true(all(result$resources$ready))
    expect_identical(result$summary$n_dictionary_assignments, 3L)
    expect_identical(result$summary$n_categories, 2L)
    expect_identical(result$summary$n_gene_sets, 2L)
    expect_identical(result$summary$n_term2gene_memberships, 4L)
    expect_identical(result$summary$n_genes, 3L)
  }
})

test_that("G2 custom dictionary validation rejects invalid shape and semantics", {
  missing_column <- g2_resource_fixture()
  dictionary <- g2_dictionary()
  dictionary$source_id <- NULL
  missing_column <- g2_replace_resource_table(
    missing_column, "lab_dictionary", dictionary
  )
  expect_error(
    g2_validate_custom_resources(missing_column),
    "LISA-CUSTOM-001.*source_id.*Repair:"
  )

  empty_id <- g2_resource_fixture()
  dictionary <- g2_dictionary()
  dictionary$gene_set_id[[1L]] <- ""
  empty_id <- g2_replace_resource_table(empty_id, "lab_dictionary", dictionary)
  expect_error(
    g2_validate_custom_resources(empty_id),
    "LISA-CUSTOM-001.*gene_set_id.*Repair:"
  )

  duplicate <- g2_resource_fixture()
  dictionary <- g2_dictionary()
  duplicate <- g2_replace_resource_table(
    duplicate, "lab_dictionary", rbind(dictionary, dictionary[1L, ])
  )
  expect_error(
    g2_validate_custom_resources(duplicate),
    "LISA-CUSTOM-002.*duplicate.*Repair:"
  )

  bad_universe <- g2_resource_fixture()
  dictionary <- g2_dictionary()
  dictionary$universe[[1L]] <- "HALLMARKS"
  bad_universe <- g2_replace_resource_table(
    bad_universe, "lab_dictionary", dictionary
  )
  expect_error(
    g2_validate_custom_resources(bad_universe),
    "LISA-CUSTOM-003.*universe.*Repair:"
  )

  # Two distinct contracts, both still enforced:
  #
  # (a) A legacy `lisa_dictionary@1` custom resource keeps the historical rule:
  #     the column may exist but must be empty. A populated value is refused.
  scored <- g2_resource_fixture()
  dictionary <- g2_dictionary(legacy_score = TRUE)
  dictionary$LISA_score <- 2
  scored <- g2_replace_resource_table(
    scored, "lab_dictionary", dictionary, schema = "lisa_dictionary@1"
  )
  expect_error(
    g2_validate_custom_resources(scored),
    "LISA-CUSTOM-003.*LISA_score.*Repair:"
  )

  # (b) A legacy custom resource with an EMPTY score column is still accepted,
  #     so the transitional reader really is transitional and not a silent
  #     rejection of every old resource.
  legacy_ok <- g2_replace_resource_table(
    g2_resource_fixture(), "lab_dictionary", g2_dictionary(legacy_score = TRUE),
    schema = "lisa_dictionary@1"
  )
  expect_silent(g2_validate_custom_resources(legacy_ok))

  # (c) A table must match the schema it declares. Declaring the canonical
  #     `lisa_dictionary@2` while carrying the legacy column is refused, so the
  #     registry can never describe a resource in terms its bytes contradict.
  expect_error(
    lisaR:::lisa_validate_custom_dictionary(
      g2_dictionary(legacy_score = TRUE), "lab_dictionary@1.0.0",
      schema = "lisa_dictionary@2"
    ),
    "LISA-CUSTOM-001.*declares schema 'lisa_dictionary@2'.*carries a LISA_score"
  )
  # ...and the converse, so neither direction can drift unnoticed.
  expect_error(
    lisaR:::lisa_validate_custom_dictionary(
      g2_dictionary(), "lab_dictionary@1.0.0", schema = "lisa_dictionary@1"
    ),
    "LISA-CUSTOM-001.*declares schema 'lisa_dictionary@1'.*omits a LISA_score"
  )

  # (d) The SCOPE requirement itself: a new custom resource is accepted with no
  #     score column at all, empty or otherwise.
  expect_silent(
    lisaR:::lisa_validate_custom_dictionary(
      g2_dictionary(), "lab_dictionary@1.0.0", schema = "lisa_dictionary@2"
    )
  )

  reserved_tier <- g2_resource_fixture()
  dictionary <- g2_dictionary()
  dictionary$tier <- "core"
  reserved_tier <- g2_replace_resource_table(
    reserved_tier, "lab_dictionary", dictionary
  )
  expect_error(
    g2_validate_custom_resources(reserved_tier),
    "LISA-CUSTOM-003.*tier.*Repair:"
  )

  unsafe_id <- g2_resource_fixture()
  dictionary <- g2_dictionary()
  dictionary$category_id[[1L]] <- "CAT A"
  unsafe_id <- g2_replace_resource_table(
    unsafe_id, "lab_dictionary", dictionary
  )
  expect_error(
    g2_validate_custom_resources(unsafe_id),
    "LISA-CUSTOM-001.*non-portable ID.*Repair:"
  )
})

test_that("G2 custom bundle validation rejects map and TERM2GENE defects", {
  duplicate_map <- g2_resource_fixture()
  category_map <- g2_category_map()
  duplicate_map <- g2_replace_resource_table(
    duplicate_map, "lab_category_map",
    rbind(category_map, category_map[1L, ])
  )
  expect_error(
    g2_validate_custom_resources(duplicate_map),
    "LISA-CUSTOM-002.*category_id.*Repair:"
  )

  duplicate_order <- g2_resource_fixture()
  category_map <- g2_category_map()
  category_map$macrogroup_id[[2L]] <- category_map$macrogroup_id[[1L]]
  category_map$macrogroup_name[[2L]] <- category_map$macrogroup_name[[1L]]
  category_map$macrogroup_order[[2L]] <- category_map$macrogroup_order[[1L]]
  category_map$category_order_within_macrogroup[[2L]] <- 1L
  duplicate_order <- g2_replace_resource_table(
    duplicate_order, "lab_category_map", category_map
  )
  expect_error(
    g2_validate_custom_resources(duplicate_order),
    "LISA-CUSTOM-004.*order.*Repair:"
  )

  inconsistent_group <- g2_resource_fixture()
  category_map <- g2_category_map()
  category_map$macrogroup_id[[2L]] <- category_map$macrogroup_id[[1L]]
  inconsistent_group <- g2_replace_resource_table(
    inconsistent_group, "lab_category_map", category_map
  )
  expect_error(
    g2_validate_custom_resources(inconsistent_group),
    "LISA-CUSTOM-004.*inconsistent name or order.*Repair:"
  )

  duplicate_group_order <- g2_resource_fixture()
  category_map <- g2_category_map()
  category_map$macrogroup_order[[2L]] <- category_map$macrogroup_order[[1L]]
  duplicate_group_order <- g2_replace_resource_table(
    duplicate_group_order, "lab_category_map", category_map
  )
  expect_error(
    g2_validate_custom_resources(duplicate_group_order),
    "LISA-CUSTOM-004.*share one macrogroup_order.*Repair:"
  )

  unsafe_map_id <- g2_resource_fixture()
  category_map <- g2_category_map()
  category_map$macrogroup_id[[1L]] <- "SUPER A"
  unsafe_map_id <- g2_replace_resource_table(
    unsafe_map_id, "lab_category_map", category_map
  )
  expect_error(
    g2_validate_custom_resources(unsafe_map_id),
    "LISA-CUSTOM-001.*non-portable ID.*Repair:"
  )

  orphan <- g2_resource_fixture()
  dictionary <- g2_dictionary()
  dictionary$category_id[[1L]] <- "ORPHAN"
  orphan <- g2_replace_resource_table(orphan, "lab_dictionary", dictionary)
  expect_error(
    g2_validate_custom_resources(orphan),
    "LISA-CUSTOM-005.*category_id.*Repair:"
  )

  label_mismatch <- g2_resource_fixture()
  category_map <- g2_category_map()
  category_map$display_name[[1L]] <- "Different label"
  label_mismatch <- g2_replace_resource_table(
    label_mismatch, "lab_category_map", category_map
  )
  expect_error(
    g2_validate_custom_resources(label_mismatch),
    "LISA-CUSTOM-005.*display.*Repair:"
  )

  absent_set <- g2_resource_fixture()
  dictionary <- g2_dictionary()
  dictionary$gene_set_id[[1L]] <- "ABSENT_SET"
  absent_set <- g2_replace_resource_table(
    absent_set, "lab_dictionary", dictionary
  )
  expect_error(
    g2_validate_custom_resources(absent_set),
    "LISA-CUSTOM-006.*gene_set_id.*TERM2GENE.*Repair:"
  )

  duplicate_membership <- g2_resource_fixture()
  term2gene <- g2_term2gene()
  duplicate_membership <- g2_replace_resource_table(
    duplicate_membership, "lab_term2gene",
    rbind(term2gene, term2gene[1L, ])
  )
  expect_error(
    g2_validate_custom_resources(duplicate_membership),
    "LISA-CUSTOM-002.*TERM2GENE.*duplicate.*Repair:"
  )

  conflicting_metadata <- g2_resource_fixture()
  term2gene <- g2_term2gene()
  conflict <- term2gene[1L, , drop = FALSE]
  conflict$gs_collection <- "C2"
  conflict$gs_subcollection <- "CP"
  conflicting_metadata <- g2_replace_resource_table(
    conflicting_metadata, "lab_term2gene", rbind(term2gene, conflict)
  )
  expect_error(
    g2_validate_custom_resources(conflicting_metadata),
    "LISA-CUSTOM-006.*inconsistent collection/provenance metadata.*Repair:"
  )
})

test_that("G2 enforces schema, species, modality and compatibility", {
  species <- g2_resource_fixture("Homo sapiens")
  expect_error(
    g2_validate_custom_resources(species, species = "Mus musculus"),
    "LISA-RESOURCE-.*species.*Repair:"
  )

  schema <- g2_resource_fixture()
  selected <- schema$registry$logical_id == "lab_dictionary"
  schema$registry$schema[selected] <- "category_map@1"
  expect_error(
    g2_validate_custom_resources(schema),
    "LISA-RESOURCE-017.*schema.*Repair:"
  )

  extra_column <- g2_resource_fixture()
  extra_column$registry$administrative_note <- "must live elsewhere"
  expect_error(
    g2_validate_custom_resources(extra_column),
    "LISA-RESOURCE-001.*unsupported column.*Repair:"
  )

  modality <- g2_resource_fixture()
  expect_error(
    g2_resolve_resources(list(profile = "targeted"), modality),
    "LISA-RESOURCE-021.*modality.*Repair:"
  )

  incompatible <- g2_resource_fixture(compatibility = "lisaR>=99.0.0")
  expect_error(
    g2_validate_custom_resources(incompatible),
    "LISA-RESOURCE-022.*compatib.*Repair:"
  )

  malformed <- g2_resource_fixture(compatibility = "future maybe")
  expect_error(
    g2_validate_custom_resources(malformed),
    "LISA-RESOURCE-022.*compatib.*Repair:"
  )
})

test_that("G2 accepts profile-independent resources with one exact modality", {
  fixture <- g2_resource_fixture()
  shared <- fixture$registry$logical_id %in% c(
    "lab_term2gene", "lab_category_map"
  )
  fixture$registry$modality[shared] <- "all"

  validated <- g2_validate_custom_resources(fixture)
  expect_true(validated$valid)
  expect_setequal(
    validated$resources$modality,
    c("transcriptomic/genomic", "all")
  )

  resolved <- g2_resolve_resources(list(
    dictionary_resource = "lab_dictionary@1.0.0",
    term2gene_resource = "lab_term2gene@1.0.0",
    category_map_resource = "lab_category_map@1.0.0"
  ), fixture)
  expect_true(resolved$ready)
  expect_identical(
    resolved$resolved$dictionary_resource$modality,
    "transcriptomic/genomic"
  )
  expect_identical(resolved$resolved$term2gene_resource$modality, "all")
  expect_identical(resolved$resolved$category_map_resource$modality, "all")

  incompatible <- fixture
  row <- incompatible$registry$logical_id == "lab_category_map"
  incompatible$registry$modality[row] <- "targeted"
  expect_error(
    g2_validate_custom_resources(incompatible),
    "LISA-RESOURCE-021.*incompatible modalities.*Repair:"
  )
})

test_that("G2 validates relations for reviewed default bundles", {
  missing_category <- g2_resource_fixture()
  category_map <- g2_category_map()
  missing_category <- g2_replace_resource_table(
    missing_category, "fixture_lisa_category_map",
    category_map[category_map$category_id != "CAT_B", , drop = FALSE]
  )
  expect_error(
    g2_resolve_resources(fixture = missing_category),
    "LISA-CUSTOM-005.*orphan dictionary category_id.*Repair:"
  )
  status <- g2_resolve_resources(
    fixture = missing_category, stop_on_error = FALSE
  )
  expect_false(status$ready)
  expect_true(all(status$status$ready == FALSE))

  missing_gene_set <- g2_resource_fixture()
  term2gene <- g2_term2gene()
  missing_gene_set <- g2_replace_resource_table(
    missing_gene_set, "fixture_msigdb_term2gene",
    term2gene[term2gene$gs_name != "GS_B", , drop = FALSE]
  )
  expect_error(
    g2_resolve_resources(fixture = missing_gene_set),
    "LISA-CUSTOM-006.*gene_set_id.*TERM2GENE.*Repair:"
  )
  status <- g2_resolve_resources(
    fixture = missing_gene_set, stop_on_error = FALSE
  )
  expect_false(status$ready)
  expect_true(all(status$status$ready == FALSE))
})

test_that("G2 configured registries augment but cannot redefine sample resources", {
  fixture <- g2_resource_fixture()
  withr::local_options(list(
    lisaR.dictionary_registry = fixture$registry,
    lisaR.dictionary_cache_root = fixture$cache,
    lisaR.shared_dictionary_root = NULL
  ))
  sample <- lisaR:::lisa_resolve_pipeline_resources(
    pipeline = list(
      profile = "transcriptomic/genomic",
      dictionary_resource = "lisa_quickstart_dictionary@2.0.0",
      term2gene_resource = "lisa_quickstart_term2gene@1.0.0",
      category_map_resource = "lisa_quickstart_category_map@1.0.0"
    ),
    species = "Homo sapiens"
  )
  expect_true(sample$ready)
  expect_true(all(vapply(
    sample$resolved, function(resource) file.exists(resource$path), logical(1)
  )))
  expect_true(all(vapply(
    sample$resolved,
    function(resource) startsWith(
      resource$path,
      normalizePath(system.file(package = "lisaR"), winslash = "/")
    ),
    logical(1)
  )))

  conflict <- g2_resource_fixture()
  conflict_path <- file.path(
    conflict$cache, "lisa_quickstart_dictionary", "2.0.0", "dictionary.tsv"
  )
  g2_write_tsv(g2_dictionary(), conflict_path)
  row <- conflict$registry[conflict$registry$logical_id == "lab_dictionary", ]
  row$logical_id <- "lisa_quickstart_dictionary"
  row$version <- "2.0.0"
  row$sha256 <- lisaR:::lisa_sha256_file(conflict_path)
  conflict$registry <- rbind(conflict$registry, row)
  withr::local_options(list(
    lisaR.dictionary_registry = conflict$registry,
    lisaR.dictionary_cache_root = conflict$cache,
    lisaR.shared_dictionary_root = NULL
  ))
  expect_error(
    lisaR:::lisa_resolve_pipeline_resources(
      pipeline = list(
        profile = "transcriptomic/genomic",
        dictionary_resource = "lisa_quickstart_dictionary@2.0.0",
        term2gene_resource = "lisa_quickstart_term2gene@1.0.0",
        category_map_resource = "lisa_quickstart_category_map@1.0.0"
      ),
      species = "Homo sapiens"
    ),
    "LISA-RESOURCE-020.*redefin.*Repair:"
  )

  distinct_modality <- conflict$registry
  distinct_modality$modality[
    distinct_modality$logical_id == "lisa_quickstart_dictionary"
  ] <- "targeted"
  expect_no_error({
    active <- lisaR:::lisa_active_resource_registry(distinct_modality)
  })
  sample_rows <- active[
    active$logical_id == "lisa_quickstart_dictionary" &
      active$version == "2.0.0" &
      active$species == "Homo sapiens",
    , drop = FALSE
  ]
  expect_setequal(
    sample_rows$modality,
    c("transcriptomic/genomic", "targeted")
  )
})
