custom_resource_example <- function() {
  root <- system.file("extdata", "custom-resources", package = "lisaR")
  stopifnot(nzchar(root))
  read_tsv <- function(name) {
    utils::read.delim(
      file.path(root, name), sep = "\t", header = TRUE, quote = "",
      comment.char = "", check.names = FALSE, stringsAsFactors = FALSE,
      colClasses = "character", na.strings = character()
    )
  }
  list(
    root = root,
    manifest = read_tsv("EXAMPLE_RESOURCE_MANIFEST.tsv"),
    dictionary = read_tsv("lisa_dictionary_custom_example_v1_0.tsv"),
    legacy_dictionary = read_tsv("example_custom_dictionary.tsv"),
    category_map = read_tsv("example_custom_category_map.tsv"),
    term2gene = read_tsv("example_custom_term2gene.tsv")
  )
}

test_that("custom-resource vignette fixtures match their reviewed manifest", {
  fixture <- custom_resource_example()
  manifest_columns <- c(
    "logical_id", "version", "species", "modality", "schema", "sha256",
    "compatibility", "approved_origin", "artifact"
  )
  expect_identical(names(fixture$manifest), manifest_columns)
  expect_identical(nrow(fixture$manifest), 3L)
  expect_true(all(grepl("example", fixture$manifest$logical_id, fixed = TRUE)))
  expect_true(all(grepl(
    "synthetic", fixture$manifest$approved_origin, ignore.case = TRUE
  )))
  expect_no_error(lisaR:::lisa_read_dictionary_registry(fixture$manifest))
  observed <- vapply(fixture$manifest$artifact, function(artifact) {
    lisaR:::lisa_sha256_file(file.path(fixture$root, artifact))
  }, character(1))
  expect_identical(unname(observed), fixture$manifest$sha256)

  expect_identical(
    names(fixture$dictionary),
    c(
      "universe", "gene_set_id", "gene_set_name", "source_id",
      "category_id", "category_display_name", "tier"
    )
  )
  expect_identical(
    names(fixture$term2gene),
    c(
      "gs_collection", "gs_subcollection", "gs_name", "gs_exact_source",
      "gene_symbol"
    )
  )
  expect_true(all(c(
    "category_id", "display_name", "macrogroup_id", "macrogroup_name",
    "macrogroup_order", "category_order_within_macrogroup"
  ) %in% names(fixture$category_map)))
  expect_setequal(
    fixture$dictionary$universe,
    c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS")
  )
  expect_false("HALLMARKS" %in% fixture$dictionary$universe)
  expect_true(all(fixture$dictionary$tier == "custom"))
  expect_false("LISA_score" %in% names(fixture$dictionary))
  expect_gt(
    max(table(fixture$dictionary$gene_set_id)),
    1L
  )
  expect_true(all(
    unique(fixture$dictionary$gene_set_id) %in% fixture$term2gene$gs_name
  ))
  expect_true(all(
    unique(fixture$dictionary$category_id) %in%
      fixture$category_map$category_id
  ))
  labels <- stats::setNames(
    fixture$category_map$display_name, fixture$category_map$category_id
  )
  expect_identical(
    fixture$dictionary$category_display_name,
    unname(labels[fixture$dictionary$category_id])
  )
})

test_that("custom-resource vignette bundle installs and validates offline", {
  fixture <- custom_resource_example()
  cache <- withr::local_tempdir(pattern = "lisa-custom-vignette-")
  registry <- file.path(cache, "resource_registry.tsv")

  installed <- lapply(seq_len(nrow(fixture$manifest)), function(i) {
    row <- fixture$manifest[i, , drop = FALSE]
    install_lisa_resource(
      source = file.path(fixture$root, row$artifact),
      logical_id = row$logical_id,
      version = row$version,
      species = row$species,
      modality = row$modality,
      schema = row$schema,
      approved_origin = row$approved_origin,
      expected_sha256 = row$sha256,
      compatibility = row$compatibility,
      cache_root = cache,
      registry_path = registry
    )
  })
  expect_true(all(vapply(installed, function(x) file.exists(x$path), logical(1))))

  withr::local_options(list(
    lisaR.dictionary_registry = registry,
    lisaR.dictionary_cache_root = cache,
    lisaR.shared_dictionary_root = NULL
  ))
  refs <- stats::setNames(
    paste0(fixture$manifest$logical_id, "@", fixture$manifest$version),
    fixture$manifest$schema
  )
  validated <- validate_lisa_custom_resources(
    dictionary_resource = refs[["lisa_dictionary@2"]],
    category_map_resource = refs[["category_map@1"]],
    term2gene_resource = refs[["term2gene@1"]],
    species = "Homo sapiens"
  )
  expect_true(validated$valid)
  expect_identical(validated$species, "Homo sapiens")
  expect_identical(validated$summary$n_dictionary_assignments, 6L)
  expect_identical(validated$summary$n_categories, 3L)
  expect_identical(validated$summary$n_gene_sets, 4L)
  expect_identical(validated$summary$n_term2gene_memberships, 12L)
  expect_identical(validated$summary$n_genes, 6L)
})

test_that("custom-resource vignette invalid cases fail with stable errors", {
  fixture <- custom_resource_example()
  ids <- c(
    dictionary_resource = "lisa_example_custom_dictionary@2.0.0",
    category_map_resource = "lisa_example_custom_category_map@1.0.0",
    term2gene_resource = "lisa_example_custom_term2gene@1.0.0"
  )

  duplicate <- rbind(fixture$dictionary, fixture$dictionary[1L, ])
  expect_error(
    lisaR:::lisa_validate_custom_dictionary(
      duplicate, ids[["dictionary_resource"]]
    ),
    "LISA-CUSTOM-002.*duplicate"
  )

  hallmark <- fixture$dictionary
  hallmark$universe[[1L]] <- "HALLMARKS"
  expect_error(
    lisaR:::lisa_validate_custom_dictionary(
      hallmark, ids[["dictionary_resource"]]
    ),
    "LISA-CUSTOM-003.*HALLMARKS"
  )

  # The retained legacy `lisa_dictionary@1` fixture still enforces the
  # historical rule, so the transitional reader is covered in both directions
  # .
  expect_silent(
    lisaR:::lisa_validate_custom_dictionary(
      fixture$legacy_dictionary, "lisa_example_custom_dictionary@1.0.0"
    )
  )
  scored <- fixture$legacy_dictionary
  scored$LISA_score[[1L]] <- "4"
  expect_error(
    lisaR:::lisa_validate_custom_dictionary(
      scored, "lisa_example_custom_dictionary@1.0.0"
    ),
    "LISA-CUSTOM-003.*LISA_score"
  )

  orphan <- fixture$dictionary
  orphan$category_id[[1L]] <- "CATEGORY_NOT_IN_MAP"
  expect_error(
    lisaR:::lisa_validate_custom_resource_relations(
      orphan, fixture$category_map, fixture$term2gene, ids
    ),
    "LISA-CUSTOM-005.*orphan"
  )

  missing_set <- fixture$dictionary
  missing_set$gene_set_id[[1L]] <- "SET_NOT_IN_TERM2GENE"
  expect_error(
    lisaR:::lisa_validate_custom_resource_relations(
      missing_set, fixture$category_map, fixture$term2gene, ids
    ),
    "LISA-CUSTOM-006.*absent from TERM2GENE"
  )

  conflicting_map <- rbind(
    fixture$category_map,
    transform(
      fixture$category_map[1L, ],
      macrogroup_id = "SECOND_GROUP",
      macrogroup_name = "Second group",
      macrogroup_order = "3"
    )
  )
  expect_error(
    lisaR:::lisa_validate_custom_category_map(
      conflicting_map, ids[["category_map_resource"]]
    ),
    "LISA-CUSTOM-004.*multiple supercategories"
  )

  duplicate_membership <- rbind(
    fixture$term2gene, fixture$term2gene[1L, ]
  )
  expect_error(
    lisaR:::lisa_validate_custom_term2gene(
      duplicate_membership, ids[["term2gene_resource"]]
    ),
    "LISA-CUSTOM-002.*duplicate exact membership"
  )
})
