test_that("semantic category guide has complete bounded verified coverage", {
  catalog_path <- system.file(
    "extdata", "semantic-catalog", "semantic_category_catalog.tsv",
    package = "lisaR"
  )
  provenance_path <- system.file(
    "extdata", "semantic-catalog", "SEMANTIC_CATALOG_PROVENANCE.tsv",
    package = "lisaR"
  )
  expect_true(file.exists(catalog_path))
  expect_true(file.exists(provenance_path))

  catalog <- utils::read.delim(
    catalog_path, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", na.strings = character()
  )
  expected_columns <- c(
    "universe", "macrogroup_id", "macrogroup_name", "macrogroup_order",
    "supercategory_description", "category_id", "category_display_name",
    "category_order_within_macrogroup", "category_description",
    "core_member_gene_sets", "example_count", "example_gene_set_1",
    "example_gene_set_2", "example_gene_set_3", "source_dictionary",
    "source_dictionary_sha256", "expanded_dictionary_sha256",
    "category_map_sha256", "verification_status"
  )
  expect_identical(names(catalog), expected_columns)
  expect_identical(nrow(catalog), 158L)
  expect_identical(
    as.integer(table(factor(catalog$universe, levels = c(
      "GOBP-C2", "GOMF", "GOCC", "PATHWAYS"
    )))),
    c(40L, 37L, 12L, 69L)
  )
  expect_identical(length(unique(catalog$macrogroup_id)), 42L)
  expect_identical(
    anyDuplicated(paste(catalog$universe, catalog$category_id)), 0L
  )
  expect_false(any(c(
    "AXL_MERTK_TAM_SIGNALING", "ERK_TRANSCRIPTIONAL_OUTPUT",
    "ROS1_SIGNALING", "OTHER_UNCLASSIFIED"
  ) %in% catalog$category_id))

  examples <- as.matrix(catalog[c(
    "example_gene_set_1", "example_gene_set_2", "example_gene_set_3"
  )])
  observed_count <- rowSums(!is.na(examples) & nzchar(examples))
  expect_true(all(observed_count == as.integer(catalog$example_count)))
  expect_true(all(observed_count >= 1L & observed_count <= 3L))
  expect_equal(sum(observed_count), 375)
  expect_true(all(vapply(seq_len(nrow(examples)), function(i) {
    selected <- examples[i, !is.na(examples[i, ]) & nzchar(examples[i, ])]
    !anyDuplicated(selected)
  }, logical(1L))))

  expect_true(all(as.integer(catalog$core_member_gene_sets) >= observed_count))
  expect_true(all(nzchar(catalog$category_description)))
  expect_true(all(nzchar(catalog$supercategory_description)))
  expect_true(all(grepl("[.]$", catalog$category_description)))
  expect_true(all(grepl("[.]$", catalog$supercategory_description)))
  # C6 rebuilt the catalog against the score-free resources, so its recorded
  # provenance now names `lisa_core@2.0.0` / `lisa_expanded@2.0.0` instead of the
  # pre-C5 `@0.1.0` labels. These digests are pinned rather than recomputed so
  # that an accidental rebuild against the wrong resource is caught here; they
  # were cross-checked against the bytes of
  # inst/extdata/dictionaries/lisa_dictionary_core_runtime_v1_0.tsv and
  # lisa_dictionary_expanded_runtime_v1_0.tsv and against
  # RESOURCE_BUNDLE_MANIFEST.tsv.
  # The catalog's editorial content is unchanged by the rebuild; that is
  # asserted separately in test-score-free-dictionaries.R.
  expect_identical(
    unique(catalog$source_dictionary), "lisa_dictionary_core@1.0.0"
  )
  expect_identical(
    unique(catalog$source_dictionary_sha256),
    "f24b5bd8d9ac66b6d64c1a91c56ab567aa64a37f05dee5ce0ed3011b2bed8994"
  )
  expect_identical(
    unique(catalog$expanded_dictionary_sha256),
    "0914c2d2e40369041a2313e913898f9c8e4065779e3973b076b70ea473045736"
  )
  expect_identical(
    unique(catalog$category_map_sha256),
    "d61fcb2e1d40d0448d459daaab952975477203bb00b76cb53145aea8f61333b1"
  )
  expect_identical(unique(catalog$verification_status),
                   "verified_exact_membership")

  provenance <- utils::read.delim(
    provenance_path, sep = "\t", quote = "", comment.char = "",
    check.names = FALSE, stringsAsFactors = FALSE,
    colClasses = "character", na.strings = character()
  )
  expect_identical(nrow(provenance), 3L)
  expect_identical(provenance$bundled, rep("TRUE", 3L))
  expect_setequal(provenance$sha256, c(
    "f24b5bd8d9ac66b6d64c1a91c56ab567aa64a37f05dee5ce0ed3011b2bed8994",
    "0914c2d2e40369041a2313e913898f9c8e4065779e3973b076b70ea473045736",
    "d61fcb2e1d40d0448d459daaab952975477203bb00b76cb53145aea8f61333b1"
  ))
  # The category map itself is untouched by C6: its digest above is unchanged.
  # Only its label moved, and that move was C5's relabel of `@0.1.1` to `@1.0.0`
  # over identical bytes (see RESOURCE_MIGRATION_MANIFEST.tsv). The catalog had
  # kept the pre-C5 label; rebuilding it recorded the current identity.
  expect_true("lisa_category_map@1.0.0" %in% provenance$source)
})
