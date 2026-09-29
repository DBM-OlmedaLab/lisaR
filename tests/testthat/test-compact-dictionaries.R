test_that("bundled Quick Start resources expose a coherent custom contract", {
  root <- system.file("extdata", "quick-start", package = "lisaR")
  dictionary_path <- file.path(root, "lisa_dictionary_quickstart_v1_0.tsv")
  category_map_path <- file.path(root, "example_category_map.tsv")
  term2gene_path <- file.path(root, "example_term2gene_v2_0.tsv")
  expect_true(all(file.exists(c(
    dictionary_path, category_map_path, term2gene_path
  ))))

  dictionary <- lisaR:::read_lisa_tsv(dictionary_path)
  category_map <- lisaR:::read_lisa_tsv(category_map_path)
  term2gene <- lisaR:::read_lisa_tsv(term2gene_path)
  expect_true(all(c(
    "universe", "gene_set_id", "gene_set_name", "source_id",
    "category_id", "category_display_name", "tier"
  ) %in% names(dictionary)))
  expect_true(all(dictionary$tier == "custom"))
  # the canonical contract has no score column at all, so a new
  # custom resource never has to carry an empty one.
  expect_false("LISA_score" %in% names(dictionary))
  expect_identical(
    unique(dictionary$universe),
    c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS")
  )
  expect_identical(
    anyDuplicated(paste(
      dictionary$universe, dictionary$gene_set_id,
      dictionary$category_id, sep = "\r"
    )),
    0L
  )
  expect_identical(anyDuplicated(dictionary$gene_set_id), 0L)

  map_match <- match(dictionary$category_id, category_map$category_id)
  expect_false(anyNA(map_match))
  expect_identical(
    dictionary$category_display_name,
    category_map$display_name[map_match]
  )
  expect_true(all(unique(dictionary$gene_set_id) %in% term2gene$gs_name))
  expect_identical(
    sort(unique(term2gene$gs_name)),
    sort(unique(dictionary$gene_set_id))
  )
  gene_set_collection <- unique(term2gene[, c("gs_name", "gs_collection")])
  expect_identical(anyDuplicated(gene_set_collection$gs_name), 0L)
  dictionary_collection <- stats::setNames(
    dictionary$universe, dictionary$gene_set_id
  )
  expect_identical(
    unname(dictionary_collection[gene_set_collection$gs_name]),
    gene_set_collection$gs_collection
  )
})

test_that("bundled scientific dictionaries and category map resolve directly", {
  for (tier in c("core", "expanded")) {
    path <- lisaR:::lisa_builtin_dictionary_path(tier)
    expect_true(file.exists(path))
    expect_identical(
      lisaR:::lisa_builtin_registered_resource(
        paste0("lisa_dictionary_", tier), "1.0.0", "Homo sapiens"
      )$path,
      path
    )
  }
  map <- lisaR:::lisa_builtin_category_map_path()
  expect_true(file.exists(map))
  expect_identical(
    lisaR:::lisa_builtin_registered_resource(
      "lisa_category_map", "1.0.0", "Homo sapiens"
    )$path,
    map
  )
})
