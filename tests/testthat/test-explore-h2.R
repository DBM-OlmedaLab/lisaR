test_that("H2 requests require exact native-map and heatmap selectors", {
  expect_error(
    lisa_figure_request(
      "single_de", "GOBP-C2", "heatmap", "CAT_A", analysis_id = "A"
    ),
    "LISA-EXPLORE-011"
  )
  expect_error(
    lisa_figure_request(
      "single_de", "GOBP-C2", "kegg_pathway_map", "CAT_A",
      analysis_id = "A"
    ),
    "LISA-EXPLORE-008"
  )
  expect_error(
    lisa_figure_request(
      "single_de", "GOBP-C2", "volcano", "CAT_A", analysis_id = "A",
      entity = "hsa04110"
    ),
    "LISA-EXPLORE-010"
  )

  raw <- lisa_figure_request(
    "single_de", "GOBP-C2", "heatmap", "CAT_A", analysis_id = "A",
    variant = "raw"
  )
  map_a <- lisa_figure_request(
    "single_de", "GOBP-C2", "kegg_pathway_map", "CAT_A",
    analysis_id = "A", entity = "hsa04110"
  )
  map_b <- lisa_figure_request(
    "single_de", "GOBP-C2", "kegg_pathway_map", "CAT_B",
    analysis_id = "A", entity = "hsa04110"
  )
  expect_identical(lisaR:::lisa_explore_selection(raw)$selection$variants, "raw")
  expect_identical(
    lisaR:::lisa_explore_selection(map_a)$selection$entities, "hsa04110"
  )
  expect_identical(
    lisaR:::lisa_explore_identity_vector(map_a),
    lisaR:::lisa_explore_identity_vector(map_b)
  )
  expect_false(identical(lisa_request_id(map_a), lisa_request_id(map_b)))
})

test_that("the exact catalogue expands variants and carries map identity", {
  catalog <- data.frame(
    unit_type = "single_de", analysis_id = "A", contrast_id = "",
    collection = "GOBP-C2", category_id = "CAT_A",
    supercategory_id = "SUPER", product = "heatmap",
    summary_path = "summary.tsv", gsea_path = "gsea.tsv",
    de_path = "de.tsv", contrast_path = "", entity = "", variant = "",
    stringsAsFactors = FALSE
  )
  maps <- data.frame(
    analysis_id = "A", collection = "GOBP-C2", category_id = "CAT_A",
    kegg_id = "hsa04110", kegg_title = "Cell cycle", rank = 7L,
    kegg_cache_root = "/cache", kegg_snapshot_id = "snapshot-a",
    kegg_species = "Homo sapiens", max_abs_log2fc = "0.5",
    color_power = "1", kegg_resource_digest = paste(rep("a", 64), collapse = ""),
    stringsAsFactors = FALSE
  )
  expanded <- lisaR:::lisa_extension_exact_catalog(
    catalog, list(heatmap_variants = c("zscore", "raw"), kegg_maps = maps)
  )
  heatmaps <- expanded[expanded$product == "heatmap", , drop = FALSE]
  pathway <- expanded[expanded$product == "kegg_pathway_map", , drop = FALSE]
  expect_setequal(heatmaps$variant, c("zscore", "raw"))
  expect_identical(nrow(pathway), 1L)
  expect_identical(pathway$entity, "hsa04110")
  expect_identical(pathway$kegg_resource_digest, maps$kegg_resource_digest)
})

test_that("catalogue matching and reuse identity include variant and resources", {
  catalog <- data.frame(
    unit_type = "single_de", analysis_id = "A", contrast_id = "",
    collection = "GOBP-C2", category_id = "CAT_A", product = "heatmap",
    entity = "", variant = c("zscore", "raw"), stringsAsFactors = FALSE
  )
  raw <- lisa_figure_request(
    "single_de", "GOBP-C2", "heatmap", "CAT_A", analysis_id = "A",
    variant = "raw"
  )
  hit <- lisaR:::lisa_explore_catalog_row(list(source_run = "unused"), raw,
                                           catalog)
  expect_identical(nrow(hit), 1L)
  expect_identical(hit$variant, "raw")

  row <- data.frame(
    kegg_species = "Homo sapiens", kegg_snapshot_id = "snapshot-a",
    kegg_cache_root = "/cache", max_abs_log2fc = "0.5", color_power = "1",
    kegg_resource_digest = paste(rep("a", 64), collapse = ""),
    stringsAsFactors = FALSE
  )
  changed <- row
  changed$kegg_resource_digest <- paste(rep("b", 64), collapse = "")
  expect_false(identical(
    lisaR:::lisa_explore_resource_digest(row),
    lisaR:::lisa_explore_resource_digest(changed)
  ))
})

test_that("prepared KEGG indexes fail closed when source or snapshot changes", {
  ws <- list(source_manifest_hash = paste(rep("1", 64), collapse = ""))
  declaration <- list(
    cache_root = "/cache", snapshot_id = "snapshot-a",
    species = "Homo sapiens", organism = "hsa"
  )
  context <- lisaR:::lisa_explore_kegg_index_context(ws, declaration)
  index <- data.frame(kegg_id = "hsa04110", stringsAsFactors = FALSE)
  for (field in names(context)) index[[field]] <- context[[field]]
  expect_true(lisaR:::lisa_explore_kegg_index_is_current(
    ws, index, declaration
  ))

  changed_source <- ws
  changed_source$source_manifest_hash <- paste(rep("2", 64), collapse = "")
  expect_false(lisaR:::lisa_explore_kegg_index_is_current(
    changed_source, index, declaration
  ))
  changed_snapshot <- declaration
  changed_snapshot$snapshot_id <- "snapshot-b"
  expect_false(lisaR:::lisa_explore_kegg_index_is_current(
    ws, index, changed_snapshot
  ))
  index$kegg_species <- NULL
  expect_false(lisaR:::lisa_explore_kegg_index_is_current(
    ws, index, declaration
  ))
})

test_that("map associations deduplicate and map files are product-selected", {
  entries <- list(list(key = "fig-one", associations = list("CAT_A")))
  rewritten <- lisaR:::lisa_explore_index_upsert(entries, list(
    key = "fig-one", state = "available", associations = list("CAT_B")
  ))
  expect_identical(rewritten[[1L]]$associations, list("CAT_A", "CAT_B"))

  files <- c(
    "07_hsa04110_Cell_cycle_kegg_base.png",
    "07_hsa04110_Cell_cycle_painted.png",
    "07_hsa04110_Cell_cycle_painted.pdf",
    "07_hsa04110_Cell_cycle_painted_source.tsv",
    "07_hsa04110_Cell_cycle_painted_recipe.R"
  )
  expect_identical(
    lisaR:::lisa_explore_pick_product_file(files, "kegg_pathway_map", "png"),
    "07_hsa04110_Cell_cycle_painted.png"
  )
  expect_identical(
    lisaR:::lisa_explore_pick_product_file(
      files, "kegg_pathway_map", "source_data"
    ),
    "07_hsa04110_Cell_cycle_painted_source.tsv"
  )
})
