# Regression contracts for installation, numeric configuration and plan parity.
# Fixtures are synthetic; the exact private map is checked in release evidence.

audit_code_config <- function() {
  list(
    pipeline = list(
      schema_version = "1.0.0", profile = "targeted",
      evidence_mode = "full_de", run_hallmarks = FALSE,
      output_dir = "results/example", dry_run = TRUE,
      gsea_padj_cutoff = 0.123456789,
      duplicate_policies = list(
        de_table_duplicate_policy = "error", matrix_duplicate_policy = "error",
        mapped_id_collision_policy = "error"
      )
    ),
    collections = list("GOBP-C2", "GOMF"),
    single_de = list(
      list(analysis_id = "A", de_path = "a.tsv", species = "Homo sapiens"),
      list(analysis_id = "B", de_path = "b.tsv", species = "Homo sapiens")
    ),
    contrasts = list(list(contrast_id = "A_vs_B", analysis_a = "A", analysis_b = "B")),
    report = list(mode = "standard", formats = list(png = TRUE, svg = FALSE, pdf = FALSE),
      source_data = TRUE, recipes = FALSE)
  )
}

audit_scoped_map <- function() {
  data.frame(
    category_id = c("CAT_A", "PATHWAY_B"),
    display_name = c("Category A", "Pathway B"),
    macrogroup_id = c("SEMANTIC_A", "PATHWAYS_B"),
    macrogroup_name = c("Semantic A", "Pathways B"),
    macrogroup_order = c(1L, 1L),
    category_order_within_macrogroup = c(1L, 1L),
    family_id = c("", "PATHWAYS_B"), pathway_id = c("", "PATHWAY_B"),
    stringsAsFactors = FALSE
  )
}

test_that("standard map installation respects scopes and preserves exact bytes", {
  root <- tempfile("audit-scoped-map-")
  dir.create(root)
  source <- file.path(root, "reviewed-map.tsv")
  utils::write.table(audit_scoped_map(), source, sep = "\t", quote = FALSE,
    row.names = FALSE, na = "")
  digest <- lisaR:::lisa_sha256_file(source)
  cache <- file.path(root, "cache")
  # Keep the built-in identity protected; use a distinct test-only version
  # for this synthetic scoped hierarchy, never redefine the reviewed map.
  expect_error(install_lisa_resource(
    source, "lisa_category_map", "1.0.0", "Homo sapiens", "all",
    "category_map@1", "Synthetic map must not replace built-in", digest,
    cache_root = file.path(root, "reserved-cache")
  ), "LISA-RESOURCE-020")
  # The pre-C5 identity stays reserved too. Before the version migration it was
  # a built-in row and LISA-RESOURCE-020 refused it; it must not become
  # installable just because the active version was renamed to 1.0.0.
  expect_error(install_lisa_resource(
    source, "lisa_category_map", "0.1.1", "Homo sapiens", "all",
    "category_map@1", "Superseded identity must stay reserved", digest,
    cache_root = file.path(root, "reserved-cache")
  ), "LISA-RESOURCE-030")
  fixture_version <- "0.0.0.9000"
  installed <- install_lisa_resource(
    source, "lisa_category_map", fixture_version, "Homo sapiens", "all",
    "category_map@1", "Synthetic independently scoped hierarchies", digest,
    cache_root = cache
  )
  expect_identical(lisaR:::lisa_sha256_file(installed$path), digest)
  resolved <- lisaR:::lisa_resolve_dictionary_resource(
    "lisa_category_map", fixture_version, "Homo sapiens", modality = "targeted",
    registry = installed$registry_path, cache_root = cache
  )
  expect_identical(lisaR:::lisa_sha256_file(resolved$path), digest)
  # Installation and bundle validation use the same dispatcher.
  expect_silent(lisaR:::lisa_validate_category_map(
    audit_scoped_map(), resolved$resource_id))
  expect_error(lisaR:::lisa_validate_category_map(audit_scoped_map(), "lab_map@1"),
    "share one macrogroup_order")
})

test_that("map ordering still rejects collisions within scope and broken identities", {
  x <- audit_scoped_map()
  x$family_id <- x$pathway_id <- c("", "")
  expect_error(lisaR:::lisa_validate_category_map(x, "lisa_category_map@0.1.1"),
    "share one macrogroup_order")
  x <- audit_scoped_map()
  x$family_id[2] <- "WRONG_GROUP"
  expect_error(lisaR:::lisa_validate_category_map(x, "lisa_category_map@0.1.1"),
    "do not match their hierarchy identities")
  x <- audit_scoped_map()
  x$pathway_id <- NULL
  expect_error(lisaR:::lisa_validate_category_map(x, "lisa_category_map@0.1.1"),
    "both family_id and pathway_id")
  x <- audit_scoped_map()
  extra <- x[1, ]
  extra$category_id <- "CAT_C"
  x <- rbind(x, extra)
  expect_error(lisaR:::lisa_validate_category_map(x, "lisa_category_map@0.1.1"),
    "duplicate within a macrogroup_id")
})

test_that("list, file and normalized configurations preserve scientific precision", {
  root <- tempfile("audit-config-precision-")
  dir.create(root)
  withr::local_dir(root)
  cfg <- audit_code_config()
  expected <- cfg$pipeline$gsea_padj_cutoff
  validated <- validate_lisa_config(cfg, check_files = FALSE)
  expect_identical(validated$gsea_padj_cutoff, expected)
  # Capture the exact config that would reach the engine after the public
  # list/file boundary, without launching analyses for this precision test.
  testthat::local_mocked_bindings(
    validate_lisa_config = function(config, check_files = TRUE) {
      raw <- lisaR:::read_lisa_pipeline_config(config)
      lisaR:::lisa_validate_pipeline_config(raw)
      list(execution_ready = TRUE)
    },
    run_lisa_pipeline_from_config = function(config) {
      raw <- lisaR:::read_lisa_pipeline_config(config)
      normalized <- tempfile("normalized-", fileext = ".json")
      lisaR:::lisa_write_normalized_config(raw, normalized)
      list(effective = raw$pipeline$gsea_padj_cutoff,
        normalized = jsonlite::fromJSON(normalized)$pipeline$gsea_padj_cutoff)
    }, .package = "lisaR"
  )
  from_list <- run_lisa(cfg)
  config_path <- file.path(root, "exact.json")
  jsonlite::write_json(cfg, config_path, auto_unbox = TRUE, digits = NA)
  from_file <- run_lisa(config_path)
  expect_identical(from_list, from_file)
  expect_identical(from_list$effective, expected)
  expect_identical(from_list$normalized, expected)
  boundary <- c(0.12345, 0.12348, 0.12351)
  expect_identical(boundary <= from_list$effective, boundary <= expected)
})

test_that("external relative indexes and aliases have the same output estimate as inline rows", {
  root <- tempfile("audit-plan-indexes-")
  dir.create(root)
  cfg <- audit_code_config()
  inline <- plan_lisa_outputs(cfg)
  for (key in c("single_de", "contrasts")) {
    table <- lisaR:::lisa_config_table(cfg[[key]])
    utils::write.table(table, file.path(root, paste0(key, ".tsv")),
      sep = "\t", quote = FALSE, row.names = FALSE)
  }
  for (top_level in c(FALSE, TRUE)) {
    external <- cfg
    external$single_de <- external$contrasts <- NULL
    if (top_level) {
      external$de_index_path <- "single_de.tsv"
      external$contrast_index_path <- "contrasts.tsv"
    } else {
      external$pipeline$de_index_path <- "single_de.tsv"
      external$pipeline$contrast_index_path <- "contrasts.tsv"
    }
    path <- file.path(root, paste0("config-", top_level, ".json"))
    jsonlite::write_json(external, path, auto_unbox = TRUE, digits = NA)
    # Caller cwd differs from the file's base; both aliases resolve there.
    expect_identical(plan_lisa_outputs(path), inline)
  }
  aliased <- cfg
  aliased$de_index <- aliased$single_de
  aliased$contrast_index <- aliased$contrasts
  aliased$single_de <- aliased$contrasts <- NULL
  expect_identical(plan_lisa_outputs(aliased), inline)
  expect_identical(inline$summary$analyses, 2L)
  expect_identical(inline$summary$contrasts, 1L)
  expect_gt(inline$summary$conservative_bytes, 0)
})

test_that("configured KEGG maps require an explicit immutable cache without disabling set figures", {
  cfg <- audit_code_config()
  cfg$pipeline$run_kegg_maps <- TRUE
  expect_error(validate_lisa_config(cfg, check_files = FALSE), "LISA-KEGG-016.*cache_only")
  cfg$pipeline$kegg_access_mode <- "cache_only"
  expect_error(validate_lisa_config(cfg, check_files = FALSE), "LISA-KEGG-017")
  cfg$pipeline$kegg_cache_root <- "recipient-kegg-cache"
  cfg$pipeline$kegg_snapshot_id <- "frozen-20260729"
  validated <- validate_lisa_config(cfg, check_files = FALSE)
  expect_true(validated$config$pipeline$run_kegg_maps)
  expect_identical(validated$config$pipeline$kegg_access_mode, "cache_only")
  expect_identical(validated$config$pipeline$kegg_snapshot_id, "frozen-20260729")
  resolved <- lisaR:::lisa_resolve_config(cfg)
  expect_identical(resolved$pipeline$kegg_access_mode, "cache_only")
  expect_identical(resolved$pipeline$kegg_snapshot_id, "frozen-20260729")
  expect_identical(basename(resolved$pipeline$kegg_cache_root), "recipient-kegg-cache")
  product_plan <- lisaR:::lisa_product_plan(
    lisaR:::lisa_evidence_contract(validated$config$pipeline$evidence_mode),
    run_gene_level = FALSE, run_kegg_maps = TRUE
  )
  expect_true(product_plan$compute[product_plan$product == "kegg_maps"])
  expect_false(product_plan$compute[product_plan$product == "gene_level"])
  bridge <- paste(deparse(body(lisaR:::run_lisa_pipeline_from_config)), collapse = "\n")
  expect_match(bridge, "run_gene_level = FALSE", fixed = TRUE)
  # KEGG maps are painted inside the same run (D5), not in a -kegg-maps sibling.
  expect_false(grepl("lisa_render_kegg_maps_extension", bridge, fixed = TRUE))
  expect_match(bridge, "kegg_maps = isTRUE(validated_config$kegg_maps$value)", fixed = TRUE)
  cfg$pipeline$run_kegg_maps <- FALSE
  cfg$report$mode <- "full"
  expect_true("kegg" %in% plan_lisa_outputs(cfg)$products$product)
  expect_false(lisaR:::lisa_validate_pipeline_config(cfg)$kegg_maps$value)
})

test_that("KEGG extension promotion rebases staged painter receipts", {
  root <- tempfile("kegg-extension-rebase-")
  dir.create(root)
  target <- normalizePath(
    file.path(tempdir(), "kegg-extension-promoted"),
    winslash = "/", mustWork = FALSE
  )
  status <- data.frame(message = file.path(root, "outputs", "map.png"),
    stringsAsFactors = FALSE)
  index_dir <- file.path(root, "outputs", "gene_level", "single_de", "A",
    "collection_GOBP-C2", "kegg_painter")
  index <- data.frame(output_png = file.path(root, "outputs", "map.png"),
    output_pdf = file.path(root, "outputs", "map.pdf"), stringsAsFactors = FALSE)
  nodes <- data.frame(output_png = index$output_png, output_pdf = index$output_pdf,
    stringsAsFactors = FALSE)
  lisaR:::write_lisa_tsv(status, file.path(root, "kegg_extension_status.tsv"))
  lisaR:::write_lisa_tsv(index, file.path(index_dir, "A_GOBP-C2_gene_level_kegg_pathway_painter_index.tsv"))
  lisaR:::write_lisa_tsv(nodes, file.path(index_dir, "A_GOBP-C2_gene_level_kegg_pathway_painter_nodes.tsv"))
  expect_true(lisaR:::lisa_rebase_kegg_extension_paths(root, target))
  expect_false(grepl(root, paste(lisaR:::read_lisa_tsv(file.path(root, "kegg_extension_status.tsv")), collapse = "\n"), fixed = TRUE))
  expect_true(grepl(target, paste(lisaR:::read_lisa_tsv(file.path(index_dir,
    "A_GOBP-C2_gene_level_kegg_pathway_painter_index.tsv")), collapse = "\n"), fixed = TRUE))
})
