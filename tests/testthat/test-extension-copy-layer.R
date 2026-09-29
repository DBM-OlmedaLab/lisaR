copy_layer_source_run <- function() {
  root <- tempfile("lisaR-copy-layer-source-")
  collection <- file.path(root, "outputs", "single_de", "A", "collection_GOBP-C2")
  for (part in c("lisa_tables", "enrichment", "inputs")) {
    dir.create(file.path(collection, part), recursive = TRUE, showWarnings = FALSE)
  }
  dir.create(file.path(root, "config"), recursive = TRUE)
  categories <- data.frame(category_id = "CAT_A", macrogroup_id = "SUPER_1",
    category_display_name = "Category A", macrogroup_name = "Super one",
    macrogroup_order = 1L, category_order_within_macrogroup = 1L,
    n_genesets = 1L, color = "#0072B2", stringsAsFactors = FALSE)
  gsea <- data.frame(pathway = "SET_A", padj = 0.01, NES = 1.5,
    category_id = "CAT_A", macrogroup_id = "SUPER_1",
    category_display_name = "Category A", macrogroup_name = "Super one",
    macrogroup_order = 1L, category_order_within_macrogroup = 1L,
    color = "#0072B2", stringsAsFactors = FALSE)
  lisaR:::write_lisa_tsv(categories,
    file.path(collection, "lisa_tables", "A_semantic_GSEA_category_summary.tsv"))
  lisaR:::write_lisa_tsv(gsea,
    file.path(collection, "enrichment", "A_GSEA_semantic_annotated.tsv"))
  lisaR:::write_lisa_tsv(data.frame(symbol = c("G1", "G2")),
    file.path(collection, "inputs", "A_standardized_DE.tsv"))
  lisaR:::write_lisa_tsv(data.frame(analysis_id = "A", species = "Homo sapiens"),
    file.path(root, "config", "de_index.tsv"))
  lisaR:::write_lisa_tsv(data.frame(key = names(lisa_test_run_contract()),
    value = unname(lisa_test_run_contract()), stringsAsFactors = FALSE),
    file.path(root, "run_contract.tsv"))
  evidence_dir <- file.path(root, "report_pages", "evidence", "A", "GOBP-C2")
  dir.create(evidence_dir, recursive = TRUE)
  payload <- list(metadata = list(analysis_id = "A", collection = "GOBP-C2",
      tier = "standard", positive_contrast = "A", gsea_padj_cutoff = 0.05,
      de_padj_cutoff = 0.05, max_sets = 10L, max_genes = 10L,
      schema_version = "1.1", universe_ledger_supplied = FALSE,
      excluded_unclassified_rows = 0L, overlap_export = "displayed pairs",
      overlap_scope = "default displayed significant sets"),
    categories = list(list(category_id = "CAT_A", category_display_name = "Category A",
      category_description = "Small browser fixture", macrogroup_id = "SUPER_1",
      macrogroup_name = "Super one", support_state = "no_significant_sets",
      mean_significant_set_NES = NULL, n_mapped_sets = 0L, n_evaluable_sets = 0L,
      n_significant_sets = 0L, n_significant_positive_sets = 0L,
      n_significant_negative_sets = 0L, n_positive_sets = 0L, n_negative_sets = 0L,
      n_unique_le_genes = 0L, n_significant_unique_le_genes = 0L,
      n_le_assignments = 0L, n_significant_le_assignments = 0L)),
    sets = list(), genes = list(), leading_edges = list(), downloads = list(),
    figure_index = list(), member_figure_index = list(), formats = list(),
    source_data = FALSE, recipes = FALSE)
  writeLines(paste0('<html><body><script type="application/json" id="evidence-data">',
    jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"), "</script></body></html>"),
    file.path(evidence_dir, "index.html"), useBytes = TRUE)
  source_manifest <- lisaR:::lisa_run_manifest(root)
  lisaR:::write_lisa_tsv(source_manifest, file.path(root, "run_manifest.tsv"))
  lisaR:::write_lisa_tsv(data.frame(event = "validated", stringsAsFactors = FALSE),
    file.path(root, "run_events.tsv"))
  root
}

test_that("saved evidence references canonical artifacts with intact manifest and cold recipe", {
  skip_if_not_installed("ggplot2")
  source <- copy_layer_source_run()
  selection <- list(selection = list(mode = "selected", categories = "CAT_A",
      products = "member_gene_sets"),
    report = list(mode = "selected", formats = list(png = TRUE, svg = FALSE, pdf = TRUE),
      source_data = TRUE, recipes = TRUE))
  output <- tempfile("lisaR-canonical-evidence-")
  render_lisa_categories(source, selection, output)
  inventory <- lisaR:::read_lisa_tsv(file.path(output, "extension_inventory.tsv"))
  manifest <- lisaR:::read_lisa_tsv(file.path(output, "extension_manifest.tsv"))
  expect_identical(nrow(inventory), 4L)
  expect_setequal(tolower(tools::file_ext(inventory$path)), c("png", "pdf", "tsv", "r"))
  manifest_artifacts <- manifest[manifest$path %in% inventory$path, , drop = FALSE]
  manifest_artifacts <- manifest_artifacts[match(inventory$path, manifest_artifacts$path), , drop = FALSE]
  expect_identical(as.character(manifest_artifacts$sha256), as.character(inventory$sha256))
  expect_identical(as.numeric(manifest_artifacts$bytes), as.numeric(inventory$bytes))
  expect_false(any(grepl("/assets/category_products/", manifest$path, fixed = TRUE)))

  page <- file.path(output, "evidence", "A", "GOBP-C2", "index.html")
  evidence <- lisaR:::lisa_extension_read_payload(page, "evidence-data")
  assets <- unlist(lapply(evidence$category_products, `[[`, "assets"), recursive = FALSE)
  expect_setequal(vapply(assets, `[[`, character(1L), "source_path"), inventory$path)
  for (asset in assets) {
    row <- inventory[inventory$path == asset$source_path, , drop = FALSE]
    expect_identical(nrow(row), 1L)
    expect_identical(asset$sha256, as.character(row$sha256))
    expect_true(file.exists(file.path(dirname(page), utils::URLdecode(asset$href))))
  }
  expect_false(dir.exists(file.path(dirname(page), "assets", "category_products")))

  recipe_asset <- assets[[which(vapply(assets, function(asset)
    identical(tolower(asset$format), "r"), logical(1L)))]]
  source_asset <- assets[[which(vapply(assets, function(asset)
    identical(tolower(asset$format), "tsv"), logical(1L)))]]
  cold_output <- tempfile(fileext = ".png")
  # The recipe is run "cold" -- a fresh `Rscript --vanilla`, as a reader
  # would run it -- but it is still a lisaR recipe: line 776 calls
  # `lisaR:::lisa_read_figure_source_tsv()`, so the child needs a library that
  # contains lisaR. When the suite runs against an isolated candidate library
  # set with `.libPaths()` *in this process only*, the child inherits nothing and
  # fails to find the package. That was a harness defect, not a defect of the
  # recipe layer and not a scientific defect: handing the child the same library
  # this process is using is what "cold" was always meant to mean.
  #
  # `--vanilla` implies `--no-environ`, which stops R reading .Renviron files; it
  # does not stop R honouring an R_LIBS that is already in the child's
  # environment, which is what `env =` sets here.
  cold <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla",
    file.path(dirname(page), utils::URLdecode(recipe_asset$href)),
    file.path(dirname(page), utils::URLdecode(source_asset$href)), cold_output),
    env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep)),
    stdout = TRUE, stderr = TRUE)
  expect_identical(attr(cold, "status"), NULL, info = paste(cold, collapse = "\n"))
  expect_true(file.exists(cold_output) && file.info(cold_output)$size > 0)

  moved <- paste0(output, "-moved")
  expect_true(file.rename(output, moved))
  expect_silent(lisaR:::lisa_extension_validate_manifest(moved))
})
