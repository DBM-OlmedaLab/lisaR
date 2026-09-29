kegg_extension_fixture <- function(root, ids = c("hsa01100", "hsa00010"),
                                   titles = c("Metabolic pathways", "Glycolysis / Gluconeogenesis")) {
  source <- file.path(root, "canonical")
  work <- file.path(root, "canonical-kegg-maps")
  painter <- file.path(work, "outputs", "gene_level", "single_de", "A",
    "collection_GOBP-C2", "kegg_painter")
  dir.create(source, recursive = TRUE)
  dir.create(painter, recursive = TRUE)
  writeLines("<!doctype html><title>Canonical</title>", file.path(source, "report_index.html"))
  stems <- sprintf("%02d_%s_%s", seq_along(ids), ids,
    gsub("[^A-Za-z0-9]+", "_", titles))
  png <- file.path(painter, paste0(stems, "_painted.png"))
  pdf <- file.path(painter, paste0(stems, "_painted.pdf"))
  source_data <- file.path(painter, paste0(stems, "_painted_source.tsv"))
  recipe <- file.path(painter, paste0(stems, "_painted_recipe.R"))
  invisible(vapply(c(png, pdf, source_data, recipe), function(path) {
    writeLines(basename(path), path, useBytes = TRUE)
    TRUE
  }, logical(1L)))
  index_path <- file.path(painter, "A_GOBP-C2_gene_level_kegg_pathway_painter_index.tsv")
  lisaR:::write_lisa_tsv(data.frame(
    rank = seq_along(ids), analysis_id = "A", universe = "GOBP-C2",
    kegg_id = ids, kegg_title = titles,
    output_png = file.path(normalizePath(work),
      vapply(png, lisaR:::lisa_extension_relative, character(1L), root = work)),
    output_pdf = file.path(normalizePath(work),
      vapply(pdf, lisaR:::lisa_extension_relative, character(1L), root = work)),
    stringsAsFactors = FALSE
  ), index_path)
  list(source = source, work = work, painter = painter, index = index_path,
    png = png, pdf = pdf, source_data = source_data, recipe = recipe)
}

test_that("KEGG sibling report uses genuine multi-map identities and portable intact links", {
  root <- tempfile("kegg-report-extension-")
  dir.create(root)
  fixture <- kegg_extension_fixture(root)
  report <- lisaR:::lisa_write_kegg_extension_report(
    fixture$work, fixture$work, fixture$source,
    selection_policy = "significant_only", max_pathways_per_collection = 20L
  )
  expect_true(file.exists(report))
  manifest <- lisaR:::read_lisa_tsv(file.path(fixture$work, "kegg_report_manifest.tsv"))
  expect_identical(manifest$kegg_id, c("hsa01100", "hsa00010"))
  expect_identical(manifest$kegg_title,
    c("Metabolic pathways", "Glycolysis / Gluconeogenesis"))
  html <- paste(readLines(report, warn = FALSE), collapse = "\n")
  expect_match(html, 'data-kegg-id="hsa01100"', fixed = TRUE)
  expect_match(html, "Metabolic pathways", fixed = TRUE)
  expect_match(html, "../canonical/report_index.html", fixed = TRUE)
  expect_false(grepl(normalizePath(root), html, fixed = TRUE))
  expect_false(grepl("file:", html, fixed = TRUE))
  audit <- lisaR:::read_lisa_tsv(file.path(fixture$work, "kegg_report_link_audit.tsv"))
  expect_true(all(as.logical(audit$target_exists)))
  expect_equal(sum(audit$target_kind == "extension_asset"), 10L)
})

test_that("KEGG report fails closed on missing cache output and fabricated map identity", {
  root <- tempfile("kegg-report-invalid-")
  dir.create(root)
  fixture <- kegg_extension_fixture(root)
  unlink(fixture$pdf[[2L]])
  expect_error(
    lisaR:::lisa_kegg_extension_inventory(fixture$work, fixture$work),
    "hsa00010.*Glycolysis.*missing local output"
  )
  fixture <- kegg_extension_fixture(tempfile("kegg-report-fake-"), ids = "KEGG_SET_X",
    titles = "Not a native map")
  expect_error(
    lisaR:::lisa_kegg_extension_inventory(fixture$work, fixture$work),
    "genuine map identifiers"
  )
})

test_that("single painter reports every absent selected cache/render output", {
  index <- data.frame(
    kegg_id = c("hsa01100", "hsa00010"),
    kegg_title = c("Metabolic pathways", "Glycolysis / Gluconeogenesis"),
    output_png = c("map-1.png", ""),
    output_pdf = c("map-1.pdf", ""),
    stringsAsFactors = FALSE
  )
  nodes <- data.frame(error = "LISA-KEGG-007 required cached resource missing: kgml/hsa00010")
  expect_error(
    lisaR:::lisa_assert_selected_kegg_painter_outputs(index, nodes),
    "hsa00010.*Glycolysis.*required cached resource missing"
  )
  index$output_png[[2L]] <- "map-2.png"
  index$output_pdf[[2L]] <- "map-2.pdf"
  expect_true(lisaR:::lisa_assert_selected_kegg_painter_outputs(index, nodes))
})

test_that("KEGG maps are rendered inside the run, never in a -kegg-maps sibling", {
  bridge <- paste(deparse(body(lisaR:::run_lisa_pipeline_from_config)), collapse = "\n")
  extension <- paste(deparse(body(lisaR:::lisa_render_kegg_maps_extension)), collapse = "\n")
  # The on-demand extension keeps its established painter limit.
  expect_identical(formals(lisaR:::lisa_render_kegg_maps_extension)$max_pathways_per_collection, 20L)
  expect_match(extension, '"--top-pathways", as.character(max_pathways_per_collection)', fixed = TRUE)
  expect_false(grepl('"--top-pathways", "1"', extension, fixed = TRUE))
  # run_lisa() no longer creates sibling folders nor returns their paths (D5/D6).
  expect_false(grepl("-kegg-maps", bridge, fixed = TRUE))
  expect_false(grepl("-full-report", bridge, fixed = TRUE))
  expect_false(grepl("lisa_render_kegg_maps_extension(", bridge, fixed = TRUE))
  expect_false(grepl("render_lisa_categories(final_output_dir", bridge, fixed = TRUE))
  expect_false(grepl("extension_output_dir =", bridge, fixed = TRUE))
  expect_match(bridge, "kegg_maps = isTRUE(validated_config$kegg_maps$value)", fixed = TRUE)
})

test_that("an unusable KEGG snapshot is reported, not fatal", {
  root <- tempfile("kegg-missing-snapshot-")
  dir.create(root)
  expect_match(lisaR:::lisa_kegg_snapshot_problem(root, "absent-snapshot", "hsa"),
    "LISA-KEGG-001", fixed = TRUE)
  expect_match(lisaR:::lisa_kegg_snapshot_problem(file.path(root, "nope"), "x", "hsa"),
    "does not exist", fixed = TRUE)
})

test_that("a genuine KEGG painter failure is fatal and retains diagnostics", {
  root <- tempfile("kegg-failure-")
  gene <- file.path(root, "outputs/gene_level/single_de/A/collection_PATHWAYS")
  dir.create(file.path(gene, "kegg_painter"), recursive = TRUE)
  writeLines("fixture", file.path(gene, "A_PATHWAYS_gene_level_gene_category_contributions.tsv"))
  log <- file.path(gene, "kegg_painter", "error.log")
  writeLines("exit_code=17", log)
  f <- lisaR:::lisa_render_run_kegg_maps
  env <- new.env(parent = asNamespace("lisaR"))
  environment(f) <- env
  env$lisa_kegg_snapshot_problem <- function(...) ""
  env$lisa_write_code_identity_ledger <- function(...) data.frame()
  env$lisa_assert_run_tree_safe <- function(...) root
  env$lisa_run_post_script <- function(...) list(status = "failed", output_path = log,
    message = "exit_code=17; simulated rendering error")
  states <- character()
  record <- function(stage, analysis_id, contrast_id, collection, status, output_path, message) {
    states <<- c(states, status)
  }
  old <- options(lisaR.run_root = root)
  on.exit(options(old), add = TRUE)
  expect_error(f(list(kegg_cache_root = "cache", kegg_snapshot_id = "v1"),
    root, file.path(root, "outputs"), data.frame(), data.frame(),
    data.frame(analysis_collection = "PATHWAYS", run_kegg_layers = TRUE),
    data.frame(analysis_id = "A", scientific_name = "Homo sapiens"),
    data.frame(analysis_id = "A", collection = "PATHWAYS"),
    data.frame(contrast_id = character(), collection = character()), find.package("lisaR"), record),
    "LISA-KEGG-032.*exit_code=17")
  expect_identical(states, "failed")
  expect_true(file.exists(log))
  status <- read.delim(file.path(root, "kegg_maps_status.tsv"))
  expect_identical(status$status, "failed")
  expect_match(status$message, "exit_code=17", fixed = TRUE)
})
