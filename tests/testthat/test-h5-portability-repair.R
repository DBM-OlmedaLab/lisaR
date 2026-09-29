h5_matrix_binding_fixture <- function() {
  root <- tempfile("lisa-h5-matrix-binding-")
  matrix_dir <- file.path(root, "effective_inputs")
  collection <- file.path(
    root, "outputs", "single_de", "A", "collection_GOBP-C2"
  )
  dir.create(matrix_dir, recursive = TRUE)
  dir.create(file.path(root, "config"), recursive = TRUE)
  dir.create(file.path(collection, "lisa_tables"), recursive = TRUE)
  dir.create(file.path(collection, "enrichment"), recursive = TRUE)
  dir.create(file.path(collection, "inputs"), recursive = TRUE)
  matrix <- file.path(matrix_dir, "A_expression_matrix_path.tsv")
  lisaR:::write_lisa_tsv(
    data.frame(symbol = c("G1", "G2"), A_1 = c(1, 2), A_2 = c(3, 4)),
    matrix
  )
  row <- data.frame(
    analysis_id = "A",
    expression_matrix_path =
      "/retired/build-host/run/effective_inputs/A_expression_matrix_path.tsv",
    expression_matrix_path_input_scale = "transformed",
    expression_matrix_path_input_scale_source = "fixture:transformed",
    stringsAsFactors = FALSE
  )
  lisaR:::write_lisa_tsv(row, file.path(root, "config", "de_index.tsv"))
  lisaR:::write_lisa_tsv(
    data.frame(
      analysis_id = "A", input_class = "expression_matrix_path",
      effective_path = row$expression_matrix_path,
      sha256 = lisaR:::lisa_sha256_file(matrix), stringsAsFactors = FALSE
    ),
    file.path(root, "effective_inputs.tsv")
  )
  lisaR:::write_lisa_tsv(
    data.frame(category_id = "CAT_A", macrogroup_id = "SUPER_A",
               n_genesets = 1L, stringsAsFactors = FALSE),
    file.path(collection, "lisa_tables", "A_GSEA_category_summary.tsv")
  )
  lisaR:::write_lisa_tsv(
    data.frame(pathway = "GO_ALPHA", padj = 0.01, NES = 1,
               leadingEdge = "G1", category_id = "CAT_A",
               macrogroup_id = "SUPER_A", stringsAsFactors = FALSE),
    file.path(collection, "enrichment", "A_GSEA_test_annotated.tsv")
  )
  lisaR:::write_lisa_tsv(
    data.frame(symbol = c("G1", "G2"), log2FoldChange = c(1, -1),
               padj = c(0.01, 0.02), stringsAsFactors = FALSE),
    file.path(collection, "inputs", "A_standardized_DE.tsv")
  )
  refresh_manifest <- function(include_matrix = TRUE) {
    relative <- c(
      "effective_inputs.tsv",
      if (include_matrix) "effective_inputs/A_expression_matrix_path.tsv"
    )
    paths <- file.path(root, relative)
    lisaR:::write_lisa_tsv(
      data.frame(
        path = relative, bytes = as.numeric(file.info(paths)$size),
        sha256 = vapply(paths, lisaR:::lisa_sha256_file, character(1L)),
        stringsAsFactors = FALSE
      ),
      file.path(root, "run_manifest.tsv")
    )
  }
  refresh_manifest()
  list(root = root, row = row, matrix = matrix,
       receipt = file.path(root, "effective_inputs.tsv"),
       refresh_manifest = refresh_manifest)
}

test_that("saved heatmaps bind one contained manifest-verified effective matrix", {
  fixture <- h5_matrix_binding_fixture()
  resolved <- lisaR:::lisa_extension_effective_matrix(fixture$root, fixture$row)
  expect_identical(
    resolved$path, "effective_inputs/A_expression_matrix_path.tsv"
  )
  expect_identical(resolved$source_column, "expression_matrix_path")

  catalog <- lisaR:::lisa_extension_discover_catalog(
    fixture$root, list(heatmap_variants = "raw")
  )
  heatmap <- catalog[catalog$product == "heatmap", , drop = FALSE]
  expect_gt(nrow(heatmap), 0L)
  expect_true(all(heatmap$matrix_path == resolved$path))
  expect_true(all(heatmap$matrix_source_column == resolved$source_column))
  expect_true(all(heatmap$variant == "raw"))
  expect_true(all(catalog$matrix_path[catalog$product != "heatmap"] == ""))
})

test_that("saved heatmap matrix resolution rejects noncanonical evidence", {
  missing <- h5_matrix_binding_fixture()
  unlink(missing$matrix)
  expect_error(
    lisaR:::lisa_extension_effective_matrix(missing$root, missing$row),
    "missing, outside, or changed"
  )

  tampered <- h5_matrix_binding_fixture()
  writeLines("tampered", tampered$matrix)
  expect_error(
    lisaR:::lisa_extension_effective_matrix(tampered$root, tampered$row),
    "missing, outside, or changed"
  )

  outside <- h5_matrix_binding_fixture()
  receipt <- lisaR:::read_lisa_tsv(outside$receipt)
  receipt$effective_path <- "/uncontained/matrix.tsv"
  lisaR:::write_lisa_tsv(receipt, outside$receipt)
  outside$refresh_manifest()
  expect_error(
    lisaR:::lisa_extension_effective_matrix(outside$root, outside$row),
    "no contained effective_inputs path"
  )

  unmanifested <- h5_matrix_binding_fixture()
  unmanifested$refresh_manifest(include_matrix = FALSE)
  expect_error(
    lisaR:::lisa_extension_effective_matrix(unmanifested$root,
                                             unmanifested$row),
    "not uniquely declared in the canonical manifest"
  )

  ambiguous <- h5_matrix_binding_fixture()
  receipt <- lisaR:::read_lisa_tsv(ambiguous$receipt)
  lisaR:::write_lisa_tsv(rbind(receipt, receipt), ambiguous$receipt)
  ambiguous$refresh_manifest()
  expect_error(
    lisaR:::lisa_extension_effective_matrix(ambiguous$root, ambiguous$row),
    "no unique effective-input receipt"
  )
})

test_that("staged explicit heatmap matrices preserve saved scale provenance", {
  expect_silent(lisaR:::lisa_validate_post_script_args(
    "build_single_de_leading_edge_gene_heatmaps.R",
    c("--expression-matrix", "/staged/matrix.tsv",
      "--expression-matrix-column", "expression_matrix_path")
  ))
  builder <- system.file(
    "scripts", "build_single_de_leading_edge_gene_heatmaps.R", package = "lisaR"
  )
  expect_true(file.exists(builder))
  env <- new.env(parent = globalenv())
  for (expr in as.list(parse(builder))) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        is.call(expr[[3L]]) &&
        identical(expr[[3L]][[1L]], as.name("function"))) {
      eval(expr, envir = env)
    }
  }
  root <- tempfile("lisa-h5-staged-scale-")
  dir.create(file.path(root, "config"), recursive = TRUE)
  staged <- file.path(root, "staged-matrix.tsv")
  lisaR:::write_lisa_tsv(
    data.frame(symbol = c("G1", "G2"), A_1 = c(1, 2), A_2 = c(3, 4)),
    staged
  )
  lisaR:::write_lisa_tsv(
    data.frame(
      analysis_id = "A",
      expression_matrix_path =
        "/retired/build/effective_inputs/A_expression_matrix_path.tsv",
      expression_matrix_path_input_scale = "transformed",
      expression_matrix_path_input_scale_source = "saved:explicit",
      stringsAsFactors = FALSE
    ),
    file.path(root, "config", "de_index.tsv")
  )
  source <- env$discover_matrix_source(
    list(de_index = file.path(root, "config", "de_index.tsv")),
    list(analysis_id = "A", expression_matrix = staged,
         expression_matrix_column = "expression_matrix_path")
  )
  expect_identical(source$path, staged)
  expect_identical(source$input_scale, "transformed")
  expect_identical(source$input_scale_source, "saved:explicit")
  expect_identical(source$expression_transform, "none")
})
