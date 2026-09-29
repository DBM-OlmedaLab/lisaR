native_heatmap_fixture <- function() {
  root <- tempfile("native-heatmap-generator-")
  analysis <- "HEATMAP_FIXTURE"; collection <- "GOCC"
  gene_dir <- file.path(root, "outputs", "gene_level", "single_de", analysis, "collection_GOCC")
  input_dir <- file.path(root, "outputs", "single_de", analysis, "collection_GOCC", "inputs")
  dir.create(gene_dir, recursive = TRUE); dir.create(input_dir, recursive = TRUE)
  dir.create(file.path(root, "config")); dir.create(file.path(root, "input"))
  symbols <- c("GENE_Z", "GENE_A", "GENE_M")
  category <- "CAT_NATIVE"
  label <- 'Long native category label with "quoted support" and a deliberate R-like string: ); stop("not code") #'
  expression <- data.frame(symbol = symbols, Visit_Z = c(1.25, 3.5, 2.75),
    Visit_A = c(2.25, 1.5, 2.75), Visit_M = c(4.25, 2.5, 2.75), check.names = FALSE)
  expression_path <- file.path(root, "input", "transformed_fixture.tsv")
  write_tsv <- function(x, file) utils::write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE)
  write_tsv(expression, expression_path)
  write_tsv(data.frame(analysis_id = analysis, label = "Native fixture",
    comparison = "Long comparison \u2014 treatment versus baseline",
    positive_direction = 'Positive means "treatment"; no hidden context changes.',
    expression_matrix_path = expression_path, expression_matrix_path_input_scale = "transformed",
    expression_matrix_path_input_scale_source = "fixture:explicit-transformed"),
    file.path(root, "config", "de_index.tsv"))
  write_tsv(data.frame(symbol = symbols, log2FoldChange = c(1.5, -3, .5),
    padj = c(.000012, .025, NA_real_)), file.path(input_dir, paste0(analysis, "_standardized_DE.tsv")))
  evidence <- data.frame(symbol = symbols, analysis_id = analysis, collection = collection,
    category_id = category, category_display_name = label, macrogroup_id = "SUPER_NATIVE",
    macrogroup_name = "Native supercategory", in_gsea_leading_edge = TRUE,
    min_source_padj = c(.01, .02, .03), stringsAsFactors = FALSE)
  # Exercise native device history: CAT_NATIVE is drawn only after a distinct
  # category has already opened/closed every requested output device.
  preceding <- evidence
  preceding$category_id <- "CAT_PRECEDING"
  preceding$category_display_name <- "Preceding category for native device-history regression"
  evidence <- rbind(preceding, evidence)
  write_tsv(evidence, file.path(gene_dir, paste0(analysis, "_GOCC_gene_level_gene_category_contributions.tsv")))
  list(root = root, analysis = analysis, category = category, label = label,
    expression = expression, output = file.path(gene_dir, "leading_edge_gene_heatmaps"))
}

native_heatmap_process <- function(script, args) {
  log <- suppressWarnings(system2(lisaR:::lisa_rscript_executable(),
    c("--vanilla", shQuote(script), vapply(args, shQuote, character(1L))), stdout = TRUE, stderr = TRUE))
  status <- attr(log, "status")
  list(status = if (is.null(status)) 0L else as.integer(status), log = log)
}

test_that("native generator emits portable exact two-argument heatmap recipes in independent processes", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("digest")
  skip_if_not_installed("png")
  builder <- system.file("scripts", "build_single_de_leading_edge_gene_heatmaps.R", package = "lisaR")
  expect_true(file.exists(builder))
  for (scale in c("zscore", "raw")) {
    fixture <- native_heatmap_fixture()
    formats <- c("png", "pdf", if (isTRUE(capabilities("cairo"))) "svg")
    result <- native_heatmap_process(builder, c("--project-dir", fixture$root,
      "--analysis-id", fixture$analysis, "--universe", "GOCC", "--scale", scale,
      "--formats", paste(formats, collapse = ","), "--top-genes", "3"))
    expect_identical(result$status, 0L, info = paste(result$log, collapse = "\n"))
    if (result$status != 0L) next
    index_file <- list.files(fixture$output, pattern = "_leading_edge_gene_heatmaps_index[.]tsv$", full.names = TRUE)
    expect_length(index_file, 1L)
    index <- lisaR:::read_lisa_tsv(index_file[[1L]])
    expect_identical(as.character(index$category_id), c("CAT_PRECEDING", fixture$category))
    preceding <- index[index$category_id == "CAT_PRECEDING", , drop = FALSE]
    expect_equal(nrow(preceding), 1L)
    for (format in formats) {
      preceding_file <- as.character(preceding[[paste0("plot_", format)]][[1L]])
      expect_true(file.exists(preceding_file))
      expect_gt(file.info(preceding_file)$size, 100)
    }
    # Replay the second category, whose native PNG followed the first
    # category's PDF (and SVG when available), in an independent cold process.
    index <- index[index$category_id == fixture$category, , drop = FALSE]
    expect_equal(nrow(index), 1L)
    expect_identical(as.character(index$category_id), fixture$category)
    source <- as.character(index$matrix_tsv[[1L]])
    recipe <- as.character(index$recipe_r[[1L]])
    original <- as.character(index$plot_png[[1L]])
    pdf <- as.character(index$plot_pdf[[1L]])
    expect_true(all(file.exists(c(source, recipe, original, pdf))))
    matrix <- lisaR:::read_lisa_tsv(source)
    expect_identical(unique(as.character(matrix$sample)), c("Visit_Z", "Visit_A", "Visit_M"))
    expect_identical(unique(as.character(matrix$symbol)), c("GENE_A", "GENE_Z", "GENE_M"))
    expect_identical(unique(as.character(matrix$category_display_name)), fixture$label)
    expect_true(all(c("log2FC", "padj", "plot_value") %in% names(matrix)))
    # This regression exercises a genuine saved matrix without adding generic
    # renderer flags. Native contract travels in the emitted recipe itself.
    expect_false("figure_type" %in% names(matrix))
    expected_raw <- vapply(seq_len(nrow(matrix)), function(i) {
      row <- match(matrix$symbol[[i]], fixture$expression$symbol)
      fixture$expression[[as.character(matrix$sample[[i]])]][[row]]
    }, numeric(1L))
    expect_equal(as.numeric(matrix$value_raw), expected_raw)
    expect_equal(as.numeric(matrix$plot_value_raw), expected_raw)
    if (scale == "raw") expect_equal(as.numeric(matrix$plot_value), expected_raw)
    if (scale == "zscore") expect_equal(as.numeric(matrix$plot_value[matrix$symbol == "GENE_M"]), rep(0, 3))
    code <- paste(readLines(recipe, warn = FALSE), collapse = "\n")
    expect_match(code, "lisa-native-leading-edge-heatmap-v1", fixed = TRUE)
    expect_match(code, digest::digest(file = source, algo = "sha256"), fixed = TRUE)
    expect_false(grepl("\\b(?:sys[.])?source\\s*\\(|lisaR::", code, perl = TRUE))
    original_hashes <- vapply(c(source, recipe, original, pdf), function(file) digest::digest(file = file, algo = "sha256"), character(1L))
    portable <- tempfile("standalone heatmap recipe "); dir.create(portable)
    copied_source <- file.path(portable, "renamed matrix.tsv")
    copied_recipe <- file.path(portable, "offline recipe.R")
    file.copy(source, copied_source); file.copy(recipe, copied_recipe)
    regenerated <- file.path(portable, "regenerated.png")
    replay <- native_heatmap_process(copied_recipe, c(copied_source, regenerated))
    expect_identical(replay$status, 0L, info = paste(replay$log, collapse = "\n"))
    expect_true(file.exists(regenerated))
    if (file.exists(regenerated)) {
      before <- png::readPNG(original); after <- png::readPNG(regenerated)
      expect_identical(dim(after), dim(before))
      expect_identical(after, before, info = "Decoded pixels must agree across cold generation and replay; no tolerance/warmup.")
    }
    for (format in setdiff(formats, "png")) {
      native_file <- as.character(index[[paste0("plot_", format)]][[1L]])
      expect_true(file.exists(native_file))
      vector_output <- file.path(portable, paste0("regenerated.", format))
      vector_replay <- native_heatmap_process(copied_recipe, c(copied_source, vector_output))
      expect_identical(vector_replay$status, 0L, info = paste(vector_replay$log, collapse = "\n"))
      expect_true(file.exists(vector_output))
      if (file.exists(vector_output)) {
        expect_gt(file.info(vector_output)$size, 100)
        if (format == "pdf") expect_identical(rawToChar(readBin(vector_output, "raw", n = 5L)), "%PDF-")
        if (format == "svg") {
          svg <- paste(readLines(vector_output, warn = FALSE), collapse = "\n")
          expect_match(svg, "<svg", fixed = TRUE)
          expect_match(svg, "<path", fixed = TRUE)
        }
      }
    }
    expect_identical(vapply(c(source, recipe, original, pdf), function(file) digest::digest(file = file, algo = "sha256"), character(1L)), original_hashes)
    # A changed matrix is not silently accepted under a recipe bound to its
    # original source; the destination must remain absent after refusal.
    changed <- file.path(portable, "changed matrix.tsv")
    tampered <- matrix; tampered$plot_value[[1L]] <- as.numeric(tampered$plot_value[[1L]]) + .125
    utils::write.table(tampered, changed, sep = "\t", quote = FALSE, row.names = FALSE)
    refused <- file.path(portable, "must-not-exist.png")
    failure <- native_heatmap_process(copied_recipe, c(changed, refused))
    expect_true(failure$status != 0L)
    expect_match(paste(failure$log, collapse = "\n"), "LISA-HEATMAP-RECIPE-002", fixed = TRUE)
    expect_false(file.exists(refused))
  }
})

test_that("heatmap executable closure binds the native helper and not the unused generic recipe", {
  dependencies <- lisaR:::lisa_post_script_dependencies()[["build_single_de_leading_edge_gene_heatmaps.R"]]
  expect_setequal(dependencies, c("lisa_plot_metadata.R", "lisa_leading_edge_heatmap_native.R"))
  root <- lisaR:::lisa_resolve_package_dir()
  identity <- lisaR:::lisa_post_script_identity(root, "build_single_de_leading_edge_gene_heatmaps.R")
  expect_null(lisaR:::lisa_post_script_renderer_sha256(identity))
  helper <- identity$files[basename(identity$files$absolute_path) == "lisa_leading_edge_heatmap_native.R", , drop = FALSE]
  expect_equal(nrow(helper), 1L)
  expect_identical(as.character(helper$sha256[[1L]]), digest::digest(file = helper$absolute_path[[1L]], algo = "sha256"))
})

test_that("native heatmap geometry preserves readable dense Riaz dimensions", {
  helper <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "lisa_leading_edge_heatmap_native.R", package = "lisaR"), envir = helper)
  dense <- helper$lisa_heatmap_size(n_genes = 30L, n_samples = 52L, n_stat_cols = 2L)
  routine <- helper$lisa_heatmap_size(n_genes = 30L, n_samples = 26L, n_stat_cols = 2L)
  expect_identical(unname(dense), c(14.46, 9.0))
  expect_identical(unname(routine), c(9.52, 9.0))
  expect_gt(dense[["width"]], dense[["height"]])
  expect_gte(dense[["width"]] / 54, 0.26)
  expect_gte(dense[["height"]] / 30, 0.30)
  code <- paste(deparse(body(helper$lisa_heatmap_plot), width.cutoff = 500L), collapse = "\n")
  expect_match(code, "angle = 45", fixed = TRUE)
  expect_match(code, "size = 6.5", fixed = TRUE)
  expect_match(code, "lineheight = 1.15", fixed = TRUE)
})


test_that("native heatmap helper mutation and substitution invalidate the executable closure", {
  package <- tempfile("native-heatmap-identity-")
  scripts <- file.path(package, "scripts")
  dir.create(scripts, recursive = TRUE)
  writeLines(c("Package: lisaR", "Version: 0.0.0"), file.path(package, "DESCRIPTION"))
  package <- lisaR:::lisa_resolve_package_dir(.test_package_dir = package)
  entrypoint <- "build_single_de_leading_edge_gene_heatmaps.R"
  files <- c(entrypoint, lisaR:::lisa_post_script_dependencies()[[entrypoint]])
  expect_true(all(file.copy(system.file("scripts", files, package = "lisaR"), scripts)))
  identity <- lisaR:::lisa_post_script_identity(package, entrypoint)
  helper <- file.path(scripts, "lisa_leading_edge_heatmap_native.R")
  writeLines(c(readLines(helper, warn = FALSE), "# altered native helper"), helper)
  changed <- tryCatch(lisaR:::lisa_verify_post_script_identity(identity), error = function(e) e)
  expect_s3_class(changed, "lisa_code_identity_error")
  expect_identical(changed$code, "LISA-CODE-IDENTITY-002")
  rebound <- lisaR:::lisa_post_script_identity(package, entrypoint)
  old <- identity$files$sha256[basename(identity$files$absolute_path) == basename(helper)]
  new <- rebound$files$sha256[basename(rebound$files$absolute_path) == basename(helper)]
  expect_false(identical(old, new))
  external <- tempfile("native-helper-external-", fileext = ".R")
  expect_true(file.copy(helper, external))
  expect_true(file.remove(helper))
  if (isTRUE(suppressWarnings(file.symlink(external, helper)))) {
    substituted <- tryCatch(lisaR:::lisa_post_script_identity(package, entrypoint), error = function(e) e)
    expect_s3_class(substituted, "lisa_code_identity_error")
    expect_identical(substituted$code, "LISA-CODE-IDENTITY-001")
  }
})
