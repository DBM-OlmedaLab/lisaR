# Regression fixtures for the independently reproduced SCI-201/202/203 cases.
# Load only function definitions from the installed script: CLI bootstrap and
# main() must not execute inside the test process. This also works in R CMD check.
audit_heatmap_functions <- function() {
  script <- system.file("scripts", "build_single_de_leading_edge_gene_heatmaps.R", package = "lisaR")
  stopifnot(file.exists(script))
  env <- new.env(parent = globalenv())
  for (expr in as.list(parse(script))) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function"))) {
      eval(expr, envir = env)
    }
  }
  env
}

audit_heatmap_fixture <- function(env, filename, column, explicit = FALSE, alias = NULL) {
  root <- tempfile("audit-matrix-")
  dir.create(root)
  matrix_path <- file.path(root, filename)
  env$write_tsv(data.frame(symbol = "A", S1 = 1, S2 = 3, S3 = 7, excluded = 99), matrix_path)
  idx <- data.frame(analysis_id = "test", sample_include_regex = "^S")
  idx[[column]] <- filename
  if (!is.null(alias)) idx[[alias]] <- filename
  index_path <- file.path(root, "de_index.tsv")
  env$write_tsv(idx, index_path)
  paths <- list(de_index = index_path, original_de = file.path(root, "test.tsv"))
  cfg <- list(analysis_id = "test", expression_matrix = if (explicit) matrix_path else "",
              top_genes = 30L, scale = "zscore")
  src <- env$discover_matrix_source(paths, cfg)
  matrix <- env$read_expression_matrix(src$path, src)
  matrix <- env$filter_expression_matrix_samples(matrix, paths, cfg)
  evidence <- data.frame(symbol = "A", category_id = "CAT", in_gsea_leading_edge_bool = TRUE,
                         min_source_padj = .01)
  de <- data.frame(symbol = "A", log2FC = 2, padj = .01)
  list(source = src, matrix = matrix,
       values = env$build_category_matrix(evidence, de, matrix, "CAT", cfg))
}

test_that("SCI-201 declared VST scale survives discovery, selection and renaming", {
  env <- audit_heatmap_functions()
  generic <- audit_heatmap_fixture(env, "expression.tsv", "vst_matrix_path")
  named <- audit_heatmap_fixture(env, "expression.vst.tsv", "vst_matrix_path")
  explicit <- audit_heatmap_fixture(env, "expression.tsv", "vst_matrix_path", explicit = TRUE)
  for (result in list(generic, named, explicit)) {
    expect_identical(result$source$input_scale, "transformed")
    expect_identical(result$source$input_scale_source, "de_index:vst_matrix_path")
    expect_identical(result$source$expression_transform, "none")
    expect_identical(attr(result$matrix, "sample_cols"), c("S1", "S2", "S3"))
    expect_equal(result$values$value_raw, c(1, 3, 7))
    expect_equal(result$values$plot_value_raw, c(1, 3, 7))
    expect_equal(result$values$plot_value, as.numeric(scale(c(1, 3, 7))))
    expect_true(all(result$values$expression_transform == "none"))
  }
  expect_equal(generic$values, named$values)
  expect_equal(generic$values, explicit$values)
})

test_that("SCI-201 typed counts override names while generic paths retain legacy behavior", {
  env <- audit_heatmap_functions()
  counts <- audit_heatmap_fixture(env, "misleading.vst.tsv", "counts_matrix_path")
  tpm <- audit_heatmap_fixture(env, "expression.tsv", "tpm_matrix_path")
  generic <- audit_heatmap_fixture(env, "expression.tsv", "expression_matrix_path")
  transformed <- audit_heatmap_fixture(env, "expression.log2.tsv", "expression_matrix_path")
  for (result in list(counts, tpm, generic)) {
    expect_equal(result$values$plot_value_raw, c(1, 2, 3))
    expect_identical(result$source$expression_transform, "log2(x+1)")
  }
  expect_identical(generic$source$input_scale_source, "legacy_basename")
  expect_identical(counts$source$input_scale_source, "de_index:counts_matrix_path")
  expect_equal(transformed$values$plot_value_raw, c(1, 3, 7))
  expect_identical(transformed$source$input_scale_source, "legacy_basename")
  # An explicit unindexed override must still work when the index has no row
  # for the requested analysis; it must not index into a zero-length column.
  path <- attr(generic$matrix, "source_path")
  cfg <- list(analysis_id = "other", expression_matrix = path)
  src <- env$discover_matrix_source(list(de_index = file.path(dirname(path), "de_index.tsv")), cfg)
  expect_identical(src$input_scale_source, "legacy_basename")
})

test_that("SCI-201 same-file aliases cannot hide or contradict a declared scale", {
  env <- audit_heatmap_functions()
  for (explicit in c(FALSE, TRUE)) {
    result <- audit_heatmap_fixture(env, "expression.tsv", "expression_matrix_path",
                                    explicit = explicit, alias = "vst_matrix_path")
    expect_identical(result$source$input_scale_source, "de_index:vst_matrix_path")
    expect_equal(result$values$plot_value_raw, c(1, 3, 7))
    expect_equal(result$values$plot_value, as.numeric(scale(c(1, 3, 7))))
    expect_error(
      audit_heatmap_fixture(env, "expression.tsv", "counts_matrix_path",
                            explicit = explicit, alias = "vst_matrix_path"),
      "Conflicting expression scales"
    )
  }
  compatible <- audit_heatmap_fixture(env, "expression.tsv", "counts_matrix_path", alias = "sample_counts_path")
  expect_equal(compatible$values$plot_value_raw, c(1, 2, 3))
  expect_identical(compatible$source$input_scale_source, "de_index:counts_matrix_path;de_index:sample_counts_path")
})

test_that("SCI-201 pipeline materialization preserves per-file expression scales", {
  root <- tempfile("audit-effective-matrices-")
  dir.create(root)
  env <- audit_heatmap_functions()
  de_path <- file.path(root, "de.tsv")
  write_lisa_tsv(data.frame(symbol = c("A", "B"), log2FC = c(2, -1),
                            padj = c(.01, .2), rank_value = c(2, -1)), de_path)
  files <- file.path(root, c("expression.tsv", "expression.log2.tsv", "counts.tsv"))
  for (path in files) {
    write_lisa_tsv(data.frame(symbol = "A", S1 = 1, S2 = 3, S3 = 7), path)
  }
  make_row <- function(id, generic, vst = "", counts = "") {
    list(analysis_id = id, de_path = de_path, species = "Homo sapiens",
         symbol_col = "symbol", logfc_col = "log2FC", padj_col = "padj",
         rank_col = "rank_value", expression_matrix_path = generic,
         vst_matrix_path = vst, counts_matrix_path = counts,
         matrix_feature_col = "symbol", sample_include_regex = "^S")
  }
  rows <- list(
    make_row("aliases", files[[1L]], vst = files[[1L]]),
    make_row("legacy", files[[2L]]),
    make_row("independent", files[[3L]], vst = files[[1L]], counts = files[[3L]]),
    make_row("no_matrix", "")
  )
  captured <- NULL
  testthat::local_mocked_bindings(
    run_lisa_pipeline = function(de_index, ...) {
      idx <- read_lisa_tsv(de_index)
      values <- lapply(idx$analysis_id[1:3], function(id) {
        cfg <- list(analysis_id = id, expression_matrix = "", top_genes = 30L, scale = "zscore")
        paths <- list(de_index = de_index, original_de = de_path)
        src <- env$discover_matrix_source(paths, cfg)
        mat <- env$filter_expression_matrix_samples(env$read_expression_matrix(src$path, src), paths, cfg)
        result <- env$build_category_matrix(
          data.frame(symbol = "A", category_id = "CAT", in_gsea_leading_edge_bool = TRUE,
                     min_source_padj = .01),
          data.frame(symbol = "A", log2FC = 2, padj = .01), mat, "CAT", cfg)
        list(source = src, values = result)
      })
      explicit <- list(analysis_id = "independent", expression_matrix = idx$vst_matrix_path[[3L]])
      captured <<- list(index = idx, values = values,
                        explicit = env$discover_matrix_source(list(de_index = de_index), explicit))
      invisible(list())
    }, .package = "lisaR"
  )
  run_fixture <- function(analyses, label) {
    cfg <- list(
      pipeline = list(
        schema_version = "1.0.0", profile = "transcriptomic/genomic", evidence_mode = "full_de",
        output_dir = file.path(root, paste0("run-", label)), workers = 1L, dry_run = TRUE,
        dictionary_resource = "lisa_quickstart_dictionary@2.0.0",
        term2gene_resource = "lisa_quickstart_term2gene@1.0.0",
        category_map_resource = "lisa_quickstart_category_map@1.0.0",
        duplicate_policies = list(de_table_duplicate_policy = "error",
                                  matrix_duplicate_policy = "error", mapped_id_collision_policy = "error")
      ),
      report = list(mode = "standard", formats = list(png = TRUE, svg = FALSE, pdf = FALSE),
                    source_data = TRUE, recipes = FALSE),
      collections = list("GOBP-C2"), single_de = analyses
    )
    config <- file.path(root, paste0(label, ".json"))
    jsonlite::write_json(cfg, config, auto_unbox = TRUE)
    suppressWarnings(run_lisa_pipeline_from_config(config))
  }
  run_fixture(rows, "valid")
  expect_equal(nrow(captured$index), 4L)
  expect_false(identical(captured$index$expression_matrix_path[[1L]], captured$index$vst_matrix_path[[1L]]))
  expect_equal(captured$index$expression_matrix_path_input_scale[1:3],
               c("transformed", "transformed", "count_like"))
  expect_identical(captured$index$vst_matrix_path_input_scale[[3L]], "transformed")
  expect_true(is.na(captured$index$expression_matrix_path_input_scale[[4L]]) ||
                captured$index$expression_matrix_path_input_scale[[4L]] == "")
  expect_equal(captured$values[[1L]]$values$plot_value_raw, c(1, 3, 7))
  expect_equal(captured$values[[2L]]$values$plot_value_raw, c(1, 3, 7))
  expect_equal(captured$values[[3L]]$values$plot_value_raw, c(1, 2, 3))
  expect_identical(captured$values[[1L]]$source$input_scale_source, "de_index:vst_matrix_path")
  expect_identical(captured$values[[2L]]$source$input_scale_source, "legacy_basename")
  expect_identical(captured$explicit$input_scale, "transformed")
  expect_identical(captured$explicit$input_scale_source, "de_index:vst_matrix_path")
  expect_error(run_fixture(list(make_row("conflict", files[[1L]], vst = files[[1L]],
                                        counts = files[[1L]])), "conflict"),
               "Conflicting expression scales")
})

test_that("SCI-202 incomplete statistics remain mapped but are not evaluable", {
  map <- data.frame(category_id = c("MIXED", "INCOMPLETE", "NOMINAL"),
                    display_name = c("Mixed", "Incomplete", "Nominal"),
                    macrogroup_id = "M", macrogroup_name = "Macrogroup", macrogroup_order = 1,
                    category_order_within_macrogroup = 1:3, color = "#000000")
  input <- data.frame(category_id = c(rep("MIXED", 4), rep("INCOMPLETE", 2), rep("NOMINAL", 2)),
                      pathway = paste0("SET", 1:8), NES = c(2, 8, NA, Inf, 5, NA, 2, -1),
                      padj = c(.01, NA, .0001, .0002, Inf, .001, .01, .5))
  summary <- lisaR:::build_lisa_summaries(input, data.frame(), map, TRUE, .05, .05)$gsea
  mixed <- summary[summary$category_id == "MIXED", ]
  expect_identical(mixed$n_genesets_mapped, 4L)
  expect_identical(mixed$n_genesets_evaluable, 1L)
  expect_identical(mixed$n_genesets_significant, 1L)
  expect_equal(mixed$mean_NES, 2)
  expect_equal(mixed$mean_NES_contextual, 2)
  expect_equal(mixed$min_padj_mapped, .01)
  incomplete <- summary[summary$category_id == "INCOMPLETE", ]
  expect_identical(incomplete$n_genesets_mapped, 2L)
  expect_identical(incomplete$n_genesets_evaluable, 0L)
  expect_identical(incomplete$n_genesets_significant, 0L)
  expect_true(is.na(incomplete$mean_NES_contextual))
  expect_true(is.na(incomplete$min_padj_mapped))
  nominal <- summary[summary$category_id == "NOMINAL", ]
  expect_identical(nominal$n_genesets_evaluable, 2L)
  expect_equal(nominal$mean_NES_contextual, .5)
  expect_equal(nominal$mean_NES, 2)
  expect_equal(nominal$min_padj_mapped, .01)
})

test_that("SCI-203 zero-weight components, including missing values, do not affect a score", {
  evidence <- data.frame(priority_id = c("A", "B", "C"), gene_set_id = c("S1", "S2", "S3"),
                         category_id = "CAT", effect = c(2, 1, NA), padj = c(.9, .01, .5))
  weights <- data.frame(scenario_id = c("effect", "confidence", "breadth"),
                        effect_weight = c(1, 0, 0), confidence_weight = c(0, 1, 0),
                        breadth_weight = c(0, 0, 1))
  result <- lisaR:::lisa_gps_sensitivity_grid(evidence, weights)$scores
  effect <- result[result$scenario_id == "effect", ]
  expect_equal(effect$lisa_gps, c(1, 0, NA))
  expect_equal(effect$rank, c(1L, 2L, NA_integer_))
  confidence <- result[result$scenario_id == "confidence", ]
  expected_confidence <- (-log10(evidence$padj) - min(-log10(evidence$padj))) /
    diff(range(-log10(evidence$padj)))
  expect_equal(confidence$lisa_gps, expected_confidence)
  breadth <- result[result$scenario_id == "breadth", ]
  expect_equal(breadth$lisa_gps, rep(1, 3))
  expect_equal(breadth$rank, rep(1L, 3))
  evidence$padj <- c(NA, NA, NA)
  missing <- lisaR:::lisa_gps_sensitivity_grid(evidence, weights)$scores
  expect_equal(missing$lisa_gps[missing$scenario_id == "effect"], c(1, 0, NA))
  expect_true(all(is.na(missing$lisa_gps[missing$scenario_id == "confidence"])))
  expect_true(all(is.na(missing$rank[missing$scenario_id == "confidence"])))
  expect_equal(missing$lisa_gps[missing$scenario_id == "breadth"], rep(1, 3))
})

test_that("SCI-203 positive-weight scores retain zero and equal-weight contracts", {
  evidence <- data.frame(priority_id = c("A", "B", "C", "D"), gene_set_id = paste0("S", 1:4),
                         category_id = "CAT", effect = c(2, 1, 3, NA), padj = c(.9, .01, .5, .1))
  base <- lisaR:::lisa_gps_prioritization(evidence)
  result <- lisaR:::lisa_gps_sensitivity_grid(evidence)$scores
  equal <- result[result$scenario_id == "equal", ]
  expect_equal(equal$lisa_gps, base$lisa_gps, tolerance = 1e-14)
  expect_equal(equal$pareto_front, base$pareto_front)
  expect_equal(equal$lisa_gps[1:2], c(0, 0))
  expect_true(is.na(equal$lisa_gps[[4]]))
  expect_true(is.na(equal$rank[[4]]))
})
