resolved_input_fixture <- function() {
  root <- tempfile("lisa-resolved-")
  dir.create(root)
  # A resolved config reports the canonical spelling of a managed directory,
  # because lisaR resolves the trusted existing prefix. `tempfile()` returns a
  # platform alias of that directory on macOS (`/var` -> `/private/var`) and
  # Windows (8.3 short names), so the fixture must declare its root in the same
  # canonical spelling the assertions compare against. The assertions
  # themselves stay exact.
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  dir.create(file.path(root, "inputs"))
  dir.create(file.path(root, "indexes"))
  writeLines("symbol\tlog2FoldChange\tpadj\nA\t1\t0.4", file.path(root, "inputs", "de.tsv"))
  threshold <- 0.12345678901234566
  cfg <- list(pipeline = list(schema_version = "1.0.0",
    profile = "transcriptomic/genomic", evidence_mode = "full_de",
    output_dir = "results", gsea_padj_cutoff = threshold,
    duplicate_policies = list(de_table_duplicate_policy = "error",
      matrix_duplicate_policy = "error", mapped_id_collision_policy = "error")),
    collections = list("GOBP-C2"),
    single_de = list(list(analysis_id = "001", de_path = "inputs/de.tsv",
      species = "Homo sapiens", label = "NA")))
  json <- file.path(root, "study.json")
  jsonlite::write_json(cfg, json, auto_unbox = TRUE, digits = 17)
  yaml <- file.path(root, "study.yaml")
  writeLines(c("pipeline:", "  schema_version: '1.0.0'",
    "  profile: transcriptomic/genomic", "  evidence_mode: full_de",
    "  output_dir: results", paste0("  gsea_padj_cutoff: ", sprintf("%.17g", threshold)),
    "  duplicate_policies:", "    de_table_duplicate_policy: error",
    "    matrix_duplicate_policy: error", "    mapped_id_collision_policy: error",
    "collections: [GOBP-C2]", "single_de:", "  - analysis_id: '001'",
    "    de_path: inputs/de.tsv", "    species: Homo sapiens", "    label: 'NA'"), yaml)
  tsv_cfg <- cfg
  tsv_cfg$pipeline$de_index_path <- "indexes/de.tsv"
  tsv_cfg$single_de <- NULL
  utils::write.table(as.data.frame(cfg$single_de[[1]]), file.path(root, "indexes", "de.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  tsv_json <- file.path(root, "indexed.json")
  jsonlite::write_json(tsv_cfg, tsv_json, auto_unbox = TRUE, digits = 17)
  list(root = root, cfg = cfg, json = json, yaml = yaml, tsv_json = tsv_json,
    threshold = threshold)
}

test_that("lists YAML JSON and TSV indexes share one lossless resolved config", {
  fixture <- resolved_input_fixture()
  reference <- lisaR:::lisa_resolve_config(fixture$cfg, config_dir = fixture$root)
  for (path in c(fixture$json, fixture$yaml, fixture$tsv_json)) {
    got <- lisaR:::lisa_resolve_config(path)
    expect_identical(got$cfg, reference$cfg)
    expect_identical(got$de_index, reference$de_index)
    expect_identical(got$contract$gsea_padj_cutoff, fixture$threshold)
    expect_identical(as.numeric(got$provenance$value[
      got$provenance$field == "pipeline.gsea_padj_cutoff"]), fixture$threshold)
    expect_identical(got$de_index$de_path, file.path(fixture$root, "inputs", "de.tsv"))
    expect_identical(got$de_index$analysis_id, "001")
    expect_identical(got$de_index$label, "NA")
  }
  explicit <- reference$provenance$source[reference$provenance$field == "pipeline.gsea_padj_cutoff"]
  default <- reference$provenance$source[reference$provenance$field == "pipeline.run_ora"]
  expect_identical(explicit, "explicit")
  expect_identical(default, "default")
  expect_false(dir.exists(file.path(fixture$root, "results")))
})

test_that("public validation exposes the resolved values and their provenance", {
  fixture <- resolved_input_fixture()
  validated <- validate_lisa_config(fixture$json, check_files = FALSE)
  expect_identical(validated$config, validated$resolved_config$cfg)
  expect_identical(validated$gsea_padj_cutoff, fixture$threshold)
  expect_identical(validated$configuration_provenance, validated$resolved_config$provenance)
  expect_true(validated$config$report$category_evidence)
  expect_false(validated$config$report$legacy_gene_products)
  expect_identical(validated$config$report$evidence_max_sets, 25L)
  expect_identical(validated$config$report$evidence_max_genes, 40L)
  expect_identical(validated$config$pipeline$cache_mode, "off")
})

test_that("index aliases and numeric cells retain exact semantics", {
  fixture <- resolved_input_fixture()
  cfg <- fixture$cfg
  cfg$single_de[[1]]$label <- fixture$threshold
  cfg$single_de <- as.data.frame(cfg$single_de[[1]])
  resolved <- lisaR:::lisa_resolve_config(cfg, config_dir = fixture$root)
  expect_identical(as.numeric(resolved$de_index$label), fixture$threshold)
  expect_identical(lisaR:::lisa_config_text(fixture$threshold), sprintf("%.17g", fixture$threshold))
  cfg$single_de <- NULL
  cfg$pipeline$de_index_path <- "indexes/de.tsv"
  writeLines("analysis_id\tde_path\tspecies\tmisspelled\n001\tinputs/de.tsv\tHomo sapiens\tx",
    file.path(fixture$root, "indexes", "de.tsv"))
  expect_error(lisaR:::lisa_resolve_config(cfg, config_dir = fixture$root), "UNKNOWN")
})

test_that("report limits and opt-in cache controls are validated once", {
  fixture <- resolved_input_fixture()
  cfg <- fixture$cfg
  cfg$report <- list(evidence_max_sets = 0)
  expect_error(lisaR:::lisa_resolve_config(cfg, config_dir = fixture$root), "INTEGER")
  cfg$report <- list(evidence_max_sets = 2.5)
  expect_error(lisaR:::lisa_resolve_config(cfg, config_dir = fixture$root), "INTEGER")
  cfg$report <- NULL
  cfg$pipeline$cache_mode <- "readwrite"
  expect_error(lisaR:::lisa_resolve_config(cfg, config_dir = fixture$root), "explicit.*cache_dir")
  cfg$pipeline$cache_dir <- "cache"
  resolved <- lisaR:::lisa_resolve_config(cfg, config_dir = fixture$root)
  expect_identical(resolved$pipeline$cache_dir, file.path(fixture$root, "cache"))
  expect_false(dir.exists(file.path(fixture$root, "cache")))
  cfg$pipeline$cache_max_bytes <- Inf
  expect_error(lisaR:::lisa_resolve_config(cfg, config_dir = fixture$root), "CACHE-003")
})

de_adapter_fixture <- function() data.frame(gene = c("A", "B", "C", "D"),
  log2FoldChange = c(1, -2, 0.1, NA), stat = c(3, -4, 0.2, NA),
  pvalue = c(0.001, 0.01, 0.8, NA), padj = c(0.01, 0.05, 0.95, NA))

test_that("DE preparation preserves non-significant and NA tested rows", {
  source <- de_adapter_fixture()
  got <- prepare_lisa_de_input(source, "DESeq2", columns = list(symbol = "gene"),
    positive_direction = "treatment versus control", matrix_scale = "transformed")
  expect_identical(nrow(got$data), 4L)
  expect_identical(got$data, got$full_data)
  expect_identical(got$data$padj, source$padj)
  expect_identical(got$data$rank_value, source$stat)
  expect_identical(got$receipt$significance_filter, "none")
  expect_false(got$receipt$producing_software_verified)
  expect_identical(got$row_audit$rank_eligible, c(TRUE, TRUE, TRUE, FALSE))
  expect_error(prepare_lisa_de_input(source, "DESeq2", list(symbol = "gene"),
    positive_direction = "treatment versus control", na_policy = "error"), "missing statistics")
})

test_that("limma orientation and edgeR signed ranking are explicit", {
  source <- data.frame(symbol = c("A", "B"), logFC = c(1, -1), t = c(-3, 2),
    P.Value = c(0.01, 0.2), adj.P.Val = c(0.02, 0.2))
  got <- prepare_lisa_de_input(source, "limma", positive_direction = "B relative to A",
    rank_direction = "opposite_to_effect")
  expect_identical(got$data$rank_value, c(3, -2))
  expect_error(prepare_lisa_de_input(source, "limma", positive_direction = "B relative to A"), "signs conflict")
  edge <- data.frame(symbol = c("A", "B"), logFC = c(1, -1), F = c(9, 4),
    PValue = c(0.01, 0.2), FDR = c(0.02, 0.2), signed_rank = c(3, -2))
  expect_error(prepare_lisa_de_input(edge, "edgeR", positive_direction = "B relative to A"), "explicit signed rank")
  expect_error(prepare_lisa_de_input(edge, "edgeR", list(rank = "F"),
    positive_direction = "B relative to A"), "unsigned")
  got <- prepare_lisa_de_input(edge, "edgeR", list(rank = "signed_rank"),
    positive_direction = "B relative to A")
  expect_identical(got$data$rank_value, edge$signed_rank)
})

test_that("ambiguous IDs directions and malformed statistics fail closed", {
  source <- de_adapter_fixture()
  expect_error(prepare_lisa_de_input(source, "DESeq2", list(symbol = "gene")), "positive_direction")
  expect_error(prepare_lisa_de_input(source, "DESeq2", list(symbol = ".rownames"),
    positive_direction = "B relative to A"), "ambiguous gene identifiers")
  source$gene[1] <- NA_character_
  expect_error(prepare_lisa_de_input(source, "DESeq2", list(symbol = "gene"),
    positive_direction = "B relative to A"), "identifiers")
  source <- de_adapter_fixture()
  source$stat <- c("3", "bad", "0.2", NA)
  expect_error(prepare_lisa_de_input(source, "DESeq2", list(symbol = "gene"),
    positive_direction = "B relative to A"), "malformed")
})

test_that("explicit duplicate selection is audited without losing full input rows", {
  source <- de_adapter_fixture()
  source$gene[2] <- "a"
  expect_error(prepare_lisa_de_input(source, "DESeq2", list(symbol = "gene"),
    positive_direction = "B relative to A"), "duplicate identifiers")
  got <- prepare_lisa_de_input(source, "DESeq2", list(symbol = "gene"),
    positive_direction = "B relative to A", duplicate_policy = list(type = "select",
      criterion = "min_padj", tie_breaker = "row_number"))
  expect_identical(nrow(got$full_data), 4L)
  expect_identical(nrow(got$data), 3L)
  expect_identical(got$full_data$padj, source$padj)
  expect_identical(got$row_audit$retained, c(TRUE, FALSE, TRUE, TRUE))
  expect_identical(got$receipt$full_rows_retained, 4L)
})

test_that("explicit ID mapping preserves rows and identifier text", {
  source <- de_adapter_fixture()
  source$gene[2] <- "a"
  path <- tempfile(fileext = ".tsv")
  writeLines("source_row\toutput_id\n1\t001\n2\t002", path)
  policy <- list(type = "mapping_file", mapping_file = path)
  got <- prepare_lisa_de_input(source, "DESeq2", list(symbol = "gene"),
    positive_direction = "B relative to A", duplicate_policy = policy)
  expect_identical(got$data$symbol, c("001", "002", "C", "D"))
  expect_identical(nrow(got$full_data), 4L)
  expect_identical(got$data$rank_value, source$stat)
  writeLines("source_row\toutput_id\n1.5\t001\n2\t002", path)
  expect_error(prepare_lisa_de_input(source, "DESeq2", list(symbol = "gene"),
    positive_direction = "B relative to A", duplicate_policy = policy), "in-range integers")
})

test_that("JSON serialization preserves extreme doubles and in-memory run thresholds", {
  fixture <- resolved_input_fixture()
  values <- c(fixture$threshold, 1.2345678901234566e-50,
    .Machine$double.xmin, .Machine$double.xmax, 1 - .Machine$double.eps)
  path <- file.path(fixture$root, "normalized.json")
  lisaR:::lisa_write_normalized_config(list(values = values), path)
  expect_identical(jsonlite::fromJSON(path)$values, values)
  resolved <- lisaR:::lisa_resolve_config(fixture$json)
  lisaR:::lisa_write_normalized_config(resolved$cfg, path)
  expect_identical(lisaR:::lisa_resolve_config(path)$cfg, resolved$cfg)
  withr::local_dir(fixture$root)
  testthat::local_mocked_bindings(
    validate_lisa_config = function(...) list(execution_ready = TRUE),
    run_lisa_pipeline_from_config = function(config_path) {
      lisaR:::read_lisa_pipeline_config(config_path)$pipeline$gsea_padj_cutoff
    }, .package = "lisaR")
  for (value in values[values <= 1]) {
    fixture$cfg$pipeline$gsea_padj_cutoff <- value
    expect_identical(run_lisa(fixture$cfg), value)
  }
})

test_that("canonical objects preserve analysis sequence and scientific vectors", {
  input <- list(z = list(list(id = "second"), list(id = "first")),
    a = list(z = 1, a = 2), vector = c(3, 1, 2))
  got <- lisaR:::lisa_config_canonical_object(input)
  expect_identical(names(got), c("a", "vector", "z"))
  expect_identical(names(got$a), c("a", "z"))
  expect_identical(got$z, input$z)
  expect_identical(got$vector, input$vector)
})

test_that("a future relative-parent cache resolves before guarded initialization", {
  fixture <- resolved_input_fixture()
  config_dir <- file.path(fixture$root, "config")
  dir.create(config_dir)
  cfg <- fixture$cfg
  cfg$pipeline$output_dir <- "../results"
  cfg$pipeline$cache_mode <- "readwrite"
  cfg$pipeline$cache_dir <- "../cache/gsea"
  cfg$single_de[[1]]$de_path <- "../inputs/de.tsv"
  path <- file.path(config_dir, "study.json")
  jsonlite::write_json(cfg, path, auto_unbox = TRUE, digits = 17)
  expected_cache <- file.path(fixture$root, "cache", "gsea")
  expect_false(dir.exists(dirname(expected_cache)))

  # Exercise the public file path with a cache that does not yet exist. Input
  # paths still follow their established resolver and remain valid/read-only.
  validation <- validate_lisa_config(path, check_files = FALSE)
  expect_identical(validation$config$pipeline$cache_dir, expected_cache)
  expect_identical(validation$de_index$de_path,
    normalizePath(file.path(fixture$root, "inputs", "de.tsv"), winslash = "/"))
  plan <- plan_lisa_outputs(path, category_counts = c(`GOBP-C2` = 2L))
  expect_identical(plan$summary$analyses, 1L)
  expect_false(dir.exists(dirname(expected_cache)))
  withr::local_dir(config_dir)
  list_validation <- validate_lisa_config(cfg, check_files = FALSE)
  expect_identical(list_validation$config$pipeline$cache_dir, expected_cache)

  options <- lisaR:::lisa_content_cache_options(
    validation$config$pipeline$cache_dir, validation$config$pipeline$cache_mode,
    validation$config$pipeline$cache_max_bytes)
  slot <- lisaR:::lisa_cache_slot(options, paste(rep("a", 64L), collapse = ""), create = TRUE)
  expect_identical(slot$root, file.path(expected_cache, "lisa-fgsea-v1"))
  expect_true(dir.exists(slot$root))
  expect_false(file.exists(slot$path))
  expect_false(dir.exists(file.path(fixture$root, "results")))

  # Do not turn this repair into an unsafe lexical collapse through a missing
  # directory. Existing managed-path and generic cache guards stay strict.
  cfg$pipeline$cache_dir <- "../never-created/../other-cache"
  expect_error(validate_lisa_config(cfg, check_files = FALSE),
    "filesystem component that does not exist")
  unsafe <- lisaR:::lisa_content_cache_options(file.path(config_dir, "..", "direct-cache"), "readwrite")
  expect_error(lisaR:::lisa_cache_slot(unsafe, paste(rep("b", 64L), collapse = ""), create = TRUE),
    "Managed destinations may not contain")
})
