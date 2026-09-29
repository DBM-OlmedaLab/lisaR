test_that("Input contract schema accepts each approved profile and evidence mode", {
  base_policy <- list(de_table_duplicate_policy = "error", matrix_duplicate_policy = "error", mapped_id_collision_policy = "error")
  for (profile in c("transcriptomic/genomic", "global proteomic", "targeted", "custom")) {
    cfg <- list(pipeline = list(schema_version = "1.0.0", profile = profile, evidence_mode = "full_de", duplicate_policies = base_policy))
    expect_equal(lisaR:::lisa_validate_pipeline_config(cfg)$profile, profile)
  }
  expect_equal(lisaR:::lisa_evidence_contract("rank_only")$products, "ranked_enrichment")
  expect_equal(lisaR:::lisa_evidence_contract("effect_only")$products, "effect_summary")
  expect_error(
    lisaR:::lisa_evidence_contract("custom"),
    "LISA-EVIDENCE-004.*field=pipeline[.]evidence[.]required_columns"
  )

  custom_cfg <- list(pipeline = list(
    schema_version = "1.0.0", profile = "custom", evidence_mode = "custom",
    evidence = list(
      required_columns = c("effect", "rank"),
      products = "ranked_enrichment"
    ),
    duplicate_policies = base_policy
  ))
  custom_contract <- lisaR:::lisa_validate_pipeline_config(custom_cfg)$evidence
  expect_identical(custom_contract$required, c("effect", "rank"))
  expect_identical(custom_contract$products, "ranked_enrichment")

  custom_preflight <- lisaR:::lisa_preflight_de_table(
    data.frame(symbol = c("A", "B"), log2FC = c(1, -1), score = c(3, -2)),
    "custom", "custom",
    columns = list(
      symbol = "symbol", logfc = "log2FC", rank = "score",
      custom = list(
        required_columns = c("effect", "rank"),
        products = "ranked_enrichment"
      )
    ),
    duplicate_policy = "error",
    warning_action = "ignore_with_reason",
    ignore_reason = "Small custom-contract fixture."
  )
  expect_identical(custom_preflight$evidence$required, c("effect", "rank"))

  expect_error(
    lisaR:::lisa_evidence_contract(
      "custom", list(required_columns = "custom_score")
    ),
    "LISA-EVIDENCE-005"
  )

  custom_cfg$pipeline$evidence$required_colums <- "typo"
  expect_error(
    lisaR:::lisa_validate_pipeline_config(custom_cfg),
    "LISA-CONFIG-UNKNOWN-001.*field=pipeline[.]evidence[.]required_colums"
  )
})

test_that("Input contract requires explicit duplicate policies and emits repair-oriented errors", {
  cfg <- list(pipeline = list(schema_version = "1.0.0", profile = "targeted", evidence_mode = "effect_only"))
  expect_error(lisaR:::lisa_validate_pipeline_config(cfg), "LISA-DUP-001")
  expect_error(lisaR:::lisa_validate_duplicate_policy("first", "de_table_duplicate_policy"), "LISA-DUP-002")
  df <- data.frame(symbol = c("A", "A"), log2FoldChange = c(1, 2), padj = c(.2, .1))
  expect_error(lisaR:::lisa_preflight_de_table(df, "targeted", "full_de", duplicate_policy = "error"), "LISA-DUP-004")
  resolved <- lisaR:::lisa_preflight_de_table(df, "targeted", "full_de", duplicate_policy = list(type = "select", criterion = "min_padj", tie_breaker = "row_number"), warning_action = "ignore_with_reason", ignore_reason = "Small targeted validation fixture.")
  expect_equal(nrow(resolved$data), 1)
  expect_equal(nrow(resolved$duplicate_resolution), 2)
  aggregate_policy <- list(type = "aggregate", numeric_method = "mean", non_numeric_method = "constant")
  expect_error(lisaR:::lisa_validate_pipeline_config(list(pipeline = list(schema_version = "1.0.0", profile = "targeted", evidence_mode = "full_de", duplicate_policies = list(de_table_duplicate_policy = aggregate_policy, matrix_duplicate_policy = "error", mapped_id_collision_policy = "error")))), "LISA-DUP-010")
  expect_error(lisaR:::lisa_validate_duplicate_policy(list(type = "aggregate", numeric_method = "mean", non_numeric_method = "first"), "matrix_duplicate_policy"), "LISA-DUP-007")
})

test_that("Input contract enforces matrix and mapped-ID policies with complete audit evidence", {
  matrix <- data.frame(feature_id = c("P1", "P1", "P2"), sample_a = c(2, 4, 8), sample_b = c(4, 6, 10))
  aggregate <- list(type = "aggregate", numeric_method = "mean", non_numeric_method = "constant")
  resolved_matrix <- lisaR:::lisa_preflight_matrix(matrix, "targeted", duplicate_policy = aggregate,
    warning_action = "ignore_with_reason", ignore_reason = "Small targeted validation fixture.")
  expect_equal(resolved_matrix$data$sample_a[resolved_matrix$data$feature_id == "P1"], 3)
  expect_equal(nrow(resolved_matrix$duplicate_resolution), 3)
  expect_true(all(nzchar(resolved_matrix$duplicate_resolution$output_values)))

  mapped <- data.frame(source_id = c("ENS1", "ENS2", "ENS3"), mapped_id = c("G1", "G1", "G2"), quality = c(1, 3, 2))
  selected <- lisaR:::lisa_resolve_mapped_id_collisions(mapped, duplicate_policy = list(type = "select", criterion = "max_quality", tie_breaker = "source_id"))
  expect_equal(selected$data$source_id[selected$data$mapped_id == "G1"], "ENS2")
  expect_equal(nrow(selected$evidence), 3)

  mapping_path <- tempfile(fileext = ".tsv")
  lisaR:::write_lisa_tsv(data.frame(source_row = c(1, 2), output_id = c("G1_A", "G1_B")), mapping_path)
  mapped_by_file <- lisaR:::lisa_resolve_duplicates(data.frame(symbol = c("G1", "G1"), value = c(1, 2)), "symbol", list(type = "mapping_file", mapping_file = mapping_path), input_class = "de_table")
  expect_equal(mapped_by_file$data$symbol, c("G1_A", "G1_B"))
  expect_equal(nrow(mapped_by_file$evidence), 2)
})

test_that("Input contract preserves an empty duplicate-audit schema", {
  empty <- lisaR:::lisa_duplicate_audit_empty()
  tagged <- lisaR:::lisa_tag_duplicate_audit(empty, "de_table")
  expect_equal(nrow(tagged), 0L)
  expect_identical(names(tagged), c("input_class", names(empty)))
})

test_that("configured DE duplicate selection is canonical, auditable, and rank-safe", {
  de <- data.frame(
    symbol = c("g1", "G1", "B"), log2FC = c(1, 1, 1),
    qvalue = c(0.01, 0.20, 0.10), score = c(2, 99, 50),
    stringsAsFactors = FALSE
  )
  columns <- list(symbol = "symbol", logfc = "log2FC", padj = "qvalue", rank = "score")
  select_min_padj <- list(type = "select", criterion = "min_padj", tie_breaker = "row_number")
  expect_error(lisaR:::lisa_preflight_de_table(de, "targeted", "full_de", columns = columns,
    duplicate_policy = "error"), "LISA-DUP-004")
  min_padj <- lisaR:::lisa_preflight_de_table(de, "targeted", "full_de", columns = columns,
    duplicate_policy = select_min_padj, warning_action = "ignore_with_reason",
    ignore_reason = "Small targeted validation fixture.")
  expect_identical(min_padj$data$symbol, c("g1", "B"))
  expect_identical(min_padj$data$score, c(2, 50))
  expect_identical(min_padj$duplicate_resolution$action, c("retained", "dropped", "retained"))
  expect_identical(min_padj$duplicate_resolution$output_id, c("g1", "g1", "B"))
  expect_true(all(grepl("g1;1;0.01;2", min_padj$duplicate_resolution$output_values[1:2], fixed = TRUE)))

  select_max_abs_rank <- list(type = "select", criterion = "max_abs_rank", tie_breaker = "row_number")
  max_abs_rank <- lisaR:::lisa_preflight_de_table(de, "targeted", "full_de", columns = columns,
    duplicate_policy = select_max_abs_rank, warning_action = "ignore_with_reason",
    ignore_reason = "Small targeted validation fixture.")
  expect_identical(max_abs_rank$data$symbol, c("G1", "B"))
  expect_identical(max_abs_rank$data$score, c(99, 50))

  ranked <- make_rank_vector(standardize_de_table(min_padj$data,
    symbol_col = "symbol", rank_col = "score", logfc_col = "log2FC",
    padj_col = "qvalue", pvalue_col = NULL, gene_id_col = NULL)$de)
  expect_identical(names(ranked), c("B", "G1"))
  expect_identical(unname(ranked), c(50, 2))
  expect_error(make_rank_vector(data.frame(symbol = c("g1", "G1"), rank_value = c(2, 99))), "LISA-DUP-018")
  expect_error(standardize_de_table(data.frame(
    symbol = c("g1", "G1"), log2FC = c(1, 1), qvalue = c(.1, .1), score = c(2, 2)
  ), symbol_col = "symbol", rank_col = "score", logfc_col = "log2FC",
  padj_col = "qvalue", pvalue_col = NULL, gene_id_col = NULL), "LISA-DUP-018")

  # lisaR's established uppercase canonical rule keeps the selected source ID
  # intact here, then uses STAT3 for downstream matching; no species remapping.
  mouse <- de[1:2, , drop = FALSE]
  mouse$symbol <- c("Stat3", "STAT3")
  mouse$qvalue <- c(0.01, 0.20)
  mouse_selected <- lisaR:::lisa_preflight_de_table(mouse, "targeted", "full_de", columns = columns,
    duplicate_policy = select_min_padj, warning_action = "ignore_with_reason",
    ignore_reason = "Small targeted validation fixture.")
  expect_identical(mouse_selected$data$symbol, "Stat3")
  mouse_rank <- make_rank_vector(standardize_de_table(mouse_selected$data,
    symbol_col = "symbol", rank_col = "score", logfc_col = "log2FC",
    padj_col = "qvalue", pvalue_col = NULL, gene_id_col = NULL)$de)
  expect_identical(names(mouse_rank), "STAT3")
})

test_that("configured duplicate policies materialize effective downstream inputs", {
  root <- tempfile("lisa-effective-inputs-")
  dir.create(root)
  de_path <- file.path(root, "de.tsv")
  matrix_path <- file.path(root, "matrix.tsv")
  mapped_path <- file.path(root, "mapped.tsv")
  de_mapping_path <- file.path(root, "de-mapping.tsv")
  write_lisa_tsv(data.frame(
    symbol = c("DUP", "DUP", "B"), log2FC = c(1, 1, 1),
    padj = c(.9, .01, .2), rank_value = c(100, 1, 50)
  ), de_path)
  write_lisa_tsv(data.frame(
    feature_id = c("P1", "P1", "P2"), sample_a = c(2, 8, 4)
  ), matrix_path)
  write_lisa_tsv(data.frame(
    source_id = c("S1", "S2", "S3"), mapped_id = c("M1", "M1", "M2"),
    quality = c(1, 2, 1)
  ), mapped_path)
  write_lisa_tsv(data.frame(
    source_row = c(1, 2), output_id = c("DUP_HIGH", "DUP_LOW")
  ), de_mapping_path)

  captured <- list()
  testthat::local_mocked_bindings(
    run_lisa_pipeline = function(de_index, ...) {
      effective_index <- read_lisa_tsv(de_index)
      captured[[length(captured) + 1L]] <<- list(
        ranks = make_rank_vector(read_lisa_tsv(effective_index$de_path[[1]])),
        matrix = read_lisa_tsv(effective_index$matrix_path[[1]]),
        mapped = read_lisa_tsv(effective_index$mapped_id_path[[1]]),
        receipt = read_lisa_tsv(file.path(dirname(dirname(de_index)), "effective_inputs.tsv"))
      )
      invisible(list())
    },
    .package = "lisaR"
  )
  run_policy <- function(de_policy, workers, label) {
    cfg <- list(
      pipeline = list(
        schema_version = "1.0.0", profile = "transcriptomic/genomic", evidence_mode = "full_de",
        output_dir = file.path(root, paste0("run-", label, "-", workers)), workers = workers,
        dry_run = TRUE, dictionary_resource = "lisa_quickstart_dictionary@2.0.0",
        term2gene_resource = "lisa_quickstart_term2gene@1.0.0",
        category_map_resource = "lisa_quickstart_category_map@1.0.0",
        duplicate_policies = list(
          de_table_duplicate_policy = de_policy,
          matrix_duplicate_policy = list(type = "aggregate", numeric_method = "mean", non_numeric_method = "constant"),
          mapped_id_collision_policy = list(type = "select", criterion = "max_quality", tie_breaker = "source_id")
        )
      ),
      report = list(mode = "standard", formats = list(png = TRUE, svg = FALSE, pdf = FALSE), source_data = TRUE, recipes = FALSE),
      collections = list("GOBP-C2"),
      single_de = list(list(
        analysis_id = "duplicate_fixture", de_path = de_path, species = "Homo sapiens",
        symbol_col = "symbol", logfc_col = "log2FC", padj_col = "padj", rank_col = "rank_value",
        matrix_path = matrix_path, matrix_id_col = "feature_id",
        mapped_id_path = mapped_path, mapped_id_col = "mapped_id", source_id_col = "source_id"
      ))
    )
    config_path <- file.path(root, paste0("config-", label, "-", workers, ".json"))
    jsonlite::write_json(cfg, config_path, auto_unbox = TRUE)
    suppressWarnings(run_lisa_pipeline_from_config(config_path))
  }

  select_policy <- list(type = "select", criterion = "min_padj", tie_breaker = "row_number")
  mapping_policy <- list(type = "mapping_file", mapping_file = de_mapping_path)
  run_policy(select_policy, 1L, "select")
  run_policy(select_policy, 2L, "select")
  run_policy(mapping_policy, 1L, "mapping")

  expect_identical(captured[[1]]$ranks, captured[[2]]$ranks)
  expect_identical(unname(captured[[1]]$ranks), c(50, 1))
  expect_identical(unname(captured[[3]]$ranks), c(100, 50, 1))
  expect_false(identical(captured[[1]]$ranks, captured[[3]]$ranks))
  expect_equal(captured[[1]]$matrix$sample_a[captured[[1]]$matrix$feature_id == "P1"], 5)
  expect_identical(captured[[1]]$mapped$source_id, c("S2", "S3"))
  expect_true(all(grepl("^[[:xdigit:]]{64}$", captured[[1]]$receipt$sha256)))
  expect_setequal(captured[[1]]$receipt$input_class, c("de_table", "matrix_path", "mapped_id"))
})

test_that("contract_manifest.tsv declares full_scope_analyses/full_scope_contrasts for full mode and leaves them empty for standard mode", {
  # Regression: report_full_scope_owner() used
  # to infer "selected" purely from artifact-directory existence, which
  # cannot distinguish a genuinely unselected owner from a selected owner
  # whose products never materialized. config_pipeline.R now records the
  # authoritative selection (every de_index/contrast_index owner, for an
  # ordinary uniform FULL run) directly in contract_manifest.tsv. The
  # contrast identity must be the composite id the report generator actually
  # matches against on disk (lisa_contrast_name(): "<contrast_id>_<output_id>",
  # with output_id overriding a bare contrast_id when declared).
  root <- tempfile("lisa-full-scope-manifest-")
  dir.create(root)
  de_path_a <- file.path(root, "de_a.tsv")
  de_path_b <- file.path(root, "de_b.tsv")
  lisaR:::write_lisa_tsv(data.frame(
    symbol = c("A", "B"), log2FC = c(1, -1), padj = c(.01, .02), rank_value = c(2, -2)
  ), de_path_a)
  lisaR:::write_lisa_tsv(data.frame(
    symbol = c("A", "B"), log2FC = c(0.5, -0.5), padj = c(.03, .04), rank_value = c(1, -1)
  ), de_path_b)

  testthat::local_mocked_bindings(
    run_lisa_pipeline = function(...) invisible(list()),
    .package = "lisaR"
  )

  run_mode <- function(mode, label) {
    output_dir <- file.path(root, paste0("run-", label))
    cfg <- list(
      pipeline = list(
        schema_version = "1.0.0", profile = "transcriptomic/genomic", evidence_mode = "full_de",
        output_dir = output_dir, workers = 1L,
        dry_run = TRUE, dictionary_resource = "lisa_quickstart_dictionary@2.0.0",
        term2gene_resource = "lisa_quickstart_term2gene@1.0.0",
        category_map_resource = "lisa_quickstart_category_map@1.0.0",
        duplicate_policies = list(
          de_table_duplicate_policy = "error",
          matrix_duplicate_policy = "error",
          mapped_id_collision_policy = "error"
        )
      ),
      report = list(mode = mode, formats = list(png = TRUE, svg = FALSE, pdf = FALSE), source_data = TRUE, recipes = FALSE),
      collections = list("GOBP-C2"),
      single_de = list(
        list(analysis_id = "analysis_one", de_path = de_path_a, species = "Homo sapiens",
          symbol_col = "symbol", logfc_col = "log2FC", padj_col = "padj", rank_col = "rank_value"),
        list(analysis_id = "analysis_two", de_path = de_path_b, species = "Homo sapiens",
          symbol_col = "symbol", logfc_col = "log2FC", padj_col = "padj", rank_col = "rank_value")
      ),
      contrasts = list(list(
        contrast_id = "resp_vs_pd", output_id = "responders_vs_progressors",
        analysis_a = "analysis_one", analysis_b = "analysis_two"
      ))
    )
    config_path <- file.path(root, paste0("config-", label, ".json"))
    jsonlite::write_json(cfg, config_path, auto_unbox = TRUE)
    suppressWarnings(lisaR:::run_lisa_pipeline_from_config(config_path))
  }

  full_result <- run_mode("full", "full")
  manifest_full <- lisaR:::read_lisa_tsv(file.path(full_result$plan_dir, "contract_manifest.tsv"))
  expect_identical(
    manifest_full$value[manifest_full$key == "full_scope_analyses"],
    "analysis_one;analysis_two"
  )
  expect_identical(
    manifest_full$value[manifest_full$key == "full_scope_contrasts"],
    "resp_vs_pd_responders_vs_progressors"
  )

  standard_result <- run_mode("standard", "standard")
  manifest_standard <- lisaR:::read_lisa_tsv(file.path(standard_result$plan_dir, "contract_manifest.tsv"))
  expect_identical(manifest_standard$value[manifest_standard$key == "full_scope_analyses"], "")
  expect_identical(manifest_standard$value[manifest_standard$key == "full_scope_contrasts"], "")
})

test_that("mouse-native MH rows satisfy the HALLMARKS contract", {
  term2gene <- data.frame(
    gs_collection = "MH", gs_subcollection = "", gs_name = "HALLMARK_TEST",
    gs_exact_source = "fixture", gene_symbol = c("Stat3", "Mdk"),
    stringsAsFactors = FALSE
  )
  result <- lisaR:::build_hallmarks_lisa_inputs(term2gene)
  expect_identical(result$lisa_dict$gene_set_id, "HALLMARK_TEST")
})

test_that("Input contract profile fixtures issue approved strong warnings", {
  fixture_dir <- test_path("..", "fixtures", "unit")
  expect_warning(lisaR:::lisa_preflight_de_table(file.path(fixture_dir, "transcriptomic_de.tsv"), "transcriptomic/genomic", "full_de", duplicate_policy = "error"), "LISA-SUFFICIENCY-002")
  expect_warning(lisaR:::lisa_preflight_de_table(file.path(fixture_dir, "proteomic_de.tsv"), "global proteomic", "full_de", duplicate_policy = "error"), "strong warning")
  expect_error(lisaR:::lisa_preflight_de_table(data.frame(symbol = "P", log2FoldChange = 1, padj = .1), "global proteomic", "full_de", duplicate_policy = "error", warning_action = "error"), "LISA-SUFFICIENCY-002")
})

test_that("Input contract does not synthesize missing full-DE evidence", {
  base <- system.file("extdata/minimal", package = "lisaR")
  de <- read_lisa_tsv(file.path(base, "de_A.tsv"))
  expect_error(build_lisa_gene_level_tables(de_table = de[, c("symbol", "padj")], gsea_table = file.path(base, "gsea_A.tsv"), output_dir = tempdir()), "LISA-EVIDENCE-002")
  expect_error(lisaR:::lisa_preflight_de_table(de[, c("symbol", "log2FC")], "targeted", "full_de", duplicate_policy = "error"), "LISA-EVIDENCE-002")
})

test_that("HALLMARKS and KEGG Maps are explicit in the plan", {
  base <- system.file("extdata/minimal", package = "lisaR")
  de <- validate_lisa_de_index(file.path(base, "de_index.tsv"), base_dir = base)
  plan_default <- lisa_pipeline_plan(de, collections = c("GOBP-C2", "HALLMARKS"))
  expect_true(any(grepl("HALLMARKS", plan_default$collections, fixed = TRUE)))
  plan_disabled <- lisa_pipeline_plan(de, collections = c("GOBP-C2", "HALLMARKS"), run_hallmarks = FALSE)
  expect_false(any(grepl("HALLMARKS", plan_disabled$collections, fixed = TRUE)))
  expect_warning(lisa_pipeline_plan(de, run_kegg = TRUE), "deprecated")
  applicability <- lisaR:::lisa_product_applicability(lisaR:::lisa_evidence_contract("rank_only"), run_kegg_maps = TRUE)
  expect_equal(applicability$status[applicability$product == "ranked_enrichment"], "enabled")
  expect_true(all(applicability$status[applicability$product != "ranked_enrichment"] == "not_applicable"))
})
