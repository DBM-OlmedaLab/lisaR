test_that("the release exposes only the investigator-facing API", {
  expect_setequal(
    getNamespaceExports("lisaR"),
    c("run_lisa", "validate_lisa_config", "run_lisa_de",
      "run_lisa_contrast", "verify_lisa_run", "plan_lisa_outputs",
      "plan_lisa_extension", "render_lisa_categories", "lisa_init_project",
      "install_lisa_resource", "prepare_lisa_msigdb_resource",
      "install_lisa_example_bundle",
      "validate_lisa_custom_resources", "prepare_lisa_de_input",
      "lisa_figure_request", "lisa_request_id", "lisa_explore_open",
      "lisa_explore_catalog", "lisa_explore_plan", "lisa_explore_status",
      "lisa_explore_submit", "lisa_explore_poll", "lisa_explore_cancel",
      "lisa_explore_artifacts", "lisa_explore_stale_artifacts",
      "lisa_explore_export", "lisa_explore_execute_job",
      "lisa_explore_pump", "lisa_explore_configure_kegg",
      "lisa_explore_prepare_kegg_index",
      "lisa_explore_prepare_contrast_products", "explore_lisa_run",
      "add_lisa_category_inference", "refresh_lisa_support_grades")
  )
})

test_that("the example acquisition entry never asks for a digest or a reviewed directory", {
  arguments <- names(formals(install_lisa_example_bundle))
  expect_identical(arguments, c("example", "destination", "source", "cache_root"))
  # No sha256/hash/checksum/bundle-directory argument may reappear: the archive
  # identity and the manifest contract are internal to lisaR.
  expect_false(any(grepl("sha|hash|checksum|digest|bundle|manifest|reviewed",
                         arguments, ignore.case = TRUE)))
})

test_that("public functions have explicit required arguments", {
  expect_true(all(c("config") %in% names(formals(run_lisa))))
  expect_true(all(c("config", "check_files", "strict") %in%
                    names(formals(validate_lisa_config))))
  expect_identical(
    names(formals(validate_lisa_custom_resources)),
    c(
      "dictionary_resource", "category_map_resource",
      "term2gene_resource", "species"
    )
  )
  expect_true(all(c(
    "source", "logical_id", "version", "species", "modality", "schema",
    "approved_origin", "expected_sha256", "cache_root", "registry_path"
  ) %in% names(formals(install_lisa_resource))))
  expect_identical(
    names(formals(prepare_lisa_msigdb_resource)),
    c("source_zip", "accept_terms", "cache_root", "registry_path")
  )
  expect_true(all(c("input", "output_dir", "species", "universes") %in%
                    names(formals(run_lisa_de))))
  expect_true(all(c("contrast_a", "contrast_b", "output_dir", "universes") %in%
                    names(formals(run_lisa_contrast))))
  expect_true("run_dir" %in% names(formals(verify_lisa_run)))
  expect_true(all(c("config", "category_counts") %in% names(formals(plan_lisa_outputs))))
  expect_true(all(c("source_run", "selection") %in% names(formals(plan_lisa_extension))))
  expect_true(all(c("source_run", "selection", "output_dir") %in% names(formals(render_lisa_categories))))
  expect_identical(
    names(formals(lisa_explore_configure_kegg)),
    c("ws", "cache_root", "snapshot_id", "species", "max_abs_log2fc",
      "color_power")
  )
  expect_identical(
    names(formals(lisa_explore_prepare_kegg_index)),
    c("ws", "analyses", "collections")
  )
  expect_identical(
    names(formals(lisa_explore_prepare_contrast_products)),
    c("ws", "contrasts", "collections")
  )
  expect_identical(names(formals(lisa_init_project)), "path")
})

test_that("validate_lisa_config rejects unknown keys", {
  cfg <- list(
    pipeline = list(
      schema_version = "1.0.0",
      profile = "transcriptomic/genomic",
      evidence_mode = "full_de",
      output_dir = tempfile("lisa-output-"),
      duplicate_policies = list(
        de_table_duplicate_policy = "error",
        matrix_duplicate_policy = "error",
        mapped_id_collision_policy = "error"
      ),
      unexpected_option = TRUE
    ),
    single_de = list(list(
      analysis_id = "treatment",
      de_path = "not-required-in-this-test.tsv",
      species = "Homo sapiens"
    ))
  )
  expect_error(
    validate_lisa_config(cfg, check_files = FALSE),
    "LISA-CONFIG-UNKNOWN-001.*field=pipeline.unexpected_option"
  )
})

test_that("run_lisa rejects ambiguous input", {
  expect_error(run_lisa(1), "YAML/JSON path or a named configuration list")
  expect_error(run_lisa(list()), "YAML/JSON path or a named configuration list")
})

test_that("report scripts are installed with the package", {
  scripts <- system.file("scripts", package = "lisaR")
  expect_true(nzchar(scripts))
  expect_true(file.exists(file.path(scripts, "build_LISA_report.R")))
  expect_true(dir.exists(lisaR:::lisa_resolve_package_dir()))
})

test_that("bundled YAML and JSON examples are equivalent and valid", {
  skip_if_not_installed("yaml")
  yaml_path <- system.file("examples", "minimal-study.yaml", package = "lisaR")
  json_path <- system.file("examples", "minimal-study.json", package = "lisaR")
  expect_true(file.exists(yaml_path))
  expect_true(file.exists(json_path))

  yaml_config <- yaml::read_yaml(yaml_path)
  json_config <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
  expect_identical(yaml_config$pipeline, json_config$pipeline)
  expect_identical(unlist(yaml_config$collections, use.names = FALSE),
                   unlist(json_config$collections, use.names = FALSE))
  expect_identical(yaml_config$single_de, json_config$single_de)
  expect_identical(yaml_config$contrasts, json_config$contrasts)

  yaml_validation <- validate_lisa_config(yaml_path)
  json_validation <- validate_lisa_config(json_path)
  expect_true(yaml_validation$valid)
  expect_true(json_validation$valid)
  expect_identical(yaml_validation$analyses, 2L)
  expect_identical(yaml_validation$contrasts, 1L)
  expect_identical(yaml_validation$de_index, json_validation$de_index)
  expect_identical(yaml_validation$contrast_index, json_validation$contrast_index)
})

test_that("bundled example creates a separate plan and then uses the free scientific destination", {
  path <- system.file("examples", "minimal-study.json", package = "lisaR")
  cfg <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  cfg$pipeline$output_dir <- tempfile("lisa-example-dry-run-")
  example_data <- system.file("extdata", "minimal", package = "lisaR")
  for (i in seq_along(cfg$single_de)) {
    cfg$single_de[[i]]$de_path <- file.path(
      example_data, basename(cfg$single_de[[i]]$de_path)
    )
  }
  result <- suppressWarnings(run_lisa(cfg))
  expect_identical(result$gate, "PLAN_PASS")
  expect_false(lisaR:::lisa_path_entry_exists(result$output_dir))
  expect_true(dir.exists(result$plan_dir))
  expect_identical(result$plan_path,
                   file.path(result$plan_dir, "lisa_pipeline_plan.tsv"))
  expect_true(file.exists(result$plan_path))
  expect_true(file.exists(file.path(result$plan_dir, "config",
                                    "pipeline_config.normalized.json")))
  status <- lisaR:::read_lisa_tsv(file.path(result$plan_dir, "lisa_plan_status.tsv"))
  expect_identical(status$gate, "PLAN_PASS")
  expect_false(status$scientific_complete)
  plan_check <- verify_lisa_run(result$plan_dir)
  expect_identical(plan_check$gate, "FAIL")
  expect_identical(plan_check$artifact_type, "plan")
  expect_true("planning_artifact_not_scientific_run" %in% plan_check$findings)
  expect_identical(lisaR:::lisa_run_state(result$plan_dir), "planned")

  engine_called <- FALSE
  testthat::local_mocked_bindings(
    run_lisa_pipeline = function(output_dir, dry_run, ...) {
      engine_called <<- TRUE
      expect_false(dry_run)
      lisaR:::write_lisa_tsv(
        data.frame(stage = "fixture", status = "completed"),
        file.path(output_dir, "fixture_scientific_status.tsv")
      )
      invisible(list(plan = data.frame(), metadata = data.frame()))
    },
    .package = "lisaR"
  )
  cfg$pipeline$dry_run <- FALSE
  completed <- suppressWarnings(run_lisa(cfg))
  expect_true(engine_called)
  expect_identical(completed$output_dir, result$output_dir)
  expect_identical(completed$gate, "PASS")
  expect_true(dir.exists(completed$output_dir))
  expect_identical(verify_lisa_run(completed$output_dir)$gate, "PASS")
  expect_identical(lisaR:::lisa_run_state(completed$output_dir), "completed")
  expect_true(dir.exists(result$plan_dir))
})

test_that("a config in a subdirectory plans into its declared parent", {
  source_config <- system.file(
    "examples", "minimal-study.json", package = "lisaR"
  )
  cfg <- jsonlite::fromJSON(source_config, simplifyVector = FALSE)
  project <- tempfile("lisa-config-parent-output-")
  config_dir <- file.path(project, "config")
  dir.create(config_dir, recursive = TRUE)
  on.exit(lisa_test_cleanup_path(project), add = TRUE)
  # `output_dir` is reported in the canonical spelling of the managed parent.
  # `tempfile()` hands back a platform alias of it on macOS and Windows, so the
  # fixture declares the canonical spelling and the assertion stays exact.
  project <- normalizePath(project, winslash = "/", mustWork = TRUE)
  config_dir <- file.path(project, "config")

  cfg$pipeline$output_dir <- file.path("..", "results", "planned-study")
  example_data <- system.file("extdata", "minimal", package = "lisaR")
  for (i in seq_along(cfg$single_de)) {
    cfg$single_de[[i]]$de_path <- file.path(
      example_data, basename(cfg$single_de[[i]]$de_path)
    )
  }
  config_path <- file.path(config_dir, "study.json")
  jsonlite::write_json(
    cfg, config_path, auto_unbox = TRUE, null = "null", digits = NA
  )

  result <- suppressWarnings(run_lisa(config_path))
  expected <- gsub(
    "\\\\", "/", file.path(project, "results", "planned-study")
  )
  expect_identical(result$gate, "PLAN_PASS")
  expect_identical(result$output_dir, expected)
  expect_false(lisaR:::lisa_path_entry_exists(expected))
  expect_true(dir.exists(result$plan_dir))
  check <- verify_lisa_run(result$plan_dir)
  expect_identical(check$gate, "FAIL")
  expect_identical(check$artifact_type, "plan")
  expect_true("planning_artifact_not_scientific_run" %in% check$findings)
})

test_that("an interrupted planning transaction is failed rather than complete", {
  final <- tempfile("lisa-plan-failure-")
  tx <- lisaR:::lisa_transaction_begin(
    final, lisa_test_run_contract(), run_id = "failedplan",
    artifact_kind = "plan"
  )
  lisaR:::lisa_transaction_event(tx, "planned")
  lisaR:::lisa_transaction_event(tx, "plan_computed")
  lisaR:::lisa_transaction_abort(tx, "simulated planning failure")
  expect_identical(lisaR:::lisa_run_state(tx$staging_dir), "failed")
  expect_false(lisaR:::lisa_path_entry_exists(final))
})
