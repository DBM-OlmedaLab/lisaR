test_that("Riaz concordance detects agreement and reversed effects", {
  example_dir <- system.file(
    "examples", "riaz-gse91061", "scripts", package = "lisaR"
  )
  expect_true(nzchar(example_dir))
  source(file.path(example_dir, "_common.R"), local = TRUE)

  reference <- data.frame(
    symbol = paste0("G", 1:5),
    log2FoldChange = c(-2, -1, 0.5, 1, 3),
    stat = c(-4, -2, 1, 2, 6)
  )
  same <- riaz_effect_concordance(reference, reference)
  expect_equal(same$lfc_spearman, 1)
  expect_equal(same$stat_spearman, 1)
  expect_equal(same$lfc_directional_concordance, 1)

  reversed <- reference
  reversed$log2FoldChange <- -reversed$log2FoldChange
  reversed$stat <- -reversed$stat
  opposite <- riaz_effect_concordance(reversed, reference)
  expect_equal(opposite$lfc_spearman, -1)
  expect_equal(opposite$stat_spearman, -1)
  expect_equal(opposite$lfc_directional_concordance, 0)
})

test_that("Riaz path guard keeps derivative outputs outside a verified run", {
  example_dir <- system.file(
    "examples", "riaz-gse91061", "scripts", package = "lisaR"
  )
  expect_true(nzchar(example_dir))
  source(file.path(example_dir, "_common.R"), local = TRUE)

  root <- tempfile("riaz-verified-run-")
  dir.create(root)
  expect_true(riaz_path_is_within(root, root))
  expect_true(riaz_path_is_within(
    riaz_new_output_path(file.path(root, "showcase")), root
  ))
  expect_false(riaz_path_is_within(paste0(root, "-sibling"), root))
  expect_error(
    riaz_new_output_path(file.path(root, "missing", "..", "showcase")),
    "existing parent directory",
    fixed = TRUE
  )
  expect_identical(
    riaz_safe_relative_path(file.path("analyses", "one", "plot.png")),
    "analyses/one/plot.png"
  )
  expect_error(
    riaz_safe_relative_path(file.path("..", "escape.tsv")),
    "unsafe or non-relative",
    fixed = TRUE
  )

  dangling <- file.path(root, "dangling-link")
  link_created <- file.symlink(file.path(root, "absent-target"), dangling)
  if (isTRUE(link_created)) {
    expect_true(riaz_path_entry_exists(dangling))
    expect_error(
      riaz_new_output_path(dangling),
      "including as a symbolic link",
      fixed = TRUE
    )
  }
})

test_that("installed Riaz SHA-256 helpers delegate to the package primitive", {
  scripts <- system.file(
    "examples", "riaz-gse91061", "scripts",
    c("_common.R", "10_build_curated_showcase.R"), package = "lisaR"
  )
  expect_true(all(file.exists(scripts)))
  code <- paste(unlist(lapply(scripts, readLines, warn = FALSE)), collapse = "\n")
  expect_match(code, "lisa_sha256_file", fixed = TRUE)
  expect_false(grepl("sha256sum|certutil", code))
})

test_that("Riaz exact counts are versioned references, not hard failures", {
  validator <- system.file(
    "examples", "riaz-gse91061", "scripts",
    "06_validate_exported_inputs.R",
    package = "lisaR"
  )
  text <- readLines(validator, warn = FALSE)
  expect_false(any(grepl("nrow\\(x\\) == 22070", text)))
  expect_false(any(grepl("raw_fdr == expected", text)))
  expect_true(any(grepl("universe_delta", text, fixed = TRUE)))
  expect_true(any(grepl("raw_model_fdr_delta", text, fixed = TRUE)))
  expect_true(any(grepl("Responders minus PD", text, fixed = TRUE)))
})

test_that("installed Riaz entry points accept stable 1.x versions and reject other series", {
  script_dir <- system.file(
    "examples", "riaz-gse91061", "scripts", package = "lisaR"
  )
  accepted <- c(
    "0.6.0" = FALSE, "0.98.999" = FALSE,
    "0.99.0" = FALSE, "0.99.1" = FALSE,
    "1.0.0" = TRUE, "1.0.1" = TRUE, "1.0.99" = TRUE,
    "1.1.0" = TRUE, "1.99.99" = TRUE,
    "2.0.0" = FALSE, "2.0.1" = FALSE, "3.0.0" = FALSE
  )
  for (script in c("09_run_lisa_example.R", "10_build_curated_showcase.R")) {
    expressions <- as.list(parse(file.path(script_dir, script)))
    guards <- Filter(function(node) {
      is.call(node) &&
        (identical(node[[1L]], as.name("riaz_assert")) ||
          identical(node[[1L]], as.name("if"))) &&
        "active_version" %in% all.vars(node[[2L]]) &&
        any(c("<", ">=") %in% all.names(node[[2L]]))
    }, expressions)
    expect_length(guards, 1L)

    # Evaluate only the entry-point guard, never the analysis or file writes.
    guard_env <- new.env(parent = baseenv())
    guard_env$riaz_assert <- function(condition, message) {
      if (!isTRUE(condition)) stop(message, call. = FALSE)
      invisible(NULL)
    }
    for (version in names(accepted)) {
      guard_env$active_version <- package_version(version)
      if (accepted[[version]]) {
        expect_no_error(eval(guards[[1L]], envir = guard_env))
      } else {
        expect_error(eval(guards[[1L]], envir = guard_env), "lisaR 1.x", fixed = TRUE)
      }
    }
  }
})

test_that("Riaz canonical configuration is bounded, standard and reproducible", {
  example_root <- system.file(
    "examples", "riaz-gse91061", package = "lisaR"
  )
  yaml_path <- file.path(example_root, "config", "riaz-gse91061.yml")
  json_path <- file.path(example_root, "config", "riaz-gse91061.json")
  runner_path <- file.path(example_root, "scripts", "09_run_lisa_example.R")

  expect_true(file.exists(yaml_path))
  expect_true(file.exists(json_path))
  expect_true(file.exists(runner_path))

  yaml_cfg <- lisaR:::read_lisa_pipeline_config(yaml_path)
  json_cfg <- lisaR:::read_lisa_pipeline_config(json_path)
  expect_identical(as.integer(yaml_cfg$pipeline$workers), 4L)
  expect_identical(as.integer(json_cfg$pipeline$workers), 4L)
  expect_identical(
    jsonlite::toJSON(yaml_cfg, auto_unbox = TRUE, null = "null", digits = NA),
    jsonlite::toJSON(json_cfg, auto_unbox = TRUE, null = "null", digits = NA)
  )
  yaml_structure <- validate_lisa_config(
    yaml_path, check_files = FALSE, strict = TRUE
  )
  json_structure <- validate_lisa_config(
    json_path, check_files = FALSE, strict = TRUE
  )
  expect_true(yaml_structure$structural_valid)
  expect_true(json_structure$structural_valid)
  expect_identical(yaml_structure$schema_version, "1.0.0")
  expect_identical(yaml_structure$analyses, 5L)
  expect_identical(yaml_structure$contrasts, 2L)
  expect_identical(yaml_cfg$pipeline$schema_version, "1.0.0")
  expect_identical(yaml_cfg$pipeline$output_dir, "../results/lisa_standard")
  expect_identical(yaml_cfg$pipeline$gsea_padj_cutoff, 0.25)
  expect_false(yaml_cfg$pipeline$run_hallmarks)
  expect_false(yaml_cfg$pipeline$run_kegg_maps)
  expect_identical(
    unlist(yaml_cfg$collections, use.names = FALSE),
    c("GOBP-C2", "GOMF", "GOCC", "PATHWAYS")
  )
  expect_length(yaml_cfg$single_de, 5L)
  expect_length(yaml_cfg$contrasts, 2L)
  expect_identical(yaml_cfg$report$mode, "standard")
  expect_identical(
    yaml_cfg$report$formats,
    list(png = TRUE, svg = FALSE, pdf = FALSE)
  )
  expect_true(isTRUE(yaml_cfg$report$source_data))
  expect_true(isTRUE(yaml_cfg$report$recipes))

  runner <- paste(readLines(runner_path, warn = FALSE), collapse = "\n")
  expect_match(runner, "4-worker parallel path", fixed = TRUE)
  expect_match(runner, "validate_lisa_config", fixed = TRUE)
  expect_match(runner, "strict = TRUE", fixed = TRUE)
  expect_match(runner, "plan_lisa_outputs", fixed = TRUE)
  expect_match(runner, "lisa_resource_resolution.tsv", fixed = TRUE)
  expect_match(runner, "lisa_standard_planning_category_counts.tsv", fixed = TRUE)
  expect_match(runner, "category_counts = category_counts", fixed = TRUE)
  expect_match(runner, "resolved installed dictionary changed", ignore.case = TRUE)
  expect_match(runner, "resolved_dictionary$path", fixed = TRUE)
  expect_match(runner, "does not match the dictionary resolved", fixed = TRUE)
  expect_match(runner, "max_planned_files <- 20000L", fixed = TRUE)
  expect_match(runner, "max_planned_bytes <- 4 * 1024^3", fixed = TRUE)
  expect_match(runner, "planned_report_artifact_files", fixed = TRUE)
  expect_match(runner, "planned_nonfigure_reserve <- 2000L", fixed = TRUE)
  expect_match(runner, "planned_files_with_reserve", fixed = TRUE)
  expect_match(runner, "exec.yml", fixed = TRUE)
  expect_match(runner, "execution_config_sha256", fixed = TRUE)
  expect_match(runner, "riaz_path_entry_exists", fixed = TRUE)
  expect_match(runner, "lisa_standard_package_identity.tsv", fixed = TRUE)
  expect_match(runner, "lisa_standard_resource_identity.tsv", fixed = TRUE)
  expect_match(runner, "execution_resource_identity", fixed = TRUE)
  expect_match(runner, "resource identity changed after planning", fixed = TRUE)
  expect_match(runner, "riaz_path_is_within(source_config_path", fixed = TRUE)
  expect_match(runner, "plan_configuration$pipeline$dry_run <- TRUE", fixed = TRUE)
  expect_match(runner, 'identical(plan_result$gate, "PLAN_PASS")', fixed = TRUE)
  expect_match(runner, 'identical(verification$gate, "PASS")', fixed = TRUE)
  expect_match(runner, 'identical(verification$artifact_type, "scientific_run")',
               fixed = TRUE)
  expect_match(runner, 'canonical_table_format = "tsv"', fixed = TRUE)
  expect_false(grepl("I_ACCEPT_RIAZ_FULL_67_GIB", runner, fixed = TRUE))
  expect_false(grepl("24480293355", runner, fixed = TRUE))
  expect_false(grepl("extension_inventory.tsv", runner, fixed = TRUE))
})

test_that("Riaz presentation fits its bounded plan and completed-output budgets", {
  runner_path <- system.file(
    "examples", "riaz-gse91061", "scripts", "09_run_lisa_example.R",
    package = "lisaR"
  )
  expressions <- as.list(parse(runner_path))
  budget_env <- new.env(parent = baseenv())
  budget_env$riaz_assert <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
    invisible(NULL)
  }
  # Evaluate only the installed limits and actual guards, never the analysis,
  # resource installation, directory scan or output writes.
  for (name in c("max_planned_files", "max_planned_bytes")) {
    assignments <- Filter(function(node) {
      is.call(node) && identical(node[[1L]], as.name("<-")) &&
        identical(node[[2L]], as.name(name))
    }, expressions)
    expect_length(assignments, 1L)
    eval(assignments[[1L]], envir = budget_env)
  }
  expect_identical(budget_env$max_planned_bytes, 4294967296)
  expect_identical(budget_env$max_planned_files, 20000L)

  guard_for <- function(variable) {
    guards <- Filter(function(node) {
      is.call(node) && identical(node[[1L]], as.name("riaz_assert")) &&
        variable %in% all.vars(node[[2L]])
    }, expressions)
    expect_length(guards, 1L)
    guards[[1L]]
  }
  # Independent regression: 2,088 static units plus 28 SVG overviews, PNG,
  # sources and recipes. The former 2 GiB limit rejected this standard plan.
  canonical_planned_bytes <- 2791800832
  for (variable in c("planned_bytes", "actual_bytes")) {
    guard <- guard_for(variable)
    assign(variable, canonical_planned_bytes, envir = budget_env)
    expect_no_error(eval(guard, envir = budget_env))
    assign(variable, 4294967296, envir = budget_env)
    expect_no_error(eval(guard, envir = budget_env))
    assign(variable, 4294967297, envir = budget_env)
    expect_error(eval(guard, envir = budget_env), "4 GiB", fixed = TRUE)
    assign(variable, Inf, envir = budget_env)
    expect_error(eval(guard, envir = budget_env), "4 GiB", fixed = TRUE)
  }
  for (variable in c("planned_files_with_reserve", "actual_files")) {
    guard <- guard_for(variable)
    assign(variable, 20000L, envir = budget_env)
    expect_no_error(eval(guard, envir = budget_env))
    assign(variable, 20001L, envir = budget_env)
    expect_error(eval(guard, envir = budget_env), "20,000-file", fixed = TRUE)
  }
})

test_that("Riaz resources use ordinary C1 resolution without a project cache", {
  example_root <- system.file(
    "examples", "riaz-gse91061", package = "lisaR"
  )
  resource_script <- paste(
    readLines(
      file.path(example_root, "scripts", "08_prepare_lisa_resources.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  expect_match(resource_script, "prepare_lisa_msigdb_resource", fixed = TRUE)
  expect_match(resource_script, "LISAR_RIAZ_MSIGDB_ZIP", fixed = TRUE)
  expect_match(resource_script, "LISAR_RIAZ_ACCEPT_MSIGDB_TERMS", fixed = TRUE)
  expect_match(resource_script, "riaz_configure_normal_resources", fixed = TRUE)
  expect_match(resource_script, "lisa_resource_resolution.tsv", fixed = TRUE)
  expect_match(
    resource_script,
    "never accepts terms silently",
    fixed = TRUE
  )
  expect_false(grepl("install_lisa_resource", resource_script, fixed = TRUE))
  expect_false(grepl("resource_registry.tsv", resource_script, fixed = TRUE))
  expect_false(grepl("LISAR_RIAZ_DICTIONARY_SHA256", resource_script, fixed = TRUE))
})

test_that("Riaz prepared installer is manifest-verified and does not run enrichment", {
  script_path <- system.file(
    "examples", "riaz-gse91061", "scripts",
    "11_prepare_prepared_project.R", package = "lisaR"
  )
  expect_true(file.exists(script_path))
  prepared <- paste(readLines(script_path, warn = FALSE), collapse = "\n")
  expect_match(prepared, "prepared_bundle_manifest.tsv", fixed = TRUE)
  expect_match(prepared, "prepared_bundle_provenance.tsv", fixed = TRUE)
  expect_match(prepared, "LISAR_RIAZ_PREPARED_BUNDLE", fixed = TRUE)
  expect_match(prepared, "RIAZ_PREPARED_PROJECT=PASS", fixed = TRUE)
  expect_match(prepared, "09_run_lisa_example.R", fixed = TRUE)
  expect_match(prepared, "file.rename(stage, project_target)", fixed = TRUE)
  expect_match(prepared, "never extracts an", fixed = TRUE)
  expect_false(grepl("run_lisa\\(", prepared))
  expect_false(grepl("DESeq2::", prepared, fixed = TRUE))
})

test_that("Riaz environment variables are namespaced", {
  scripts_dir <- system.file(
    "examples", "riaz-gse91061", "scripts", package = "lisaR"
  )
  scripts <- list.files(scripts_dir, pattern = "[.]R$", full.names = TRUE)
  code <- paste(unlist(lapply(scripts, readLines, warn = FALSE)), collapse = "\n")
  expect_false(grepl('Sys.getenv\\("RIAZ_', code))
  expect_match(code, "LISAR_RIAZ_PROJECT_DIR", fixed = TRUE)
  expect_match(code, "LISAR_RIAZ_INPUT_CACHE", fixed = TRUE)
  expect_match(code, "LISAR_RIAZ_BENCHMARK_DIR", fixed = TRUE)
})

test_that("all installed Riaz stage scripts parse", {
  scripts_dir <- system.file(
    "examples", "riaz-gse91061", "scripts", package = "lisaR"
  )
  scripts <- sort(list.files(
    scripts_dir, pattern = "[.]R$", full.names = TRUE
  ))
  expect_gte(length(scripts), 14L)
  for (script in scripts) {
    expect_no_error(parse(script))
  }
})

test_that("Riaz showcase discovers standard-run outputs through indexes", {
  showcase_path <- system.file(
    "examples", "riaz-gse91061", "scripts",
    "10_build_curated_showcase.R", package = "lisaR"
  )
  showcase <- paste(readLines(showcase_path, warn = FALSE), collapse = "\n")
  expect_match(showcase, "verify_lisa_run", fixed = TRUE)
  expect_match(showcase, "scientific_run", fixed = TRUE)
  expect_match(showcase, '"config", "de_index.tsv"', fixed = TRUE)
  expect_match(showcase, '"config", "contrast_index.tsv"', fixed = TRUE)
  expect_match(showcase, "contrast_status.tsv", fixed = TRUE)
  expect_match(showcase, "n_genesets_mapped", fixed = TRUE)
  expect_match(showcase, "n_genesets_evaluable", fixed = TRUE)
  expect_match(showcase, "n_genesets_significant", fixed = TRUE)
  expect_match(showcase, "RIAZ_SHOWCASE_GATE=PASS", fixed = TRUE)
  expect_match(showcase, "riaz_path_is_within(output_root, source_root)", fixed = TRUE)
  expect_match(showcase, "outside the verified source run", fixed = TRUE)
  expect_match(showcase, "selected showcase source escaped", ignore.case = TRUE)
  expect_match(showcase, "riaz_safe_relative_path", fixed = TRUE)
  expect_match(showcase, "appeared before promotion", fixed = TRUE)
  expect_false(grepl("response_dynamics_response_dynamics", showcase, fixed = TRUE))
  expect_false(grepl("_RIAZ_GSE91061_GSEA", showcase, fixed = TRUE))
  expect_false(grepl("gene_cards|volcano|heatmap", showcase))
  expect_no_error(parse(showcase_path))
})

test_that("Riaz configuration and documentation share one output contract", {
  example_root <- system.file(
    "examples", "riaz-gse91061", package = "lisaR"
  )
  readme_path <- file.path(example_root, "README.md")
  vignette_path <- vignette_source_path("riaz-gse91061-worked-example")

  expect_true(file.exists(readme_path))
  expect_true(file.exists(vignette_path))
  readme <- paste(readLines(readme_path, warn = FALSE), collapse = "\n")
  vignette <- paste(readLines(vignette_path, warn = FALSE), collapse = "\n")
  documentation <- paste(readme, vignette, sep = "\n")

  expect_match(readme, "standard", ignore.case = TRUE)
  # Editorial wording is not an API. Check the online guide route and the
  # executable workflow contract, rather than superseded prose fragments.
  expect_match(readme, "https://olmedalab.org/lisaR/reader/riaz-gse91061-worked-example.html", fixed = TRUE)
  expect_match(vignette, "validate_lisa_config", fixed = TRUE)
  expect_match(vignette, "plan_lisa_outputs", fixed = TRUE)
  expect_match(vignette, "09_run_lisa_example.R", fixed = TRUE)
  expect_match(vignette, "build_selected_LISA_report.R", fixed = TRUE)
  expect_match(vignette, "scientific-report.html#assemble-an-integrated-selected-full-report", fixed=TRUE)
  expect_match(vignette, "scientific_run", fixed = TRUE)
  expect_false(grepl("build_lisa_report_package", documentation, fixed = TRUE))
  expect_false(grepl("A zero may mean", documentation, fixed = TRUE))
  expect_false(grepl("gzip-compressed", documentation, fixed = TRUE))
  expect_false(grepl("TSV and XLSX tables", documentation, fixed = TRUE))
})

test_that("parallel-execution vignette records the user-facing safety contract", {
  vignette_path <- system.file(
    "doc", "parallel-execution.Rmd", package = "lisaR"
  )
  source_path <- system.file(
    "", package = "lisaR"
  )
  candidates <- c(
    vignette_path,
    file.path(source_path, "vignettes", "parallel-execution.Rmd"),
    file.path(getwd(), "vignettes", "parallel-execution.Rmd"),
    file.path(getwd(), "..", "vignettes", "parallel-execution.Rmd"),
    file.path(getwd(), "..", "..", "vignettes", "parallel-execution.Rmd")
  )
  candidates <- candidates[file.exists(candidates)]
  skip_if(length(candidates) == 0L, "source vignette is unavailable in this test installation")

  text <- readLines(candidates[[1]], warn = FALSE)
  expect_true(any(grepl("pipeline:", text, fixed = TRUE)))
  expect_true(any(grepl("workers: 4", text, fixed = TRUE)))
  expect_true(any(grepl("parallel_task_seeds.tsv", text, fixed = TRUE)))
  expect_true(any(grepl("Do not rename the staging tree as a result", text, fixed = TRUE)))
  expect_false(any(grepl("For maintainers and auditors", text, fixed = TRUE)))
})


test_that("Riaz resource receipts preserve a single path with spaces", {
  source(system.file("examples/riaz-gse91061/scripts/_common.R",package="lisaR"),local=TRUE)
  path <- tempfile(fileext=".tsv")
  on.exit(unlink(path))
  writeLines(c("resource_cache", "D:/lisa test/cache"),path)
  receipt <- riaz_read_tsv(path)
  expect_named(receipt,"resource_cache")
  expect_identical(receipt$resource_cache,"D:/lisa test/cache")
})
