strict_config_fixture <- function() {
  list(
    pipeline = list(
      schema_version = "1.0.0",
      profile = "targeted",
      evidence_mode = "full_de",
      output_dir = tempfile("lisa-strict-config-"),
      workers = 1L,
      dry_run = TRUE,
      run_ora = FALSE,
      run_hallmarks = TRUE,
      run_kegg_maps = FALSE,
      duplicate_policies = list(
        de_table_duplicate_policy = "error",
        matrix_duplicate_policy = "error",
        mapped_id_collision_policy = "error"
      )
    ),
    report = list(
      mode = "standard",
      formats = list(png = TRUE, svg = FALSE, pdf = FALSE),
      source_data = TRUE,
      recipes = FALSE
    ),
    collections = list("GOBP-C2"),
    single_de = list(list(
      analysis_id = "analysis_a",
      de_path = "analysis-a.tsv",
      species = "Homo sapiens",
      symbol_col = "symbol",
      logfc_col = "log2FC",
      padj_col = "padj",
      expression_matrix_path = "matrix.tsv",
      matrix_feature_col = "symbol",
      sample_include_regex = "^A_",
      label = "Analysis A",
      comparison = "A minus reference",
      positive_direction = "Positive values are higher in A.",
      model_note = "Paired fixture."
    )),
    contrasts = list(list(
      contrast_id = "a_vs_b",
      output_id = "profile_a_vs_b",
      analysis_a = "analysis_a",
      analysis_b = "analysis_a",
      contrast_label = "A versus A",
      comparison = "Profile A minus profile A",
      positive_direction = "Positive values favour the first profile."
    ))
  )
}

test_that("pipeline.package_dir can never select executable code", {
  for (strict in c(TRUE, FALSE)) {
    cfg <- strict_config_fixture()
    cfg$pipeline$package_dir <- tempfile("untrusted-lisa-tree-")
    error <- tryCatch(
      lisaR:::lisa_validate_pipeline_config(cfg, strict = strict),
      error = identity
    )
    expect_s3_class(error, "lisa_code_selection_error")
    expect_identical(error$code, "LISA-CONFIG-CODE-001")
    expect_identical(error$field, "pipeline.package_dir")
    expect_match(conditionMessage(error), "remove pipeline.package_dir", fixed = TRUE)
    expect_match(conditionMessage(error), "install the intended lisaR version", fixed = TRUE)
  }
})

test_that("configuration booleans are exact scalar JSON/YAML booleans", {
  invalid <- list("treu", 0L, 1L, 2L, c(TRUE, FALSE), NA)
  for (value in invalid) {
    cfg <- strict_config_fixture()
    cfg$pipeline$run_ora <- value
    expect_error(
      lisaR:::lisa_validate_pipeline_config(cfg),
      "LISA-CONFIG-BOOL-001.*field=pipeline[.]run_ora"
    )
  }

  cfg <- strict_config_fixture()
  cfg$report$formats$png <- c(TRUE, FALSE)
  expect_error(
    lisaR:::lisa_validate_pipeline_config(cfg),
    "LISA-CONFIG-BOOL-001.*field=report[.]formats[.]png"
  )
  cfg <- strict_config_fixture()
  cfg$report$source_data <- 1L
  expect_error(
    lisaR:::lisa_validate_pipeline_config(cfg),
    "LISA-CONFIG-BOOL-001.*field=report[.]source_data"
  )
})

test_that("path-bearing pipeline identifiers fail portable checks during prevalidation", {
  invalid <- list(
    run_id = c("../escape", "CON", "caf\u00e9", "trailing."),
    file_label_prefix = c("A B", "LPT1.txt", "x|y", "trailing ")
  )
  for (field in names(invalid)) {
    for (value in invalid[[field]]) {
      cfg <- strict_config_fixture()
      cfg$pipeline[[field]] <- value
      expect_error(
        lisaR:::lisa_validate_pipeline_config(cfg),
        paste0("pipeline[.]", field),
        info = paste(field, value)
      )
    }
  }

  cfg <- strict_config_fixture()
  cfg$pipeline$run_id <- "run-20260830_A"
  cfg$pipeline$file_label_prefix <- "RIAZ_GSE91061"
  expect_silent(lisaR:::lisa_validate_pipeline_config(cfg))

  cfg$pipeline$run_id <- "../prevalidation-escape"
  expect_error(
    validate_lisa_config(cfg, check_files = FALSE),
    "pipeline[.]run_id"
  )
})

test_that("unknown nested configuration keys fail with their complete field path", {
  mutations <- list(
    list(path = "pipeline.positive_direciton", mutate = function(x) {
      x$pipeline$positive_direciton <- "typo"; x
    }),
    list(path = "single_de[1].positive_direciton", mutate = function(x) {
      x$single_de[[1]]$positive_direciton <- "typo"; x
    }),
    list(path = "single_de[1].logfc_clo", mutate = function(x) {
      x$single_de[[1]]$logfc_clo <- "log2FC"; x
    }),
    list(path = "contrasts[1].analysis_c", mutate = function(x) {
      x$contrasts[[1]]$analysis_c <- "analysis_a"; x
    }),
    list(
      path = "pipeline.duplicate_policies.de_table_duplicate_policy.criteron",
      mutate = function(x) {
        x$pipeline$duplicate_policies$de_table_duplicate_policy <- list(
          type = "select", criterion = "min_padj", tie_breaker = "row_number",
          criteron = "min_padj"
        ); x
      }
    ),
    list(
      path = "pipeline.duplicate_policies.de_table_duplicate_policy.tie_breaker.directon",
      mutate = function(x) {
        x$pipeline$duplicate_policies$de_table_duplicate_policy <- list(
          type = "select", criterion = "min_padj",
          tie_breaker = list(column = "row_number", directon = "ascending")
        ); x
      }
    ),
    list(path = "report.source_date", mutate = function(x) {
      x$report$source_date <- TRUE; x
    }),
    list(path = "source_data[1].source_rol", mutate = function(x) {
      x$source_data <- list(list(
        path = "source.tsv", role = "input", target_subdir = "registered",
        allowlisted_paths = "source.tsv", source_rol = "typo"
      )); x
    })
  )
  for (mutation in mutations) {
    error <- tryCatch(
      lisaR:::lisa_validate_pipeline_config(mutation$mutate(strict_config_fixture())),
      error = identity
    )
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error), "LISA-CONFIG-UNKNOWN-001", fixed = TRUE)
    expect_match(conditionMessage(error), paste0("field=", mutation$path), fixed = TRUE)
  }
})

test_that("strict FALSE ignores unknown keys but still validates known fields", {
  cfg <- strict_config_fixture()
  cfg$future_top_level <- list(note = "ignored by compatibility validation")
  cfg$pipeline$future_pipeline_key <- "ignored"
  cfg$single_de[[1]]$future_analysis_key <- "ignored"
  cfg$contrasts[[1]]$future_contrast_key <- "ignored"
  cfg$report$future_report_key <- "ignored"
  cfg$report$formats$jpeg <- TRUE
  expect_silent(
    validated <- lisaR:::lisa_validate_pipeline_config(cfg, strict = FALSE)
  )
  expect_identical(validated$report$formats[["png"]], TRUE)

  cfg$report$formats$png <- "treu"
  expect_error(
    lisaR:::lisa_validate_pipeline_config(cfg, strict = FALSE),
    "LISA-CONFIG-BOOL-001.*field=report[.]formats[.]png"
  )
})

test_that("strict is one explicit non-missing logical value", {
  for (value in list("false", 0L, NA, c(TRUE, FALSE))) {
    expect_error(
      lisaR:::lisa_validate_pipeline_config(
        strict_config_fixture(), strict = value
      ),
      "LISA-CONFIG-BOOL-001.*field=strict"
    )
  }
})

test_that("the public API rejects nested typos before table flattening", {
  cfg <- strict_config_fixture()
  cfg$single_de[[1]]$logfc_clo <- "log2FC"
  expect_error(
    validate_lisa_config(cfg, check_files = FALSE),
    "LISA-CONFIG-UNKNOWN-001.*field=single_de\\[1\\][.]logfc_clo"
  )
})

test_that("YAML, JSON, and in-memory run input retain strict failures", {
  skip_if_not_installed("yaml")
  root <- tempfile("lisa-strict-files-")
  dir.create(root)
  cases <- list(
    json = function(x, path) jsonlite::write_json(
      x, path, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"
    ),
    yaml = function(x, path) yaml::write_yaml(x, path)
  )
  for (extension in names(cases)) {
    cfg <- strict_config_fixture()
    cfg$single_de[[1]]$logfc_clo <- "log2FC"
    path <- file.path(root, paste0("typo.", extension))
    cases[[extension]](cfg, path)
    expect_error(
      validate_lisa_config(path, check_files = FALSE),
      "LISA-CONFIG-UNKNOWN-001.*field=single_de\\[1\\][.]logfc_clo"
    )

    cfg <- strict_config_fixture()
    cfg$pipeline$dry_run <- "treu"
    path <- file.path(root, paste0("boolean.", extension))
    cases[[extension]](cfg, path)
    expect_error(
      validate_lisa_config(path, check_files = FALSE),
      "LISA-CONFIG-BOOL-001.*field=pipeline[.]dry_run"
    )
  }

  cfg <- strict_config_fixture()
  cfg$pipeline$run_ora <- NA
  expect_error(
    run_lisa(cfg),
    "LISA-CONFIG-BOOL-001.*field=pipeline[.]run_ora"
  )
})

test_that("R allowlists and JSON schema agree for strict nested objects", {
  schema <- jsonlite::fromJSON(
    system.file("schema", "lisa-config.schema.json", package = "lisaR"),
    simplifyVector = FALSE
  )
  expect_false(schema$`$defs`$singleAnalysis$additionalProperties)
  expect_false(schema$`$defs`$contrast$additionalProperties)
  expect_false(schema$`$defs`$duplicatePolicy$oneOf[[2]]$additionalProperties)
  expect_setequal(
    names(schema$`$defs`$singleAnalysis$properties),
    lisaR:::lisa_config_allowed_keys("analysis")
  )
  expect_setequal(
    names(schema$`$defs`$contrast$properties),
    lisaR:::lisa_config_allowed_keys("contrast")
  )
  expect_setequal(
    names(schema$`$defs`$duplicatePolicy$oneOf[[2]]$properties),
    lisaR:::lisa_config_allowed_keys("duplicate_policy")
  )
  expect_false("package_dir" %in% names(schema$properties$pipeline$properties))
  expect_false("package_dir" %in% lisaR:::lisa_config_allowed_keys("pipeline"))
  expect_false(schema$properties$report$additionalProperties)
  expect_false(schema$properties$report$properties$formats$additionalProperties)
  expect_false(schema$properties$pipeline$properties$evidence$additionalProperties)
  expect_identical(
    schema$properties$pipeline$properties$evidence$required,
    list("required_columns")
  )
  evidence_rule <- Filter(function(rule)
    identical(rule$`if`$properties$evidence_mode$const, "custom"),
    schema$properties$pipeline$allOf)
  expect_length(evidence_rule, 1L)
  expect_identical(evidence_rule[[1L]]$then$required, list("evidence"))
  expect_false(schema$properties$source_data$items$additionalProperties)
  expect_setequal(
    names(schema$properties$source_data$items$properties),
    lisaR:::lisa_config_allowed_keys("source_data")
  )
  for (key in c("run_kegg_maps", "run_ora", "run_hallmarks", "dry_run", "resume")) {
    expect_identical(schema$properties$pipeline$properties[[key]]$type, "boolean")
  }
  for (key in c("source_data", "recipes")) {
    expect_identical(schema$properties$report$properties[[key]]$type, "boolean")
  }
})

test_that("R and JSON schema share the portable technical-ID contract", {
  schema <- jsonlite::fromJSON(
    system.file("schema", "lisa-config.schema.json", package = "lisaR"),
    simplifyVector = FALSE
  )
  portable <- schema$`$defs`$portableId
  expect_identical(portable$pattern, lisaR:::lisa_portable_id_pattern())

  expected_ref <- "#/$defs/portableId"
  references <- c(
    schema$properties$pipeline$properties$file_label_prefix$`$ref`,
    schema$properties$pipeline$properties$run_id$`$ref`,
    schema$properties$source_data$items$properties$target_subdir$`$ref`,
    schema$`$defs`$singleAnalysis$properties$analysis_id$`$ref`,
    schema$`$defs`$contrast$properties$contrast_id$`$ref`,
    schema$`$defs`$contrast$properties$output_id$`$ref`
  )
  expect_identical(unname(unlist(references)), rep(expected_ref, 6L))

  cases <- data.frame(
    value = c(
      "A", "sample-A_1", "GOBP-C2", "snapshot.2026-08-30", "COM10",
      "CON", "con.txt", "LPT9.csv", "trailing.", "trailing ",
      ".hidden", "A:B", "A/B", "..", "caf\u00e9", "nai\u0308ve", "\u03b1"
    ),
    accepted = c(rep(TRUE, 5L), rep(FALSE, 12L)),
    stringsAsFactors = FALSE
  )
  schema_accepts <- function(value) {
    grepl(portable$pattern, value, perl = TRUE, useBytes = TRUE) &&
      !grepl(portable$not$pattern, value, perl = TRUE, useBytes = TRUE)
  }
  r_accepts <- function(value) {
    isTRUE(tryCatch({
      lisaR:::lisa_safe_id(value)
      TRUE
    }, error = function(e) FALSE))
  }
  observed_schema <- vapply(cases$value, schema_accepts, logical(1))
  observed_r <- vapply(cases$value, r_accepts, logical(1))
  expect_identical(unname(observed_schema), cases$accepted)
  expect_identical(unname(observed_r), cases$accepted)
  expect_identical(observed_schema, observed_r)
})

test_that("Unicode remains available in human-facing labels and titles", {
  cfg <- strict_config_fixture()
  cfg$pipeline$report_title <- "Respuesta inmunitaria: caf\u00e9 y se\u00f1al \u03b1"
  cfg$single_de[[1]]$label <- "Respondedores en tratamiento"
  cfg$single_de[[1]]$comparison <- "Despu\u00e9s frente a antes"
  cfg$contrasts[[1]]$contrast_label <- "Din\u00e1mica de la respuesta"
  expect_silent(lisaR:::lisa_validate_pipeline_config(cfg))
})

test_that("all bundled study configurations remain valid under strict key checks", {
  paths <- c(
    system.file("examples", "minimal-study.yaml", package = "lisaR"),
    system.file("examples", "minimal-study.json", package = "lisaR"),
    system.file("examples", "quick-start", "study.yml", package = "lisaR"),
    system.file("examples", "riaz-gse91061", "config", "riaz-gse91061.yml", package = "lisaR"),
    system.file("examples", "riaz-gse91061", "config", "riaz-gse91061.json", package = "lisaR")
  )
  expect_true(all(nzchar(paths)))
  for (path in paths) {
    cfg <- lisaR:::read_lisa_pipeline_config(path)
    expect_true(lisaR:::lisa_validate_raw_config(cfg, strict = TRUE), info = path)
    technical_ids <- c(
      unlist(cfg$pipeline[c("run_id", "file_label_prefix")], use.names = FALSE),
      unlist(lapply(cfg$single_de, function(row) row$analysis_id), use.names = FALSE),
      unlist(lapply(cfg$contrasts, function(row) {
        c(row$contrast_id, row$output_id)
      }), use.names = FALSE),
      unlist(lapply(cfg$source_data, function(row) row$target_subdir),
             use.names = FALSE)
    )
    technical_ids <- as.character(technical_ids)
    technical_ids <- technical_ids[!is.na(technical_ids) & nzchar(technical_ids)]
    expect_silent(
      vapply(technical_ids, lisaR:::lisa_safe_id, character(1))
    )
  }
})
