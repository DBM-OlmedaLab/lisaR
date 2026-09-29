schema_one_fixture <- function() {
  list(pipeline = list(schema_version = "1.0.0", profile = "targeted",
    evidence_mode = "rank_only", output_dir = "results/schema-one",
    duplicate_policies = list(de_table_duplicate_policy = "error",
      matrix_duplicate_policy = "error", mapped_id_collision_policy = "error")),
    collections = list("GOBP-C2"),
    single_de = list(list(analysis_id = "001", de_path = "unread-input.tsv",
      species = "Homo sapiens", symbol_col = "symbol", rank_col = "stat")))
}

test_that("schema 1.0 and NES selection survive public YAML and JSON validation", {
  root <- tempfile("lisa-schema-one-")
  dir.create(root)
  selections <- list(c("clean", "percentages", "direction", "dispersion"), "clean", "percentages", "direction",
    "dispersion", c("dispersion", "direction"))
  for (selected in selections) {
    cfg <- schema_one_fixture()
    cfg$report <- list(category_nes_variants = as.list(selected))
    json_path <- file.path(root, "study.json")
    yaml_path <- file.path(root, "study.yaml")
    jsonlite::write_json(cfg, json_path, auto_unbox = TRUE, null = "null", digits = 17)
    yaml::write_yaml(cfg, yaml_path)
    for (path in c(json_path, yaml_path)) {
      got <- validate_lisa_config(path, check_files = FALSE)
      expect_true(got$structural_valid)
      expect_identical(got$schema_version, "1.0.0")
      expect_identical(got$config$report$category_nes_variants, as.list(selected))
      expect_identical(got$de_index$analysis_id, "001")
      row <- got$configuration_provenance[
        got$configuration_provenance$field == "report.category_nes_variants", ]
      expect_identical(row$source, "explicit")
    }
  }
  expect_false(dir.exists(file.path(root, "results")))
})

test_that("public validation rejects empty variants without restoring defaults", {
  root <- tempfile("lisa-schema-invalid-")
  dir.create(root)
  for (invalid in list(list(), NULL, FALSE, c("clean", "clean"), "unknown")) {
    cfg <- schema_one_fixture()
    cfg$report <- list(category_nes_variants = invalid)
    path <- file.path(root, "study.json")
    jsonlite::write_json(cfg, path, auto_unbox = TRUE, null = "null")
    expect_error(validate_lisa_config(path, check_files = FALSE), "LISA-CONFIG-NES-VARIANTS-001")
    path <- file.path(root, "study.yaml")
    yaml::write_yaml(cfg, path)
    expect_error(validate_lisa_config(path, check_files = FALSE), "LISA-CONFIG-NES-VARIANTS-001")
  }
  cfg <- schema_one_fixture()
  omitted <- validate_lisa_config(cfg, check_files = FALSE)
  expect_identical(omitted$config$report$category_nes_variants,
    as.list(c("clean", "percentages", "direction", "dispersion")))
  expect_false(dir.exists(file.path(root, "results")))
})

test_that("variant selection and the explicit table-only mode stay independent", {
  cfg <- schema_one_fixture()
  cfg$report <- list(category_nes_variants = list("clean"),
    formats = list(png = FALSE, svg = FALSE, pdf = FALSE))
  got <- validate_lisa_config(cfg, check_files = FALSE)
  expect_true(got$structural_valid)
  expect_identical(got$config$report$category_nes_variants, list("clean"))
  expect_false(any(unlist(got$config$report$formats)))
  cfg$report$category_nes_variants <- list()
  expect_error(validate_lisa_config(cfg, check_files = FALSE), "NES-VARIANTS-001")
})

test_that("all shipped study templates explicitly expose the current presentation options", {
  relative <- c("minimal-study.json", "minimal-study.yaml", "quick-start/study.yml",
    "riaz-gse91061/config/riaz-gse91061.json", "riaz-gse91061/config/riaz-gse91061.yml")
  for (file in relative) {
    path <- system.file("examples", file, package = "lisaR")
    cfg <- if (grepl("[.]json$", file)) jsonlite::read_json(path) else yaml::read_yaml(path)
    expect_identical(cfg$pipeline$schema_version, "1.0.0")
    expect_identical(unlist(cfg$report$category_nes_variants, use.names = FALSE),
      c("clean", "percentages", "direction", "dispersion"))
    expect_identical(cfg$report$mode, "standard")
    expect_true(cfg$report$category_evidence)
    expect_false(cfg$report$legacy_gene_products)
    expect_true(validate_lisa_config(path, check_files = FALSE)$structural_valid)
  }
  schema <- jsonlite::read_json(system.file("schema", "lisa-config.schema.json", package = "lisaR"))
  expect_identical(schema$properties$pipeline$properties$schema_version$const,
    lisaR:::lisa_pipeline_schema_version())
  variants <- schema$properties$report$properties$category_nes_variants
  expect_equal(variants$minItems, 1)
  expect_equal(variants$maxItems, 4)
  expect_true(variants$uniqueItems)
})
