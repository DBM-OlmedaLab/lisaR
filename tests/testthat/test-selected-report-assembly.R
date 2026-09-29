test_that("selected report contracts reject implicit or unknown analyses", {
  env <- new.env(parent = asNamespace("lisaR"))
  sys.source(system.file("scripts", "build_selected_LISA_report.R", package = "lisaR"), env)
  root <- tempfile("selected-contract-")
  dir.create(file.path(root, "config"), recursive = TRUE)
  env$selected_write(data.frame(analysis_id = c("A", "B", "EXCLUDED")), file.path(root, "config/de_index.tsv"))
  env$selected_write(data.frame(contrast_id = c("A_B", "A_E"), output_id = c("A_B", "A_E"),
    contrast_a = "A", contrast_b = c("B", "EXCLUDED")), file.path(root, "config/contrast_index.tsv"))
  contract <- env$selected_contract(root, c("B", "A"), "A_B")
  expect_identical(contract$de$analysis_id, c("B", "A"))
  expect_identical(contract$contrasts$contrast_a, "A")
  expect_identical(contract$contrasts$contrast_b, "B")
  expect_error(env$selected_contract(root, c("A", "B"), "A_E"), "only to selected")
  expect_error(env$selected_contract(root, c("A", "UNKNOWN"), "A_B"), "absent")
  statuses <- data.frame(analysis_id = c("A", "B", "EXCLUDED", ""),
    contrast_id = c("", "", "", "A_E"), status = "completed")
  expect_identical(env$selected_rows(statuses, contract)$analysis_id, c("A", "B"))
  expect_error(env$selected_ids("A,A"), "unique")
  expect_error(env$selected_ids("A,../B"), "unique")
  expect_error(env$parse_args(c("--title", "a", "--title", "b")), "repeated")
})

test_that("declared FULL scope does not silently fall back to folder discovery", {
  env <- new.env(parent = asNamespace("lisaR"))
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), env)
  root <- tempfile("full-selected-"); dir.create(root)
  utils::write.table(data.frame(key = c("analyses", "contrast_owners"), value = c("A;B", "C_C")),
    file.path(root, "report_selection.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  expect_true(env$report_full_scope_owner(root, "A"))
  expect_true(env$report_full_scope_owner(root, "C_C", TRUE))
  expect_false(env$report_full_scope_owner(root, "EXCLUDED"))
  dir.create(file.path(root, "artifacts", "EXCLUDED"), recursive = TRUE)
  expect_false(env$report_full_scope_owner(root, "EXCLUDED"))
})

test_that("the same selected generator supports a one-DE report without contrasts", {
  env <- new.env(parent = asNamespace("lisaR"))
  sys.source(system.file("scripts", "build_selected_LISA_report.R", package = "lisaR"), env)
  root <- tempfile("single-only-")
  dir.create(file.path(root, "config"), recursive = TRUE)
  env$selected_write(data.frame(analysis_id = "protein"), file.path(root, "config/de_index.tsv"))
  env$selected_write(data.frame(contrast_id = character(), output_id = character(),
    contrast_a = character(), contrast_b = character()), file.path(root, "config/contrast_index.tsv"))
  expect_identical(env$selected_ids("none"), character())
  contract <- env$selected_contract(root, "protein", env$selected_ids("none"))
  expect_identical(contract$de$analysis_id, "protein")
  expect_equal(nrow(contract$contrasts), 0L)
  expect_identical(contract$owners, character())
  cfg <- env$parse_args(c("--source-dir", root, "--output-dir", "new",
    "--analyses", "protein", "--contrasts", "none", "--complete-missing", "true"))
  expect_identical(cfg$artifacts_dir, "")
  expect_identical(cfg$complete_missing, "true")
  writeLines("", file.path(root, "config/contrast_index.tsv"))
  blank <- env$selected_contract(root, "protein", character())
  expect_identical(blank$owners, character())
  expect_identical(names(blank$contrasts), c("contrast_id", "output_id", "contrast_a", "contrast_b"))
  expect_error(env$selected_contract(root, "protein", "missing"), "absent")
  writeLines("invalid_header\ncorrupt", file.path(root, "config/contrast_index.tsv"))
  expect_error(env$selected_contract(root, "protein", character()), "incomplete")
})


test_that("saved evidence script anchors are normalized without touching other content", {
  env <- new.env(parent=asNamespace("lisaR"))
  sys.source(system.file("scripts", "build_LISA_report.R", package="lisaR"), env)
  markup <- '<a href="reproduce_contrast_evidence.R" download>Offline base-R figure recipe</a><a href="data.tsv">R values</a>'
  expected <- '<a href="reproduce_contrast_evidence.R" download>R script</a><a href="data.tsv">R values</a>'
  expect_identical(env$report_script_link_labels(markup), expected)
  expect_identical(env$report_script_link_labels(expected), expected)
})


test_that("figure captions use configured example labels instead of storage identifiers", {
  env <- new.env(parent=asNamespace("lisaR"))
  sys.source(system.file("scripts", "build_LISA_report.R", package="lisaR"), env)
  env$report_display_labels <- c(technical_id="Tumor vs adjacent non-tumoral tissue")
  x <- env$pretty_label("technical_id_PATHWAYS_gene_level_top_recurrent_genes.png")
  expect_match(x, "Tumor vs adjacent non-tumoral tissue", fixed=TRUE)
  expect_false(grepl("technical",x))
})
