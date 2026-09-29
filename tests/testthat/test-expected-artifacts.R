lisa_expected_artifact_fixture <- function(root, requested_report_mode = "standard",
                                            run_ora = FALSE) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  de_index <- data.frame(
    analysis_id = c("A", "B"),
    species = "Homo sapiens",
    stringsAsFactors = FALSE
  )
  contrast_index <- data.frame(
    contrast_id = "A_vs_B",
    analysis_a = "A",
    analysis_b = "B",
    stringsAsFactors = FALSE
  )
  registry <- lisaR:::lisa_collection_registry("GOBP-C2")
  registry$run_ora <- registry$run_ora & isTRUE(run_ora)
  expected <- lisaR:::lisa_expected_artifacts(
    de_index = de_index,
    contrast_index = contrast_index,
    registry = registry,
    output_dir = root,
    file_label_prefix = "semantic",
    run_reports = TRUE,
    run_ora = run_ora,
    plot_formats = "png",
    export_formats = "tsv",
    source_data = TRUE,
    recipes = FALSE,
    requested_report_mode = requested_report_mode
  )
  list(
    expected = expected,
    single_status = data.frame(
      analysis_id = c("A", "B"), collection = "GOBP-C2",
      status = "completed", stringsAsFactors = FALSE
    ),
    contrast_status = data.frame(
      contrast_id = "A_vs_B", collection = "GOBP-C2",
      status = "completed", stringsAsFactors = FALSE
    ),
    post_lisa_status = data.frame(
      stage = c("root_html_report", rep("single_de_category_inference", 2)),
      analysis_id = c("", "A", "B"), contrast_id = "",
      collection = "GOBP-C2", status = "completed",
      stringsAsFactors = FALSE
    )
  )
}

lisa_materialize_expected_witnesses <- function(expected, root) {
  required <- expected[expected$expectation == "required", , drop = FALSE]
  for (encoded in required$witnesses) {
    if (!nzchar(encoded)) next
    witnesses <- strsplit(encoded, "|", fixed = TRUE)[[1L]]
    for (witness in witnesses) {
      if (startsWith(witness, "dir:")) {
        dir.create(
          file.path(root, substring(witness, 5L)),
          recursive = TRUE, showWarnings = FALSE
        )
      } else if (startsWith(witness, "file:")) {
        path <- file.path(root, substring(witness, 6L))
        dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
        writeLines(c("field", "material"), path)
      } else {
        stop("unknown test witness: ", witness)
      }
    }
  }
  invisible(root)
}

test_that("expected artifact ledger classifies canonical and external scopes", {
  root <- tempfile("lisaR-expected-artifacts-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  standard <- lisa_expected_artifact_fixture(root)$expected

  expect_identical(
    standard$expectation[standard$stage == "single_de_lisa"],
    rep("required", 2L)
  )
  single_witnesses <- standard$witnesses[
    standard$stage == "single_de_lisa"
  ]
  expect_true(all(grepl("_GSEA_universe_ledger.tsv", single_witnesses,
                        fixed = TRUE)))
  expect_true(all(grepl("_GSEA_OTHER_UNCLASSIFIED.tsv", single_witnesses,
                        fixed = TRUE)))
  expect_identical(
    standard$expectation[standard$stage == "category_contrasts"],
    "required"
  )
  expect_identical(
    standard$expectation[standard$stage == "single_de_ora"],
    rep("optional", 2L)
  )
  expect_identical(
    standard$expectation[standard$stage == "root_html_report"],
    "required"
  )
  expect_identical(
    standard$expectation[
      standard$stage == "full_report_extension_contract"
    ],
    "optional"
  )

  full <- lisa_expected_artifact_fixture(
    root, requested_report_mode = "full"
  )$expected
  extension <- full[full$stage == "full_report_extension_contract", ]
  expect_identical(extension$expectation, "required")
  expect_identical(extension$validation_scope, "external_extension")
  expect_match(extension$external_contract, "extension_receipt.tsv", fixed = TRUE)
  expect_match(extension$external_contract, "extension_inventory.tsv", fixed = TRUE)
  expect_match(extension$external_contract, "extension_manifest.tsv", fixed = TRUE)

  pathways <- lisaR:::lisa_expected_artifacts(
    de_index = data.frame(analysis_id = "A", species = "Homo sapiens"),
    contrast_index = NULL,
    registry = lisaR:::lisa_collection_registry("PATHWAYS"),
    output_dir = root,
    file_label_prefix = "semantic",
    run_reports = FALSE
  )
  pathway_witness <- pathways$witnesses[
    pathways$artifact_id == "single_de_lisa:A:PATHWAYS"
  ]
  expect_match(pathway_witness, "_GSEA_universe_ledger.tsv", fixed = TRUE)
  expect_false(grepl("_GSEA_OTHER_UNCLASSIFIED.tsv", pathway_witness,
                     fixed = TRUE))
})

test_that("required canonical artifacts need completed status and material witnesses", {
  root <- tempfile("lisaR-artifact-validation-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  fixture <- lisa_expected_artifact_fixture(root)
  lisa_materialize_expected_witnesses(fixture$expected, root)

  validation <- lisaR:::lisa_reconcile_expected_artifacts(
    fixture$expected,
    single_status = fixture$single_status,
    contrast_status = fixture$contrast_status,
    post_lisa_status = fixture$post_lisa_status,
    output_dir = root
  )
  required <- validation[validation$expectation == "required", ]
  expect_true(all(required$validation_status == "validated"))
  expect_true(lisaR:::lisa_assert_required_artifacts(validation))

  omitted <- fixture$single_status[fixture$single_status$analysis_id != "A", ]
  missing_status <- lisaR:::lisa_reconcile_expected_artifacts(
    fixture$expected,
    single_status = omitted,
    contrast_status = fixture$contrast_status,
    post_lisa_status = fixture$post_lisa_status,
    output_dir = root
  )
  a_row <- missing_status[
    missing_status$artifact_id == "single_de_lisa:A:GOBP-C2", , drop = FALSE
  ]
  expect_identical(a_row$validation_status, "missing_status")
  expect_error(
    lisaR:::lisa_assert_required_artifacts(missing_status),
    "LISA-ARTIFACT-001.*single_de_lisa:A:GOBP-C2",
    class = "lisa_artifact_validation_error"
  )
})

test_that("deleted and zero-byte canonical witnesses fail closed", {
  root <- tempfile("lisaR-artifact-witness-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  fixture <- lisa_expected_artifact_fixture(root)
  lisa_materialize_expected_witnesses(fixture$expected, root)
  expected_a <- fixture$expected[
    fixture$expected$artifact_id == "single_de_lisa:A:GOBP-C2",
    , drop = FALSE
  ]
  witnesses <- strsplit(expected_a$witnesses, "|", fixed = TRUE)[[1L]]
  gsea <- witnesses[grepl("_GSEA_semantic_annotated[.]tsv$", witnesses)]
  summary <- witnesses[grepl("_GSEA_category_summary[.]tsv$", witnesses)]
  expect_length(gsea, 1L)
  expect_length(summary, 1L)

  unlink(file.path(root, substring(gsea, 6L)))
  deleted <- lisaR:::lisa_reconcile_expected_artifacts(
    fixture$expected, fixture$single_status, fixture$contrast_status,
    fixture$post_lisa_status, root
  )
  deleted_a <- deleted[
    deleted$artifact_id == "single_de_lisa:A:GOBP-C2", , drop = FALSE
  ]
  expect_identical(deleted_a$validation_status, "missing_witness")
  expect_match(deleted_a$missing_witnesses, "_GSEA_semantic_annotated.tsv", fixed = TRUE)
  expect_error(
    lisaR:::lisa_assert_required_artifacts(deleted),
    class = "lisa_artifact_validation_error"
  )

  lisa_materialize_expected_witnesses(fixture$expected, root)
  summary_path <- file.path(root, substring(summary, 6L))
  writeBin(raw(), summary_path)
  expect_identical(as.numeric(file.info(summary_path)$size), 0)
  zero <- lisaR:::lisa_reconcile_expected_artifacts(
    fixture$expected, fixture$single_status, fixture$contrast_status,
    fixture$post_lisa_status, root
  )
  zero_a <- zero[
    zero$artifact_id == "single_de_lisa:A:GOBP-C2", , drop = FALSE
  ]
  expect_identical(zero_a$validation_status, "missing_witness")
  expect_match(zero_a$missing_witnesses, "_GSEA_category_summary.tsv", fixed = TRUE)
})

test_that("full mode validates only its explicit external-extension contract", {
  root <- tempfile("lisaR-full-extension-contract-")
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  fixture <- lisa_expected_artifact_fixture(
    root, requested_report_mode = "full"
  )
  lisa_materialize_expected_witnesses(fixture$expected, root)
  lisaR:::write_lisa_tsv(
    data.frame(key = "report_mode", value = "full"),
    file.path(root, "contract_manifest.tsv")
  )
  lisaR:::write_lisa_tsv(
    data.frame(product = "gene_level", status = "disabled"),
    file.path(root, "product_plan.tsv")
  )

  validation <- lisaR:::lisa_reconcile_expected_artifacts(
    fixture$expected, fixture$single_status, fixture$contrast_status,
    fixture$post_lisa_status, root
  )
  extension <- validation[
    validation$stage == "full_report_extension_contract", , drop = FALSE
  ]
  expect_identical(extension$validation_scope, "external_extension")
  expect_identical(
    extension$validation_status, "external_contract_validated"
  )
  expect_false(file.exists(file.path(root, "extension_receipt.tsv")))
  expect_true(lisaR:::lisa_assert_required_artifacts(validation))
})

test_that("dry pipeline materializes expectations before any computation", {
  base <- system.file("extdata", "minimal", package = "lisaR")
  skip_if_not(nzchar(base), "minimal installed fixture is unavailable")
  output <- tempfile("lisaR-expected-plan-")
  on.exit(unlink(output, recursive = TRUE, force = TRUE), add = TRUE)

  result <- lisaR:::run_lisa_pipeline(
    de_index = file.path(base, "de_index.tsv"),
    contrast_index = file.path(base, "contrast_index.tsv"),
    dictionary_dir = "unused",
    term2gene = "unused",
    output_dir = output,
    collections = "GOBP-C2",
    base_dir = base,
    run_reports = TRUE,
    workers = 1L,
    dry_run = TRUE
  )
  path <- file.path(output, "expected_artifacts.tsv")
  expect_true(file.exists(path))
  expected <- utils::read.delim(
    path, check.names = FALSE, stringsAsFactors = FALSE
  )
  expect_gt(nrow(expected), 0L)
  expect_true(all(c(
    "single_de_lisa", "category_contrasts", "root_html_report"
  ) %in% expected$stage))
  expect_false(file.exists(file.path(output, "artifact_validation.tsv")))
  expect_gt(nrow(result$expected_artifacts), 0L)
  expect_identical(nrow(result$artifact_validation), 0L)
})
