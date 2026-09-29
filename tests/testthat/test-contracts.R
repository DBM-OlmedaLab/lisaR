test_that("stable package identity agrees with DESCRIPTION and keeps schema contracts", {
  versions <- lisa_contract_versions()
  description <- read.dcf(system.file("DESCRIPTION", package = "lisaR"))

  expect_identical(versions$package_development_version, "1.0.0")
  expect_identical(versions$package_r_version, "1.0.0")
  expect_identical(versions$package_r_version, unname(description[1L, "Version"]))
  expect_identical(unname(description[1L, "Config/lisaR/release-stage"]), "stable")
  expect_identical(unname(description[1L, "Config/lisaR/release-label"]), "1.0.0")
  expect_identical(unname(description[1L, "Config/lisaR/pipeline-production-label"]), "lisaR-1.0.0")
  expect_identical(unname(description[1L, "Config/lisaR/pipeline-production-version"]), "1")
  expect_identical(versions$pipeline_schema_version, "1.0.0")
  expect_identical(versions$report_schema_version, "7.0.0")
  expect_identical(versions$contract_manifest_version, "1.0.0")
})

test_that("minimal contract manifest preserves version identity", {
  manifest <- new_lisa_contract_manifest(
    source_tree = "test-tree",
    configuration_schema_version = "pending"
  )

  expect_s3_class(manifest, "lisa_contract_manifest")
  expect_identical(manifest$source_tree, "test-tree")
  expect_identical(manifest$configuration_schema_version, "pending")
  expect_identical(manifest$report_schema_version, "7.0.0")
})

test_that("mixed-species KEGG contrasts are represented as explicit omissions", {
  code <- paste(deparse(body(lisaR:::run_lisa_pipeline_post_lisa)), collapse = "\n")
  expect_match(code, "skipped_species_mismatch", fixed = TRUE)
  expect_match(code, "shared organism map would be misleading", fixed = TRUE)
  expect_false(grepl("LISA-SPECIES-007", code, fixed = TRUE))
})

test_that("contrast KEGG species inherit from the public analysis columns", {
  contrast <- data.frame(
    contrast_id = "a_minus_b",
    analysis_a = "analysis_a",
    analysis_b = "analysis_b",
    stringsAsFactors = FALSE
  )
  human <- data.frame(
    analysis_id = c("analysis_a", "analysis_b"),
    scientific_name = c("Homo sapiens", "Homo sapiens"),
    stringsAsFactors = FALSE
  )
  mixed <- human
  mixed$scientific_name[[2]] <- "Mus musculus"

  expect_identical(lisaR:::lisa_contrast_species(contrast, human), "Homo sapiens")
  expect_identical(
    lisaR:::lisa_contrast_species(contrast, mixed),
    c("Homo sapiens", "Mus musculus")
  )
  expect_error(
    lisaR:::lisa_contrast_species(contrast, human[1, , drop = FALSE]),
    "LISA-SPECIES-007"
  )
})

test_that("registered category-map order is authoritative", {
  category_map <- data.frame(
    category_id = c("SYN_UNUSED", "SYN_Z", "SYN_A"),
    macrogroup_order = c(1L, 2L, 1L),
    category_order_within_macrogroup = c(1L, 1L, 2L),
    stringsAsFactors = FALSE
  )
  dictionary <- data.frame(
    universe = c("PATHWAYS", "PATHWAYS"),
    category_id = c("SYN_Z", "SYN_A"),
    stringsAsFactors = FALSE
  )

  observed <- lisaR:::prepare_registered_category_map(
    category_map, dictionary
  )
  expect_identical(observed$category_id, c("SYN_A", "SYN_Z"))
  expect_false("SYN_UNUSED" %in% observed$category_id)
})

test_that("fgsea warnings are preserved as notes instead of flooding the console", {
  pathways <- list(set_1 = c("A", "B", "C"))
  ranks <- c(A = 2, B = 2, C = -1)

  expect_warning(
    observed <- lisaR:::run_fgsea_lisa(
      pathways, ranks, 1, 10, 1, 100, 0
    ),
    NA
  )
  expect_match(
    paste(attr(observed, "lisa_notes"), collapse = "\n"),
    "ties in the preranked stats"
  )
})

test_that("canonical pipeline never silently ignores inline gene-level work", {
  blocked_output <- tempfile("lisaR-inline-gene-level-")
  error <- expect_error(
    lisaR:::run_lisa_pipeline(
      de_index = data.frame(),
      dictionary_dir = "unused",
      term2gene = "unused",
      output_dir = blocked_output,
      run_gene_level = TRUE
    ),
    "LISA-GENE-LEVEL-001",
    fixed = TRUE
  )
  expect_match(conditionMessage(error), "plan_lisa_extension()", fixed = TRUE)
  expect_match(conditionMessage(error), "render_lisa_categories()", fixed = TRUE)
  expect_false(lisaR:::lisa_path_entry_exists(blocked_output))

  for (invalid in list("yes", 1, NA, c(FALSE, TRUE))) {
    expect_error(
      lisaR:::run_lisa_pipeline(
        de_index = data.frame(),
        dictionary_dir = "unused",
        term2gene = "unused",
        output_dir = tempfile("lisaR-invalid-gene-level-"),
        run_gene_level = invalid
      ),
      "LISA-CONFIG-BOOL-001.*field=run_gene_level"
    )
  }

  regular_output <- tempfile("lisaR-canonical-no-gene-level-")
  on.exit(unlink(regular_output, recursive = TRUE, force = TRUE), add = TRUE)
  result <- lisaR:::run_lisa_pipeline(
    de_index = data.frame(
      analysis_id = "analysis_a",
      de_path = "unused.tsv",
      species = "Homo sapiens",
      stringsAsFactors = FALSE
    ),
    dictionary_dir = "unused",
    term2gene = "unused",
    output_dir = regular_output,
    collections = "GOBP-C2",
    plot_formats = "png",
    export_formats = "tsv",
    run_gene_level = FALSE,
    run_reports = FALSE,
    source_data = FALSE,
    workers = 1L,
    dry_run = TRUE
  )
  expect_true(dir.exists(regular_output))
  expect_false(any(grepl("gene_level", result$plan$stage, fixed = TRUE)))
})
