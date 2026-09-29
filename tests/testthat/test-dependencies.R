test_that("optional dependency diagnostics name the package and repair action", {
  message <- tryCatch(
    lisaR:::lisa_require_optional("definitely-not-installed", feature = "the test feature"),
    error = conditionMessage
  )

  expect_match(message, "definitely-not-installed", fixed = TRUE)
  expect_match(message, "Install it with", fixed = TRUE)
  expect_match(message, "install.packages", fixed = TRUE)
  expect_match(message, "--profile=full", fixed = TRUE)
})

test_that("declared optional package gates report availability", {
  status <- lisa_optional_dependency_status(c("jsonlite", "definitely-not-installed"))

  expect_equal(status$package, c("jsonlite", "definitely-not-installed"))
  expect_type(status$available, "logical")
  expect_false(status$available[[2]])
})

test_that("configuration dependencies follow the requested optional outputs", {
  none_available <- function(package) FALSE

  core <- lisaR:::lisa_config_dependency_requirements(
    pipeline = list(
      run_gene_level = FALSE,
      run_reports = FALSE,
      run_kegg_maps = FALSE,
      export_formats = "tsv"
    ),
    species = "Homo sapiens",
    available = none_available
  )
  expect_setequal(
    core$package,
    c("fgsea", "ggplot2", "hommel", "jsonlite", "openxlsx", "yaml")
  )

  full <- lisaR:::lisa_config_dependency_requirements(
    pipeline = list(
      run_gene_level = TRUE,
      run_reports = TRUE,
      run_kegg_maps = TRUE,
      export_formats = c("tsv", "xlsx")
    ),
    species = "Homo sapiens",
    report = list(mode = "full"),
    available = none_available
  )
  retired_layout_package <- paste0("cow", "plot")
  expect_true(all(c(
    "ggrepel", "gridExtra", "patchwork",
    "AnnotationDbi", "KEGGREST", "png", "org.Hs.eg.db"
  ) %in% full$package))
  expect_false(retired_layout_package %in% full$package)
  expect_false("org.Mm.eg.db" %in% full$package)
})

test_that("missing configured dependencies fail before analysis with repair commands", {
  message <- tryCatch(
    lisaR:::lisa_assert_config_dependencies(
      pipeline = list(
        run_gene_level = FALSE,
        run_reports = FALSE,
        run_kegg_maps = TRUE,
        export_formats = "tsv"
      ),
      species = "Homo sapiens",
      report = list(mode = "standard"),
      available = function(package) FALSE
    ),
    error = conditionMessage
  )

  expect_match(message, "LISA-DEPENDENCY-001", fixed = TRUE)
  expect_match(message, "install.packages", fixed = TRUE)
  expect_match(message, "BiocManager::install", fixed = TRUE)
  expect_match(message, "--profile=full", fixed = TRUE)
})
