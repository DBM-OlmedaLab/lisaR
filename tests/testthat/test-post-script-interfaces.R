extract_parse_args <- function(path) {
  expressions <- as.list(parse(path))
  is_parse_args <- vapply(expressions, function(expression) {
    is.call(expression) &&
      identical(expression[[1L]], as.name("<-")) &&
      identical(expression[[2L]], as.name("parse_args"))
  }, logical(1))
  stopifnot(sum(is_parse_args) == 1L)

  environment <- new.env(parent = baseenv())
  environment$default_project_dir <- NA_character_
  environment$kegg_species_contract <- function(species) {
    list(scientific_name = species)
  }
  eval(expressions[[which(is_parse_args)]], envir = environment)
  environment$parse_args
}

interface_value <- function(flag) {
  if (flag == "species") return("Homo sapiens")
  if (flag == "kegg-access-mode") return("cache_only")
  if (flag == "scale") return("zscore")
  if (flag == "formats") return("png")
  if (flag == "render-graphs") return("true")
  if (flag %in% c("source-data", "recipes")) return("true")
  # build_selected_LISA_report.R mode switches are strict true/false.
  if (flag %in% c(
    "plan-only", "complete-missing", "assemble", "kegg-maps"
  )) return("true")
  if (flag == "gsea-padj-cutoff") return("0.25")
  if (flag == "overlap-export") return("displayed")
  if (flag == "index-only") return("false")
  if (flag == "kegg-id") return("hsa04110")
  if (flag == "tables-only") return("false")
  if (flag == "exact-product") return("contrast_gene_card")
  if (flag == "variant") return("all_categories_top7")
  if (grepl(
    "^(top|max|min|lfc|color|de-padj|label)",
    flag
  )) return("1")
  "test-value"
}

test_that("every executable installed script has one registered interface", {
  script_dir <- system.file("scripts", package = "lisaR")
  scripts <- list.files(
    script_dir,
    pattern = "^build_.*[.]R$",
    full.names = FALSE
  )
  expect_setequal(names(lisaR:::lisa_post_script_interfaces()), scripts)
})

test_that("registered arguments are accepted by each script parser", {
  script_dir <- system.file("scripts", package = "lisaR")
  interfaces <- lisaR:::lisa_post_script_interfaces()

  for (script_name in names(interfaces)) {
    parser <- extract_parse_args(file.path(script_dir, script_name))
    flags <- interfaces[[script_name]]
    values <- vapply(
      flags,
      interface_value,
      character(1)
    )
    args <- as.vector(rbind(paste0("--", flags), values))
    expect_no_error(parser(args))
  }
})

test_that("the orchestrator rejects unknown or malformed script arguments", {
  expect_error(
    lisaR:::lisa_validate_post_script_args(
      "build_single_de_kegg_pathway_painter.R",
      c("--project-dir", "x", "--unknown-option", "y")
    ),
    "Unsupported argument"
  )
  expect_error(
    lisaR:::lisa_validate_post_script_args(
      "build_LISA_report.R",
      c("--project-dir")
    ),
    "flag/value pairs"
  )
  expect_error(
    lisaR:::lisa_validate_post_script_args(
      "build_LISA_report.R",
      c("--project-dir", "safe", "--project-dir", "other")
    ),
    "Duplicate post-processing argument.*--project-dir"
  )
  expect_false(any(vapply(
    lisaR:::lisa_post_script_interfaces(),
    function(flags) "lisa-internal-renderer-sha256" %in% flags,
    logical(1)
  )))
  expect_error(
    lisaR:::lisa_validate_post_script_args(
      "build_LISA_report.R",
      c(
        "--project-dir", "x",
        "--lisa-internal-renderer-sha256", paste(rep("a", 64L), collapse = "")
      )
    ),
    "Unsupported argument"
  )
})

test_that("the root HTML report constrains desktop overflow and renders model notes", {
  script <- system.file("scripts", "build_LISA_report.R", package = "lisaR")
  code <- paste(readLines(script, warn = FALSE), collapse = "\n")
  expect_false(grepl("grid-template-columns:380px", code, fixed = TRUE))
  expect_match(code, "max-width:1580px;width:100%;min-width:0", fixed = TRUE)
  expect_match(code, "lisa_report_shell", fixed = TRUE)
  expect_match(code, "de_index$model_note", fixed = TRUE)
})
