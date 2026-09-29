runner_r_path_literal <- function(path) {
  encodeString(
    normalizePath(path, winslash = "/", mustWork = FALSE),
    quote = '"'
  )
}

runner_script_candidates <- function() {
  c(
    file.path("tools", "maintenance", "run_lisa_prepared_project.R"),
    testthat::test_path("..", "..", "tools", "maintenance",
                        "run_lisa_prepared_project.R")
  )
}

# The runner refuses to start unless its reviewed inventory equals R/ exactly.
# A fixture that hardcodes its own copy of that list therefore rots silently
# every time an R/ file is added, and the failure surfaces as an unrelated
# "inventory does not match" error instead of the boundary being tested. Read
# the list out of the runner itself so there is one source of truth. Only the
# `source_files <- c(...)` expression is evaluated; nothing else in the runner
# is executed.
runner_reviewed_sources <- function(runner) {
  env <- new.env(parent = baseenv())
  found <- FALSE
  for (expression in as.list(parse(runner))) {
    if (is.call(expression) &&
        identical(expression[[1L]], as.name("<-")) &&
        identical(expression[[2L]], as.name("source_files"))) {
      eval(expression, envir = env)
      found <- TRUE
    }
  }
  stopifnot(
    found,
    is.character(env$source_files),
    length(env$source_files) > 0L,
    !anyDuplicated(env$source_files)
  )
  env$source_files
}

test_that("the maintenance runner's reviewed inventory equals the package R/ directory", {
  runner <- runner_script_candidates()[file.exists(runner_script_candidates())][1L]
  skip_if(is.na(runner) || !nzchar(runner),
          "maintenance tools are excluded from the source tarball")
  runner <- normalizePath(runner, winslash = "/", mustWork = TRUE)
  package_dir <- normalizePath(
    file.path(dirname(runner), "..", ".."), winslash = "/", mustWork = TRUE
  )
  observed <- sort(list.files(file.path(package_dir, "R"), pattern = "[.]R$"))
  expect_gt(length(observed), 0L)
  # Exact equality, in both directions: the runner's own check is a hard stop,
  # so a stale inventory disables the tool completely.
  expect_identical(sort(runner_reviewed_sources(runner)), observed)
})

test_that("prepared-project runner rejects unreviewed source files before sourcing", {
  candidates <- c(
    file.path("tools", "maintenance", "run_lisa_prepared_project.R"),
    testthat::test_path("..", "..", "tools", "maintenance", "run_lisa_prepared_project.R")
  )
  runner <- candidates[file.exists(candidates)][1L]
  skip_if(is.na(runner) || !nzchar(runner), "maintenance tools are excluded from the source tarball")
  runner <- normalizePath(runner, mustWork = TRUE)

  fixture <- tempfile("lisa-runner-shadow-")
  shadow_package <- file.path(fixture, "shadow-lisaR")
  project_dir <- file.path(fixture, "project")
  dictionary_dir <- file.path(fixture, "dictionaries")
  project_root <- file.path(fixture, "lisa-project")
  dir.create(file.path(shadow_package, "R"), recursive = TRUE)
  dir.create(project_dir, recursive = TRUE)
  dir.create(dictionary_dir, recursive = TRUE)
  dir.create(project_root, recursive = TRUE)
  dictionary <- file.path(dictionary_dir, "dictionary.tsv")
  category_map <- file.path(dictionary_dir, "category-map.tsv")
  term2gene <- file.path(fixture, "term2gene.tsv")
  writeLines("fixture", dictionary)
  writeLines("fixture", category_map)
  writeLines("fixture", term2gene)
  writeLines(
    c("Package: lisaR", "Version: 0.6.0", "Title: shadow fixture"),
    file.path(shadow_package, "DESCRIPTION")
  )
  sentinel <- file.path(fixture, "SENTINEL_EXECUTED")
  writeLines(
    sprintf("writeLines('executed', %s)", runner_r_path_literal(sentinel)),
    file.path(shadow_package, "R", "zzz_sentinel.R")
  )

  output <- suppressWarnings(system2(
    lisaR:::lisa_rscript_executable(),
    c(
      "--vanilla", shQuote(runner),
      "--package-dir", shQuote(shadow_package),
      "--project-dir", shQuote(project_dir),
      "--dictionary", shQuote(dictionary),
      "--category-map", shQuote(category_map),
      "--term2gene", shQuote(term2gene),
      "--lisa-project-root", shQuote(project_root),
      "--report-title", shQuote("shadow boundary fixture"),
      "--collections", "HALLMARKS"
    ),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_false(identical(status, 0L))
  expect_match(paste(output, collapse = "\n"), "reviewed lisaR source inventory")
  expect_false(file.exists(sentinel))
})

test_that("prepared-project runner rejects symlinked reviewed sources", {
  candidates <- c(
    file.path("tools", "maintenance", "run_lisa_prepared_project.R"),
    testthat::test_path("..", "..", "tools", "maintenance", "run_lisa_prepared_project.R")
  )
  runner <- candidates[file.exists(candidates)][1L]
  skip_if(is.na(runner) || !nzchar(runner),
          "maintenance tools are excluded from the source tarball")
  runner <- normalizePath(runner, mustWork = TRUE)
  # Derived from the runner, never duplicated here.
  reviewed_sources <- runner_reviewed_sources(runner)
  fixture <- tempfile("lisa-runner-source-link-")
  shadow_package <- file.path(fixture, "shadow-lisaR")
  source_dir <- file.path(shadow_package, "R")
  project_dir <- file.path(fixture, "project")
  dictionary_dir <- file.path(fixture, "dictionaries")
  project_root <- file.path(fixture, "lisa-project")
  dir.create(source_dir, recursive = TRUE)
  dir.create(project_dir, recursive = TRUE)
  dir.create(dictionary_dir, recursive = TRUE)
  dir.create(project_root, recursive = TRUE)
  dictionary <- file.path(dictionary_dir, "dictionary.tsv")
  category_map <- file.path(dictionary_dir, "category-map.tsv")
  term2gene <- file.path(fixture, "term2gene.tsv")
  writeLines("fixture", dictionary)
  writeLines("fixture", category_map)
  writeLines("fixture", term2gene)
  writeLines(
    c("Package: lisaR", "Version: 0.6.0", "Title: shadow fixture"),
    file.path(shadow_package, "DESCRIPTION")
  )
  expect_true(all(file.create(file.path(source_dir, reviewed_sources))))
  linked_source <- file.path(source_dir, "trusted_rds.R")
  expect_true(file.remove(linked_source))
  sentinel <- file.path(fixture, "SYMLINK_SOURCE_EXECUTED")
  external <- file.path(fixture, "external-source.R")
  writeLines(
    sprintf("writeLines('executed', %s)", runner_r_path_literal(sentinel)),
    external
  )
  if (!isTRUE(suppressWarnings(file.symlink(external, linked_source)))) {
    skip("symbolic links are unavailable on this platform")
  }
  output <- suppressWarnings(system2(
    lisaR:::lisa_rscript_executable(),
    c(
      "--vanilla", shQuote(runner),
      "--package-dir", shQuote(shadow_package),
      "--project-dir", shQuote(project_dir),
      "--dictionary", shQuote(dictionary),
      "--category-map", shQuote(category_map),
      "--term2gene", shQuote(term2gene),
      "--lisa-project-root", shQuote(project_root),
      "--report-title", shQuote("symlink boundary fixture")
    ),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_false(identical(status, 0L))
  expect_match(paste(output, collapse = "\n"), "symbolic link")
  expect_false(file.exists(sentinel))
})

test_that("prepared-project runner rejects symlinked explicit resources", {
  candidates <- c(
    file.path("tools", "maintenance", "run_lisa_prepared_project.R"),
    testthat::test_path("..", "..", "tools", "maintenance", "run_lisa_prepared_project.R")
  )
  runner <- candidates[file.exists(candidates)][1L]
  skip_if(is.na(runner) || !nzchar(runner),
          "maintenance tools are excluded from the source tarball")
  runner <- normalizePath(runner, winslash = "/", mustWork = TRUE)
  package_dir <- normalizePath(
    file.path(dirname(runner), "..", ".."), winslash = "/", mustWork = TRUE
  )

  fixture <- tempfile("lisa-runner-resource-link-")
  project_dir <- file.path(fixture, "project")
  project_root <- file.path(fixture, "lisa-project")
  dir.create(project_dir, recursive = TRUE)
  dir.create(project_root, recursive = TRUE)
  dictionary_target <- file.path(fixture, "dictionary-target.tsv")
  dictionary_link <- file.path(fixture, "dictionary-link.tsv")
  category_map <- file.path(fixture, "category-map.tsv")
  term2gene <- file.path(fixture, "term2gene.tsv")
  writeLines("fixture", dictionary_target)
  writeLines("fixture", category_map)
  writeLines("fixture", term2gene)
  if (!isTRUE(suppressWarnings(file.symlink(
      dictionary_target, dictionary_link
  )))) {
    skip("symbolic links are unavailable on this platform")
  }

  output <- suppressWarnings(system2(
    lisaR:::lisa_rscript_executable(),
    c(
      "--vanilla", shQuote(runner),
      "--package-dir", shQuote(package_dir),
      "--project-dir", shQuote(project_dir),
      "--dictionary", shQuote(dictionary_link),
      "--category-map", shQuote(category_map),
      "--term2gene", shQuote(term2gene),
      "--lisa-project-root", shQuote(project_root),
      "--report-title", shQuote("resource symlink boundary fixture")
    ),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_false(identical(status, 0L))
  expect_match(paste(output, collapse = "\n"), "rejects symbolic-link")
})

test_that("prepared-project runner accepts dot syntax for its reviewed source", {
  candidates <- c(
    file.path("tools", "maintenance", "run_lisa_prepared_project.R"),
    testthat::test_path("..", "..", "tools", "maintenance",
                        "run_lisa_prepared_project.R")
  )
  runner <- candidates[file.exists(candidates)][1L]
  skip_if(is.na(runner) || !nzchar(runner),
          "maintenance tools are excluded from the source tarball")
  runner <- normalizePath(runner, winslash = "/", mustWork = TRUE)
  package_dir <- normalizePath(
    file.path(dirname(runner), "..", ".."), winslash = "/", mustWork = TRUE
  )

  fixture <- tempfile("lisa-runner-dot-path-")
  project_dir <- file.path(fixture, "project")
  project_root <- file.path(fixture, "lisa-project")
  dir.create(project_dir, recursive = TRUE)
  dir.create(project_root, recursive = TRUE)
  on.exit(unlink(fixture, recursive = TRUE, force = TRUE), add = TRUE)
  dictionary <- file.path(fixture, "dictionary.tsv")
  category_map <- file.path(fixture, "category-map.tsv")
  term2gene <- file.path(fixture, "term2gene.tsv")
  writeLines("fixture", dictionary)
  writeLines("fixture", category_map)
  writeLines("fixture", term2gene)

  previous <- setwd(package_dir)
  on.exit(setwd(previous), add = TRUE)
  output <- suppressWarnings(system2(
    lisaR:::lisa_rscript_executable(),
    c(
      "--vanilla", shQuote(runner),
      "--package-dir", ".",
      "--project-dir", shQuote(project_dir),
      "--dictionary", shQuote(dictionary),
      "--category-map", shQuote(category_map),
      "--term2gene", shQuote(term2gene),
      "--lisa-project-root", shQuote(project_root),
      "--report-title", shQuote("dot-path boundary fixture")
    ),
    stdout = TRUE, stderr = TRUE
  ))
  setwd(previous)
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  message <- paste(output, collapse = "\n")
  expect_false(identical(status, 0L))
  expect_match(message, "Missing de_index")
  expect_false(grepl("rejects symbolic-link", message, fixed = TRUE))
})

test_that("prepared-project runner loads trusted-RDS and lollipop helpers before run_LISA_DE", {
  candidates <- c(
    file.path("tools", "maintenance", "run_lisa_prepared_project.R"),
    testthat::test_path("..", "..", "tools", "maintenance", "run_lisa_prepared_project.R")
  )
  runner <- candidates[file.exists(candidates)][1L]
  skip_if(is.na(runner) || !nzchar(runner), "maintenance tools are excluded from the source tarball")
  runner <- normalizePath(runner, mustWork = TRUE)

  # Derived from the runner, never duplicated here.
  reviewed_sources <- runner_reviewed_sources(runner)
  fixture <- tempfile("lisa-runner-trusted-rds-")
  shadow_package <- file.path(fixture, "shadow-lisaR")
  project_dir <- file.path(fixture, "project")
  dictionary_dir <- file.path(fixture, "dictionaries")
  project_root <- file.path(fixture, "lisa-project")
  source_dir <- file.path(shadow_package, "R")
  dir.create(source_dir, recursive = TRUE)
  dir.create(file.path(project_dir, "config"), recursive = TRUE)
  dir.create(file.path(project_dir, "input"), recursive = TRUE)
  dir.create(file.path(project_dir, "input_matrices_vst"), recursive = TRUE)
  dir.create(dictionary_dir, recursive = TRUE)
  dir.create(project_root, recursive = TRUE)
  dictionary <- file.path(dictionary_dir, "dictionary.tsv")
  category_map <- file.path(dictionary_dir, "category-map.tsv")
  term2gene <- file.path(fixture, "term2gene.tsv")
  writeLines("fixture", dictionary)
  writeLines("fixture", category_map)
  writeLines("fixture", term2gene)
  writeLines(
    c("Package: lisaR", "Version: 0.6.0", "Title: shadow fixture"),
    file.path(shadow_package, "DESCRIPTION")
  )
  expect_true(all(file.create(file.path(source_dir, reviewed_sources))))

  writeLines(
    "lisa_sha256_file <- function(path) paste(rep('0', 64L), collapse = '')",
    file.path(source_dir, "security.R")
  )

  sentinel <- file.path(fixture, "TRUSTED_RDS_HELPER_AVAILABLE")
  writeLines(
    "lisa_prepare_de_input <- function(...) 'trusted-rds-helper-loaded'",
    file.path(source_dir, "trusted_rds.R")
  )
  writeLines(
    "lisa_write_lollipop_figure_source <- function(...) 'lollipop-helper-loaded'",
    file.path(source_dir, "lollipop_reproduction.R")
  )
  writeLines(
    c(
      "run_LISA_DE <- function(...) {",
      "  if (!exists('lisa_prepare_de_input', mode = 'function', inherits = TRUE)) stop('trusted RDS helper missing')",
      "  if (!identical(lisa_prepare_de_input(), 'trusted-rds-helper-loaded')) stop('wrong trusted RDS helper')",
      "  if (!exists('lisa_write_lollipop_figure_source', mode = 'function', inherits = TRUE)) stop('lollipop helper missing')",
      "  if (!identical(lisa_write_lollipop_figure_source(), 'lollipop-helper-loaded')) stop('wrong lollipop helper')",
      sprintf("  writeLines('ready', %s)", runner_r_path_literal(sentinel)),
      "  invisible(TRUE)",
      "}"
    ),
    file.path(source_dir, "run_LISA_DE.R")
  )
  writeLines(
    c(
      "run_lisa_pipeline <- function(...) {",
      "  dots <- list(...)",
      "  source_runtime <- getOption('lisaR.source_runtime_functions', character())",
      "  if (!'lisa_prepare_de_input' %in% source_runtime) stop('source runtime inventory missing trusted RDS helper')",
      "  if (!'lisa_write_lollipop_figure_source' %in% source_runtime) stop('source runtime inventory missing lollipop helper')",
      sprintf("  if (!identical(dots$dictionary_path, %s)) stop('wrong dictionary path')", runner_r_path_literal(dictionary)),
      sprintf("  if (!identical(dots$category_map_path, %s)) stop('wrong category-map path')", runner_r_path_literal(category_map)),
      sprintf("  if (!identical(dots$term2gene, %s)) stop('wrong TERM2GENE path')", runner_r_path_literal(term2gene)),
      "  run_LISA_DE()",
      "}"
    ),
    file.path(source_dir, "lisa_pipeline.R")
  )
  writeLines(
    c(
      "read_lisa_tsv <- function(path) utils::read.delim(path, sep = '\\t', check.names = FALSE, stringsAsFactors = FALSE)",
      "write_lisa_tsv <- function(x, path) {",
      "  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)",
      "  utils::write.table(x, path, sep = '\\t', quote = FALSE, row.names = FALSE, na = '')",
      "  invisible(path)",
      "}"
    ),
    file.path(source_dir, "utils_io.R")
  )

  de_path <- file.path(project_dir, "input", "de.tsv")
  matrix_path <- file.path(project_dir, "input_matrices_vst", "matrix.tsv")
  writeLines("symbol\tlogFC\nG1\t1", de_path)
  writeLines("symbol\tS1\nG1\t1", matrix_path)
  utils::write.table(
    data.frame(
      analysis_id = "fixture", de_path = de_path,
      counts_matrix_path = matrix_path, stringsAsFactors = FALSE
    ),
    file.path(project_dir, "config", "de_index.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  writeLines("contrast_id", file.path(project_dir, "config", "contrast_index.tsv"))
  output <- suppressWarnings(system2(
    lisaR:::lisa_rscript_executable(),
    c(
      "--vanilla", shQuote(runner),
      "--package-dir", shQuote(shadow_package),
      "--project-dir", shQuote(project_dir),
      "--dictionary", shQuote(dictionary),
      "--category-map", shQuote(category_map),
      "--term2gene", shQuote(term2gene),
      "--lisa-project-root", shQuote(project_root),
      "--report-title", shQuote("trusted RDS source fixture"),
      "--collections", "HALLMARKS"
    ),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_identical(status, 0L, info = paste(output, collapse = "\n"))
  expect_true(file.exists(sentinel))
  expect_identical(readLines(sentinel, warn = FALSE), "ready")
  receipt <- utils::read.delim(
    file.path(project_dir, "manifests", "development_runner_inputs.tsv"),
    sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
  )
  expect_setequal(
    receipt$key,
    c(
      "package_dir", "project_dir", "dictionary", "dictionary_sha256",
      "category_map", "category_map_sha256", "term2gene",
      "term2gene_sha256", "lisa_project_root", "report_title", "collections"
    )
  )
  expect_true(all(grepl(
    "^[0-9a-f]{64}$", receipt$value[grepl("_sha256$", receipt$key)]
  )))
})
