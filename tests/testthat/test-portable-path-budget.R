test_that("Windows path budget counts UTF-16 and identifies the destination", {
  expect_equal(lisaR:::lisa_path_utf16_length(c("a", "\U0001F9EC")), c(1, 2))
  expect_invisible(lisaR:::lisa_assert_portable_path("D:/lisa/a.tsv", windows = TRUE))
  long <- paste0("D:/", strrep("a", 240), ".tsv")
  expect_error(lisaR:::lisa_assert_portable_path(long, windows = TRUE), "LISA-PATH-001.*247.*D:/")
  expect_invisible(lisaR:::lisa_assert_portable_path(long, windows = FALSE))
  for (name in c("CON.tsv", "nul", "LPT1.png", "name.", "name "))
    expect_error(lisaR:::lisa_assert_portable_path(paste0("D:/lisa/", name), windows=TRUE), "LISA-PATH-002")
})

test_that("short staging retains identity without repeating the run basename", {
  a <- lisaR:::lisa_short_staging_path("D:/lisa/riaz/out", "run-123")
  expect_equal(nchar(basename(a)), 16L)
  expect_identical(dirname(a), "D:/lisa/riaz")
  expect_identical(a, lisaR:::lisa_short_staging_path("D:/lisa/riaz/out", "run-123"))
  expect_false(identical(a, lisaR:::lisa_short_staging_path("D:/lisa/riaz/other", "run-123")))
  expect_false(identical(a, lisaR:::lisa_short_staging_path("D:/lisa/riaz/out", "run-124")))
})

test_that("contrast file keys separate labels without truncation collisions", {
  a <- lisaR:::lisa_contrast_plot_stem("response_separation_by_time", "RIAZ", "opposite_direction")
  b <- lisaR:::lisa_contrast_plot_stem("response_separation_by_time_2", "RIAZ", "opposite_direction")
  expect_lt(nchar(a), 24)
  expect_false(identical(a,b))
  expect_match(a, "^cx_[a-f0-9]{12}_opp$")
})

test_that("short project map has distinct portable destinations", {
  root <- tempfile("short-project-"); dir.create(root)
  on.exit(lisa_test_cleanup_path(root), add = TRUE)
  template <- system.file("examples", "project_paths.tsv", package="lisaR")
  mapping <- read.delim(template, check.names=FALSE)
  expect_equal(sum(mapping$example=="riaz"),17L)
  expect_equal(sum(mapping$example=="cptac"),6L)
  for (example in unique(mapping$example)) {
    m <- mapping[mapping$example==example,]
    expect_false(anyDuplicated(tolower(m$path)) > 0)
    expect_true(all(nchar(m$path) < 45))
  }
  expect_false(any(grepl("[.][.]|^/|^[A-Za-z]:", mapping$path)))
})


test_that("analysis and contrast paths cannot collide only by letter case", {
  de <- data.frame(analysis_id=c("Case", "case"), de_path=c("a.tsv", "b.tsv"))
  expect_error(lisaR:::validate_lisa_de_index(de, require_files=FALSE), "case-insensitive")
  cx <- data.frame(contrast_id=c("Case", "case"), analysis_a="a", analysis_b="b")
  expect_error(lisaR:::validate_lisa_contrast_index(cx), "case-insensitive")
})


test_that("full extensions reuse shared saved inputs without overwriting conflicts", {
  root <- withr::local_tempdir("shared-input-")
  source <- file.path(root,"source"); work <- file.path(root,"stage","work")
  dir.create(source); dir.create(work,recursive=TRUE)
  withr::local_options(lisaR.run_root=file.path(root,"stage"))
  writeLines("unchanged input",file.path(source,"matrix.tsv"))
  stage <- function() lisaR:::lisa_extension_stage_files(source,work,"matrix.tsv")
  expect_invisible(stage()); expect_invisible(stage())
  expect_identical(readLines(file.path(work,"matrix.tsv")),"unchanged input")
  writeLines("conflicting input",file.path(work,"matrix.tsv"))
  expect_error(stage(),"staged input differs")
  expect_identical(readLines(file.path(work,"matrix.tsv")),"conflicting input")
})


test_that("the active sample map changes only its explanatory note", {
  root <- system.file("extdata", "quick-start", package="lisaR")
  old <- read.delim(file.path(root,"example_category_map.tsv"))
  new <- read.delim(file.path(root,"example_category_map_v1_1.tsv"))
  expect_identical(old[setdiff(names(old),"notes")], new[setdiff(names(new),"notes")])
  expect_false(any(grepl("teaching", unlist(new), ignore.case=TRUE)))
  expect_identical(lisaR:::lisa_sha256_file(file.path(root,"example_category_map_v1_1.tsv")), "d79a663529678b0728176382d75d853bb94f31bb47ada0f04e7564ce37433b0e")
})


test_that("current sample resources contain no retired descriptive labels", {
  root <- system.file("extdata", "quick-start", package="lisaR")
  cfg <- yaml::read_yaml(system.file("examples/quick-start/study.yml",package="lisaR"))
  expect_identical(cfg$pipeline$dictionary_resource,"lisa_dictionary_quickstart@1.1.0")
  expect_identical(cfg$pipeline$category_map_resource,"lisa_quickstart_category_map@1.1.0")
  expect_identical(cfg$pipeline$term2gene_resource,"lisa_quickstart_term2gene@2.1.0")
  for (name in c("lisa_dictionary_quickstart_v1_1.tsv","example_category_map_v1_1.tsv","example_term2gene_v2_1.tsv")) {
    expect_false(any(grepl("teaching",readLines(file.path(root,name)),ignore.case=TRUE)),info=name)
  }
  old <- read.delim(file.path(root,"example_term2gene_v2_0.tsv"))
  new <- read.delim(file.path(root,"example_term2gene_v2_1.tsv"))
  expect_identical(old[setdiff(names(old),"gs_exact_source")],new[setdiff(names(new),"gs_exact_source")])
  expect_identical(lisaR:::lisa_sha256_file(file.path(root,"example_term2gene_v2_1.tsv")),"1e1434af1ec53c931004eae23f2e8cb2c163314b2424789bdf9fcc41164699ec")
})
