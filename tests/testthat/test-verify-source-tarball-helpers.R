test_that("retired-layout scan exempts NEWS.md and THIRD_PARTY_NOTICES.md but not code", {
  candidates <- c(
    file.path("..", "verify_source_tarball_helpers.R"),
    testthat::test_path("..", "verify_source_tarball_helpers.R")
  )
  helper <- candidates[file.exists(candidates)][1L]
  skip_if(is.na(helper) || !nzchar(helper), "verify_source_tarball_helpers.R is excluded from the source tarball")
  source(normalizePath(helper, mustWork = TRUE), local = TRUE)

  retired_package <- paste0("cow", "plot")
  fixture <- tempfile("lisa-retired-layout-")
  dir.create(file.path(fixture, "R"), recursive = TRUE)
  writeLines(
    sprintf("- inst/THIRD_PARTY_NOTICES.md records the transitive-only status of `%s`.", retired_package),
    file.path(fixture, "NEWS.md")
  )
  writeLines(
    sprintf("`%s` (GPL-2 only) reaches an installation only as a transitive dependency of `fgsea`.", retired_package),
    file.path(fixture, "THIRD_PARTY_NOTICES.md")
  )
  writeLines(
    sprintf("layout <- %s::plot_grid(a, b)", retired_package),
    file.path(fixture, "R", "regressed_layout.R")
  )

  files <- list.files(fixture, recursive = TRUE, full.names = TRUE)
  hits <- lisa_retired_layout_hits(files, retired_layout_package = retired_package)

  expect_true(any(grepl("regressed_layout[.]R$", hits)))
  expect_false(any(grepl("NEWS[.]md$", hits)))
  expect_false(any(grepl("THIRD_PARTY_NOTICES[.]md$", hits)))
})
