test_that("theme paths match existing assets to not-yet-written report pages", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), env)
  root <- tempfile("theme-")
  dir.create(file.path(root, "report_assets"), recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  css <- file.path(root, "report_assets", "lisa_pastel.css")
  writeLines("body { color: black; }", css)
  canonical <- normalizePath(root, winslash = "/", mustWork = TRUE)
  # On macOS tempfile() may spell this /var while canonical spells it
  # /private/var. Also covers the Windows short/long-name alias boundary.
  page <- file.path(root, "report_pages", "evidence", "A", "GOCC", "index.html")
  expect_false(file.exists(page))
  expect_identical(env$report_relative_path(page, canonical),
    "report_pages/evidence/A/GOCC/index.html")
  expect_identical(env$report_relative_path(css, dirname(page)),
    "../../../../report_assets/lisa_pastel.css")
  html <- env$report_apply_theme("<html><head></head><body></body></html>", page, canonical)
  expect_match(html, 'data-style-scope="analysis"', fixed = TRUE)
  expect_false(grepl(canonical, html, fixed = TRUE))
  if (!identical(.Platform$OS.type, "windows")) {
    future <- paste0("/", basename(tempfile("lisa-absent-")), "/planned/index.html")
    expect_false(file.exists(dirname(dirname(future))))
    expect_identical(env$report_path_components(future)$root, "/")
    # Use the actual filesystem root; /tmp is a macOS alias for /private/tmp.
    expect_identical(env$report_relative_path(future, "/"),
      sub("^/", "", future))
    foreign <- "//__lisar_absent_unc__/share/report/page.html"
    expect_identical(env$report_path_components(foreign)$root,
      "//__lisar_absent_unc__/share")
    expect_identical(env$report_relative_path(
      foreign, "//__lisar_absent_unc__/share/assets"), "../report/page.html")
    expect_error(env$report_relative_path(foreign, "/"),
      "across filesystem roots", fixed = TRUE)
    expect_error(env$report_relative_path(
      foreign, "//__lisar_absent_unc__/other/assets"),
      "across filesystem roots", fixed = TRUE)
  }
})

# Reproduced on Linux before the correction; this is not native macOS/Windows QA.
test_that("planned report pages resolve existing directory aliases without escaping the portable root", {
  env <- new.env(parent=globalenv())
  sys.source(system.file("scripts/build_LISA_report.R",package="lisaR"),env)
  root <- withr::local_tempdir("planned-report-alias-")
  real <- file.path(root,"real-report")
  dir.create(file.path(real,"report_assets"),recursive=TRUE)
  asset <- file.path(real,"report_assets/theme.css")
  writeLines("body{}",asset)
  alias <- file.path(root,"alias-report")
  linked <- suppressWarnings(file.symlink(real,alias))
  if (!isTRUE(linked)) skip("directory symlink capability unavailable")
  page <- file.path(alias,"report_pages/evidence/index.html")
  before <- tools::md5sum(asset)
  expect_false(file.exists(page))
  expect_identical(env$report_relative_path(page,real),"report_pages/evidence/index.html")
  expect_identical(env$report_relative_path(asset,dirname(page)),"../../report_assets/theme.css")
  expect_identical(tools::md5sum(asset),before)
  expect_false(file.exists(dirname(page)))
})
