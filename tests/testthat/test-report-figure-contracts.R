test_that("report relative paths do not depend on GNU realpath", {
  script <- system.file("scripts", "build_LISA_report.R", package = "lisaR")
  code <- paste(readLines(script, warn = FALSE), collapse = "\n")
  expect_false(grepl("realpath", code, fixed = TRUE))

  env <- new.env(parent = globalenv())
  sys.source(script, envir = env)
  project <- tempfile("portable-report-paths-")
  page_dir <- file.path(project, "report_pages")
  figure_dir <- file.path(project, "outputs")
  dir.create(page_dir, recursive = TRUE)
  dir.create(figure_dir, recursive = TRUE)
  page <- file.path(page_dir, "single_de.html")
  figure <- file.path(figure_dir, "figure.png")
  file.create(page, figure)

  empty_path <- tempfile("no-realpath-bin-")
  dir.create(empty_path)
  old_path <- Sys.getenv("PATH")
  on.exit(Sys.setenv(PATH = old_path), add = TRUE)
  Sys.setenv(PATH = empty_path)
  expect_identical(env$rel_path(figure, page), "../outputs/figure.png")
  expect_match(env$figure_id_for_path(figure, project), "^F_[0-9a-f]{16}$")
})

test_that("report path handling rejects Windows device namespaces", {
  env <- new.env(parent = globalenv())
  sys.source(
    system.file("scripts", "build_LISA_report.R", package = "lisaR"),
    envir = env
  )
  project <- tempfile("report-device-path-")
  dir.create(project)
  expect_false(env$report_is_absolute_path("C: reactive protein"))
  expect_true(env$report_is_absolute_path("///srv/private/input.tsv"))
  expect_error(env$report_is_absolute_path("//server"), "incomplete UNC")
  expect_error(
    env$report_portable_text(
      "C:relative.tsv", project, force_path = TRUE
    ),
    "drive-relative"
  )
  expect_error(
    env$report_portable_path("\\\\?\\C:\\private\\input.tsv", project),
    "device-path"
  )
  expect_error(
    env$report_portable_path("\\\\.\\C:\\private\\input.tsv", project),
    "device-path"
  )
  expect_error(
    env$report_portable_text("file://?/C:/private/input.tsv", project),
    "device-path"
  )
  expect_error(
    env$report_path_components("\\\\?\\C:\\private\\input.tsv"),
    "device-path"
  )
})

test_that("report navigation links only to QC sections that exist", {
  env <- new.env(parent = globalenv())
  sys.source(
    system.file("scripts", "build_LISA_report.R", package = "lisaR"),
    envir = env
  )
  project <- tempfile("report-qc-navigation-")
  dir.create(file.path(project, "config"), recursive = TRUE)
  page <- file.path(project, "report_index.html")

  without_post_lisa <- env$page_shell(
    "Fixture", "overview", "", project, page, root = TRUE
  )
  expect_false(grepl("qc-post-lisa-status-tsv", without_post_lisa, fixed = TRUE))

  utils::write.table(
    data.frame(stage = "fixture", status = "completed"),
    file.path(project, "post_lisa_status.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  with_post_lisa <- env$page_shell(
    "Fixture", "overview", "", project, page, root = TRUE
  )
  expect_true(grepl("qc-post-lisa-status-tsv", with_post_lisa, fixed = TRUE))
})

test_that("report source resolver accepts exact figure data and rejects generic indexes", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  d <- file.path(tempdir(), ".", "figure_source_resolver")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  png <- file.path(d, "01_IMMUNE_category_gene_card.png")
  source <- file.path(d, "01_IMMUNE_category_gene_card_source.tsv")
  index <- file.path(d, "cards_index.tsv")
  file.create(png)
  utils::write.table(data.frame(figure_id = "card", figure_type = "gene_card", symbol = "A"),
    source, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(png_path = png), index, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_identical(
    env$source_data_for_figure(png),
    normalizePath(source, winslash = "/", mustWork = TRUE)
  )

  unlink(source)
  expect_identical(env$source_data_for_figure(png), "")

  candidate_root <- file.path(tempdir(), ".", "figure_source_candidate")
  plot_dir <- file.path(candidate_root, "plots")
  table_dir <- file.path(candidate_root, "lisa_tables")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  candidate_png <- file.path(plot_dir, "A_GSEA_lollipop.png")
  candidate_tsv <- file.path(table_dir, "A_GSEA_category_summary.tsv")
  file.create(candidate_png, candidate_tsv)
  expect_identical(
    env$source_data_for_figure(candidate_png),
    normalizePath(candidate_tsv, winslash = "/", mustWork = TRUE)
  )
})

test_that("KEGG report state rebases stale index paths to the local painter directory", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  project <- file.path(tempdir(), "kegg_report_state")
  painter <- file.path(project, "outputs", "gene_level", "single_de", "C1", "collection_GOBP-C2", "kegg_painter")
  dir.create(painter, recursive = TRUE, showWarnings = FALSE)
  map <- file.path(painter, "01_mmu00000_fixture_painted.png")
  file.create(map)
  utils::write.table(data.frame(stage = "single_de_kegg_pathway_painter", analysis_id = "C1",
      contrast_id = "", collection = "GOBP-C2", status = "completed", message = ""),
    file.path(project, "post_lisa_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(output_png = "/obsolete/cluster/path/01_mmu00000_fixture_painted.png"),
    file.path(painter, "C1_GOBP-C2_kegg_pathway_painter_index.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(kegg_id = "mmu00000", x = 1, y = 1),
    file.path(painter, "C1_GOBP-C2_kegg_pathway_painter_nodes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  state <- env$kegg_report_state(project, "single", "C1", "GOBP-C2", painter, map)
  expect_identical(state$state, "rendered")
})

test_that("category member plots are recursively discoverable for report integration", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  collection <- file.path(tempdir(), "collection_GOBP-C2")
  member_dir <- file.path(collection, "plots", "C1_GSEA_category_pathways")
  dir.create(member_dir, recursive = TRUE, showWarnings = FALSE)
  member <- file.path(member_dir, "01_IMMUNE.png")
  summary <- file.path(collection, "plots", "C1_GSEA_lollipop.png")
  file.create(member, summary)
  images <- env$pngs(collection)
  expect_true(member %in% images)
  expect_equal(sum(grepl("GSEA_category_pathways", images)), 1L)
})

test_that("figure-specific source contracts are gzip-compressed and readable", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  project <- file.path(tempdir(), "compressed_figure_contract")
  figure_dir <- file.path(project, "outputs", "gene_level", "single_de", "C1")
  page_dir <- file.path(project, "report_pages")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(page_dir, recursive = TRUE, showWarnings = FALSE)
  png <- file.path(figure_dir, "01_IMMUNE_category_gene_card.png")
  source <- file.path(figure_dir, "01_IMMUNE_category_gene_card_source.tsv")
  page <- file.path(page_dir, "single_de.html")
  file.create(png, page)
  utils::write.table(data.frame(
    symbol = c("A", "B"), log2FC = c(1, -1), padj = c(0.01, 0.02),
    gene_contribution_score = c(2, 1), plotted_order = 1:2
  ), source, sep = "\t", quote = FALSE, row.names = FALSE)

  old_package_dir <- get0("report_package_dir", envir = env, inherits = FALSE)
  assign("report_package_dir", system.file(package = "lisaR"), envir = env)
  on.exit(assign("report_package_dir", old_package_dir, envir = env), add = TRUE)
  contract <- env$figure_source_contract(png, page)
  source_path <- file.path(dirname(page), contract$source_tsv)

  expect_match(source_path, "[.]tsv[.]gz$")
  expect_true(file.exists(source_path))
  expect_identical(readBin(source_path, "raw", n = 2), as.raw(c(0x1f, 0x8b)))
  x <- env$read_tsv(source_path)
  expect_equal(nrow(x), 2L)
  expect_identical(unique(x$figure_id), contract$figure_id)
  expect_silent(env$validate_figure_source_contracts(project))
})

test_that("portable report copies compress downloadable TSV tables", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  project <- file.path(tempdir(), "compressed_report_files")
  page_dir <- file.path(project, "report_pages")
  dir.create(page_dir, recursive = TRUE, showWarnings = FALSE)
  source <- file.path(project, "large_table.tsv")
  page <- file.path(page_dir, "downloads.html")
  file.create(page)
  utils::write.table(data.frame(symbol = c("A", "B"), value = c(1, 2)),
    source, sep = "\t", quote = FALSE, row.names = FALSE)

  href <- env$short_file_copy(source, page)
  portable <- file.path(dirname(page), href)
  expect_match(portable, "[.]tsv[.]gz$")
  expect_identical(readBin(portable, "raw", n = 2), as.raw(c(0x1f, 0x8b)))
  expect_equal(nrow(env$read_tsv(portable)), 2L)

  empty_source <- file.path(project, "empty_table.tsv")
  file.create(empty_source)
  empty_href <- env$short_file_copy(empty_source, page)
  empty_portable <- normalizePath(
    file.path(dirname(page), empty_href), mustWork = TRUE
  )
  expect_identical(readBin(empty_portable, "raw", n = 2), as.raw(c(0x1f, 0x8b)))
  empty_connection <- gzfile(empty_portable, open = "rt")
  empty_lines <- readLines(empty_connection, warn = FALSE)
  close(empty_connection)
  expect_length(empty_lines, 0L)
})

test_that("portable TSV copies hide Unix and Windows paths without changing science", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  project <- tempfile("private-run-")
  page_dir <- file.path(project, "report_pages")
  dir.create(file.path(project, "inputs"), recursive = TRUE)
  dir.create(page_dir, recursive = TRUE)
  page <- file.path(page_dir, "downloads.html")
  file.create(page)
  source <- file.path(project, "status.tsv")
  internal <- file.path(project, "inputs", "study.tsv")
  private_unix <- "/srv/private-study/receipts/external.tsv"
  private_windows <- "C:\\Users\\fixture-user\\private-study\\windows.csv"
  private_file_uri <- "file:///home/fixture-user/private-study/uri.tsv"
  private_windows_uri <- "file:///C:/Users/fixture-user/private-study/uri-windows.csv"
  original <- data.frame(
    symbol = c("TP53", "MYC", "STAT1", "JUN", "FOS"),
    score = c("1.2300", "-2.500", "0.000", "0.750", "-0.250"),
    input_path = c(
      internal, private_unix, private_windows, private_file_uri,
      private_windows_uri
    ),
    role = c(
      "internal_input", rep("external_receipt", 4L)
    ),
    message = c(
      "scientific row retained",
      paste("read", private_unix),
      paste("read", private_windows),
      paste("read", private_file_uri),
      paste("read", private_windows_uri)
    ),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    original, source, sep = "\t", quote = FALSE, row.names = FALSE
  )
  canonical_before <- readBin(source, "raw", n = file.info(source)$size)

  assign("report_project_dir", project, envir = env)
  html <- env$table_html(original, NULL, "privacy", source, 20)
  href <- env$short_file_copy(source, page)
  portable <- normalizePath(file.path(dirname(page), href), mustWork = TRUE)
  copied <- env$report_read_delimited_text(portable)

  expect_false(grepl(project, html, fixed = TRUE))
  expect_false(grepl("/srv/private-study", html, fixed = TRUE))
  expect_false(grepl("C:\\Users\\fixture-user", html, fixed = TRUE))
  expect_false(grepl("file:", html, fixed = TRUE))
  portable_text <- env$report_read_portable_text_file(portable)
  expect_false(grepl(project, portable_text, fixed = TRUE))
  expect_false(grepl("/srv/private-study", portable_text, fixed = TRUE))
  expect_false(grepl("C:/Users/fixture-user", portable_text, fixed = TRUE))
  expect_false(grepl("file:", portable_text, fixed = TRUE))
  expect_identical(copied$symbol, original$symbol)
  expect_identical(copied$score, original$score)
  expect_identical(copied$role, original$role)
  expect_identical(
    copied$input_path,
    c(
      "inputs/study.tsv", "external.tsv", "windows.csv", "uri.tsv",
      "uri-windows.csv"
    )
  )
  expect_match(copied$message[[2L]], "external.tsv", fixed = TRUE)
  expect_match(copied$message[[3L]], "windows.csv", fixed = TRUE)
  expect_match(copied$message[[4L]], "uri.tsv", fixed = TRUE)
  expect_match(copied$message[[5L]], "uri-windows.csv", fixed = TRUE)
  expect_identical(
    readBin(source, "raw", n = file.info(source)$size), canonical_before
  )

  csv_source <- file.path(project, "status.csv")
  utils::write.csv(original, csv_source, quote = TRUE, row.names = FALSE)
  csv_href <- env$short_file_copy(csv_source, page)
  csv_portable <- normalizePath(
    file.path(dirname(page), csv_href), mustWork = TRUE
  )
  csv_copied <- env$report_read_delimited_text(csv_portable)
  expect_identical(csv_copied$symbol, original$symbol)
  expect_identical(csv_copied$score, original$score)
  expect_identical(csv_copied$role, original$role)
  expect_identical(
    csv_copied$input_path,
    c(
      "inputs/study.tsv", "external.tsv", "windows.csv", "uri.tsv",
      "uri-windows.csv"
    )
  )
  expect_false(grepl(
    "C:/Users/fixture-user|/srv/private-study|file:",
    env$report_read_portable_text_file(csv_portable), perl = TRUE
  ))
})

test_that("portable paths resolve an existing symlinked ancestor", {
  env <- new.env(parent = globalenv())
  sys.source(
    system.file("scripts", "build_LISA_report.R", package = "lisaR"),
    envir = env
  )
  base <- tempfile("portable-path-alias-")
  dir.create(base)
  on.exit(lisa_test_cleanup_path(base), add = TRUE)
  project <- file.path(base, "project")
  alias <- file.path(base, "alias")
  dir.create(file.path(project, "inputs"), recursive = TRUE)
  if (!isTRUE(suppressWarnings(file.symlink(project, alias)))) {
    skip("file.symlink() is unavailable in this test environment")
  }

  observed <- env$report_path_within_root(
    file.path(alias, "inputs", "missing.tsv"), project
  )
  expect_identical(observed, "inputs/missing.tsv")
})

test_that("shareable privacy gate rejects residual file URIs", {
  env <- new.env(parent = globalenv())
  sys.source(
    system.file("scripts", "build_LISA_report.R", package = "lisaR"),
    envir = env
  )
  project <- tempfile("file-uri-report-")
  dir.create(file.path(project, "report_pages"), recursive = TRUE)
  page <- file.path(project, "report_pages", "unsafe.html")
  writeLines('<html><a href="file:///home/fixture-user/private.tsv">x</a></html>', page)
  expect_error(
    env$report_assert_no_private_paths(project, page),
    "LISA-REPORT-PRIVACY-004"
  )
})

test_that("file URI sanitation preserves scientific profile labels and web URLs", {
  env <- new.env(parent = globalenv())
  sys.source(
    system.file("scripts", "build_LISA_report.R", package = "lisaR"),
    envir = env
  )
  project <- tempfile("file-uri-boundary-")
  dir.create(project)
  retained <- c(
    "profile:immune",
    "molecular_profile:MYC",
    "https://example.org/profile:data",
    "my+file:data",
    "urn:file:data"
  )
  expect_identical(env$report_portable_text(retained, project), retained)
  expect_identical(
    env$report_portable_text(
      "profile: file:///srv/private/a.tsv", project
    ),
    "profile: a.tsv"
  )
  expect_identical(
    env$report_portable_text(
      "source_file:///srv/private/a.tsv", project
    ),
    "source_a.tsv"
  )
})

test_that("shareable manifest is exact, relocatable and has one SELF exception", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  project <- tempfile("shareable-report-")
  relocate <- tempfile("relocated-report-")
  for (directory in c(
    "report_pages", "report_assets", "report_media", "report_files",
    "report_figure_data"
  )) {
    dir.create(file.path(project, directory), recursive = TRUE)
  }
  writeLines(
    '<html><a href="report_pages/page.html">page</a></html>',
    file.path(project, "report_index.html")
  )
  writeLines("body{}", file.path(project, "report_assets", "style.css"))
  writeBin(as.raw(c(1, 2, 3)), file.path(project, "report_media", "figure.png"))
  table_connection <- gzfile(
    file.path(project, "report_files", "table.tsv.gz"), "wt"
  )
  writeLines(c("symbol\tscore", "TP53\t1.2"), table_connection)
  close(table_connection)
  source_connection <- gzfile(
    file.path(project, "report_figure_data", "F_0000000000000001.tsv.gz"),
    "wt"
  )
  writeLines(c("symbol\tfigure_id", "TP53\tF_0000000000000001"), source_connection)
  close(source_connection)
  writeLines(
    "invisible(TRUE)",
    file.path(project, "report_figure_data", "F_0000000000000001_recipe.R")
  )
  writeLines(c(
    '<html><head><link href="../report_assets/style.css"></head><body>',
    '<img src="../report_media/figure.png">',
    '<a href="../report_files/table.tsv.gz">table</a>',
    '<a href="../report_figure_data/F_0000000000000001.tsv.gz">source</a>',
    '<a href="../report_figure_data/F_0000000000000001_recipe.R">recipe</a>',
    '<a href="../shareable_report_manifest.tsv">manifest</a>',
    '</body></html>'
  ), file.path(project, "report_pages", "page.html"))

  env$write_shareable_report_manifest(project)
  manifest <- env$validate_shareable_report_manifest(project)
  self <- manifest$relative_path == "shareable_report_manifest.tsv"
  expect_equal(sum(self), 1L)
  expect_identical(manifest$sha256[self], "SELF")
  expect_true(all(grepl("^[0-9a-f]{64}$", manifest$sha256[!self])))

  dir.create(relocate, recursive = TRUE)
  for (relative in manifest$relative_path) {
    target <- file.path(relocate, relative)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    expect_true(file.copy(file.path(project, relative), target))
  }
  expect_silent(env$validate_shareable_report_manifest(relocate))
})

test_that("shareable report discovery rejects lexical file and directory links", {
  env <- new.env(parent = globalenv())
  sys.source(
    system.file("scripts", "build_LISA_report.R", package = "lisaR"),
    envir = env
  )
  project <- tempfile("shareable-link-report-")
  pages <- file.path(project, "report_pages")
  dir.create(pages, recursive = TRUE)
  on.exit(lisa_test_cleanup_path(project), add = TRUE)

  target <- file.path(pages, "real.html")
  link <- file.path(pages, "alias.html")
  writeLines("<html></html>", target)
  if (!isTRUE(suppressWarnings(file.symlink(target, link)))) {
    skip("file.symlink() is unavailable in this test environment")
  }
  on.exit(if (lisaR:::lisa_path_is_link(link)) {
    lisaR:::lisa_link_delete(link)
  }, add = TRUE)
  expect_true(lisaR:::lisa_path_is_link(link))
  expect_error(env$report_shareable_files(project), "LISA-REPORT-SHARE-001")
  lisaR:::lisa_link_delete(link)

  lisa_test_cleanup_path(pages)
  real_pages <- file.path(project, "real-pages")
  dir.create(real_pages)
  writeLines("<html></html>", file.path(real_pages, "page.html"))
  if (!isTRUE(suppressWarnings(file.symlink(real_pages, pages)))) {
    skip("directory symbolic links are unavailable in this test environment")
  }
  on.exit(if (lisaR:::lisa_path_is_link(pages)) {
    lisaR:::lisa_link_delete(pages)
  }, add = TRUE)
  expect_true(lisaR:::lisa_path_is_link(pages))
  expect_error(env$report_shareable_files(project), "LISA-REPORT-SHARE-001")
})

test_that("shareable manifest write and read reject a lexical link", {
  env <- new.env(parent = globalenv())
  sys.source(
    system.file("scripts", "build_LISA_report.R", package = "lisaR"),
    envir = env
  )
  project <- tempfile("shareable-manifest-link-")
  dir.create(project)
  on.exit(lisa_test_cleanup_path(project), add = TRUE)
  writeLines("<html></html>", file.path(project, "report_index.html"))
  outside <- tempfile("shareable-manifest-sentinel-")
  writeLines("outside-original", outside, useBytes = TRUE)
  on.exit(unlink(outside, force = TRUE), add = TRUE)
  outside_before <- readBin(outside, "raw", n = file.info(outside)$size)
  manifest_link <- file.path(project, "shareable_report_manifest.tsv")
  if (!isTRUE(suppressWarnings(file.symlink(outside, manifest_link)))) {
    skip("file.symlink() is unavailable in this test environment")
  }
  on.exit(if (lisaR:::lisa_path_is_link(manifest_link)) {
    lisaR:::lisa_link_delete(manifest_link)
  }, add = TRUE)
  expect_true(lisaR:::lisa_path_is_link(manifest_link))
  expect_error(
    env$write_shareable_report_manifest(project),
    "LISA-REPORT-SHARE-001"
  )
  expect_error(
    env$read_shareable_report_manifest(project),
    "LISA-REPORT-SHARE-001"
  )
  expect_identical(
    readBin(outside, "raw", n = file.info(outside)$size),
    outside_before
  )
})

test_that("shareable receipts reject a file mutated while hashing", {
  env <- new.env(parent = globalenv())
  sys.source(
    system.file("scripts", "build_LISA_report.R", package = "lisaR"),
    envir = env
  )
  project <- tempfile("shareable-receipt-swap-")
  dir.create(project)
  on.exit(lisa_test_cleanup_path(project), add = TRUE)
  component <- file.path(project, "report_index.html")
  writeLines("original", component, useBytes = TRUE)
  real_sha256 <- lisaR:::lisa_sha256_file
  mutated <- FALSE
  testthat::local_mocked_bindings(
    lisa_sha256_file = function(path) {
      digest <- real_sha256(path)
      if (!mutated && identical(
        lisaR:::lisa_path_key(path), lisaR:::lisa_path_key(component)
      )) {
        mutated <<- TRUE
        writeLines("changed-after-hash-with-a-different-size", path,
                   useBytes = TRUE)
      }
      digest
    },
    .package = "lisaR"
  )
  expect_error(
    env$report_shareable_file_receipt(component, project),
    "LISA-REPORT-SHARE-012"
  )
  expect_true(mutated)
})

test_that("shareable report discovery rejects non-regular entries", {
  skip_on_os("windows")
  mkfifo <- Sys.which("mkfifo")
  skip_if(!nzchar(mkfifo), "mkfifo is unavailable")
  env <- new.env(parent = globalenv())
  sys.source(
    system.file("scripts", "build_LISA_report.R", package = "lisaR"),
    envir = env
  )
  project <- tempfile("shareable-fifo-report-")
  pages <- file.path(project, "report_pages")
  dir.create(pages, recursive = TRUE)
  on.exit(lisa_test_cleanup_path(project), add = TRUE)
  fifo_path <- file.path(pages, "pipe.tsv")
  status <- system2(mkfifo, shQuote(fifo_path))
  expect_identical(status, 0L)
  on.exit(unlink(fifo_path, force = TRUE), add = TRUE)
  expect_identical(
    as.character(fs::file_info(fifo_path, follow = FALSE)$type[[1L]]),
    "FIFO"
  )
  expect_error(env$report_shareable_files(project), "LISA-REPORT-SHARE-011")
})

test_that("report pages use one media copy rather than a duplicated subtree", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  project <- tempfile("single-report-media-")
  page <- file.path(project, "report_pages", "single_de.html")
  figure <- file.path(project, "outputs", "figure.png")
  dir.create(dirname(page), recursive = TRUE)
  dir.create(dirname(figure), recursive = TRUE)
  file.create(page)
  writeBin(as.raw(c(1, 2, 3)), figure)

  href <- env$short_media_copy(figure, page)
  expect_match(href, "^[.][.]/report_media/", perl = TRUE)
  expect_false(dir.exists(file.path(project, "report_pages", "report_media")))
  expect_true(file.exists(file.path(dirname(page), href)))
})

test_that("figure recipe regenerates paired A/B contrast heatmaps", {
  skip_if_not_installed("ggplot2")
  project <- file.path(tempdir(), "paired_contrast_heatmap_recipe")
  dir.create(project, recursive = TRUE, showWarnings = FALSE)
  source <- file.path(project, "contrast_heatmap.tsv.gz")
  output <- file.path(project, "contrast_heatmap.png")
  con <- gzfile(source, open = "wt")
  utils::write.table(data.frame(
    symbol = c("A", "B"), log2FC_A = c(1.2, -0.7),
    log2FC_B = c(-0.5, 0.9), label_A = "Treatment A",
    label_B = "Treatment B", lfc_cap = 1.5,
    figure_type = "heatmap", selected_for_plot = TRUE
  ), con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  recipe <- system.file("scripts", "reproduce_lisa_figure.R", package = "lisaR")
  log <- system2(lisaR:::lisa_rscript_executable(),
    c(shQuote(recipe), shQuote(source), shQuote(output)),
    env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep)),
    stdout = TRUE, stderr = TRUE)
  status <- attr(log, "status")
  if (is.null(status)) status <- 0L
  expect_equal(status, 0L, info = paste(log, collapse = "\n"))
  expect_true(file.exists(output))
  expect_gt(file.info(output)$size, 0)
})

test_that("known unsupported contrast layers render explicit inapplicable state", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  html <- env$layer_section(
    "C5 - GOBP-C2 - Volcano overlays", character(), "contrasts.html",
    "Contrast-specific volcano overlays are not generated by the current contrast engine.",
    empty_state = "inapplicable"
  )
  expect_match(html, "not generated by this engine", fixed = TRUE)
  expect_false(grepl("missing expected layer", html, fixed = TRUE))
})

test_that("missing required layers preserve their actual family and context", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  html <- env$missing_layer_section(
    "Treatment dynamics - GOBP-C2 - Heatmaps", "contrasts.html"
  )
  expect_match(html, 'data-lisa-layer-family="Heatmaps"', fixed = TRUE)
  expect_match(html, 'data-lisa-layer-context="Treatment dynamics - GOBP-C2"', fixed = TRUE)
  expect_identical(env$report_expected_layer_diagnostics(html),
    "family=Heatmaps; context=Treatment dynamics - GOBP-C2")
})

test_that("report headings collapse only index-confirmed repeated contrast identifiers", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  expect_identical(
    env$short_title("response_dynamics_response_dynamics"),
    "response dynamics response dynamics"
  )
  expect_identical(env$short_title("single_single"), "single single")
  expect_identical(env$short_title("A_vs_B_profile"), "A vs B profile")
  index <- data.frame(
    contrast_id = c("response_dynamics", "A_vs_B"),
    output_id = c("response_dynamics", "profile"),
    stringsAsFactors = FALSE
  )
  expect_identical(
    env$contrast_output_id(index, "response_dynamics_response_dynamics"),
    "response_dynamics"
  )
  expect_identical(env$contrast_output_id(index, "A_vs_B_profile"), "profile")
  expect_identical(
    env$contrast_short_title(index, "response_dynamics_response_dynamics"),
    "response dynamics"
  )
  expect_identical(
    env$contrast_short_title(index, "A_vs_B_profile"),
    "A vs B profile"
  )
  expect_identical(
    env$contrast_short_title(index, "single_single"),
    "single single"
  )
})

test_that("HALLMARKS-only inapplicable category layers do not masquerade as missing", {
  env <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = env)
  html <- env$layer_section(
    "C1 - HALLMARKS - Gene sets by LISA category", character(), "single_de.html",
    "Not applicable to HALLMARKS.", empty_state = "inapplicable"
  )
  expect_match(html, "not generated by this engine", fixed = TRUE)
  expect_false(grepl("missing expected layer", html, fixed = TRUE))
})
