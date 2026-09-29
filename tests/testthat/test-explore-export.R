
test_that("export copies the whole STANDARD plus available figures, rendering nothing", {
  built <- explore_shared_bundle()
  shared <- built$shared
  destination <- built$dir
  bundle <- built$bundle

  # Export assembles; it never generates.
  expect_identical(length(built$trace$post), 0L)
  expect_identical(built$trace$pathway, 0L)
  expect_identical(bundle$extensions, 1L)

  # The complete STANDARD is present. No category disappears for not having been
  # visited. Pages the export deliberately annotates are compared separately
  # below, so this check stays exact for everything else.
  standard <- lisaR:::lisa_explore_relative_paths(shared$run)
  expect_gt(length(standard), 1000L)
  expect_identical(bundle$standard_files, length(standard))
  missing <- standard[!file.exists(file.path(destination, standard))]
  expect_identical(missing, character(0))

  untouched <- setdiff(standard, bundle$attached_pages)
  expect_identical(
    unname(tools::md5sum(file.path(destination, untouched))),
    unname(tools::md5sum(file.path(shared$run, untouched))))

  # The report entry point and its navigation pages came across.
  expect_true(file.exists(file.path(destination, "report_index.html")))
  expect_gt(length(list.files(file.path(destination, "report_pages"),
                              recursive = TRUE)), 0L)

  # The source run is untouched by exporting.
  expect_identical(explore_tree_digest(shared$run), shared$before)
})

test_that("the bundle records where each figure attaches in the report", {
  built <- explore_shared_bundle()
  destination <- built$dir

  presentation <- read.delim(file.path(destination, "report_presentation.tsv"),
                             sep = "\t", stringsAsFactors = FALSE)
  expect_identical(nrow(presentation), 1L)
  expect_identical(presentation$product, "volcano")
  expect_identical(presentation$category_id, "SYN_SIGNAL")

  # The figure attaches at its own category route in the existing report; the
  # export does not invent a separate gallery.
  # the route is the category evidence sheet the reader actually opens,
  # carrying the category as the query route the report's own navigator emits --
  # not the per-collection navigator page the figure used to be stacked under.
  expect_identical(presentation$attach_route, explore_attach_route())
  expect_identical(presentation$attach_anchor, "category-SYN_SIGNAL")
  expect_identical(presentation$attach_query, "category=SYN_SIGNAL")
  expect_identical(presentation$attach_surface, "evidence")
  expect_true(file.exists(file.path(built$shared$run, presentation$attach_route)))

  # All four files of the one figure are addressable from the bundle.
  for (column in c("png", "pdf", "source_data", "recipe")) {
    expect_true(nzchar(presentation[[column]]), info = column)
    expect_true(file.exists(file.path(destination, presentation[[column]])),
                info = column)
  }

  # The retired "KEGG gene sets" gallery is absent and the legitimate
  # "member gene sets" content is untouched.
  expect_false(any(grepl("kegg_gene_sets", presentation$png)))
  expect_identical(built$bundle$extensions, 1L)
})

# --- finding 1: the figure must be VISIBLE in the exported report ------------

test_that("the exported category page itself shows the figure and its downloads", {
  built <- explore_shared_bundle()
  destination <- built$dir
  route <- explore_attach_route()

  expect_identical(built$bundle$attached_pages, route)
  page <- file.path(destination, route)
  html <- paste(readLines(page, warn = FALSE), collapse = "\n")

  # A route table alone is not enough: the page must carry the figure.
  expect_match(html, 'data-lisa-explore-attachment="true"', fixed = TRUE)
  expect_match(html, "SYN_SIGNAL", fixed = TRUE)

  # Exactly one block, so repeated exports cannot stack duplicates.
  expect_identical(
    length(gregexpr('data-lisa-explore-attachment="true"', html,
                    fixed = TRUE)[[1L]]), 1L)

  # The image and all four downloads resolve relative to the page.
  img <- regmatches(html, gregexpr(
    '(?<=<img class="lisa-explore-attachment-image" src=")[^"]+', html,
    perl = TRUE))[[1L]]
  downloads <- regmatches(html, gregexpr('(?<=<a download href=")[^"]+', html,
                                         perl = TRUE))[[1L]]
  expect_identical(length(img), 1L)
  expect_identical(length(downloads), 4L)
  for (href in c(img, downloads)) {
    expect_true(file.exists(file.path(dirname(page), href)), info = href)
  }
  expect_identical(sum(grepl("[.]png$", downloads)), 1L)
  expect_identical(sum(grepl("[.]pdf$", downloads)), 1L)
  expect_identical(sum(grepl("_source[.]tsv$", downloads)), 1L)
  expect_identical(sum(grepl("_recipe[.]R$", downloads)), 1L)

  # Styled by a contained, relatively linked stylesheet.
  expect_true(file.exists(file.path(destination, "report_assets/lisa-explore.css")))
  expect_match(html, 'data-lisa-explore-style="true"', fixed = TRUE)

  # No alternative gallery replaced the source report.
  expect_identical(length(list.files(destination, pattern = "gallery",
                                     recursive = TRUE, ignore.case = TRUE)), 0L)

  # The attachment is written before hashing, so the manifest describes the page
  # the user actually opens.
  manifest <- read.delim(file.path(destination, "explore_export_manifest.tsv"),
                         sep = "\t", stringsAsFactors = FALSE)
  row <- manifest[manifest$path == route, , drop = FALSE]
  expect_identical(nrow(row), 1L)
  expect_identical(row$sha256[[1L]], lisaR:::lisa_sha256_file(page))

  # And the SOURCE page was not modified to achieve any of this.
  expect_identical(
    unname(tools::md5sum(file.path(built$shared$run, route))),
    unname(tools::md5sum(file.path(built$shared$run, route))))
  expect_identical(explore_tree_digest(built$shared$run), built$shared$before)
})

test_that("attaching an already attached page is a no-op", {
  built <- explore_shared_bundle()
  route <- explore_attach_route()
  html <- paste(readLines(file.path(built$dir, route), warn = FALSE),
                collapse = "\n")
  expect_identical(
    lisaR:::lisa_explore_presentation_attach_html(html, built$bundle$presentation,
                                                  route),
    html)
})

test_that("the bundle opens from a new path with R and Shiny closed", {
  built <- explore_shared_bundle()

  # Move a copy to a directory whose name contains spaces, as a user would.
  moved_parent <- file.path(tempdir(), "moved bundle dir")
  unlink(moved_parent, recursive = TRUE, force = TRUE)
  on.exit(unlink(moved_parent, recursive = TRUE, force = TRUE), add = TRUE)
  moved <- file.path(moved_parent, "lisa report copy")
  explore_copy_bundle(moved)

  manifest <- read.delim(file.path(moved, "explore_export_manifest.tsv"),
                         sep = "\t", stringsAsFactors = FALSE)
  expect_gt(nrow(manifest), 1000L)
  expect_true(all(file.exists(file.path(moved, manifest$path))))

  # Content survives relocation: re-hash a sample across the bundle, including
  # the figure files, rather than trusting the move.
  sample_rows <- unique(c(
    which(manifest$path == "report_index.html"),
    grep("[.](png|pdf)$", manifest$path)[seq_len(min(20L,
      length(grep("[.](png|pdf)$", manifest$path))))],
    grep("^explore_extensions/", manifest$path)))
  for (row in sample_rows) {
    path <- file.path(moved, manifest$path[[row]])
    expect_identical(lisaR:::lisa_sha256_file(path), manifest$sha256[[row]],
                     info = manifest$path[[row]])
  }

  # The attached figure still resolves after the move. This is the check that
  # makes relocation meaningful rather than a file-count.
  route <- explore_attach_route()
  page <- file.path(moved, route)
  html <- paste(readLines(page, warn = FALSE), collapse = "\n")
  links <- c(
    regmatches(html, gregexpr(
      '(?<=<img class="lisa-explore-attachment-image" src=")[^"]+', html,
      perl = TRUE))[[1L]],
    regmatches(html, gregexpr('(?<=<a download href=")[^"]+', html,
                              perl = TRUE))[[1L]])
  expect_identical(length(links), 5L)
  for (href in links) {
    expect_true(file.exists(file.path(dirname(page), href)), info = href)
  }

  # Nothing in the bundle links back to a working directory or a local server.
  receipt <- read.delim(file.path(moved, "explore_export_receipt.tsv"),
                        sep = "\t", stringsAsFactors = FALSE)
  expect_identical(as.integer(receipt$generators_invoked), 0L)
  expect_identical(receipt$complete_missing, "false")
  expect_identical(as.integer(receipt$html_files_with_absolute_references), 0L)
  expect_identical(as.integer(receipt$symlinks), 0L)
  # read.delim coerces the literal "TRUE" to logical, so compare as character.
  expect_identical(as.character(receipt$relocatable), "TRUE")
  expect_identical(as.integer(receipt$attached_pages), 1L)
  expect_identical(length(list.files(moved, recursive = TRUE,
                                     all.files = TRUE, pattern = "[.]Rproj$")),
                   0L)
})

test_that("export refuses to write inside the source run or over an existing bundle", {
  built <- explore_shared_bundle()
  expect_error(lisa_explore_export(built$shared$ws,
                                   file.path(built$shared$run, "bundle")),
               "LISA-EXPLORE-030")
  # The already-built shared bundle is an existing destination. This refusal
  # happens before anything is copied, so it costs nothing.
  expect_error(lisa_explore_export(built$shared$ws, built$dir),
               "LISA-EXPLORE-031")
})

test_that("an export taken while nothing is available still carries the STANDARD", {
  run <- explore_fixture_run()
  root <- file.path(tempdir(), "lisa-explore-empty-ws")
  unlink(root, recursive = TRUE, force = TRUE)
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  ws <- lisa_explore_open(run, root)

  destination <- file.path(tempdir(), "lisa-explore-bundle-empty")
  unlink(destination, recursive = TRUE, force = TRUE)
  on.exit(unlink(destination, recursive = TRUE, force = TRUE), add = TRUE)

  traced <- explore_with_generator_trace(lisa_explore_export(ws, destination))
  expect_identical(length(traced$post), 0L)
  expect_identical(traced$value$extensions, 0L)
  expect_gt(traced$value$standard_files, 1000L)
  expect_true(file.exists(file.path(destination, "report_index.html")))
  # An empty presentation table, not a missing one.
  expect_true(file.exists(file.path(destination, "report_presentation.tsv")))
  # Nothing to attach means no page was touched and no stylesheet was added.
  expect_identical(traced$value$attached_pages, character(0))
  expect_false(file.exists(file.path(destination, "report_assets/lisa-explore.css")))
})
