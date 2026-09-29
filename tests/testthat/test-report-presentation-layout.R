presentation_report_environment <- function() {
  e <- new.env(parent = globalenv())
  sys.source(system.file("scripts", "build_LISA_report.R", package = "lisaR"), envir = e)
  e$report_package_dir <- system.file(package = "lisaR")
  e
}

test_that("assembled report shares one portable header and navigation inventory", {
  e <- presentation_report_environment()
  root <- tempfile("report presentation with spaces ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "config"), recursive = TRUE)
  dir.create(file.path(root, "report_pages", "gene_evidence"), recursive = TRUE)
  lisaR:::write_lisa_tsv(data.frame(stage = c("single_de_lisa", "root_html_report", "single_de_gene_evidence")),
    file.path(root, "lisa_pipeline_plan.tsv"))
  utils::write.table(data.frame(analysis_id = c("A1", "A2"),
    label = c("Responders", "Nonresponders")), file.path(root, "config", "de_index.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  section <- e$report_navigation_section("A1 - GOCC - LISA categories", "LISA categories")
  contexts <- list(list(id = "A1", label = "Responders", kind = "analysis",
    collections = list(list(id = "GOCC", label = "GOCC", sections = list(section)))))
  html <- e$page_shell("Study", "single_de", "<h1>Analyses</h1>", root,
    file.path(root, "report_pages", "single_de.html"), navigation = contexts)
  expect_match(html, 'data-lisa-route="analyses" href="single_de.html" aria-current="page"', fixed = TRUE)
  expect_equal(length(regmatches(html, gregexpr('aria-current="page"', html, fixed = TRUE))[[1L]]), 1L)
  expect_match(html, 'id="lisa-report-navigation"', fixed = TRUE)
  expect_match(html, '"label":"Responders"', fixed = TRUE)
  expect_match(html, '"kind":"lisa-categories"', fixed = TRUE)
  expect_match(html, 'href="../report_index.html"', fixed = TRUE)
  expect_match(html, 'href="gene_evidence/index.html"', fixed = TRUE)
  expect_match(html, 'id="lisa-main"', fixed = TRUE)
  expect_equal(length(regmatches(html, gregexpr("<main ", html, fixed = TRUE))[[1L]]), 1L)
  expect_false(grepl(root, html, fixed = TRUE))
  expect_true(file.exists(file.path(root, "report_assets", "lisa-shell", "lisa_shell.css")))
})

test_that("navigation inventory preserves sections, context and safe JSON", {
  e <- presentation_report_environment()
  root <- tempfile("navigation inventory ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "report_pages"), recursive = TRUE)
  contexts <- list(list(id = "A_001", scientific_id = "A_001", label = "Response </script> & baseline",
    kind = "analysis", collections = list(list(id = "GOCC", label = "GOCC", sections = list(
      e$report_navigation_section("A - GOCC - LISA categories", "LISA categories"),
      e$report_navigation_section("A - GOCC - Heatmaps", "Heatmaps"))))))
  html <- e$report_navigation_json("single_de", contexts)
  expect_false(grepl("Response </script>", html, fixed = TRUE))
  text <- sub('^.*?id="lisa-report-navigation">', "", html)
  text <- sub('</script>$', "", text)
  decoded <- jsonlite::fromJSON(text, simplifyVector = FALSE)
  expect_identical(decoded$contexts[[1L]]$label, "Response </script> & baseline")
  expect_length(decoded$contexts[[1L]]$collections[[1L]]$sections, 2L)
  target <- e$report_write_navigation_inventory(root, contexts, list())
  persisted <- jsonlite::read_json(target, simplifyVector = FALSE)
  expect_identical(persisted$contexts[[1L]]$scientific_id, "A_001")
  expect_identical(persisted$contexts[[1L]]$collections[[1L]]$sections[[2L]]$href,
    "report_pages/single_de.html#a-gocc-heatmaps")
  expect_false(grepl(root, paste(readLines(target), collapse = ""), fixed = TRUE))
})

test_that("evidence static assets refresh atomically even when hard-linked from another project tree", {
  e <- presentation_report_environment()
  root <- tempfile("evidence asset refresh ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  origin <- tempfile("evidence asset refresh origin ")
  on.exit(unlink(origin, recursive = TRUE), add = TRUE)

  for (kind in c("category-evidence", "contrast-evidence")) {
    page_root <- if (identical(kind, "category-evidence")) "evidence" else "contrast_evidence"
    evidence_dir <- file.path(root, "report_pages", page_root, "A", "GOCC")
    asset_dir <- file.path(evidence_dir, "assets")
    dir.create(asset_dir, recursive = TRUE)
    writeLines("<html><body>Evidence</body></html>", file.path(evidence_dir, "index.html"))

    origin_dir <- file.path(origin, kind)
    dir.create(origin_dir, recursive = TRUE)

    stale <- list()
    for (extension in c("css", "js")) {
      origin_file <- file.path(origin_dir, paste0(kind, ".", extension))
      writeLines(paste0("/* stale ", kind, " ", extension, " inherited from another project tree */"), origin_file)
      destination <- file.path(asset_dir, paste0(kind, ".", extension))
      # Simulate the real incident: the evidence directory's static assets
      # were inherited via a hard link from another, older project tree, so
      # destination and origin_file share one inode before the refresh.
      expect_true(file.link(origin_file, destination))
      stale[[extension]] <- readLines(origin_file)
    }

    result <- e$report_refresh_evidence_static_assets(evidence_dir, kind)
    expect_identical(result, asset_dir)

    installed_root <- system.file(kind, package = "lisaR")
    for (extension in c("css", "js")) {
      destination <- file.path(asset_dir, paste0(kind, ".", extension))
      origin_file <- file.path(origin_dir, paste0(kind, ".", extension))
      installed <- readLines(file.path(installed_root, paste0("viewer.", extension)))
      # Refreshed destination now matches the currently installed package,
      # not the stale content that was there before.
      expect_identical(readLines(destination), installed)
      # The separate hard-linked origin copy is completely untouched --
      # the shared inode was replaced at the destination, never mutated.
      expect_identical(readLines(origin_file), stale[[extension]])
      if (.Platform$OS.type != "windows") {
        # `stat -c` is GNU coreutils only; BSD/macOS stat spells the same
        # request `-f`. `ls -di` is POSIX and prints the inode on both.
        inode <- function(path) {
          line <- suppressWarnings(
            system2("ls", c("-di", shQuote(path)), stdout = TRUE, stderr = FALSE)
          )
          sub("^\\s*([0-9]+)\\s.*$", "\\1", paste(line, collapse = " "))
        }
        # Fail loudly if the probe itself stops working: without this, a broken
        # probe would return the same non-answer for both paths and the
        # distinct-inode assertion below would pass vacuously.
        expect_match(inode(destination), "^[0-9]+$")
        expect_false(identical(inode(destination), inode(origin_file)))
      }
    }

    # Running the refresh again against the now-independent destination is a
    # stable, idempotent no-op: no error, identical resulting content.
    expect_no_error(e$report_refresh_evidence_static_assets(evidence_dir, kind))
    for (extension in c("css", "js")) {
      destination <- file.path(asset_dir, paste0(kind, ".", extension))
      installed <- readLines(file.path(installed_root, paste0("viewer.", extension)))
      expect_identical(readLines(destination), installed)
    }
  }
})

test_that("inline navigator rebases exact category links and has no nested document", {
  e <- presentation_report_environment()
  root <- tempfile("inline category report ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  nav <- file.path(root, "report_pages", "category_navigation", "A", "GOCC")
  evidence <- file.path(root, "report_pages", "evidence", "A", "GOCC", "index.html")
  dir.create(nav, recursive = TRUE)
  dir.create(dirname(evidence), recursive = TRUE)
  writeLines("<html><body>Evidence</body></html>", evidence)
  svg <- paste0('<svg xmlns="http://www.w3.org/2000/svg" role="group" aria-labelledby="navigation-title navigation-description">',
    '<title id="navigation-title">Category evidence</title><desc id="navigation-description">Scope A GOCC</desc>',
    '<a class="category-row" href="../../../evidence/A/GOCC/index.html?category=001&amp;analysis_id=A&amp;collection=GOCC&amp;tier=core" ',
    'target="_top" tabindex="0"><text>Gene set</text></a></svg>')
  writeLines(svg, file.path(nav, "category_navigation.svg"))
  page <- file.path(root, "report_pages", "single_de.html")
  inline <- e$report_inline_category_navigation(nav, page)
  expect_match(inline, 'href="evidence/A/GOCC/index.html?category=001&amp;analysis_id=A&amp;collection=GOCC&amp;tier=core"', fixed = TRUE)
  expect_match(inline, 'tabindex="0"', fixed = TRUE)
  expect_false(grepl('id="navigation-title"', inline, fixed = TRUE))
  expect_match(inline, 'id="nav-[0-9a-f]+-navigation-title"')
  cover <- e$report_category_evidence_section("A - GOCC - Category evidence", nav, evidence, page)
  expect_match(cover, "Category evidence</h2>", fixed = TRUE)
  expect_match(cover, "<svg", fixed = TRUE)
  expect_false(grepl("iframe|Accessible category table|<!doctype", cover))
  expect_match(cover, "Explore a category", fixed = TRUE)
  writeLines(gsub("../../../evidence/A/GOCC/index.html", "https://example.invalid/index.html", svg, fixed = TRUE),
    file.path(nav, "category_navigation.svg"))
  expect_error(e$report_inline_category_navigation(nav, page), "report-relative")
})

test_that("normal FULL assembly references canonical exact-ID products once", {
  e <- presentation_report_environment()
  root <- tempfile("canonical category products ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  evidence_dir <- file.path(root, "report_pages", "evidence", "A", "GOCC")
  product_dir <- file.path(root, "outputs", "single_de", "A", "collection_GOCC", "product")
  dir.create(file.path(evidence_dir, "tables"), recursive = TRUE)
  dir.create(product_dir, recursive = TRUE)
  writeLines('<script type="application/json" id="evidence-data">{"categories":[]}</script>',
    file.path(evidence_dir, "index.html"))
  utils::write.table(data.frame(key = c("analysis_id", "collection"), value = c("A", "GOCC")),
    file.path(evidence_dir, "tables", "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(category_id = "CAT_A"), file.path(evidence_dir, "tables", "categories.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(c("category_id", "CAT_A"), file.path(product_dir, "CAT_A_source.tsv"))
  writeBin(as.raw(c(137, 80, 78, 71)), file.path(product_dir, "CAT_A.png"))
  writeLines("not categorised", file.path(product_dir, "orphan.txt"))
  attached <- e$report_attach_category_products(file.path(evidence_dir, "index.html"),
    list(list(product = "volcano", label = "Volcano overlays", path = product_dir)),
    "A", "GOCC", project_dir = root)
  expect_identical(attached$category_products[[1L]]$category_id, "CAT_A")
  expect_identical(length(attached$unclassified), 1L)
  html <- paste(readLines(file.path(evidence_dir, "index.html"), warn = FALSE), collapse = "\n")
  payload <- sub('^.*id="evidence-data">', "", html)
  payload <- sub('</script>.*$', "", payload)
  decoded <- jsonlite::fromJSON(payload, simplifyVector = FALSE)
  hrefs <- vapply(decoded$category_products[[1L]]$assets, `[[`, character(1L), "href")
  expect_true(all(startsWith(hrefs, "../../../../outputs/")))
  expect_false(dir.exists(file.path(evidence_dir, "assets", "category_products")))
  expect_true(all(file.exists(file.path(evidence_dir, utils::URLdecode(hrefs)))))
  for (asset in decoded$category_products[[1L]]$assets) {
    expect_true(startsWith(asset$source_path, "outputs/"))
    expect_identical(asset$sha256, digest::digest(
      file = file.path(root, asset$source_path), algo = "sha256"
    ))
    expect_true(file.exists(file.path(product_dir, asset$source_name)))
  }
  e$write_shareable_report_manifest(root)
  manifest <- e$validate_shareable_report_manifest(root)
  expect_true(all(vapply(decoded$category_products[[1L]]$assets, `[[`,
    character(1L), "source_path") %in% manifest$relative_path))
  expect_false(any(grepl("assets/category_products", manifest$relative_path,
    fixed = TRUE)))
  expect_identical(length(unique(e$report_category_product_family(c(
    "SYN_A_leading_edge_gene_heatmap_matrix.tsv",
    "SYN_A_leading_edge_gene_heatmap_zscore.png"
  )))), 1L)
  # A page missing its exact owner key is not a valid attachment target: the
  # saved source/image pair remains unclassified instead of crashing or being
  # assigned from its filename.
  utils::write.table(data.frame(key = "collection", value = "GOCC"),
    file.path(evidence_dir, "tables", "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  missing_owner <- e$report_attach_category_products(file.path(evidence_dir, "index.html"),
    list(list(product = "volcano", label = "Volcano overlays", path = product_dir)),
    "A", "GOCC", project_dir = root)
  expect_length(missing_owner$category_products, 0L)
  expect_true(all(c(file.path(product_dir, "CAT_A_source.tsv"),
    file.path(product_dir, "CAT_A.png")) %in% missing_owner$unclassified))

  # The standard member-gene-set chart is already part of an evidence sheet;
  # a FULL inventory must not stage it as a duplicate saved-product card.
  utils::write.table(data.frame(key = c("analysis_id", "collection"), value = c("A", "GOCC")),
    file.path(evidence_dir, "tables", "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  excluded <- e$report_attach_category_products(file.path(evidence_dir, "index.html"),
    list(list(product = "member_gene_sets", label = "Member gene sets", path = product_dir)),
    "A", "GOCC", exclude_products = "member_gene_sets",
    project_dir = root)
  expect_length(excluded$category_products, 0L)
})

test_that("normal FULL contrast assembly binds canonical artifact JSON", {
  e <- presentation_report_environment()
  root <- tempfile("canonical contrast products ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  evidence_dir <- file.path(
    root, "report_pages", "contrast_evidence", "A_vs_B", "GOCC"
  )
  product_dir <- file.path(
    root, "artifacts", "contrasts", "A_vs_B", "GOCC",
    "contrast_heatmap"
  )
  dir.create(file.path(evidence_dir, "tables"), recursive = TRUE)
  dir.create(product_dir, recursive = TRUE)
  writeLines(
    '<script type="application/json" id="contrast-evidence-data">{"categories":[]}</script>',
    file.path(evidence_dir, "index.html")
  )
  utils::write.table(
    data.frame(
      key = c("contrast_id", "collection"),
      value = c("A_vs_B", "GOCC")
    ),
    file.path(evidence_dir, "tables", "metadata.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    data.frame(category_id = "CAT_A"),
    file.path(evidence_dir, "tables", "categories.tsv"), sep = "\t",
    quote = FALSE, row.names = FALSE
  )
  utils::write.table(
    data.frame(category_id = "CAT_A", value = 1),
    file.path(product_dir, "CAT_A_leading_edge_gene_heatmap_matrix.tsv"),
    sep = "\t",
    quote = FALSE, row.names = FALSE
  )
  writeBin(
    as.raw(c(137, 80, 78, 71)),
    file.path(product_dir, "CAT_A_leading_edge_gene_heatmap_zscore.png")
  )
  attached <- e$report_attach_category_products(
    file.path(evidence_dir, "index.html"),
    list(list(
      product = "contrast_heatmap", label = "Category heatmaps",
      path = product_dir
    )),
    "A_vs_B", "GOCC", contrast = TRUE, project_dir = root
  )
  expect_length(attached$category_products, 1L)
  assets <- attached$category_products[[1L]]$assets
  expect_true(all(vapply(
    assets, function(asset) startsWith(asset$href, "../../../../artifacts/"),
    logical(1L)
  )))
  expect_true(all(vapply(
    assets, function(asset) startsWith(asset$source_path, "artifacts/"),
    logical(1L)
  )))
  expect_false(dir.exists(file.path(
    evidence_dir, "assets", "category_products"
  )))
  e$write_shareable_report_manifest(root)
  expect_no_error(e$validate_shareable_report_manifest(root))
})

test_that("complete contrast heatmap cards replace only the duplicate gallery expectation", {
  e <- presentation_report_environment()
  complete <- list(category_products = list(
    list(category_id = "CAT_A", product = "contrast_heatmap", assets = list(
      list(format = "png"), list(format = "pdf"), list(format = "tsv"), list(format = "r")
    )),
    list(category_id = "CAT_B", product = "contrast_heatmap", assets = list(
      list(format = "png"), list(format = "tsv"), list(format = "r")
    ))
  ))
  expect_true(e$report_has_complete_category_product(complete, "contrast_heatmap"))
  expect_false(e$report_has_complete_category_product(complete, "contrast_profile"))
  incomplete <- complete
  incomplete$category_products[[2L]]$assets <- incomplete$category_products[[2L]]$assets[1:2]
  expect_false(e$report_has_complete_category_product(incomplete, "contrast_heatmap"))
  expect_false(e$report_has_complete_category_product(list(category_products = list()), "contrast_heatmap"))
})

test_that("FULL attachments follow standard evidence in installed viewer templates", {
  category_template <- paste(readLines(system.file("category-evidence", "viewer.html", package = "lisaR"), warn = FALSE), collapse = "\n")
  contrast_template <- paste(readLines(system.file("contrast-evidence", "viewer.html", package = "lisaR"), warn = FALSE), collapse = "\n")
  expect_lt(regexpr('id="member-chart-panel"', category_template, fixed = TRUE)[[1L]],
    regexpr('id="category-products"', category_template, fixed = TRUE)[[1L]])
  expect_lt(regexpr("Gene-set enrichment in A and B", contrast_template, fixed = TRUE)[[1L]],
    regexpr('id="category-products"', contrast_template, fixed = TRUE)[[1L]])
})

test_that("report assembly relocates legacy saved-product panels after standard evidence", {
  e <- presentation_report_environment()
  legacy <- paste0(
    '<main><section id="category-products" class="panel">products</section>',
    '<section id="member-chart-panel">standard evidence</section></main>'
  )
  assembled <- e$report_move_category_products_to_end(legacy)
  expect_lt(regexpr('id="member-chart-panel"', assembled, fixed = TRUE)[[1L]],
    regexpr('id="category-products"', assembled, fixed = TRUE)[[1L]])
})

# ---------------------------------------------------------------------------
# Regression for the 2026-09-09 delivered-FULL navigation defect. Every prior
# check verified attachment, hrefs, gating and hashes -- all of which were and
# remain correct -- and all of them passed on a report that was rejected on
# sight, because a selected owner's visible section set was IDENTICAL to an
# unselected owner's: the extras were subtracted from the main pages and never
# surfaced anywhere a reader could see them. These tests assert the property
# nobody had asserted: that a selected owner gains a visible, navigable,
# distinctly-kinded section, and that an unselected owner gains nothing.
# ---------------------------------------------------------------------------

full_products_fixture <- function(root, kind = "category") {
  contrast <- identical(kind, "contrast")
  page_root <- if (contrast) "contrast_evidence" else "evidence"
  owner <- if (contrast) "A_vs_B" else "A"
  evidence_dir <- file.path(root, "report_pages", page_root, owner, "GOCC")
  dir.create(file.path(evidence_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  writeLines(paste0('<script type="application/json" id="',
    if (contrast) "contrast-evidence-data" else "evidence-data", '">{"categories":[]}</script>'),
    file.path(evidence_dir, "index.html"))
  metadata <- if (contrast) {
    data.frame(key = c("contrast_id", "analysis_a", "analysis_b", "collection", "tier"),
      value = c("A_vs_B", "A", "B", "GOCC", "tier1"))
  } else {
    data.frame(key = c("analysis_id", "collection", "tier"), value = c("A", "GOCC", "tier1"))
  }
  utils::write.table(metadata, file.path(evidence_dir, "tables", "metadata.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(data.frame(category_id = c("CAT A/1", "CAT_B")),
    file.path(evidence_dir, "tables", "categories.tsv"), sep = "\t", quote = FALSE,
    row.names = FALSE)
  evidence_dir
}

# A payload shaped exactly like report_attach_category_products() returns it,
# with the artifacts actually present on disk so href resolution is real.
full_products_payload <- function(root, evidence_dir, specs) {
  products <- list()
  for (spec in specs) {
    source_path <- file.path("artifacts", spec$dir, paste0(spec$stem, ".png"))
    target <- file.path(root, source_path)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    writeBin(as.raw(c(137, 80, 78, 71)), target)
    table_path <- file.path("artifacts", spec$dir, paste0(spec$stem, "_source.tsv"))
    writeLines(c("category_id", spec$category), file.path(root, table_path))
    products[[length(products) + 1L]] <- list(category_id = spec$category,
      product = spec$product, label = spec$label, assets = list(
        list(format = "png", href = "unused", source_path = source_path,
          source_name = basename(source_path)),
        list(format = "tsv", href = "unused", source_path = table_path,
          source_name = basename(table_path))))
  }
  list(category_products = products, unclassified = character())
}

test_that("a selected owner gains a FULL figures section and an unselected owner gains none", {
  e <- presentation_report_environment()
  root <- tempfile("full products section ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  evidence_dir <- full_products_fixture(root)
  page_file <- file.path(root, "report_pages", "single_de.html")
  dir.create(dirname(page_file), recursive = TRUE, showWarnings = FALSE)
  attachment <- full_products_payload(root, evidence_dir, list(
    list(dir = "A/GOCC/volcano", stem = "CAT_A_volcano", category = "CAT A/1",
      product = "volcano", label = "Volcano overlays"),
    list(dir = "A/GOCC/volcano", stem = "CAT_B_volcano", category = "CAT_B",
      product = "volcano", label = "Volcano overlays"),
    list(dir = "A/GOCC/heatmap", stem = "CAT_A_heatmap", category = "CAT A/1",
      product = "heatmap", label = "Gene heatmaps")))
  title <- "A - GOCC - FULL figures"
  section <- e$report_full_products_section(title, attachment,
    file.path(evidence_dir, "index.html"), page_file, "category", "A", "GOCC",
    tier = "tier1", project_dir = root)

  expect_match(section, '<h2>FULL figures</h2>', fixed = TRUE)
  expect_match(section, paste0('id="', e$slug(title), '"'), fixed = TRUE)
  # Per-family counts are by distinct category, using each product's own label.
  expect_match(section, "Volcano overlays: 2 categories", fixed = TRUE)
  expect_match(section, "Gene heatmaps: 1 category", fixed = TRUE)

  # An unselected owner has an empty payload and therefore no section at all:
  # this is the assertion that fails on the delivered attempt1 tree, where the
  # selected and unselected label sets were equal.
  expect_identical(e$report_full_products_section(title,
    list(category_products = list(), unclassified = character()),
    file.path(evidence_dir, "index.html"), page_file, "category", "A", "GOCC",
    project_dir = root), "")
  expect_identical(e$report_full_products_section(title, NULL,
    file.path(evidence_dir, "index.html"), page_file, "category", "A", "GOCC",
    project_dir = root), "")
})

test_that("the FULL figures section reaches the dropdown with its own kind", {
  e <- presentation_report_environment()
  root <- tempfile("full products navigation ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dir.create(file.path(root, "report_pages"), recursive = TRUE)
  section <- e$report_navigation_section("A - GOCC - FULL figures", "FULL figures")
  expect_identical(section$kind, "full-products")
  # Not the slug(label) fallback that the missing switch case produced.
  expect_false(identical(section$kind, e$slug("FULL figures")))
  expect_identical(section$id, "a-gocc-full-figures")
  contexts <- list(list(id = "A", scientific_id = "A", label = "A", kind = "analysis",
    collections = list(list(id = "GOCC", label = "GOCC", sections = list(
      e$report_navigation_section("A - GOCC - Category evidence", "Category evidence"),
      section)))))
  target <- e$report_write_navigation_inventory(root, contexts, list())
  persisted <- jsonlite::read_json(target, simplifyVector = FALSE)
  entry <- persisted$contexts[[1L]]$collections[[1L]]$sections[[2L]]
  expect_identical(entry$kind, "full-products")
  expect_identical(entry$href, "report_pages/single_de.html#a-gocc-full-figures")
})

test_that("the FULL figures section references canonical paths and creates no file", {
  e <- presentation_report_environment()
  root <- tempfile("full products canonical ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  evidence_dir <- full_products_fixture(root)
  page_file <- file.path(root, "report_pages", "single_de.html")
  dir.create(dirname(page_file), recursive = TRUE, showWarnings = FALSE)
  attachment <- full_products_payload(root, evidence_dir, list(
    list(dir = "A/GOCC/volcano", stem = "CAT_A_volcano", category = "CAT A/1",
      product = "volcano", label = "Volcano overlays")))
  # An unreferenced sibling file: it must not appear in the section, which
  # proves the section is built from the payload and not from a second scan.
  writeBin(as.raw(c(137, 80, 78, 71)),
    file.path(root, "artifacts", "A", "GOCC", "volcano", "NOT_IN_PAYLOAD.png"))
  snapshot <- function() {
    files <- list.files(file.path(root, c("artifacts", "outputs")), recursive = TRUE,
      full.names = TRUE, all.files = TRUE, no.. = TRUE)
    data.frame(path = files, mtime = file.info(files)$mtime, stringsAsFactors = FALSE)
  }
  before <- snapshot()
  section <- e$report_full_products_section("A - GOCC - FULL figures", attachment,
    file.path(evidence_dir, "index.html"), page_file, "category", "A", "GOCC",
    tier = "tier1", project_dir = root)
  expect_identical(snapshot(), before)
  expect_false(grepl("NOT_IN_PAYLOAD", section, fixed = TRUE))

  srcs <- regmatches(section, gregexpr('src="[^"]+"', section, perl = TRUE))[[1L]]
  hrefs <- regmatches(section, gregexpr('href="[^"]+"', section, perl = TRUE))[[1L]]
  extract <- function(x) sub('^[^=]+="', "", sub('"$', "", x))
  expect_gt(length(srcs), 0L)
  for (src in extract(srcs)) {
    # Every displayed image is the canonical artifact, referenced in place.
    expect_true(startsWith(src, "../artifacts/") || startsWith(src, "../outputs/"))
    expect_true(file.exists(file.path(dirname(page_file), src)))
  }
  downloads <- Filter(function(x) startsWith(x, "../artifacts/") || startsWith(x, "../outputs/"),
    extract(hrefs))
  expect_gt(length(downloads), 0L)
  for (href in downloads) expect_true(file.exists(file.path(dirname(page_file), href)))
})

test_that("FULL figures labels are inherited verbatim and never synthesised", {
  e <- presentation_report_environment()
  root <- tempfile("full products kegg labels ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  evidence_dir <- full_products_fixture(root)
  page_file <- file.path(root, "report_pages", "single_de.html")
  dir.create(dirname(page_file), recursive = TRUE, showWarnings = FALSE)
  # Both KEGG families in one fixture so a shared code path cannot conflate
  # them: artifacts/<a>/<col>/kegg holds KEGG *gene-set* plots (hsa IDs are the
  # set's source identifier, not a painted diagram), while kegg_painter holds
  # the genuine painted maps. KEGG gene-set presentation is excluded from the
  # FULL section by explicit product decision, so the section must show the
  # painted maps only and must NOT relabel the excluded family into them.
  attachment <- full_products_payload(root, evidence_dir, list(
    list(dir = "A/GOCC/kegg", stem = "CAT_A_kegg", category = "CAT A/1",
      product = "kegg", label = "KEGG gene sets"),
    list(dir = "A/GOCC/kegg_painter", stem = "CAT_B_painted", category = "CAT_B",
      product = "kegg_painter", label = "KEGG maps")))
  section <- e$report_full_products_section("A - GOCC - FULL figures", attachment,
    file.path(evidence_dir, "index.html"), page_file, "category", "A", "GOCC",
    tier = "tier1", project_dir = root)
  # The painted-map family survives, counted under its own recorded label.
  expect_match(section, "KEGG maps: 1 category", fixed = TRUE)
  # The excluded gene-set family is absent as a family and as a category: its
  # count is not silently folded into "KEGG maps".
  expect_false(grepl("KEGG gene sets", section, fixed = TRUE))
  expect_false(grepl("CAT A/1", section, fixed = TRUE))
  expect_false(grepl("KEGG maps: 2 categories", section, fixed = TRUE))
  # No family name is invented anywhere in the prose.
  expect_false(grepl("painted pathway", section, ignore.case = TRUE))
  expect_false(grepl("pathway map", section, ignore.case = TRUE))

  # A payload of excluded gene-set products alone yields no section at all,
  # rather than an empty or mislabelled one.
  kegg_only <- full_products_payload(root, evidence_dir, list(
    list(dir = "A/GOCC/kegg", stem = "CAT_A_kegg", category = "CAT A/1",
      product = "kegg", label = "KEGG gene sets")))
  expect_identical(
    e$report_full_products_section("A - GOCC - FULL figures", kegg_only,
      file.path(evidence_dir, "index.html"), page_file, "category", "A", "GOCC",
      tier = "tier1", project_dir = root),
    ""
  )
})

test_that("evidence deep links use the exact URL form each viewer parses", {
  e <- presentation_report_environment()
  root <- tempfile("evidence deeplinks ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  page_file <- file.path(root, "report_pages", "single_de.html")
  dir.create(dirname(page_file), recursive = TRUE, showWarnings = FALSE)
  category_index <- file.path(full_products_fixture(root, "category"), "index.html")
  contrast_index <- file.path(full_products_fixture(root, "contrast"), "index.html")

  category_link <- e$report_evidence_deeplink(category_index, page_file, "category",
    "A", "GOCC", "tier1", "CAT A/1")
  expect_identical(category_link,
    "evidence/A/GOCC/index.html?category=CAT%20A%2F1&analysis_id=A&collection=GOCC&tier=tier1")

  contrast_link <- e$report_evidence_deeplink(contrast_index, page_file, "contrast",
    "A_vs_B", "GOCC", "tier1", "CAT A/1")
  # A hash, and the key is category_id. A bare `category=` here would load the
  # page and silently select the FIRST category instead of the requested one.
  expect_match(contrast_link, "#contrast_id=A_vs_B&", fixed = TRUE)
  expect_match(contrast_link, "category_id=CAT%20A%2F1", fixed = TRUE)
  expect_match(contrast_link, "analysis_a=A&analysis_b=B", fixed = TRUE)
  expect_false(grepl("[?&#]category=", contrast_link))
  expect_false(grepl("?", contrast_link, fixed = TRUE))
  # Scope keys are taken from the page's own metadata, so a link can never
  # disagree with the page it points at.
  expect_match(e$report_evidence_deeplink(contrast_index, page_file, "contrast"),
    "analysis_a=A&analysis_b=B", fixed = TRUE)
})

test_that("the evidence products anchor is injected idempotently and only when products exist", {
  e <- presentation_report_environment()
  root <- tempfile("evidence products anchor ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  for (kind in c("category-evidence", "contrast-evidence")) {
    contrast <- identical(kind, "contrast-evidence")
    data_id <- if (contrast) "contrast-evidence-data" else "evidence-data"
    template <- paste(readLines(system.file(kind, "viewer.html", package = "lisaR"),
      warn = FALSE), collapse = "\n")
    page <- sub(if (contrast) "<!-- CONTRAST_DATA -->" else "<!-- EVIDENCE_DATA -->",
      paste0('<script type="application/json" id="', data_id, '">{"categories":[]}</script>'),
      template, fixed = TRUE)
    evidence_dir <- file.path(root, kind)
    dir.create(evidence_dir, recursive = TRUE, showWarnings = FALSE)
    index <- file.path(evidence_dir, "index.html")
    writeLines(page, index)
    original <- paste(readLines(index, warn = FALSE), collapse = "\n")

    products <- list(list(category_id = "CAT_A", product = "volcano",
      label = "Volcano overlays", assets = list(list(format = "png", href = "x.png"))))
    e$report_replace_evidence_payload(index, products, data_id)
    injected <- paste(readLines(index, warn = FALSE), collapse = "\n")
    expect_match(injected, 'id="evidence-products-anchor"', fixed = TRUE)
    expect_equal(length(gregexpr('id="evidence-products-anchor"', injected,
      fixed = TRUE)[[1L]]), 1L)
    # Plain HTML with a real fragment link: present in --dump-dom and usable
    # with JavaScript disabled.
    expect_match(injected, 'href="#category-products"', fixed = TRUE)
    expect_match(injected, "1 saved figure product across 1 category.", fixed = TRUE)
    # Placed after the page intro and before the standard evidence, while
    # #category-products itself stays where the layout tests pin it.
    anchor_at <- regexpr('id="evidence-products-anchor"', injected, fixed = TRUE)[[1L]]
    intro_at <- regexpr('class="page-intro"', injected, fixed = TRUE)[[1L]]
    products_at <- regexpr('id="category-products"', injected, fixed = TRUE)[[1L]]
    expect_lt(intro_at, anchor_at)
    expect_lt(anchor_at, products_at)
    if (!contrast) {
      expect_lt(regexpr('id="member-chart-panel"', injected, fixed = TRUE)[[1L]], products_at)
    }

    # Idempotent: a second assembly replaces the block, never appends one.
    e$report_replace_evidence_payload(index, products, data_id)
    expect_identical(paste(readLines(index, warn = FALSE), collapse = "\n"), injected)

    # An empty payload injects nothing and removes any earlier block.
    e$report_replace_evidence_payload(index, list(), data_id)
    emptied <- paste(readLines(index, warn = FALSE), collapse = "\n")
    expect_false(grepl('id="evidence-products-anchor"', emptied, fixed = TRUE))

    # Byte-identity for the empty case, asserted on the injector itself: the
    # surrounding payload rewrite legitimately changes the embedded JSON, so it
    # is the injector that must be a strict no-op here.
    expect_identical(e$report_inject_evidence_products_anchor(original, list()), original)
    expect_identical(e$report_inject_evidence_products_anchor(original, NULL), original)
    once <- e$report_inject_evidence_products_anchor(original, products)
    expect_identical(e$report_inject_evidence_products_anchor(once, products), once)
    expect_identical(e$report_inject_evidence_products_anchor(once, list()), original)
  }
})

test_that("both evidence viewers reserve thumb space and settle a jump until stable", {
  for (kind in c("category-evidence", "contrast-evidence")) {
    css <- paste(readLines(system.file(kind, "viewer.css", package = "lisaR"), warn = FALSE),
      collapse = "\n")
    js <- paste(readLines(system.file(kind, "viewer.js", package = "lisaR"), warn = FALSE),
      collapse = "\n")
    # Reserved box before the lazy image loads.
    expect_match(css, ".category-product-thumb{", fixed = TRUE)
    expect_match(css, "aspect-ratio", fixed = TRUE)
    expect_match(js, "category-product-thumb", fixed = TRUE)
    # The settle-until-stable scroll ported from lisa_shell.js:99-160.
    expect_match(js, "scrollIntoView", fixed = TRUE)
    expect_match(js, "ResizeObserver", fixed = TRUE)
    expect_match(js, "requestAnimationFrame", fixed = TRUE)
    expect_match(js, "document.body", fixed = TRUE)
    for (type in c("wheel", "keydown", "touchstart")) expect_match(js, type, fixed = TRUE)
    # Bounded deadline, so the guard can never hold the page hostage.
    expect_match(js, "2000", fixed = TRUE)
    # The injected anchor is wired to the per-category count and the settle.
    expect_match(js, "evidence-products-anchor", fixed = TRUE)
    expect_match(js, "evidence-products-jump", fixed = TRUE)
  }
})

# Regression for the 2026-09-09 delivered-FULL defect: a scoped FULL was built
# by calling the extension API (render_lisa_categories()/plan_lisa_extension())
# on its own, which produces a separate "LISA derived report" gallery instead
# of the normal standard shell with extras attached in place. The normal
# assembly contract is: every analysis/contrast keeps its standard evidence
# page; only the explicitly selected owners gain attached FULL products
# (referenced from `artifacts/`, never copied); unselected owners are left
# exactly as their standard page, with no "missing"/error marker and no
# separate gallery file anywhere in the tree.
test_that("selected-owner FULL gating attaches extras only to the chosen owners across a 5-analysis/2-contrast base", {
  e <- presentation_report_environment()
  root <- tempfile("selected full gating ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)

  analyses <- c(
    "responders_on_vs_pre", "pd_on_vs_pre", "responders_vs_pd_pre",
    "responders_vs_pd_on", "differential_longitudinal_response"
  )
  contrasts <- c("response_dynamics", "response_separation_by_time")
  selected_analyses <- c("responders_vs_pd_on", "responders_vs_pd_pre")
  selected_contrast <- "response_separation_by_time"

  make_evidence_page <- function(owner, contrast) {
    kind <- if (contrast) "contrast_evidence" else "evidence"
    owner_key <- if (contrast) "contrast_id" else "analysis_id"
    data_id <- if (contrast) "contrast-evidence-data" else "evidence-data"
    evidence_dir <- file.path(root, "report_pages", kind, owner, "GOCC")
    dir.create(file.path(evidence_dir, "tables"), recursive = TRUE)
    writeLines(paste0('<script type="application/json" id="', data_id, '">{"categories":[]}</script>'),
      file.path(evidence_dir, "index.html"))
    utils::write.table(data.frame(key = c(owner_key, "collection"), value = c(owner, "GOCC")),
      file.path(evidence_dir, "tables", "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    utils::write.table(data.frame(category_id = "CAT_A"), file.path(evidence_dir, "tables", "categories.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE)
    evidence_dir
  }
  make_artifacts <- function(owner, contrast) {
    product_dir <- if (contrast) {
      file.path(root, "artifacts", "contrasts", owner, "GOCC", "contrast_heatmap")
    } else {
      file.path(root, "artifacts", owner, "GOCC", "volcano")
    }
    dir.create(product_dir, recursive = TRUE)
    writeLines(c("category_id", "CAT_A"), file.path(product_dir, "CAT_A_source.tsv"))
    writeBin(as.raw(c(137, 80, 78, 71)), file.path(product_dir, "CAT_A.png"))
    product_dir
  }

  evidence_dirs <- stats::setNames(vapply(analyses, make_evidence_page, character(1L), contrast = FALSE), analyses)
  contrast_dirs <- stats::setNames(vapply(contrasts, make_evidence_page, character(1L), contrast = TRUE), contrasts)
  for (analysis in selected_analyses) make_artifacts(analysis, contrast = FALSE)
  make_artifacts(selected_contrast, contrast = TRUE)

  attach_for <- function(owner, evidence_dir, contrast) {
    product <- if (contrast) "contrast_heatmap" else "volcano"
    product_dir <- if (contrast) {
      file.path(root, "artifacts", "contrasts", owner, "GOCC", "contrast_heatmap")
    } else {
      file.path(root, "artifacts", owner, "GOCC", "volcano")
    }
    e$report_attach_category_products(file.path(evidence_dir, "index.html"),
      list(list(product = product, label = "extra", path = product_dir)),
      owner, "GOCC", contrast = contrast, project_dir = root)
  }

  analysis_results <- Map(attach_for, analyses, evidence_dirs, MoreArgs = list(contrast = FALSE))
  contrast_results <- Map(attach_for, contrasts, contrast_dirs, MoreArgs = list(contrast = TRUE))

  for (analysis in selected_analyses) {
    expect_length(analysis_results[[analysis]]$category_products, 1L)
    href <- analysis_results[[analysis]]$category_products[[1L]]$assets[[1L]]$href
    expect_true(startsWith(href, "../../../../artifacts/"))
  }
  for (analysis in setdiff(analyses, selected_analyses)) {
    expect_length(analysis_results[[analysis]]$category_products, 0L)
    html <- paste(readLines(file.path(evidence_dirs[[analysis]], "index.html"), warn = FALSE), collapse = "\n")
    expect_match(html, '"categories":\\[\\]', perl = TRUE)
  }
  expect_length(contrast_results[[selected_contrast]]$category_products, 1L)
  for (contrast_name in setdiff(contrasts, selected_contrast)) {
    expect_length(contrast_results[[contrast_name]]$category_products, 0L)
  }

  # FULL is the standard shell plus attachments, never a second root gallery:
  # every evidence page from the 5+2 base still exists, and nothing produced a
  # separate top-level report/index outside report_pages/.
  all_index_pages <- list.files(root, pattern = "^index[.]html$", recursive = TRUE, full.names = TRUE)
  expect_length(all_index_pages, length(analyses) + length(contrasts))
  expect_false(file.exists(file.path(root, "index.html")))
  expect_false(any(vapply(all_index_pages, function(p) grepl("LISA derived report",
    paste(readLines(p, warn = FALSE), collapse = "\n"), fixed = TRUE), logical(1L))))

  e$write_shareable_report_manifest(root)
  expect_no_error(e$validate_shareable_report_manifest(root))
})

test_that("NES report selector respects subsets and keeps exact recipe inputs", {
  e <- presentation_report_environment()
  root <- tempfile("nes report ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  collection <- file.path(root, "outputs", "single_de", "A", "collection_GOCC")
  dir.create(file.path(collection, "qc"), recursive = TRUE)
  folder <- file.path(collection, "plots", "category_nes")
  dir.create(folder, recursive = TRUE)
  dir.create(file.path(root, "report_pages"))
  source <- file.path(folder, "category_nes_source.tsv")
  writeLines(c('"category_id"\t"mean_NES"', '"001"\t"1.2345678901234567"'), source)
  checksum <- digest::digest(file = source, algo = "sha256")
  # A stale run-local recipe must never become the code copied to the report.
  writeLines("stop('untrusted fixture recipe')", file.path(folder, "category_nes_recipe.R"))
  variants <- c("direction", "dispersion")
  rows <- lapply(variants, function(v) {
    stem <- paste0("category_nes_", v)
    writeBin(as.raw(c(137, 80, 78, 71)), file.path(folder, paste0(stem, ".png")))
    jsonlite::write_json(list(variant = v, source_sha256 = checksum,
      x_limits = c(-3, 3), schema_version = "1.0"), file.path(folder, paste0(stem, "_settings.json")), auto_unbox = TRUE)
    data.frame(variant = v, status = "rendered", format = "png",
      path = file.path("plots", "category_nes", paste0(stem, ".png")),
      source = "plots/category_nes/category_nes_source.tsv",
      settings = file.path("plots", "category_nes", paste0(stem, "_settings.json")),
      recipe = "plots/category_nes/category_nes_recipe.R")
  })
  manifest <- file.path(collection, "qc", "A_category_nes_variants.tsv")
  utils::write.table(do.call(rbind, rows), manifest, sep = "\t", quote = FALSE, row.names = FALSE)
  html <- e$report_nes_variants_section(collection, file.path(root, "report_pages", "single_de.html"), "A GOCC")
  expect_match(html, '<option value="direction">', fixed = TRUE)
  expect_match(html, '<option value="dispersion">', fixed = TRUE)
  expect_false(grepl('<option value="clean">', html, fixed = TRUE))
  expect_match(html, 'data-nes-panel="direction">', fixed = TRUE)
  expect_match(html, 'data-nes-panel="dispersion" hidden', fixed = TRUE)
  expect_match(html, '>Mean + median</option>', fixed = TRUE)
  expect_match(html, '>Percentiles</option>', fixed = TRUE)
  expect_false(grepl("Category NES views", html, fixed = TRUE))
  expect_match(html, "Download this view", fixed = TRUE)
  expect_match(html, "not a confidence interval", fixed = TRUE)
  copies <- list.files(file.path(root, "report_figure_data"), recursive = TRUE, full.names = TRUE)
  copied_source <- copies[basename(copies) == basename(source)]
  expect_length(copied_source, 1L)
  expect_identical(digest::digest(file = copied_source, algo = "sha256"), checksum)
  recipe <- copies[basename(copies) == "reproduce_category_nes.R"]
  expect_length(recipe, 1L)
  expect_false(any(grepl("untrusted fixture", readLines(recipe), fixed = TRUE)))
  expect_true(any(grepl("lisa_reproduce_category_nes_variant", readLines(recipe), fixed = TRUE)))
  settings <- file.path(folder, "category_nes_direction_settings.json")
  jsonlite::write_json(list(variant = "direction", source_sha256 = "incorrect"), settings, auto_unbox = TRUE)
  expect_error(e$report_nes_variants_section(collection, file.path(root, "report_pages", "single_de.html"), "A GOCC"),
    "settings disagree")
})

test_that("four NES views and contrast subsets keep distinct artifacts", {
  e <- presentation_report_environment()
  root <- tempfile("four views ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  collection <- file.path(root, "outputs", "contrast", "A_B", "collection_GOCC")
  folder <- file.path(collection, "plots", "category_nes")
  dir.create(folder, recursive = TRUE)
  dir.create(file.path(collection, "qc"), recursive = TRUE)
  dir.create(file.path(root, "report_pages"))
  variants <- c("clean", "percentages", "direction", "dispersion")
  subsets <- c("all", "same_direction", "opposite_direction")
  rows <- list()
  for (subset in subsets) for (variant in variants) {
    stem <- paste("category", subset, variant, sep = "_")
    source <- file.path(folder, paste0(stem, "_source.tsv"))
    writeLines(c("category_id\tmacrogroup_name\tmean_NES", "001\tCell biology\t1.5"), source)
    settings <- file.path(folder, paste0(stem, "_settings.json"))
    jsonlite::write_json(list(variant = variant, plot_set = subset,
      source_sha256 = digest::digest(file = source, algo = "sha256")), settings, auto_unbox = TRUE)
    figure <- file.path(folder, paste0(stem, ".png"))
    writeBin(as.raw(c(137, 80, 78, 71)), figure)
    recipe <- file.path(folder, "category_recipe.R")
    writeLines("stop('do not use runtime recipe')", recipe)
    rows[[length(rows) + 1L]] <- data.frame(variant = variant, plot_set = subset,
      status = "rendered", format = "png", path = paste0("plots/category_nes/", basename(figure)),
      source = paste0("plots/category_nes/", basename(source)),
      settings = paste0("plots/category_nes/", basename(settings)),
      recipe = "plots/category_nes/category_recipe.R")
  }
  manifest <- file.path(collection, "qc", "A_B_category_nes_variants.tsv")
  utils::write.table(do.call(rbind, rows), manifest, sep = "\t", quote = FALSE, row.names = FALSE)
  page <- file.path(root, "report_pages", "contrasts.html")
  html <- e$report_nes_variants_section(collection, page, "A-B GOCC")
  for (label in c("Clean", "Percentages", "Mean + median", "Percentiles", "All", "Same direction", "Opposite direction"))
    expect_match(html, paste0(">", label, "</option>"), fixed = TRUE)
  expect_equal(length(regmatches(html, gregexpr('class="nes-variant-panel"', html, fixed = TRUE))[[1L]]), 12L)
  expect_equal(length(regmatches(html, gregexpr(' hidden>', html, fixed = TRUE))[[1L]]), 11L)
  expect_match(html, 'data-nes-plot-set="opposite_direction" data-nes-panel="percentages" hidden', fixed = TRUE)
  copies <- list.files(file.path(root, "report_figure_data"), recursive = TRUE, full.names = TRUE)
  expect_length(copies[grepl("_source.tsv$", copies)], 12L)
  expect_length(copies[grepl("_settings.json$", copies)], 12L)
  broken <- do.call(rbind, rows); broken$plot_set[[1L]] <- "unknown"
  utils::write.table(broken, manifest, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_error(e$report_nes_variants_section(collection, page, "A-B GOCC"), "identities")
})

test_that("unified category section replaces duplicate summaries but retains gene-set overview", {
  e <- presentation_report_environment()
  e$report_nes_variants_section <- function(...) '<div data-nes-variants="one">Four views</div>'
  e$figure_card <- function(path, ...) paste0('<figure>', basename(path), '</figure>')
  html <- e$report_lisa_categories_section("A - GOCC - LISA categories", "collection",
    c("A_GSEA_lollipop.png", "A_GSEA_lollipop_direction_stats.png", "A_GSEA_pathway_dotplot.png", "A_ORA_barplot.png"),
    "single_de.html")
  expect_match(html, '<h2>LISA categories</h2>', fixed = TRUE)
  expect_match(html, "Supercategory bands", fixed = TRUE)
  expect_false(grepl("lollipop|Category NES views", html))
  expect_match(html, "A_GSEA_pathway_dotplot.png", fixed = TRUE)
  expect_match(html, "A_ORA_barplot.png", fixed = TRUE)
  expect_match(html, "Gene-set overview", fixed = TRUE)
  e$report_nes_variants_section <- function(...) ""
  e$report_output_policy <- data.frame(product = c("png", "svg", "pdf"), requested = FALSE)
  no_figures <- e$report_lisa_categories_section("A GOCC", "collection", character(), "single_de.html")
  expect_match(no_figures, "not_requested", fixed = TRUE)
  expect_false(grepl("<img|data-nes-select|missing expected", no_figures))
})

test_that("scope context uses readable labels and preserves exact identifiers", {
  e <- presentation_report_environment()
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path), add = TRUE)
  utils::write.table(data.frame(key = c("tier", "positive_contrast", "gsea_padj_cutoff"),
    value = c("core", "Treatment versus baseline", "0.123456789")),
    path, sep = "\t", quote = FALSE, row.names = FALSE)
  context <- e$report_scope_context("001", "Responders", "GOCC", path)
  text <- gsub("&quot;", '"', context, fixed = TRUE)
  parsed <- jsonlite::fromJSON(text)
  expect_identical(parsed$analysis, "Responders")
  expect_identical(parsed$exact_ids$analysis_id, "001")
  expect_identical(parsed$direction, "Treatment versus baseline")
  expect_match(parsed$cutoff, "0.123456789", fixed = TRUE)
  expect_false(grepl(path, context, fixed = TRUE))
})


test_that("figure cards reject missing requested SVG and omit unrequested SVG", {
  e <- presentation_report_environment()
  root <- tempfile("optional-svg-"); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  png <- file.path(root, "figure.png"); writeBin(as.raw(c(137, 80, 78, 71)), png)
  e$short_media_copy <- function(path, page_file) basename(path)
  e$figure_source_contract <- function(path, page_file) list(source_tsv = "", recipe_r = "")
  e$report_output_policy <- data.frame(product = c("svg", "source_data", "recipes"), requested = FALSE)
  page <- file.path(root, "index.html")
  expect_false(grepl("SVG missing", e$figure_card(png, page), fixed = TRUE))
  e$report_output_policy$requested[[1L]] <- TRUE
  expect_error(e$figure_card(png, page), "Requested figure format missing:.*svg")
  writeLines('<svg xmlns="http://www.w3.org/2000/svg"/>', file.path(root, "figure.svg"))
  card <- e$figure_card(png, page)
  expect_match(card, 'href="figure.svg">SVG</a>', fixed = TRUE)
  expect_false(grepl("SVG missing", card, fixed = TRUE))
})

test_that("report_full_scope_owner detects a selected owner from either its saved artifacts or its full-pipeline outputs", {
  e <- presentation_report_environment()
  root <- tempfile("full scope owner ")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  expect_false(e$report_full_scope_owner(root, "unselected_analysis"))
  expect_false(e$report_full_scope_owner(root, "unselected_contrast", contrast = TRUE))
  dir.create(file.path(root, "artifacts", "selected_via_artifacts"), recursive = TRUE)
  expect_true(e$report_full_scope_owner(root, "selected_via_artifacts"))
  dir.create(file.path(root, "outputs", "gene_level", "single_de", "selected_via_pipeline"), recursive = TRUE)
  expect_true(e$report_full_scope_owner(root, "selected_via_pipeline"))
  dir.create(file.path(root, "artifacts", "contrasts", "selected_contrast"), recursive = TRUE)
  expect_true(e$report_full_scope_owner(root, "selected_contrast", contrast = TRUE))
  expect_false(e$report_full_scope_owner(root, "selected_via_artifacts", contrast = TRUE))
})

# Regression: directory existence
# alone cannot distinguish "never selected" from "selected but its expected
# products never materialized" -- both look identical on disk (no
# directory). When contract_manifest.tsv declares full_scope_analyses/
# full_scope_contrasts, that declared list must be authoritative and must
# NOT additionally require the directory to exist: a declared-selected
# owner with a totally absent directory must still return TRUE, so its
# layers fail closed as "missing" rather than silently downgrading to
# "not_requested".
test_that("report_full_scope_owner treats a declared contract_manifest.tsv selection as authoritative even when the owner's directory is absent", {
  e <- presentation_report_environment()
  root <- tempfile("declared full scope ")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  utils::write.table(
    data.frame(
      key = c("report_mode", "full_scope_analyses", "full_scope_contrasts"),
      value = c("full", "selected_but_missing_analysis", "selected_but_missing_contrast")
    ),
    file.path(root, "contract_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
  )
  # Neither owner has ANY directory anywhere (artifacts/, outputs/gene_level/...):
  # this is the exact "reassembly declared intent but never produced the
  # product" / partial-pipeline-failure scenario the fix must catch.
  expect_true(e$report_full_scope_owner(root, "selected_but_missing_analysis"))
  expect_true(e$report_full_scope_owner(root, "selected_but_missing_contrast", contrast = TRUE))
  # An owner not named in the declared list stays genuinely unselected, even
  # though its directory is equally absent -- must remain eligible for
  # "not_requested", never silently promoted to "missing".
  expect_false(e$report_full_scope_owner(root, "genuinely_unselected_analysis"))
  expect_false(e$report_full_scope_owner(root, "genuinely_unselected_contrast", contrast = TRUE))
  # A declared-but-empty key (ordinary STANDARD run: nothing in full scope)
  # must resolve to an empty list, not fall back to the directory heuristic.
  utils::write.table(
    data.frame(key = c("report_mode", "full_scope_analyses", "full_scope_contrasts"),
      value = c("standard", "", "")),
    file.path(root, "contract_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
  )
  dir.create(file.path(root, "artifacts", "would_have_matched_old_heuristic"), recursive = TRUE)
  expect_false(e$report_full_scope_owner(root, "would_have_matched_old_heuristic"))
})

# End-to-end trace of the fail-closed guarantee itself (LISA-REPORT-LAYER-004):
# report_full_scope_owner() feeding the exact empty_state ternary used at the
# real call sites, through layer_section()/missing_layer_section(), into
# report_expected_layer_diagnostics() -- the same scan the top-level report
# builder runs over the assembled HTML before it will ship a report. This
# proves a selected-but-absent owner cannot silently render as
# "not_requested" without invoking the full multi-thousand-line pipeline
# fixture; report_full_scope_owner() is exercised for real (not stubbed),
# and the downstream gate logic is exercised for real on its actual inputs.
test_that("a contract-declared selected owner with no materialized products fails closed through the real LISA-REPORT-LAYER-004 gate logic", {
  e <- presentation_report_environment()
  root <- tempfile("fail closed gate ")
  dir.create(root, recursive = TRUE)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  utils::write.table(
    data.frame(key = c("report_mode", "full_scope_analyses", "full_scope_contrasts"),
      value = c("full", "responders_on_vs_pre", "")),
    file.path(root, "contract_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
  )
  # No artifacts/, outputs/gene_level/... directory exists anywhere for this
  # (or any) owner: the "upstream failure wiped everything" scenario.
  analysis_full_scope <- e$report_full_scope_owner(root, "responders_on_vs_pre", contrast = FALSE)
  expect_true(analysis_full_scope)
  gene_collection_dir <- file.path(root, "outputs", "gene_level", "single_de", "responders_on_vs_pre", "collection_GOCC")
  spec <- list(path = file.path(gene_collection_dir, "category_gene_cards"),
    files = e$pngs(file.path(gene_collection_dir, "category_gene_cards"), "\\.png$"),
    note = "Category GeneCards.",
    empty_state = if (!analysis_full_scope) "not_requested" else "missing")
  expect_identical(spec$empty_state, "missing")
  expect_length(spec$files, 0L)
  section_html <- e$layer_section("responders_on_vs_pre - GOCC - GeneCards", spec$files,
    file.path(root, "report_pages", "single_de.html"), spec$note, empty_state = spec$empty_state)
  expect_match(section_html, "missing expected layer", fixed = TRUE)
  expect_false(grepl("not_requested", section_html, fixed = TRUE))
  diagnostics <- e$report_expected_layer_diagnostics(section_html)
  expect_length(diagnostics, 1L)
  expect_match(diagnostics, "family=GeneCards", fixed = TRUE)

  # Contrast to the "genuinely unselected" owner: same total absence of any
  # directory, but not named in full_scope_analyses -- must resolve to
  # "not_requested" and must NOT trip the fail-closed gate.
  unselected_full_scope <- e$report_full_scope_owner(root, "differential_longitudinal_response", contrast = FALSE)
  expect_false(unselected_full_scope)
  unselected_empty_state <- if (!unselected_full_scope) "not_requested" else "missing"
  unselected_html <- e$layer_section("differential_longitudinal_response - GOCC - GeneCards", character(),
    file.path(root, "report_pages", "single_de.html"), "Category GeneCards.", empty_state = unselected_empty_state)
  expect_match(unselected_html, "not_requested", fixed = TRUE)
  expect_length(e$report_expected_layer_diagnostics(unselected_html), 0L)
})
