extension_navigation_units <- function() {
  data.frame(
    unit_type = c("single_de", "single_de", "single_de", "contrast"),
    analysis_id = c("Analysis A", "Analysis A", "Analysis B", ""),
    contrast_id = c("", "", "", "A_vs_B_pair"),
    collection = c("GOCC", "GOCC", "GOMF", "GOCC"),
    category_id = c("CAT_A", "CAT_A", "CAT_B", "CAT_A"),
    supercategory_id = c("SUPER_1", "SUPER_1", "SUPER_2", "SUPER_1"),
    product = c("member_gene_sets", "gene_cards", "volcano", "contrast_profile"),
    stringsAsFactors = FALSE
  )
}

extension_navigation_html <- function(paths, units = extension_navigation_units(), source_run = tempfile("immutable-source-")) {
  root <- tempfile("lisa-extension-navigation-")
  dir.create(root)
  # Index attachment reads the saved source/matrix metadata; filenames alone
  # deliberately cannot classify an artifact.
  for (path in paths) {
    target <- file.path(root, path)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    prefixes <- if (nrow(units)) vapply(seq_len(nrow(units)), function(i) {
      contrast <- identical(as.character(units$unit_type[[i]]), "contrast")
      paste(c("artifacts", if (contrast) "contrasts", if (contrast) units$contrast_id[[i]] else units$analysis_id[[i]],
        units$collection[[i]], units$product[[i]]), collapse = "/")
    }, character(1L)) else character()
    matched <- which(startsWith(path, paste0(prefixes, "/")))
    category <- if (length(matched) == 1L) as.character(units$category_id[[matched]]) else ""
    if (grepl("_(source|matrix)[.]tsv$", path, ignore.case = TRUE)) {
      writeLines(c("category_id", category), target, useBytes = TRUE)
    } else writeBin(as.raw(0L), target)
  }
  inventory <- data.frame(path = paths,
    bytes = as.numeric(file.info(file.path(root, paths))$size),
    sha256 = vapply(file.path(root, paths), lisaR:::lisa_sha256_file, character(1L)),
    stringsAsFactors = FALSE)
  lisaR:::lisa_extension_write_index(root,
    list(mode = "full", source_run = source_run, units = units), inventory)
  list(root = root, html = paste(readLines(file.path(root, "index.html"), warn = FALSE), collapse = "\n"),
    nav = jsonlite::read_json(file.path(root, "navigation_inventory.json"), simplifyVector = FALSE))
}

extension_saved_evidence <- function(owner = "Analysis A", collection = "GOCC", contrast = FALSE,
    exact_owner = owner, categories = "CAT_A") {
  root <- tempfile("lisa-extension-source-"); dir.create(root)
  page <- file.path(root, "report_pages", if (contrast) "contrast_evidence" else "evidence", owner, collection, "index.html")
  dir.create(dirname(page), recursive = TRUE)
  metadata <- stats::setNames(list(exact_owner, collection), c(if (contrast) "contrast_id" else "analysis_id", "collection"))
  payload <- list(metadata = metadata,
    categories = lapply(categories, function(id) list(category_id = id)), sets = list(), genes = list(), leading_edges = list())
  writeLines(paste0('<html><body><script type="application/json" id="', if (contrast) "contrast-evidence-data" else "evidence-data", '">',
    jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null"), "</script></body></html>"), page, useBytes = TRUE)
  root
}

test_that("full extension has real fixed contextual inventory and one shared full-logo shell", {
  paths <- c(
    "artifacts/Analysis A/GOCC/member_gene_sets/CAT_A.png",
    "artifacts/Analysis A/GOCC/member_gene_sets/CAT_A_source.tsv",
    "artifacts/Analysis A/GOCC/gene_cards/CAT_A_category_gene_card.png",
    "artifacts/Analysis A/GOCC/gene_cards/CAT_A_category_gene_card_source.tsv",
    "artifacts/Analysis B/GOMF/volcano/CAT_B_category_volcano_overlay.svg",
    "artifacts/Analysis B/GOMF/volcano/CAT_B_category_volcano_overlay_source.tsv",
    "artifacts/contrasts/A_vs_B_pair/GOCC/contrast_profile/CAT_A.png",
    "artifacts/contrasts/A_vs_B_pair/GOCC/contrast_profile/CAT_A_source.tsv"
  )
  result <- extension_navigation_html(paths)
  expect_identical(length(result$nav$contexts), 3L)
  expect_identical(vapply(result$nav$contexts, `[[`, character(1L), "kind"), c("analysis", "analysis", "contrast"))
  expect_identical(vapply(result$nav$contexts, `[[`, character(1L), "scientific_id"),
    c("Analysis A", "Analysis B", "A_vs_B_pair"))
  expect_identical(vapply(result$nav$contexts[[1L]]$collections[[1L]]$sections, `[[`, character(1L), "label"),
    "Saved category products")
  # Analysis A has two product families for CAT_A but one category card.
  expect_identical(length(gregexpr('data-extension-card', result$html, fixed = TRUE)[[1L]]), 3L)
  expect_match(result$html, 'data-lisa-shell', fixed = TRUE)
  expect_match(result$html, 'data-lisa-nav-select="context"', fixed = TRUE)
  expect_match(result$html, 'data-lisa-nav-select="collection"', fixed = TRUE)
  expect_match(result$html, 'data-lisa-nav-select="section"', fixed = TRUE)
  expect_match(result$html, 'LISA_logo_A1_muted_red_S_automated_annotation_final.svg', fixed = TRUE)
  expect_true(file.exists(file.path(result$root, "report_assets/lisa-shell/lisa_shell.js")))
  for (context in result$nav$contexts) {
    expect_match(result$html, paste0('data-lisa-nav-context="', context$id, '"'), fixed = TRUE)
    for (collection in context$collections) {
      expect_match(result$html, paste0('data-lisa-nav-collection="', collection$id, '"'), fixed = TRUE)
      for (section in collection$sections) {
        expect_identical(section$href, paste0("#", section$id))
        expect_match(result$html, paste0('id="', section$id, '"'), fixed = TRUE)
      }
    }
  }
  expect_identical(length(gregexpr('data-lisa-nav-context="[^"]+" open', result$html, perl = TRUE)[[1L]]), 1L)
  expect_identical(length(gregexpr('data-lisa-nav-collection="[^"]+" data-lisa-context-json="[^"]+" open',
    result$html, perl = TRUE)[[1L]]), 1L)
})

test_that("full extension offers one figure card with matching available formats and reproducibility files", {
  base <- "artifacts/Analysis A/GOCC/member_gene_sets/CAT_A"
  files <- paste0(base, c(".png", ".svg", ".pdf", "_source.tsv", "_settings.json", "_recipe.R"))
  result <- extension_navigation_html(files, extension_navigation_units()[1L, , drop = FALSE])
  cards <- gregexpr('class="lisa-extension-card"', result$html, fixed = TRUE)[[1L]]
  expect_identical(length(cards), 1L)
  expect_match(result$html, 'src="artifacts/Analysis%20A/GOCC/member_gene_sets/CAT_A.png"', fixed = TRUE)
  expect_false(grepl('<img[^>]*[.]pdf', result$html, perl = TRUE))
  for (file in files) expect_match(result$html,
    paste0('download href="', lisaR:::lisa_extension_url_path(file), '"'), fixed = TRUE)
  expect_match(result$html, "Source TSV", fixed = TRUE)
  expect_match(result$html, "Settings JSON", fixed = TRUE)
  expect_match(result$html, "Recipe R", fixed = TRUE)
  expect_match(result$html, 'data-extension-super', fixed = TRUE)
  expect_match(result$html, 'data-supercategory="SUPER_1"', fixed = TRUE)
  expect_match(result$html, 'value="SUPER_1"', fixed = TRUE)
  section <- result$nav$contexts[[1L]]$collections[[1L]]$sections[[1L]]
  expect_identical(section$category_products[[1L]]$category_id, "CAT_A")
  expect_identical(section$category_products[[1L]]$product, "member_gene_sets")
  expect_identical(length(section$global_products), 0L)
})

test_that("extension reuses saved evidence with canonical artifact references instead of copies", {
  source_run <- extension_saved_evidence()
  paths <- c("artifacts/Analysis A/GOCC/gene_cards/CAT_A_category_gene_card.png",
    "artifacts/Analysis A/GOCC/gene_cards/CAT_A_category_gene_card_source.tsv")
  result <- extension_navigation_html(paths, extension_navigation_units()[2L, , drop = FALSE], source_run)
  section <- result$nav$contexts[[1L]]$collections[[1L]]$sections[[1L]]
  expect_identical(section$kind, "category-evidence")
  expect_identical(section$label, "Evidence")
  expect_match(section$href, "^evidence/Analysis A/GOCC/index[.]html$")
  evidence <- file.path(result$root, section$href)
  expect_true(file.exists(evidence))
  html <- paste(readLines(evidence, warn = FALSE), collapse = "\n")
  expect_match(html, 'id="evidence-data"', fixed = TRUE)
  payload <- lisaR:::lisa_extension_read_payload(evidence, "evidence-data")
  hrefs <- unlist(lapply(payload$category_products, function(product)
    vapply(product$assets, `[[`, character(1L), "href")), use.names = FALSE)
  expect_true(all(startsWith(hrefs, "../../../artifacts/")))
  expect_true(all(file.exists(file.path(dirname(evidence), utils::URLdecode(hrefs)))))
  expect_match(html, 'id="lisa-report-navigation"', fixed = TRUE)
  expect_false(grepl('<div data-extension-gallery', result$html, fixed = TRUE))
  expect_false(dir.exists(file.path(dirname(evidence), "assets", "category_products")))
})

test_that("saved evidence selects the union of product categories and preserves recipe bytes", {
  source_run <- extension_saved_evidence(categories = c("CAT_A", "CAT_B", "CAT_UNUSED"))
  source_dir <- file.path(source_run, "report_pages", "evidence", "Analysis A", "GOCC")
  dir.create(file.path(source_dir, "tables"))
  lisaR:::write_lisa_tsv(data.frame(category_id = c("CAT_A", "CAT_B", "CAT_UNUSED"), value = c(3, 7, 11)),
    file.path(source_dir, "tables", "categories.tsv"))
  writeLines("# Saved offline fixture recipe: copy exactly; never execute.", file.path(source_dir, "reproduce_category_evidence.R"))
  units <- extension_navigation_units()[1:2, , drop = FALSE]
  units$category_id[[2L]] <- "CAT_B"
  paths <- c("artifacts/Analysis A/GOCC/member_gene_sets/CAT_A_source.tsv",
    "artifacts/Analysis A/GOCC/gene_cards/CAT_B_source.tsv")
  result <- extension_navigation_html(paths, units, source_run)
  section <- result$nav$contexts[[1L]]$collections[[1L]]$sections[[1L]]
  page <- file.path(result$root, section$href)
  payload <- lisaR:::lisa_extension_read_payload(page, "evidence-data")
  expect_identical(vapply(payload$categories, `[[`, character(1L), "category_id"), c("CAT_A", "CAT_B"))
  expect_identical(vapply(payload$category_products, `[[`, character(1L), "category_id"), c("CAT_A", "CAT_B"))
  table <- lisaR:::read_lisa_tsv(file.path(dirname(page), "tables", "categories.tsv"))
  expect_identical(as.character(table$category_id), c("CAT_A", "CAT_B"))
  expect_identical(lisaR:::lisa_sha256_file(file.path(source_dir, "reproduce_category_evidence.R")),
    lisaR:::lisa_sha256_file(file.path(dirname(page), "reproduce_category_evidence.R")))
  for (ext in c("css", "js")) expect_identical(
    lisaR:::lisa_sha256_file(system.file("category-evidence", paste0("viewer.", ext), package = "lisaR")),
    lisaR:::lisa_sha256_file(file.path(dirname(page), "assets", paste0("category-evidence.", ext))))
})

test_that("contrast evidence resolves only the exact saved index identity and copies its recipe", {
  source_run <- extension_saved_evidence("A_vs_B_pair", contrast = TRUE, exact_owner = "A_vs_B")
  dir.create(file.path(source_run, "config"))
  index <- data.frame(contrast_id = "A_vs_B", output_id = "pair")
  lisaR:::write_lisa_tsv(index, file.path(source_run, "config", "contrast_index.tsv"))
  source_dir <- file.path(source_run, "report_pages", "contrast_evidence", "A_vs_B_pair", "GOCC")
  recipe <- file.path(source_dir, "reproduce_contrast_evidence.R")
  writeLines("# Saved paired evidence fixture recipe", recipe)
  paths <- "artifacts/contrasts/A_vs_B_pair/GOCC/contrast_profile/CAT_A_source.tsv"
  unit <- extension_navigation_units()[4L, , drop = FALSE]
  result <- extension_navigation_html(paths, unit, source_run)
  section <- result$nav$contexts[[1L]]$collections[[1L]]$sections[[1L]]
  expect_identical(section$kind, "contrast-evidence")
  page <- file.path(result$root, section$href)
  payload <- lisaR:::lisa_extension_read_payload(page, "contrast-evidence-data")
  expect_identical(payload$metadata$contrast_id, "A_vs_B")
  href <- payload$category_products[[1L]]$assets[[1L]]$href
  expect_match(href, "^../../../artifacts/contrasts/A_vs_B_pair/GOCC/contrast_profile/")
  expect_true(file.exists(file.path(dirname(page), utils::URLdecode(href))))
  expect_false(dir.exists(file.path(dirname(page), "assets", "category_products")))
  expect_identical(lisaR:::lisa_sha256_file(recipe),
    lisaR:::lisa_sha256_file(file.path(dirname(page), basename(recipe))))
  for (ext in c("css", "js")) expect_identical(
    lisaR:::lisa_sha256_file(system.file("contrast-evidence", paste0("viewer.", ext), package = "lisaR")),
    lisaR:::lisa_sha256_file(file.path(dirname(page), "assets", paste0("contrast-evidence.", ext))))
  index$output_id <- "other"
  lisaR:::write_lisa_tsv(index, file.path(source_run, "config", "contrast_index.tsv"))
  mismatch <- extension_navigation_html(paths, unit, source_run)
  expect_identical(mismatch$nav$contexts[[1L]]$collections[[1L]]$sections[[1L]]$kind, "category-products")
})

test_that("nested evidence navigation resolves saved-gallery anchors and product downloads", {
  source_run <- extension_saved_evidence()
  units <- extension_navigation_units()[c(1L, 3L), , drop = FALSE]
  paths <- c("artifacts/Analysis A/GOCC/member_gene_sets/CAT_A_source.tsv",
    "artifacts/Analysis B/GOMF/volcano/CAT_B_source.tsv")
  result <- extension_navigation_html(paths, units, source_run)
  section <- result$nav$contexts[[1L]]$collections[[1L]]$sections[[1L]]
  moved <- paste0(result$root, "-moved")
  expect_true(file.rename(result$root, moved))
  page <- file.path(moved, section$href)
  html <- paste(readLines(page, warn = FALSE), collapse = "\n")
  match <- regmatches(html, regexec('<script id="lisa-report-navigation" type="application/json">(.*?)</script>', html, perl = TRUE))[[1L]]
  nav <- jsonlite::fromJSON(match[[2L]], simplifyVector = FALSE)
  expect_identical(nav$mode, "external")
  expect_identical(nav$current$context, result$nav$contexts[[1L]]$id)
  expect_identical(nav$current$collection, "GOCC")
  expect_identical(nav$current$section, section$id)
  expect_identical(nav$contexts[[1L]]$collections[[1L]]$sections[[1L]]$href, "index.html")
  fallback <- nav$contexts[[2L]]$collections[[1L]]$sections[[1L]]
  expect_match(fallback$href, "^../../../index[.]html#extension-section-")
  for (context in nav$contexts) for (collection in context$collections) for (s in collection$sections) {
    target <- strsplit(s$href, "#", fixed = TRUE)[[1L]][[1L]]
    expect_true(file.exists(file.path(dirname(page), utils::URLdecode(target))))
    for (product in s$category_products) for (asset in product$assets)
      expect_true(file.exists(file.path(dirname(page), utils::URLdecode(asset$href))))
  }
})

test_that("referenced product names are URL-safe with exact source provenance", {
  unit <- extension_navigation_units()[1L, , drop = FALSE]
  paths <- paste0("artifacts/Analysis A/GOCC/member_gene_sets/biological label café #", strrep("long_", 12L),
    c("_source.tsv", ".png"))
  result <- extension_navigation_html(paths, unit, extension_saved_evidence())
  page <- file.path(result$root, result$nav$contexts[[1L]]$collections[[1L]]$sections[[1L]]$href)
  payload <- lisaR:::lisa_extension_read_payload(page, "evidence-data")
  assets <- payload$category_products[[1L]]$assets
  expect_setequal(vapply(assets, `[[`, character(1L), "source_path"), paths)
  for (asset in assets) {
    expect_match(asset$href, "^../../../artifacts/Analysis%20A/GOCC/member_gene_sets/")
    expect_match(asset$href, "caf%C3%A9%20%23", fixed = TRUE)
    expect_identical(asset$source_name, basename(asset$source_path))
    expect_identical(asset$sha256, lisaR:::lisa_sha256_file(file.path(result$root, asset$source_path)))
    expect_identical(asset$sha256, lisaR:::lisa_sha256_file(file.path(dirname(page), utils::URLdecode(asset$href))))
  }
  expect_false(dir.exists(file.path(dirname(page), "assets", "category_products")))
})

test_that("serialized category-product references fail closed after tampering", {
  source_run <- extension_saved_evidence()
  paths <- c("artifacts/Analysis A/GOCC/member_gene_sets/CAT_A.png",
    "artifacts/Analysis A/GOCC/member_gene_sets/CAT_A_source.tsv")
  result <- extension_navigation_html(paths, extension_navigation_units()[1L, , drop = FALSE], source_run)
  page <- file.path(result$root, "evidence", "Analysis A", "GOCC", "index.html")
  html <- paste(readLines(page, warn = FALSE), collapse = "\n")
  writeLines(sub("../../../artifacts/", "../../../other/", html, fixed = TRUE), page, useBytes = TRUE)
  inventory <- data.frame(path = paths,
    sha256 = vapply(file.path(result$root, paths), lisaR:::lisa_sha256_file, character(1L)),
    stringsAsFactors = FALSE)
  expected <- stats::setNames(list(paths), "evidence/Analysis A/GOCC/index.html")
  expect_error(lisaR:::lisa_extension_validate_product_references(result$root, expected, inventory),
    "failed reconciliation")
})

test_that("installed evidence viewers allow only legacy-local or canonical artifact paths", {
  node <- Sys.which("node")
  skip_if(!nzchar(node), "Node is optional; the viewer-path regression needs its JavaScript runtime")
  scripts <- c(system.file("category-evidence", "viewer.js", package = "lisaR"),
    system.file("contrast-evidence", "viewer.js", package = "lisaR"))
  skip_if(any(!nzchar(scripts)))
  probe <- tempfile(fileext = ".js")
  on.exit(unlink(probe), add = TRUE)
  writeLines(c(
    "const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');",
    "for(const file of process.argv.slice(2)){",
    " const source=fs.readFileSync(file,'utf8');",
    " const start=source.indexOf('  function safeCategoryAsset(path) {');",
    " const end=source.indexOf('\\n  function ',start+4);",
    " assert.ok(start>=0 && end>start,`safeCategoryAsset not found in ${file}`);",
    " const safe=vm.runInNewContext(`(()=>{${source.slice(start,end)};return safeCategoryAsset;})()`);",
    " for(const value of ['assets/category_products/legacy.png','../../../artifacts/Analysis%20A/GOCC/plot%20caf%C3%A9%20%23.png','../../../artifacts/a/%252F.png','../../../../outputs/gene_level/A/plot.png','../../../../artifacts/A/GOCC/plot.png']) assert.equal(safe(value),true,value);",
    " for(const value of ['../../../../../artifacts/a.png','../../../../other/a.png','../../../other/a.png','../../../artifacts/../outside.png','../../../artifacts/%2e%2e/outside.png','../../../artifacts/a%2Fb.png','assets/%2e%2e/outside.png','https://example.invalid/a.png','/artifacts/a.png','..\\\\..\\\\..\\\\artifacts\\\\a.png','../../../artifacts/%','../../../artifacts/ab.png']) assert.equal(safe(value),false,value);",
    "}",
    "process.stdout.write('PASS canonical category-product path gates\\n');"
  ), probe, useBytes = TRUE)
  result <- system2(node, c(probe, scripts), stdout = TRUE, stderr = TRUE)
  expect_identical(attr(result, "status"), NULL, info = paste(result, collapse = "\n"))
  expect_match(paste(result, collapse = "\n"), "PASS canonical category-product path gates", fixed = TRUE)
})

test_that("KEGG catalog discovery uses the saved source significance cutoff", {
  root <- tempfile("extension-kegg-cutoff-")
  collection <- file.path(root, "outputs", "single_de", "A", "collection_PATHWAYS")
  for (part in c("lisa_tables", "enrichment", "inputs")) dir.create(file.path(collection, part), recursive = TRUE)
  dir.create(file.path(root, "config"))
  lisaR:::write_lisa_tsv(data.frame(analysis_id = "A"), file.path(root, "config", "de_index.tsv"))
  categories <- data.frame(category_id = c("CAT_LOW", "CAT_MID", "CAT_HIGH"), macrogroup_id = "SUPER", n_genesets = 1L)
  lisaR:::write_lisa_tsv(categories, file.path(collection, "lisa_tables", "A_GSEA_category_summary.tsv"))
  gsea <- data.frame(pathway = c("KEGG_LOW", "KEGG_MID", "KEGG_HIGH"), padj = c(.01, .20, .30),
    category_id = categories$category_id, macrogroup_id = "SUPER")
  lisaR:::write_lisa_tsv(gsea, file.path(collection, "enrichment", "A_GSEA_semantic_annotated.tsv"))
  lisaR:::write_lisa_tsv(data.frame(symbol = "GENE"), file.path(collection, "inputs", "A_standardized_DE.tsv"))
  for (cutoff in c(.25, .05)) {
    lisaR:::write_lisa_tsv(data.frame(key = "gsea_padj_cutoff", value = cutoff), file.path(root, "contract_manifest.tsv"))
    catalog <- lisaR:::lisa_extension_discover_catalog(root)
    expect_setequal(catalog$category_id[catalog$product == "kegg"],
      if (cutoff == .25) c("CAT_LOW", "CAT_MID") else "CAT_LOW")
    expect_setequal(catalog$category_id[catalog$product == "member_gene_sets"], categories$category_id)
  }
})

test_that("heatmap display-scale suffix maps to its exact matrix and recipe family", {
  unit <- extension_navigation_units()[1L, , drop = FALSE]
  unit$product <- "heatmap"
  base <- "artifacts/Analysis A/GOCC/heatmap/A_semantic_CAT_A_leading_edge_gene_heatmap"
  files <- paste0(base, c("_zscore.png", "_zscore.svg", "_matrix.tsv", "_recipe.R"))
  result <- extension_navigation_html(files, unit)
  expect_identical(length(gregexpr('class="lisa-extension-card"', result$html, fixed = TRUE)[[1L]]), 1L)
  expect_match(result$html, 'data-category="CAT_A" data-supercategory="SUPER_1"', fixed = TRUE)
  expect_match(result$html, "Matrix TSV", fixed = TRUE)
  expect_match(result$html, "Recipe R", fixed = TRUE)
  expect_identical(length(unique(lisaR:::lisa_extension_figure_family(files))), 1L)
})

test_that("table-only and PDF-only full extensions never invent inline images", {
  unit <- extension_navigation_units()[1L, , drop = FALSE]
  base <- "artifacts/Analysis A/GOCC/member_gene_sets/CAT_A"
  table <- extension_navigation_html(paste0(base, "_source.tsv"), unit)
  expect_false(grepl('<article[^>]*>.*?<img', table$html, perl = TRUE))
  expect_match(table$html, "No inline image was requested", fixed = TRUE)
  expect_match(table$html, "Source TSV", fixed = TRUE)
  pdf <- extension_navigation_html(paste0(base, ".pdf"), unit)
  expect_false(grepl('<img[^>]*[.]pdf', pdf$html, perl = TRUE))
  expect_match(pdf$html, 'CAT_A.pdf">CAT_A.pdf</a>', fixed = TRUE)
  misc <- extension_navigation_html("artifacts/Analysis A/GOCC/member_gene_sets/extra.tsv", unit)
  expect_match(misc$html, "Unclassified saved files (1)", fixed = TRUE)
  expect_match(misc$html, "No category product had an exact source metadata attachment", fixed = TRUE)
})

test_that("extension index never assigns a category from a filename without exact source metadata", {
  unit <- extension_navigation_units()[1L, , drop = FALSE]
  result <- extension_navigation_html(
    "artifacts/Analysis A/GOCC/member_gene_sets/CAT_A.png", unit
  )
  expect_false(grepl('data-category="CAT_A"', result$html, fixed = TRUE))
  expect_match(result$html, "Unclassified saved files (1)", fixed = TRUE)
})

test_that("extension galleries are paginated and navigation labels cannot terminate the data script", {
  unit <- extension_navigation_units()[1L, , drop = FALSE]
  unit$analysis_id <- 'A </script><script>alert(1)</script>'
  paths <- character()
  result <- extension_navigation_html(paths, unit)
  expect_match(result$html, "data-extension-search", fixed = TRUE)
  expect_match(result$html, "data-extension-prev", fixed = TRUE)
  expect_match(result$html, "data-extension-next", fixed = TRUE)
  expect_match(result$html, "size=12", fixed = TRUE)
  nav_text <- paste(readLines(file.path(result$root, "navigation_inventory.json"), warn = FALSE), collapse = "")
  expect_false(grepl("</script>", nav_text, fixed = TRUE))
  expect_identical(result$nav$contexts[[1L]]$scientific_id, unit$analysis_id)
})

test_that("extension index rejects traversal and preserves unassigned saved artifacts", {
  expect_error(extension_navigation_html("artifacts/../private.png"), "contained relative")
  expect_error(extension_navigation_html("/outside.png"), "contained relative")
  expect_error(extension_navigation_html(c("artifacts/x.png", "artifacts/x.png")), "contained relative")
  result <- extension_navigation_html("artifacts/legacy/plot #1.png", data.frame())
  expect_identical(result$nav$contexts[[1L]]$label, "Saved results")
  expect_match(result$html, 'href="artifacts/legacy/plot%20%231.png"', fixed = TRUE)
  expect_false(grepl("%2F", result$html, fixed = TRUE))
})
