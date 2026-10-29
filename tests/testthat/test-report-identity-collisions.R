identity_report_environment <- function() {
  e <- new.env(parent = globalenv())
  script <- Sys.getenv("LISAR_IDENTITY_TEST_SCRIPT", system.file("scripts/build_LISA_report.R", package = "lisaR"))
  sys.source(script, e)
  e$report_package_dir <- system.file(package = "lisaR")
  e
}

identity_real_paths <- function() {
  paste0("outputs/gene_level/category_contrasts/CELLS_03_B16_cells_LoF_sh2_vs_sh3/collection_PATHWAYS/contrast_category_cards/",
    c("20_PAK_SIGNALING", "51_ALK_SIGNALING"), "_contrast_category_card")
}

test_that("report identity separates known linear-checksum collisions deterministically", {
  e <- identity_report_environment()
  pairs <- list(c("CAT_ACA", "CAT_BAB"), paste0(identity_real_paths(), "_recipe.R"),
    paste0(identity_real_paths(), ".png"))
  for (pair in pairs) {
    expect_identical(e$report_legacy_text_hash(pair[[1L]]), e$report_legacy_text_hash(pair[[2L]]))
    hashes <- vapply(pair, e$stable_text_hash, character(1L))
    expect_equal(length(unique(hashes)), 2L)
    expect_true(all(grepl("^[0-9a-f]{64}$", hashes)))
    expect_identical(hashes, vapply(pair, e$stable_text_hash, character(1L)))
    expect_identical(unname(hashes), vapply(pair, digest::digest, character(1L), algo = "sha256", serialize = FALSE, USE.NAMES = FALSE))
  }
  root_a <- tempfile(); root_b <- tempfile()
  dir.create(root_a); dir.create(root_b)
  on.exit(unlink(c(root_a, root_b), recursive = TRUE), add = TRUE)
  expect_identical(e$figure_id_for_path(file.path(root_a, "outputs", "caf\u00e9.png"), root_a),
    e$figure_id_for_path(file.path(root_b, "outputs", "caf\u00e9.png"), root_b))
  expect_identical(e$stable_text_hash(""), digest::digest("", algo = "sha256", serialize = FALSE))
  expect_error(e$stable_text_hash(c("a", "b")), "LISA-REPORT-ID-001")
  expect_error(e$stable_text_hash(NA_character_), "LISA-REPORT-ID-001")
})

test_that("an identity collision cannot reuse an existing cache or contract", {
  e <- identity_report_environment()
  e$report_text_sha256 <- function(text) paste(rep("a", 64L), collapse = "")
  expect_match(e$stable_text_hash("first"), "^[a-f0-9]{64}$")
  expect_identical(e$stable_text_hash("first"), e$stable_text_hash("first"))
  expect_error(e$stable_text_hash("different"), "LISA-REPORT-ID-002")
})

test_that("colliding native contrast names retain distinct recipes and scientific source rows", {
  e <- identity_report_environment()
  root <- tempfile("report-native-collision-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  owner <- "CELLS_03_B16_cells_LoF_sh2_vs_sh3"
  evidence <- file.path(root, "report_pages/contrast_evidence", owner, "PATHWAYS")
  dir.create(file.path(evidence, "tables"), recursive = TRUE)
  page <- file.path(evidence, "index.html")
  writeLines('<script type="application/json" id="contrast-evidence-data">{"categories":[]}</script>', page)
  lisaR:::write_lisa_tsv(data.frame(key = c("contrast_id", "collection"), value = c("C", "PATHWAYS")), file.path(evidence, "tables/metadata.tsv"))
  ids <- c("PAK_SIGNALING", "ALK_SIGNALING")
  lisaR:::write_lisa_tsv(data.frame(category_id = ids), file.path(evidence, "tables/categories.tsv"))
  stems <- file.path(root, identity_real_paths())
  dir.create(dirname(stems[[1L]]), recursive = TRUE)
  for (i in seq_along(stems)) {
    writeBin(as.raw(c(137, 80, 78, 71)), paste0(stems[[i]], ".png"))
    lisaR:::write_lisa_tsv(data.frame(category_id = ids[[i]], symbol = paste0("GENE_", i),
      figure_type = "contrast_gene_card", log2FC_A = i, log2FC_B = -i, plotted_order = 1L), paste0(stems[[i]], "_source.tsv"))
    writeLines("# immutable original recipe", paste0(stems[[i]], "_recipe.R"))
  }
  originals <- list.files(dirname(stems[[1L]]), full.names = TRUE)
  before <- vapply(originals, lisaR:::lisa_sha256_file, character(1L))
  attachment <- e$report_attach_category_products(page,
    list(list(product = "contrast_gene_cards", label = "Prioritized genes (paired)", path = dirname(stems[[1L]]))),
    "C", "PATHWAYS", contrast = TRUE, project_dir = root)
  refs <- e$report_category_product_reference_rows(root)
  expect_equal(nrow(refs), 6L)
  recipes <- refs[refs$format == "r", , drop = FALSE]
  expect_equal(length(unique(recipes$source_path)), 2L)
  expect_true(all(grepl("/F_[0-9a-f]{64}_recipe[.]R$", recipes$source_path)))
  expect_identical(vapply(originals, lisaR:::lisa_sha256_file, character(1L)), before)
  # The same collision also affected figure_source_contract's existence
  # shortcut, silently attaching the first category's TSV to the second one.
  main_page <- file.path(root, "report_pages/contrasts.html")
  contracts <- lapply(paste0(stems, ".png"), e$figure_source_contract, page_file = main_page)
  expect_equal(length(unique(vapply(contracts, `[[`, character(1L), "source_tsv"))), 2L)
  for (i in seq_along(contracts)) {
    data <- read.delim(file.path(dirname(main_page), contracts[[i]]$source_tsv), check.names = FALSE)
    expect_identical(data$category_id, ids[[i]])
    expect_identical(data$symbol, paste0("GENE_", i))
    expect_equal(data$log2FC_A, i)
    expect_equal(data$log2FC_B, -i)
  }
  expect_equal(nrow(e$validate_figure_source_contracts(root)), 2L)
  # Preserve the original fail-closed duplicate rule, including the same
  # asset repeated within one product and a product duplicated in the page.
  original_payload <- e$report_embedded_json_payload(page, "contrast-evidence-data")
  bad <- original_payload
  bad$category_products[[1L]]$assets <- c(bad$category_products[[1L]]$assets, bad$category_products[[1L]]$assets[1L])
  e$report_replace_evidence_payload(page, bad$category_products, "contrast-evidence-data")
  expect_error(e$report_category_product_reference_rows(root), "LISA-REPORT-PRODUCT-004")
  e$report_replace_evidence_payload(page, c(original_payload$category_products, original_payload$category_products[1L]), "contrast-evidence-data")
  expect_error(e$report_category_product_reference_rows(root), "LISA-REPORT-PRODUCT-004")
  e$report_replace_evidence_payload(page, original_payload$category_products, "contrast-evidence-data")
  manifest <- data.frame(relative_path = refs$source_path, sha256 = refs$sha256)
  expect_silent(e$report_validate_category_product_references(root, manifest))
  bad <- original_payload
  bad$category_products[[1L]]$assets[[1L]]$sha256 <- paste(rep("0", 64L), collapse = "")
  e$report_replace_evidence_payload(page, bad$category_products, "contrast-evidence-data")
  expect_error(e$report_validate_category_product_references(root, manifest), "LISA-REPORT-PRODUCT-005")
  e$report_replace_evidence_payload(page, original_payload$category_products, "contrast-evidence-data")
  writeLines("corrupted after manifest", file.path(root, refs$source_path[[1L]]))
  expect_error(e$report_validate_category_product_references(root, manifest), "LISA-REPORT-PRODUCT-005")
})

test_that("legacy NES markup refresh binds new SHA-256 bundles without changing old assets", {
  e <- identity_report_environment()
  root <- tempfile("legacy-identity-refresh-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  collection <- file.path(root, "outputs/category_contrasts/C_O/collection_GOCC")
  folder <- file.path(collection, "plots/category_nes")
  dir.create(folder, recursive = TRUE)
  dir.create(file.path(collection, "qc"))
  dir.create(file.path(root, "config"), recursive = TRUE)
  evidence <- file.path(root, "report_pages/contrast_evidence/C_O/GOCC/index.html")
  dir.create(dirname(evidence), recursive = TRUE)
  writeLines("<html><body>Saved evidence</body></html>", evidence)
  lisaR:::write_lisa_tsv(data.frame(contrast_id = "C", output_id = "O"), file.path(root, "config/contrast_index.tsv"))
  source <- file.path(folder, "nes_source.tsv")
  writeLines(c("category_id\tmean_NES", "CAT_A\t1.5"), source)
  figure <- file.path(folder, "nes.png")
  writeBin(as.raw(c(137, 80, 78, 71)), figure)
  settings <- file.path(folder, "nes_settings.json")
  jsonlite::write_json(list(variant = "clean", source_sha256 = digest::digest(file = source, algo = "sha256")), settings, auto_unbox = TRUE)
  recipe <- file.path(folder, "nes_recipe.R")
  writeLines("# untrusted runtime recipe, never used", recipe)
  manifest <- file.path(collection, "qc/C_category_nes_variants.tsv")
  lisaR:::write_lisa_tsv(data.frame(variant = "clean", status = "rendered", format = "png",
    path = "plots/category_nes/nes.png", source = "plots/category_nes/nes_source.tsv",
    settings = "plots/category_nes/nes_settings.json", recipe = "plots/category_nes/nes_recipe.R"), manifest)
  e$report_output_policy <- data.frame(product = c("category_evidence"), requested = FALSE)
  page <- file.path(root, "report_pages/contrasts.html")
  current_hash <- e$stable_text_hash
  e$stable_text_hash <- e$report_legacy_text_hash
  legacy <- e$report_nes_variants_section(collection, page, "Legacy")
  expect_match(legacy, 'data-nes-variants="N_[0-9a-f]{16}"')
  writeLines(paste0('<html><body>', legacy,
    '<section class="panel contrast-evidence-cover" id="cover"><a href="', e$rel_path(evidence,page), '">Evidence</a></section></body></html>'), page)
  writeLines('<html><body><a href="report_pages/contrasts.html">Contrasts</a></body></html>', file.path(root, "report_index.html"))
  e$write_shareable_report_manifest(root)
  expect_silent(e$validate_shareable_report_manifest(root))
  old_assets <- c(list.files(file.path(root, "outputs"), recursive = TRUE, full.names = TRUE),
    list.files(file.path(root, "report_figure_data/nes_variants"), recursive = TRUE, full.names = TRUE))
  before <- vapply(old_assets, lisaR:::lisa_sha256_file, character(1L))
  e$stable_text_hash <- current_hash
  expect_message(e$report_refresh_category_links(root), "CATEGORY_LINK_REFRESH_PASS")
  expect_silent(e$validate_shareable_report_manifest(root))
  updated <- paste(readLines(page, warn = FALSE), collapse = "\n")
  expect_match(updated, 'data-nes-variants="N_[0-9a-f]{64}"')
  expect_identical(vapply(old_assets, lisaR:::lisa_sha256_file, character(1L)), before)
  rows <- e$read_shareable_report_manifest(root)
  bundles <- rows$relative_path[grepl("^report_figure_data/nes_variants/N_[0-9a-f]{64}/", rows$relative_path)]
  expect_gte(length(bundles), 3L)
  expect_true(all(file.exists(file.path(root, bundles))))
})
