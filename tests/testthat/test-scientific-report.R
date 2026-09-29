scientific_report_fixture <- function(root) {
  unlink(root, recursive = TRUE)
  dir.create(file.path(root, "d"), recursive = TRUE)
  dir.create(file.path(root, "g"))
  dir.create(file.path(root, "m"))
  tab <- data.frame(category_id = c("P01", "P02"), symbol = c("STAT3", "COL1A1"), score = c(2, 1))
  write.table(tab, file.path(root, "d", "T0001.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  png <- as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
  writeBin(png, file.path(root, "g", "F0001.png"))
  writeLines('<svg xmlns="http://www.w3.org/2000/svg"><rect width="10" height="10"/></svg>', file.path(root, "g", "F0002.svg"))
  roles <- c(
    "outputs/single_de/A/collection_PATHWAYS/lisa_tables/A_semantic_GSEA_category_summary.tsv",
    "outputs/single_de/A/collection_PATHWAYS/plots/A_semantic_GSEA_category_card.png",
    "outputs/single_de/A/collection_PATHWAYS/plots/A_semantic_GSEA_category_card.svg"
  )
  paths <- c("d/T0001.tsv", "g/F0001.png", "g/F0002.svg")
  reg <- data.frame(table_id = c("T0001", "F0001", "F0002"), path = paths,
    artifact_type = c("table", "figure", "figure"), role = roles,
    analysis_id = "A", contrast_id = "", collection = "PATHWAYS",
    sha256 = vapply(file.path(root, paths), lisaR:::lisa_report_package_sha256, character(1)),
    source_count = 1L, stringsAsFactors = FALSE)
  write.table(reg, file.path(root, "m", "registry.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(analysis_id = "A", label = "Secretome GoF"), file.path(root, "m", "de.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(contrast_id = character(), contrast_label = character()), file.path(root, "m", "cx.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(analysis_collection = "PATHWAYS"), file.path(root, "m", "col.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(source_path = roles, role = c("table", "figure", "figure"), source_sha256 = reg$sha256), file.path(root, "m", "prov.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  files <- sort(list.files(root, recursive = TRUE, full.names = FALSE))
  write.table(data.frame(path = files, sha256 = vapply(file.path(root, files), lisaR:::lisa_report_package_sha256, character(1))), file.path(root, "m", "sha256.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  root
}

test_that("scientific report renders a clean production-style multipage bundle", {
  fx <- scientific_report_fixture(file.path(tempdir(), "scientific_report_fixture"))
  out <- file.path(tempdir(), "scientific_report_output"); unlink(out, recursive = TRUE)
  result <- build_lisa_scientific_report(fx, out, "Treatment pilot", "Treatment/Control")
  expect_true(file.exists(result$index))
  expect_true(lisaR:::lisa_scientific_report_root_layout_ok(out))
  expect_equal(result$cards, 1L)
  expect_equal(result$figures, 2L)
  expect_equal(result$broken_links, 0L)
  expect_lte(result$max_path, 120L)
  expect_lte(result$max_html, 1024^2)
  html <- paste(readLines(list.files(file.path(out, "pages"), pattern = "^s.*html$", full.names = TRUE)[[1]], warn = FALSE), collapse = "\n")
  expect_match(html, "fig-card", fixed = TRUE)
  expect_match(html, "Source data", fixed = TRUE)
  expect_match(html, "assets/report.css", fixed = TRUE)
  expect_match(html, "source_data/single_de/a/pathways/T0001_semantic_GSEA_category_summary.tsv", fixed = TRUE)
  expect_false(grepl("Other visuals|Offline portable report|R-native renderer|Checksummed source fixture", html))
  expect_true(file.exists(file.path(out, "README.md")))
  expect_true(file.exists(file.path(out, "PROJECT.md")))
  expect_true(file.exists(file.path(out, "assets", "logo.svg")))
  expect_true(file.exists(file.path(out, "assets", "lisa-logo-rectangular.svg")))
  index_html <- paste(readLines(file.path(out, "index.html"), warn = FALSE), collapse = "\n")
  expect_match(index_html, 'class="hero-logo"', fixed = TRUE)
  expect_match(index_html, 'src="assets/lisa-logo-rectangular.svg"', fixed = TRUE)
  expect_false(grepl("lisa-logo-rectangular.svg", html, fixed = TRUE))
  css <- paste(readLines(file.path(out, "assets", "report.css"), warn = FALSE), collapse = "\n")
  expect_match(css, ".brand img{width:84px;height:84px}", fixed = TRUE)
  expect_match(css, ".hero-logo{", fixed = TRUE)
  expect_true(file.exists(file.path(out, "source_data", "INDEX.tsv")))
  expect_length(lisaR:::lisa_scientific_report_verify_links(out), 0L)
  cards <- read_lisa_tsv(file.path(out, "metadata", "figure_cards.tsv"))
  expect_equal(cards$source_table_id, "T0001")
  expect_true(nzchar(cards$png_path) && nzchar(cards$svg_path))
  manifest <- read_lisa_tsv(file.path(out, "metadata", "manifest.tsv"))
  expect_true(all(vapply(seq_len(nrow(manifest)), function(i) {
    identical(lisaR:::lisa_report_package_sha256(file.path(out, manifest$path[[i]])), manifest$sha256[[i]])
  }, logical(1))))
})

test_that("scientific report packaged logos are the approved compact and rectangular assets", {
  compact <- system.file("report_assets", "LISA_logo_C_compact_icon_muted_red_S.svg", package = "lisaR")
  rectangular <- system.file("report_assets", "LISA_logo_A1_muted_red_S_automated_annotation_final.svg", package = "lisaR")
  expect_true(file.exists(compact))
  expect_true(file.exists(rectangular))
  expect_match(paste(readLines(compact, warn = FALSE), collapse = "\n"), "viewBox=\"0 0 800 800\"", fixed = TRUE)
  # The approved rectangular asset is the CROPPED artwork: viewBox "40 70 1090
  # 340". The old "0 0 1600 500" expectation described the pre-crop draft and
  # was never updated. Provenance: this file is byte-unchanged since the
  # approved commit 0035cc3 (verified with `git diff 0035cc3 -- <path>`), so the
  # test is pinned to the intended asset, not to a fresh visual approval.
  expect_match(paste(readLines(rectangular, warn = FALSE), collapse = "\n"), "viewBox=\"40 70 1090 340\"", fixed = TRUE)
  # Pin both assets by content so a silent substitution fails here even if a
  # replacement happens to keep the same viewBox attribute.
  expect_identical(
    lisaR:::lisa_sha256_file(compact),
    "9fb546d321916b5ef9fc82bd59fb58e500d08455c9d9a2d9f8c378c1794b58d2"
  )
  expect_identical(
    lisaR:::lisa_sha256_file(rectangular),
    "8f9e41768418db601a4c2184312466ced6f508d0b737be52c63ef538a58c167e"
  )
})

test_that("leading-edge heatmaps preserve input sample-column order without clustering", {
  path <- system.file("scripts", "build_single_de_leading_edge_gene_heatmaps.R", package = "lisaR")
  expect_true(file.exists(path))
  code <- paste(readLines(path, warn = FALSE), collapse = "\n")
  expect_match(code, "lapply(sample_cols, function(sample)", fixed = TRUE)
  expect_match(code, "factor(long$sample, levels = sample_cols)", fixed = TRUE)
  expect_match(code, "match(long$sample, sample_cols)", fixed = TRUE)
  expect_false(grepl("cluster_cols", code, fixed = TRUE))
  expect_false(grepl("hclust(", code, fixed = TRUE))
  expect_match(code, "displayed in their input-matrix order", fixed = TRUE)
})

test_that("scientific report fails closed for incomplete fixtures and embeds PDF-only cards", {
  bad <- tempfile("scientific_report_invalid-"); dir.create(file.path(bad, "m"), recursive = TRUE)
  expect_error(build_lisa_scientific_report(bad, tempfile("scientific_report_output-")), "incomplete fixture")
  fx <- scientific_report_fixture(file.path(tempdir(), "scientific_report_fixture-pdf"))
  reg <- read_lisa_tsv(file.path(fx, "m", "registry.tsv"))
  file.rename(file.path(fx, "g", "F0001.png"), file.path(fx, "g", "F0001.pdf"))
  reg$path[reg$table_id == "F0001"] <- "g/F0001.pdf"
  reg$role[reg$table_id == "F0001"] <- sub("png$", "pdf", reg$role[reg$table_id == "F0001"])
  reg$sha256[reg$table_id == "F0001"] <- lisaR:::lisa_report_package_sha256(file.path(fx, "g", "F0001.pdf"))
  write.table(reg, file.path(fx, "m", "registry.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  unlink(file.path(fx, "g", "F0002.svg")); reg <- reg[reg$table_id != "F0002", , drop = FALSE]
  write.table(reg, file.path(fx, "m", "registry.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  files <- sort(setdiff(list.files(fx, recursive = TRUE, full.names = FALSE), "m/sha256.tsv"))
  write.table(data.frame(path = files, sha256 = vapply(file.path(fx, files), lisaR:::lisa_report_package_sha256, character(1))), file.path(fx, "m", "sha256.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  out <- tempfile("scientific_report_output-")
  result <- build_lisa_scientific_report(fx, out)
  expect_equal(result$cards, 1L)
  html <- paste(readLines(list.files(file.path(out, "pages"), pattern = "^s.*html$", full.names = TRUE)[[1]], warn = FALSE), collapse = "\n")
  expect_match(html, "application/pdf", fixed = TRUE)
})

test_that("scientific report source-data priorities remain biologically aligned", {
  roles <- c("semantic_GSEA_category_summary.tsv", "kegg_pathway_painter_index.tsv",
             "gene_category_contributions.tsv", "category_gene_cards_index.tsv",
             "recurrent_gene_screen.tsv")
  gene <- lisaR:::lisa_scientific_report_source_score(roles, "gene_prioritization")
  expect_gt(gene[[3]], gene[[4]])
  recurrent <- lisaR:::lisa_scientific_report_source_score(roles, "recurrent_genes")
  expect_gt(recurrent[[5]], recurrent[[3]])
  kegg <- lisaR:::lisa_scientific_report_source_score(roles, "kegg_painted_maps")
  expect_gt(kegg[[2]], kegg[[1]])
  heat <- lisaR:::lisa_scientific_report_source_score(roles, "supporting_gene_heatmaps")
  expect_gt(heat[[1]], heat[[5]])
})

test_that("scientific report.1 assigns every production figure family explicitly and has no Other bucket", {
  cases <- c(
    "x_GSEA_lollipop.png" = "lisa_summary",
    "x_GSEA_lollipop_direction_stats.png" = "lisa_summary",
    "x_GSEA_pathway_dotplot.png" = "lisa_summary",
    "x_ORA_barplot.png" = "lisa_summary",
    "plots/x_GSEA_category_pathways/CAT.png" = "lisa_category_gene_sets",
    "category_gene_cards/01_CAT_category_card.png" = "gene_prioritization",
    "category_volcano_overlays/01_CAT_volcano.png" = "volcano_overlays",
    "recurrent_gene_screen/top_recurrent_genes.png" = "recurrent_genes",
    "leading_edge_gene_heatmaps/CAT.png" = "supporting_gene_heatmaps",
    "kegg_painter/01_mmu00010_x_painted.png" = "kegg_painted_maps",
    "contrast_dumbbell/x.png" = "category_shifts",
    "contrast_category_cards/01_CAT_contrast_category_card.png" = "contrast_gene_cards",
    "paired_gene_heatmap.png" = "paired_gene_heatmaps",
    "gene_category_network/network.png" = "gene_category_networks",
    "contrast_kegg_pathway_painter/01_mmu00010_x_contrast_painted.png" = "contrast_kegg_painted_maps"
  )
  got <- vapply(names(cases), function(x) lisaR:::lisa_scientific_report_classify(x, "figure", x), character(1))
  expect_identical(unname(got), unname(cases))
  expect_false(any(grepl("other", unname(got), ignore.case = TRUE)))
  expect_true(all(lisa_builtin_collection_registry()$run_enrichmentmap == FALSE))
})

test_that("scientific report.2 separates within-category prioritization from recurrence", {
  expect_identical(lisaR:::lisa_scientific_report_layer_label("gene_prioritization"), "Gene prioritization")
  expect_identical(lisaR:::lisa_scientific_report_layer_label("recurrent_genes"), "Cross-category recurrent genes")
  expect_match(lisaR:::lisa_scientific_report_layer_note("gene_prioritization"),
               "sqrt\\(LISA support\\).*log2FC.*DE FDR")
  expect_match(lisaR:::lisa_scientific_report_layer_note("recurrent_genes"),
               "separately from within-category gene prioritization")
})

test_that("scientific report omits redundant contrast gene prioritization pages", {
  expect_false("contrast_gene_prioritization" %in% c(
    "category_shifts", "contrast_gene_cards", "paired_gene_heatmaps",
    "gene_category_networks", "contrast_kegg_painted_maps"
  ))
  renderer <- paste(c(
    deparse(body(lisaR:::build_lisa_scientific_report)),
    deparse(body(lisaR:::lisa_scientific_report_layer_label)),
    deparse(body(lisaR:::lisa_scientific_report_make_nav)),
    deparse(body(lisaR:::lisa_scientific_report_source_score))
  ), collapse = "\n")
  expect_false(grepl("Contrast gene prioritization", renderer, fixed = TRUE))
  expect_false(grepl('add_table_page("contrast", id, collection, "contrast_gene_prioritization"', renderer, fixed = TRUE))
  expect_match(renderer, "contrast_gene_cards = c(paired_gene_evidence", fixed = TRUE)
})
