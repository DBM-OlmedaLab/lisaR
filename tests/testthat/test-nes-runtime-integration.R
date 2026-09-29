test_that("report NES variants default to four and reject invalid selections", {
  all_variants <- c("clean", "percentages", "direction", "dispersion")
  expect_identical(lisaR:::lisa_validate_report_config(list())$category_nes_variants, all_variants)
  for (value in list("percentages", "dispersion", c("direction", "clean"), list("direction", "clean"))) {
    got <- lisaR:::lisa_validate_report_config(list(report = list(category_nes_variants = value)))
    expect_identical(got$category_nes_variants, unlist(value, use.names = FALSE))
  }
  invalid <- list(NULL, list(), character(), c("clean", "clean"),
    c("clean", NA_character_), "all", "", 1, TRUE, list(clean = "clean"), list(c("clean", "direction")))
  for (value in invalid) {
    expect_error(lisaR:::lisa_validate_report_config(list(report = list(category_nes_variants = value))),
      "LISA-CONFIG-NES-VARIANTS-001")
  }
})

test_that("deep Riaz view downloads stay portable without truncating scientific identity", {
  stem <- lisaR:::lisa_nes_runtime_stem
  id <- "response_separation_by_time"
  keys <- c(stem(id, "RIAZ_GSE91061", TRUE), stem(paste0(id, "_other"), "RIAZ_GSE91061", TRUE),
    stem(id, "other prefix", TRUE), stem(id, "RIAZ_GSE91061", FALSE))
  expect_identical(anyDuplicated(keys), 0L)
  expect_true(all(nchar(keys) <= 37L))
  folder <- paste0("D:/sample-inbox/lisaR_1234567/evaluation/R/outputs/category_contrasts/",
    id, "_", id, "/collection_PATHWAYS")
  paths <- sub("^file:", "", lisaR:::lisa_nes_expected_witnesses(folder, id, "RIAZ_GSE91061",
    c("clean", "percentages", "direction", "dispersion"), c("png", "svg", "pdf"), TRUE))
  expect_true(all(nchar(paths, type = "chars") < 260L))
  expect_true(any(endsWith(paths, "_opposite_direction_dispersion_settings.json")))
})

test_that("one NES variant survives normalized JSON as an array", {
  cfg <- list(pipeline = list(schema_version = "1.0.0", profile = "targeted",
    evidence_mode = "full_de", duplicate_policies = list(de_table_duplicate_policy = "error",
      matrix_duplicate_policy = "error", mapped_id_collision_policy = "error")),
    single_de = list(list(analysis_id = "A", de_path = "a.tsv", species = "Homo sapiens")),
    collections = list("GOBP-C2"), report = list(category_nes_variants = list("dispersion")))
  resolved <- lisaR:::lisa_resolve_config(cfg)
  expect_identical(resolved$contract$report$category_nes_variants, "dispersion")
  expect_identical(resolved$cfg$report$category_nes_variants, list("dispersion"))
  path <- tempfile(fileext = ".json")
  lisaR:::lisa_write_normalized_config(resolved$cfg, path)
  again <- jsonlite::read_json(path, simplifyVector = FALSE)
  expect_identical(again$report$category_nes_variants, list("dispersion"))
  row <- resolved$provenance[resolved$provenance$field == "report.category_nes_variants", ]
  expect_identical(row$value, "dispersion")
  expect_identical(row$source, "explicit")
})

nes_runtime_fixture <- function() {
  root <- tempfile("nes-runtime-"); dir.create(root)
  map <- data.frame(category_id = c("SIGNAL", "CONTEXT", "ONLY_B", "EMPTY"),
    display_name = c("Signal", "Context", "Only B", "Empty"),
    macrogroup_id = "SUPER", macrogroup_name = "Shared supercategory",
    macrogroup_order = 1, category_order_within_macrogroup = 1:4,
    color = "#2166ac", stringsAsFactors = FALSE)
  # ONLY_B has significant support only in B, not an absent/unmeasured A.
  # The existing contrast contract requires a finite contextual counterpart
  # whenever either side is significant; it never substitutes zero for absence.
  a <- data.frame(category_id = c("SIGNAL", "SIGNAL", "CONTEXT", "ONLY_B"),
    pathway = c("001", "s2", "cx", "only"), NES = c(2, -1, .5, -.5),
    padj = c(.01, .02, .7, .7), universe = "GOBP-C2", stringsAsFactors = FALSE)
  b <- data.frame(category_id = c("SIGNAL", "CONTEXT", "ONLY_B"),
    pathway = c("001", "cx", "only"), NES = c(-2, 2, 1),
    padj = c(.01, .01, .01), universe = "GOBP-C2", stringsAsFactors = FALSE)
  write_side <- function(id, members, keep_map) {
    dirs <- lisaR:::make_output_dirs(file.path(root, id))
    summary <- lisaR:::build_lisa_summaries(members, data.frame(), keep_map,
      include_empty_categories = TRUE, gsea_padj_cutoff = .05, ora_padj_cutoff = .05)$gsea
    path <- file.path(dirs$lisa_tables, paste0(id, "_semantic_GSEA_category_summary.tsv"))
    lisaR:::write_tsv_local(summary, path)
    lisaR:::write_tsv_local(members, file.path(dirs$enrichment, paste0(id, "_GSEA_semantic_annotated.tsv")))
    list(summary = summary, path = path, members = members)
  }
  list(root = root, a = write_side("A", a, map),
    b = write_side("B", b, map))
}

test_that("exact member discovery rejects ambiguous provenance and preserves string IDs", {
  fixture <- nes_runtime_fixture()
  got <- lisaR:::lisa_nes_members_for_summary(fixture$a$path, "GOBP-C2")
  expect_identical(got$pathway[[1L]], "001")
  expect_error(lisaR:::lisa_nes_members_for_summary(fixture$a$path, "GOMF"), "collection")
  standalone <- file.path(fixture$root, "summary.tsv")
  file.copy(fixture$a$path, standalone)
  expect_null(lisaR:::lisa_nes_members_for_summary(standalone, "GOBP-C2"))
})

test_that("single runtime selection writes only chosen variants without changing summaries", {
  skip_if_not_installed("ggplot2")
  fixture <- nes_runtime_fixture()
  before <- fixture$a$summary
  prepared <- lisaR:::lisa_prepare_category_nes_variants(before, fixture$a$members, .05)
  dirs <- lisaR:::make_output_dirs(file.path(fixture$root, "single-views"))
  got <- lisaR:::lisa_write_runtime_nes_variants(prepared, dirs, "A", "GOBP-C2", "semantic",
    "direction", "png")
  expect_identical(got$variant, "direction")
  expect_identical(got$status, "rendered")
  canonical_recipe <- lisaR:::lisa_post_script_path(lisaR:::lisa_resolve_package_dir(),
    "reproduce_lisa_category_nes.R")
  expect_identical(lisaR:::lisa_sha256_file(file.path(dirname(dirs$plots), got$recipe)),
    lisaR:::lisa_sha256_file(canonical_recipe))
  expect_true(all(file.exists(file.path(dirname(dirs$plots), unlist(got[c("path", "source", "settings", "recipe")])))))
  expect_false(any(grepl("_(clean|dispersion)\\.", list.files(file.path(dirs$plots, "category_nes")))))
  expect_identical(fixture$a$summary, before)
  settings <- jsonlite::read_json(file.path(dirname(dirs$plots), got$settings), simplifyVector = TRUE)
  expect_identical(settings$variant, "direction")
  expect_equal(settings$quantiles, prepared$metadata$quantiles)
  expect_identical(lisaR:::lisa_write_runtime_nes_variants(prepared, dirs,
    "A", "GOBP-C2", "semantic", "clean", "tiff")$status, "not_requested_supported_format")
})

test_that("contrast runtime keeps original means and delta with separate optional distributions", {
  skip_if_not_installed("ggplot2")
  fixture <- nes_runtime_fixture()
  original_a <- readBin(fixture$a$path, "raw", n = file.info(fixture$a$path)$size)
  original_b <- readBin(fixture$b$path, "raw", n = file.info(fixture$b$path)$size)
  # The runtime consumes persisted summaries, whose blank string metadata uses
  # the historical TSV round-trip representation. Compare that same input
  # boundary rather than conflating NA labels with numerical differences.
  expected <- lisaR:::build_lisa_contrast_summary(
    lisaR:::read_lisa_contrast_summary(fixture$a$path, "A", .05),
    lisaR:::read_lisa_contrast_summary(fixture$b$path, "B", .05),
    "A", "B", "GOBP-C2", include_missing_categories = TRUE,
    direction_epsilon = 0, gsea_padj_cutoff = .05)
  expected <- lisaR:::order_lisa_contrast_summary(expected, "supracategory")
  out <- file.path(fixture$root, "contrast")
  got <- lisaR:::run_LISA_contrast(fixture$a$path, fixture$b$path, out,
    comparison_name = "A_vs_B", contrast_a_label = "A", contrast_b_label = "B",
    universes = "GOBP-C2", export_formats = "tsv", plot_formats = "png",
    plot_sets = "all", annotate_gene_sets = FALSE, gsea_padj_cutoff = .05,
    verbose = FALSE)
  expect_equal(got$summary, expected)
  expect_identical(got$manifest$annotation_variant, "plain")
  expect_identical(got$nes_variant_manifest$variant, c("clean", "percentages", "direction", "dispersion"))
  expect_true(all(got$nes_variant_manifest$status == "rendered"))
  src <- utils::read.delim(file.path(out, got$nes_variant_manifest$source[[1L]]),
    quote = '"', check.names = FALSE)
  expect_equal(src$delta_mean_NES[src$side == "A"], expected$delta_mean_NES)
  context_a <- src[src$side == "A" & src$category_id == "CONTEXT", ]
  expect_equal(context_a$display_mean_NES, .5)
  expect_equal(context_a$n_genesets_significant, 0)
  expect_true(all(is.na(context_a[c("median_NES", "p25_NES", "p75_NES")])) )
  only_b_a <- src[src$side == "A" & src$category_id == "ONLY_B", ]
  expect_identical(only_b_a$display_mean_NES, -.5)
  expect_equal(only_b_a$n_genesets_significant, 0)
  expect_identical(only_b_a$endpoint_state, "contextual_not_significant")
  expect_true(all(is.na(only_b_a[c("mean_NES", "median_NES", "p25_NES", "p75_NES")])))
  expect_equal(readBin(fixture$a$path, "raw", length(original_a)), original_a)
  expect_equal(readBin(fixture$b$path, "raw", length(original_b)), original_b)
})

test_that("NES fixture does not relax the canonical missing-counterpart contract", {
  fixture <- nes_runtime_fixture()
  absent_a <- fixture$a$summary[fixture$a$summary$category_id != "ONLY_B", , drop = FALSE]
  expect_error(lisaR:::build_lisa_contrast_summary(absent_a, fixture$b$summary,
    "A", "B", "GOBP-C2", include_missing_categories = TRUE,
    direction_epsilon = 0, gsea_padj_cutoff = .05), "LISA-CONTRAST-SUPPORT-002")
})

test_that("summary-only contrasts record missing distributions without inventing support", {
  skip_if_not_installed("ggplot2")
  fixture <- nes_runtime_fixture()
  a <- file.path(fixture$root, "standalone_a.tsv"); file.copy(fixture$a$path, a)
  b <- file.path(fixture$root, "standalone_b.tsv"); file.copy(fixture$b$path, b)
  out <- file.path(fixture$root, "summary-only")
  got <- lisaR:::run_LISA_contrast(a, b, out, comparison_name = "legacy",
    universes = "GOBP-C2", export_formats = "tsv", plot_formats = "png",
    plot_sets = "all", annotate_gene_sets = FALSE, gsea_padj_cutoff = .05, verbose = FALSE)
  expect_true(all(got$nes_variant_manifest$status == "unavailable_missing_member_table"))
  expect_false(dir.exists(file.path(out, "plots", "category_nes")))
  expect_true(file.exists(paste0(got$manifest$plot_stem[[1L]], ".png")))
})
