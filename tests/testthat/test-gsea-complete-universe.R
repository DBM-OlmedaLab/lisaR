g11_write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = TRUE, na = ""
  )
  path
}

g11_run_fixture <- function(universe) {
  root <- tempfile(paste0("lisa-g11-", tolower(universe), "-"))
  dir.create(root)
  genes <- paste0("GENE", sprintf("%02d", seq_len(30L)))
  gene_sets <- c("GS_CLASSIFIED", "GS_RESIDUAL_A", "GS_RESIDUAL_B")
  term2gene <- do.call(rbind, lapply(seq_along(gene_sets), function(index) {
    data.frame(
      gs_collection = universe,
      gs_subcollection = "SYNTHETIC",
      gs_name = gene_sets[[index]],
      gs_exact_source = "G11_TEST",
      gene_symbol = genes[seq.int((index - 1L) * 10L + 1L, index * 10L)],
      stringsAsFactors = FALSE
    )
  }))
  dictionary <- data.frame(
    universe = universe,
    gene_set_id = "GS_CLASSIFIED",
    gene_set_name = "Classified test set",
    source_id = "G11_TEST",
    category_id = "CAT_CLASSIFIED",
    category_display_name = "Classified category",
    tier = "core",
    stringsAsFactors = FALSE
  )
  category_map <- data.frame(
    category_id = "CAT_CLASSIFIED",
    display_name = "Classified category",
    macrogroup_id = "SUPER_CLASSIFIED",
    macrogroup_name = "Classified supercategory",
    macrogroup_order = 1L,
    category_order_within_macrogroup = 1L,
    notes = "G11 test fixture",
    stringsAsFactors = FALSE
  )
  de <- data.frame(
    symbol = genes,
    stat = seq(3, -3, length.out = length(genes)),
    log2FoldChange = seq(2, -2, length.out = length(genes)),
    padj = seq(0.001, 0.9, length.out = length(genes)),
    stringsAsFactors = FALSE
  )
  paths <- list(
    de = g11_write_tsv(de, file.path(root, "de.tsv")),
    term2gene = g11_write_tsv(term2gene, file.path(root, "term2gene.tsv")),
    dictionary = g11_write_tsv(dictionary, file.path(root, "dictionary.tsv")),
    category_map = g11_write_tsv(category_map, file.path(root, "category_map.tsv")),
    output = file.path(root, "output")
  )
  list(root = root, paths = paths, gene_sets = gene_sets)
}

test_that("semantic universe selection is independent of dictionary coverage", {
  term2gene <- data.frame(
    gs_collection = c("C5", "C2", "C5", "C5"),
    gs_subcollection = c("GO:BP", "CP:BIOCARTA", "GO:MF", "GO:CC"),
    gs_name = c("GOBP_TEST", "BIOCARTA_TEST", "GOMF_TEST", "GOCC_TEST"),
    gs_exact_source = "G11_TEST",
    gene_symbol = c("A", "B", "C", "D"),
    stringsAsFactors = FALSE
  )
  c2_sources <- c("BIOCARTA", "KEGG_LEGACY", "KEGG_MEDICUS", "PID",
                  "REACTOME", "WIKIPATHWAYS", "CP_OTHER")
  gobp <- lisaR:::select_lisa_term2gene_universe(
    term2gene, "GOBP-C2", c2_sources, dictionary_gene_sets = "GOBP_TEST"
  )
  pathways <- lisaR:::select_lisa_term2gene_universe(
    term2gene, "PATHWAYS", c2_sources, dictionary_gene_sets = character()
  )
  gomf <- lisaR:::select_lisa_term2gene_universe(
    term2gene, "GOMF", c2_sources, dictionary_gene_sets = character()
  )
  gocc <- lisaR:::select_lisa_term2gene_universe(
    term2gene, "GOCC", c2_sources, dictionary_gene_sets = character()
  )
  expect_setequal(unique(gobp$gs_name), c("GOBP_TEST", "BIOCARTA_TEST"))
  expect_setequal(unique(pathways$gs_name), c("GOBP_TEST", "BIOCARTA_TEST"))
  expect_identical(unique(gomf$gs_name), "GOMF_TEST")
  expect_identical(unique(gocc$gs_name), "GOCC_TEST")
})

test_that("GSEA ledger accounts for every source set and fails closed", {
  fixture <- g11_run_fixture("GOBP-C2")
  on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
  term2gene <- utils::read.delim(fixture$paths$term2gene,
                                stringsAsFactors = FALSE)
  dictionary <- utils::read.delim(fixture$paths$dictionary,
                                  stringsAsFactors = FALSE)
  ranks <- stats::setNames(seq(3, -3, length.out = 30L),
                           paste0("GENE", sprintf("%02d", seq_len(30L))))
  ledger <- lisaR:::build_lisa_gsea_ledger(
    term2gene, dictionary, ranks, min_gs_size = 3L, max_gs_size = 20L,
    universe = "GOBP-C2"
  )
  expect_identical(nrow(ledger), 3L)
  expect_true(all(ledger$included_in_gsea_universe))
  expect_true(all(ledger$eligible_for_gsea))
  expect_identical(sum(ledger$classification_status == "classified"), 1L)
  expect_identical(sum(ledger$classification_status == "unclassified"), 2L)
  expect_true(all(
    ledger$classification_bucket[ledger$classification_status == "unclassified"] ==
      "OTHER_UNCLASSIFIED"
  ))

  result <- data.frame(
    pathway = ledger$pathway,
    pval = 0.5, padj = 0.5, ES = 0, NES = 0, size = 10L,
    leadingEdge = "", stringsAsFactors = FALSE
  )
  completed <- lisaR:::complete_lisa_gsea_results(result, ledger)
  expect_identical(completed$gsea_result_status, rep("tested", 3L))
  expect_error(
    lisaR:::complete_lisa_gsea_results(result[-1, ], ledger),
    "LISA-GSEA-UNIVERSE-004"
  )
})

test_that("run_LISA_DE writes residual tables except for PATHWAYS", {
  skip_if_not_installed("fgsea")
  for (universe in c("GOBP-C2", "PATHWAYS")) {
    fixture <- g11_run_fixture(universe)
    on.exit(unlink(fixture$root, recursive = TRUE, force = TRUE), add = TRUE)
    result <- suppressWarnings(lisaR:::run_LISA_DE(
      input = fixture$paths$de,
      input_type = "de_table",
      output_dir = fixture$paths$output,
      comparison_name = "g11_fixture",
      species = "Homo sapiens",
      lisa_dictionary = "fixture",
      lisa_dictionary_path = fixture$paths$dictionary,
      term2gene_path = fixture$paths$term2gene,
      category_map_path = fixture$paths$category_map,
      universes = universe,
      rank_col = "stat",
      run_gsea = TRUE,
      run_ora = FALSE,
      min_gs_size = 3L,
      max_gs_size = 20L,
      n_threads = 1L,
      fgsea_nperm = 100L,
      gsea_padj_cutoff = 1,
      outputs = if (universe == "GOBP-C2") {
        c("tables", "qc", "plots")
      } else {
        c("tables", "qc")
      },
      plots = if (universe == "GOBP-C2") {
        c("lollipop", "direction_lollipop", "category_pathways", "dotplot")
      } else {
        character()
      },
      export_formats = "tsv",
      plot_formats = "png",
      figure_recipes = FALSE,
      file_label_prefix = "semantic",
      lisa_project_root = fixture$root,
      verbose = FALSE
    ))
    expect_identical(result$nes_variant_manifest$variant, c("clean", "percentages", "direction", "dispersion"))
    if (universe == "GOBP-C2") {
      expect_true(all(result$nes_variant_manifest$status == "rendered"))
      expect_true(all(file.exists(file.path(fixture$paths$output, result$nes_variant_manifest$path))))
      expect_false(any(grepl("OTHER_UNCLASSIFIED",
        readLines(file.path(fixture$paths$output, result$nes_variant_manifest$source[[1L]])))))
    } else {
      expect_true(all(result$nes_variant_manifest$status == "not_requested_plot_output"))
    }
    expect_setequal(result$gsea_ledger$pathway, fixture$gene_sets)
    expect_identical(nrow(result$gsea_unclassified), 2L)
    expect_setequal(
      result$gsea_unclassified$pathway,
      c("GS_RESIDUAL_A", "GS_RESIDUAL_B")
    )
    expect_true(file.exists(file.path(
      fixture$paths$output, "qc", "g11_fixture_GSEA_universe_ledger.tsv"
    )))
    residual_path <- file.path(
      fixture$paths$output, "lisa_tables",
      "g11_fixture_semantic_GSEA_OTHER_UNCLASSIFIED.tsv"
    )
    expect_identical(file.exists(residual_path), universe != "PATHWAYS")
    expect_false(any(
      is.na(result$summaries$gsea$category_id) |
        result$summaries$gsea$category_id == "OTHER_UNCLASSIFIED"
    ))
    report <- readLines(file.path(
      fixture$paths$output, "g11_fixture_semantic_RUN_REPORT.txt"
    ), warn = FALSE)
    expect_true(any(report == "GSEA universe gene sets: 3"))
    expect_true(any(report == "GSEA unclassified gene sets: 2"))
    if (universe == "GOBP-C2") {
      plot_files <- list.files(
        file.path(fixture$paths$output, "plots"), recursive = TRUE,
        full.names = TRUE
      )
      expect_true(any(grepl("[.]png$", plot_files)))
      tabular_plot_sources <- plot_files[grepl("[.]tsv$", plot_files)]
      expect_gt(length(tabular_plot_sources), 0L)
      plot_text <- unlist(lapply(
        tabular_plot_sources, readLines, warn = FALSE
      ), use.names = FALSE)
      expect_false(any(grepl("GS_RESIDUAL", plot_text, fixed = TRUE)))
      expect_false(any(grepl("OTHER_UNCLASSIFIED", plot_text,
                             fixed = TRUE)))
    }
  }
})

test_that("unclassified rows are rejected by every plot boundary", {
  unclassified <- data.frame(
    pathway = "GS_RESIDUAL",
    category_id = NA_character_,
    classification_status = "unclassified",
    classification_bucket = "OTHER_UNCLASSIFIED",
    stringsAsFactors = FALSE
  )
  expect_error(
    lisaR:::lisa_assert_classified_plot_rows(unclassified, "G11 test"),
    "LISA-PLOT-UNCLASSIFIED-001"
  )
  expect_identical(nrow(lisaR:::lisa_classified_gsea_rows(unclassified)), 0L)

  script <- system.file(
    "scripts", "build_single_de_enrichmentmap.R", package = "lisaR"
  )
  if (!nzchar(script) || !file.exists(script)) {
    script <- file.path(
      testthat::test_path("..", ".."), "inst", "scripts",
      "build_single_de_enrichmentmap.R"
    )
  }
  expect_true(file.exists(script))
  code <- paste(readLines(script, warn = FALSE), collapse = "\n")
  expect_match(code, "gsea <- classified_rows(read_tsv", fixed = TRUE)
  expect_match(code, "OTHER_UNCLASSIFIED", fixed = TRUE)
})

test_that("each run drops map categories and supercategories absent from its tier", {
  dictionary <- data.frame(
    category_id = "CAT_ACTIVE", category_display_name = "Active",
    stringsAsFactors = FALSE
  )
  category_map <- data.frame(
    category_id = c("CAT_ACTIVE", "CAT_OTHER_TIER", "CAT_EMPTY"),
    display_name = c("Active", "Other tier", "Empty"),
    macrogroup_id = c("SUPER_ACTIVE", "SUPER_OTHER", "SUPER_EMPTY"),
    macrogroup_name = c("Active", "Other", "Empty"),
    macrogroup_order = 1:3,
    category_order_within_macrogroup = 1L,
    stringsAsFactors = FALSE
  )
  registered <- lisaR:::prepare_registered_category_map(
    category_map, dictionary
  )
  low_level <- lisaR:::augment_category_map(category_map, dictionary)
  expect_identical(registered$category_id, "CAT_ACTIVE")
  expect_identical(low_level$category_id, "CAT_ACTIVE")
  expect_identical(unique(registered$macrogroup_id), "SUPER_ACTIVE")
  expect_identical(unique(low_level$macrogroup_id), "SUPER_ACTIVE")
})
