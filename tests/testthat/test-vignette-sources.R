test_that("vignette sources and figure semantics are discoverable", {
  retained <- c("configuration-reference", "function-reference", "getting-started",
    "parallel-execution", "riaz-gse91061-worked-example", "scientific-report",
    "quick-start", "how-lisar-works", "resources-and-provenance",
    "troubleshooting", "custom-dictionaries-and-supercategories",
    "semantic-category-guide", "installation",
    "preparing-de-inputs", "cptac-ccrcc-proteomics",
    "dictionary-construction", "licence-and-citation", "shiny-app")
  for (name in retained) {
    expect_true(file.exists(vignette_source_path(name)))
    expect_match(read_vignette_source(name), "VignetteEngine{knitr::rmarkdown}",
                 fixed = TRUE)
  }

  semantics_path <- system.file("semantics", "figure-semantics.tsv", package = "lisaR")
  expect_true(file.exists(semantics_path))
  semantics <- utils::read.delim(semantics_path, sep = "\t", check.names = FALSE)
  expect_true(all(c("figure_type", "x_encoding", "y_encoding",
    "colour_encoding", "source_fields", "interpretation_limit") %in% names(semantics)))
  expect_setequal(semantics$figure_type, c("standard_lisa_overview",
    "lisa_category_gene_sets", "ora_category_bar", "lisa_dumbbell",
    "volcano_overlay", "gene_card", "contrast_gene_card",
    "leading_edge_heatmap", "contrast_heatmap", "network", "kegg_single",
    "kegg_contrast"))
  expect_false(any(!nzchar(semantics$interpretation_limit)))
})
