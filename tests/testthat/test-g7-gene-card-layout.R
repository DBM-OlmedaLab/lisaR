g7_package_path <- function(...) {
  normalizePath(
    system.file(..., package = "lisaR"),
    mustWork = TRUE
  )
}

g7_builder_path <- function() {
  g7_package_path(
    "scripts", "build_single_de_category_gene_cards.R"
  )
}

g7_gene_card_evidence <- function(n = 20L, include_lisa_support = TRUE) {
  index <- seq_len(n)
  evidence <- data.frame(
    category_id = "G7_LAYOUT",
    gene_contribution_score = rev(index) + 0.25,
    padj = pmin(0.9, index / 100),
    log2FC = rep(c(-2.1, 1.8, 0.7, -0.4), length.out = n),
    symbol = sprintf("LONG_GENE_SYMBOL_%02d", index),
    rank_value = rev(index),
    in_gsea_leading_edge = rep(c(TRUE, TRUE, FALSE), length.out = n),
    in_ora_overlap = rep(c(FALSE, TRUE, TRUE), length.out = n),
    category_display_name = "A deliberately long category title for clipping regression",
    macrogroup_name = "A deliberately long supercategory label",
    color = "#4C6A92",
    de_is_significant = rep(c(TRUE, FALSE), length.out = n),
    category_recurrence = rep(1:4, length.out = n),
    source_geneset_count = rep(1:3, length.out = n),
    stringsAsFactors = FALSE
  )
  if (include_lisa_support) {
    evidence$lisa_support_score <- rev(index) / 3
  }
  evidence
}

test_that("full dependency contract excludes the retired layout package", {
  retired_layout_package <- paste0("cow", "plot")
  full <- lisaR:::lisa_config_dependency_requirements(
    pipeline = list(run_kegg_maps = FALSE),
    species = "Homo sapiens",
    report = list(mode = "full"),
    available = function(package) FALSE
  )

  expect_true(all(c("ggrepel", "gridExtra", "patchwork") %in% full$package))
  expect_false(retired_layout_package %in% full$package)
})

test_that("distributable metadata and GeneCard builder use grid and patchwork", {
  retired_layout_package <- paste0("cow", "plot")
  description <- readLines(g7_package_path("DESCRIPTION"), warn = FALSE)
  profiles <- utils::read.delim(
    g7_package_path("standalone", "dependency-profiles.tsv"),
    sep = "\t", check.names = FALSE, stringsAsFactors = FALSE
  )
  builder <- readLines(g7_builder_path(), warn = FALSE)

  expect_false(any(grepl(retired_layout_package, description, fixed = TRUE)))
  expect_false(retired_layout_package %in% profiles$package)
  expect_false(any(grepl(retired_layout_package, builder, fixed = TRUE)))
  expect_true(any(grepl("grid::textGrob", builder, fixed = TRUE)))
  expect_true(any(grepl("patchwork::wrap_elements", builder, fixed = TRUE)))
  expect_true(any(grepl("patchwork::wrap_plots", builder, fixed = TRUE)))
})

test_that("GeneCard builder renders both table branches without the retired package", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("gridExtra")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("png")
  fixture_dir <- tempfile("g7-gene-card-layout-")
  dir.create(fixture_dir)
  input_paths <- file.path(fixture_dir, c("with_support.tsv", "without_support.tsv"))
  for (index in seq_along(input_paths)) {
    include_lisa_support <- index == 1L
    evidence <- g7_gene_card_evidence(
      n = if (include_lisa_support) 20L else 3L,
      include_lisa_support = include_lisa_support
    )
    if (include_lisa_support) evidence$gene_contribution_score[[1L]] <- 145.63
    if (!include_lisa_support) evidence$log2FC <- 1.2
    utils::write.table(
      evidence, input_paths[[index]], sep = "\t", quote = FALSE,
      row.names = FALSE
    )
  }

  runner_path <- file.path(fixture_dir, "render-fixtures.R")
  writeLines(c(
    "args <- commandArgs(trailingOnly = TRUE)",
    "environment <- new.env(parent = globalenv())",
    "sys.source(args[[1]], envir = environment)",
    "stopifnot(isTRUE(all.equal(environment$gene_card_lfc_breaks(c(-2.1, -0.4, 0.7, 1.8)), c(-2.1, 0, 1.8))))",
    "stopifnot(isTRUE(all.equal(environment$gene_card_lfc_breaks(c(0.2, 0.8, 1.4)), c(0.2, 0.8, 1.4))))",
    "stopifnot(isTRUE(all.equal(environment$gene_card_lfc_breaks(rep(1.2, 4)), 1.2)))",
    "find_score_plot <- function(plot) {",
    "  if (identical(plot$labels$y, 'DE-weighted contribution')) return(plot)",
    "  children <- plot$patches$plots",
    "  if (length(children) == 0L) return(NULL)",
    "  for (child in children) {",
    "    found <- Recall(child)",
    "    if (!is.null(found)) return(found)",
    "  }",
    "  NULL",
    "}",
    "for (index in 2:3) {",
    "  evidence <- utils::read.delim(args[[index]], sep = '\\t', check.names = FALSE)",
    "  card <- environment$build_card_plot(evidence, 'G7_LAYOUT', top_genes = 15L)",
    "  stopifnot(inherits(card, 'patchwork'))",
    "  score_plot <- find_score_plot(card)",
    "  stopifnot(!is.null(score_plot))",
    "  score_build <- ggplot2::ggplot_build(score_plot)",
    "  text_layers <- which(vapply(score_build$data, function(layer) 'label' %in% names(layer), logical(1)))",
    "  stopifnot(length(text_layers) == 1L)",
    "  score_text <- score_build$data[[text_layers[[1L]]]]",
    "  score_max <- max(score_text$y[is.finite(score_text$y)])",
    "  score_upper <- score_build$layout$panel_params[[1L]]$x.range[[2L]]",
    "  stopifnot(score_upper / score_max >= 1.20 - 1e-8)",
    "  if (index == 2L) stopifnot('145.63' %in% score_text$label)",
    "  grDevices::pdf(NULL)",
    "  grob <- tryCatch(patchwork::patchworkGrob(card), finally = grDevices::dev.off())",
    "  stopifnot(inherits(grob, 'gtable'))",
    "  stem <- sub('[.]tsv$', '', args[[index]])",
    "  png_path <- paste0(stem, '.png')",
    "  pdf_path <- paste0(stem, '.pdf')",
    "  ggplot2::ggsave(png_path, card, width = 8.9, height = 10.4, dpi = 72, bg = 'white', limitsize = FALSE)",
    "  ggplot2::ggsave(pdf_path, card, width = 8.9, height = 10.4, bg = 'white', limitsize = FALSE)",
    "  stopifnot(file.info(png_path)$size > 1000, file.info(pdf_path)$size > 1000)",
    "  image <- png::readPNG(png_path)",
    "  stopifnot(identical(dim(image)[1:2], c(748L, 640L)))",
    "  bottom_rgb <- image[(dim(image)[[1]] - 2L):dim(image)[[1]], , 1:3, drop = FALSE]",
    "  stopifnot(all(bottom_rgb > 0.98))",
    "}"
  ), runner_path, useBytes = TRUE)

  # `system2(env = )` is not portable: on Windows the NAME=value tokens are
  # appended to the child's command line (see ?system2), where Rscript reads
  # the first one as the script to run. Export the one variable that matters
  # instead -- `--vanilla` already covers R_PROFILE_USER/R_ENVIRON_USER, and
  # `/dev/null` is not a Windows path. R_TESTS must be cleared because
  # R CMD check points it at a startup file that breaks nested R processes.
  output <- lisaR:::lisa_with_child_environment(
    c(R_TESTS = ""),
    suppressWarnings(system2(
      file.path(R.home("bin"), "Rscript"),
      c("--vanilla", runner_path, g7_builder_path(), input_paths),
      stdout = TRUE,
      stderr = TRUE
    ))
  )
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  expect_identical(as.integer(status), 0L, info = paste(output, collapse = "\n"))
})
