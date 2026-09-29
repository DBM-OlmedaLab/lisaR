test_that("single-analysis plots expose the configured comparison and direction", {
  skip_if_not_installed("ggplot2")
  base <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  plot <- lisaR:::lisa_apply_plot_context(
    base,
    display_title = "Responders: ON vs PRE",
    comparison_subtitle = "CR/PR ON-treatment minus CR/PR PRE-treatment",
    positive_direction = "Positive values indicate an on-treatment increase in responders."
  )
  expect_identical(plot$labels$title, "Responders: ON vs PRE")
  expect_match(
    plot$labels$subtitle,
    "CR/PR ON-treatment minus CR/PR PRE-treatment",
    fixed = TRUE
  )
  expect_match(
    plot$labels$subtitle,
    "Positive values indicate an on-treatment increase in responders.",
    fixed = TRUE
  )
})

test_that("configured analysis labels replace internal analysis IDs", {
  index <- data.frame(
    analysis_id = "responders_on_vs_pre",
    label = "Responders: ON vs PRE",
    stringsAsFactors = FALSE
  )
  expect_identical(
    lisaR:::lisa_de_label(index, "responders_on_vs_pre"),
    "Responders: ON vs PRE"
  )
})

test_that("the contrast runner receives the index used for human labels", {
  expect_true(
    "de_index" %in% names(formals(lisaR:::run_lisa_pipeline_contrasts))
  )
})

test_that("production code resolution ignores a source checkout in the working directory", {
  root <- tempfile("lisa-script-priority-")
  dir.create(file.path(root, "inst", "scripts"), recursive = TRUE)
  writeLines(
    c(
      "Package: lisaR",
      "Version: 0.0.0",
      "Title: Source Checkout Fixture",
      "Description: Minimal package metadata for testing source layout detection.",
      "License: MIT"
    ),
    file.path(root, "DESCRIPTION")
  )
  canonical <- file.path(root, "inst", "scripts", "worker.R")
  writeLines("canonical", canonical)
  old <- getOption("lisaR.test_package_dir", NULL)
  old_wd <- getwd()
  on.exit({
    options(lisaR.test_package_dir = old)
    setwd(old_wd)
  }, add = TRUE)
  options(lisaR.test_package_dir = root)
  setwd(root)
  resolved <- lisaR:::lisa_resolve_package_dir()
  expected_package_root <- system.file(package = "lisaR")
  if (!file.exists(file.path(expected_package_root, "DESCRIPTION"))) {
    # pkgload maps system.file(package = "lisaR") to source/inst; production
    # resolution intentionally returns the containing source package root.
    expected_package_root <- dirname(expected_package_root)
  }
  expect_identical(
    as.character(resolved),
    normalizePath(expected_package_root, winslash = "/")
  )
  expect_false(identical(as.character(resolved), normalizePath(root, winslash = "/")))

  injected <- lisaR:::lisa_resolve_package_dir(.test_package_dir = root)
  expect_true(isTRUE(attr(injected, "lisaR.test_injection")))
  expect_identical(
    lisaR:::lisa_post_script_path(injected, "worker.R"),
    normalizePath(canonical, winslash = "/")
  )

  writeLines(sub("^Package: lisaR$", "Package: notlisaR", readLines(
    file.path(root, "DESCRIPTION"), warn = FALSE
  )), file.path(root, "DESCRIPTION"))
  expect_error(
    lisaR:::lisa_resolve_package_dir(.test_package_dir = root),
    "not a lisaR source"
  )
})

test_that("standalone postprocessors build comparison-first subtitles", {
  helper <- system.file("scripts", "lisa_plot_metadata.R", package = "lisaR")
  expect_true(file.exists(helper))
  environment <- new.env(parent = baseenv())
  sys.source(helper, envir = environment)
  metadata <- list(
    label = "Baseline: responders vs PD",
    comparison = "CR/PR PRE-treatment minus PD PRE-treatment",
    positive_direction = "Positive values are higher in future responders at baseline.",
    model_note = "Adjusted for prior ipilimumab."
  )
  subtitle <- environment$lisa_plot_subtitle(
    metadata, prefix = "T-cell receptor signalling", suffix = "PATHWAYS"
  )
  expect_match(subtitle, metadata$comparison, fixed = TRUE)
  expect_match(subtitle, metadata$positive_direction, fixed = TRUE)
  expect_false(grepl("^responders_vs_pd_pre", subtitle))
})

test_that("standalone plot subtitles are portable across shipped graphics formats", {
  skip_if_not_installed("ggplot2")
  helper <- system.file("scripts", "lisa_plot_metadata.R", package = "lisaR")
  expect_true(file.exists(helper))
  environment <- new.env(parent = baseenv())
  sys.source(helper, envir = environment)
  metadata <- list(
    comparison = "(CR/PR ON \u2212 PRE) minus (PD ON \u2212 PRE)",
    positive_direction = "en dash \u2013 em dash \u2014",
    label = "Longitudinal interaction",
    model_note = ""
  )
  subtitle <- environment$lisa_plot_subtitle(metadata)
  expect_identical(
    subtitle,
    "(CR/PR ON - PRE) minus (PD ON - PRE)\nen dash - em dash -"
  )
  expect_false(grepl("[\u2013\u2014\u2212]", subtitle))

  plot <- ggplot2::ggplot(
    data.frame(x = 1, y = 1),
    ggplot2::aes(x, y)
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(subtitle = subtitle)
  devices <- list(
    png = list(device = "png"),
    pdf = list(device = "pdf", useDingbats = FALSE),
    svg = list(device = grDevices::svg)
  )
  for (extension in names(devices)) {
    path <- tempfile("lisa-device-text-", fileext = paste0(".", extension))
    arguments <- c(
      list(filename = path, plot = plot, width = 4, height = 3),
      devices[[extension]]
    )
    expect_warning(do.call(ggplot2::ggsave, arguments), NA)
    expect_true(file.exists(path))
    expect_gt(file.info(path)$size, 0)
  }
})

test_that("contrast subtitles are device-safe without changing the input manifest", {
  skip_if_not_installed("ggplot2")
  make_summary <- function(path, mean_nes) {
    summary <- data.frame(
      category_id = "SYN_CATEGORY_A",
      category_display_name = "T-cell receptor signaling",
      macrogroup_id = "IMMUNE",
      macrogroup_name = "Immune signaling",
      macrogroup_order = 1,
      category_order_within_macrogroup = 1,
      color = "#336699",
      n_genesets = 4,
      mean_NES = mean_nes,
      median_NES = mean_nes,
      consistency = 1,
      same_direction_pct = 100,
      min_padj = 0.01,
      stringsAsFactors = FALSE
    )
    write.table(summary, path, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  root <- tempfile("lisa-contrast-device-text-")
  dir.create(root)
  summary_a <- file.path(root, "a_GSEA_category_summary.tsv")
  summary_b <- file.path(root, "b_GSEA_category_summary.tsv")
  make_summary(summary_a, 1.5)
  make_summary(summary_b, 0.5)
  subtitles <- c(
    response_dynamics =
      "LISA profile(Responders ON \u2212 PRE) minus LISA profile(PD ON \u2212 PRE)",
    response_separation_by_time =
      "LISA profile(CR/PR \u2212 PD at ON) minus LISA profile(CR/PR \u2212 PD at PRE)"
  )

  for (case_id in names(subtitles)) {
    output_dir <- file.path(root, case_id)
    expect_warning(
      lisaR:::run_LISA_contrast(
        contrast_a = summary_a,
        contrast_b = summary_b,
        output_dir = output_dir,
        comparison_name = case_id,
        contrast_a_label = "Contrast A",
        contrast_b_label = "Contrast B",
        comparison_subtitle = subtitles[[case_id]],
        universes = "PATHWAYS",
        plot_sets = "all",
        export_formats = "tsv",
        plot_formats = c("png", "pdf", "svg"),
        annotate_gene_sets = FALSE,
        verbose = FALSE
      ),
      NA
    )
    input_manifest <- read.delim(
      file.path(output_dir, "qc", paste0(case_id, "_input_manifest.tsv")),
      check.names = FALSE
    )
    expect_identical(input_manifest$comparison_subtitle[[1]], subtitles[[case_id]])
    plot_stem <- file.path(
      output_dir, "plots",
      lisaR:::lisa_contrast_plot_stem(case_id, "LISA", "all")
    )
    for (extension in c("png", "pdf", "svg")) {
      path <- paste0(plot_stem, ".", extension)
      expect_true(file.exists(path))
      expect_gt(file.info(path)$size, 0)
    }
  }
})
