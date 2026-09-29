test_that("the run output policy is the one source of plot formats", {
  root <- tempfile("lisa-policy-")
  dir.create(root)
  expect_identical(lisaR:::lisa_requested_plot_formats(root), "png")
  lisaR:::write_lisa_tsv(
    lisaR:::lisa_output_policy_table(c("png", "svg"), TRUE, TRUE),
    file.path(root, "report_output_policy.tsv"))
  expect_identical(lisaR:::lisa_requested_plot_formats(root), c("png", "svg"))
  lisaR:::write_lisa_tsv(
    lisaR:::lisa_output_policy_table("png", TRUE, TRUE),
    file.path(root, "report_output_policy.tsv"))
  expect_identical(lisaR:::lisa_requested_plot_formats(root), "png")
})

test_that("missing-product completion no longer hard-codes png,pdf (C3)", {
  body_text <- paste(deparse(body(lisaR:::lisa_complete_full_products)), collapse = "\n")
  expect_false(grepl('"--formats", "png,pdf"', body_text, fixed = TRUE))
  expect_match(body_text, "lisa_requested_plot_formats(root)", fixed = TRUE)
  expect_true("families" %in% names(formals(lisaR:::lisa_complete_full_products)))
})

test_that("the recurrent-gene collector keeps the SVG its generator wrote (C1)", {
  root <- tempfile("lisa-collect-")
  src <- file.path(root, "src"); dest <- file.path(root, "dest")
  dir.create(src, recursive = TRUE)
  stem <- file.path(src, "A_top_recurrent_genes")
  for (ext in c(".png", ".svg", ".pdf", "_source.tsv", "_recipe.R")) {
    writeLines("x", paste0(stem, ext))
  }
  index <- file.path(src, "A_recurrent_gene_screen_index.tsv")
  lisaR:::write_lisa_tsv(data.frame(plot_png = paste0(stem, ".png"),
    plot_pdf = paste0(stem, ".pdf"), plot_svg = paste0(stem, ".svg")), index)
  old <- options(lisaR.run_root = root)
  on.exit(options(old), add = TRUE)
  collected <- lisaR:::lisa_extension_collect_inventory(index, src, dest,
    columns = c("plot_png", "plot_pdf", "plot_svg"),
    stem_pattern = lisaR:::lisa_format_pattern(),
    companions = c(".png", ".pdf", ".svg", "_source.tsv", "_recipe.R"))
  expect_true("A_top_recurrent_genes.svg" %in% collected)
  expect_true(file.exists(file.path(dest, "A_top_recurrent_genes.svg")))
  arch <- paste(deparse(body(lisaR:::lisa_extension_render_unit)), collapse = "\n")
  expect_match(arch, '"plot_svg"', fixed = TRUE)
})

test_that("product cards flag a missing SVG only when SVG was requested (C4)", {
  pill <- lisaR:::lisa_svg_status_pill
  expect_identical(pill(c("png", "tsv"), "gene_cards", "a.png", FALSE), "")
  expect_identical(pill(c("png", "svg"), "gene_cards", c("a.png", "a.svg"), TRUE), "")
  expect_match(pill(c("png", "tsv"), "gene_cards", "a.png", TRUE), "SVG missing", fixed = TRUE)
  expect_match(pill("png", "kegg_pathway_map", "x_painted.png", TRUE),
    "SVG not available for this product type", fixed = TRUE)
  viewer <- readLines(system.file("category-evidence", "viewer.js", package = "lisaR"))
  expect_true(any(grepl("svgStatus", viewer, fixed = TRUE)))
  viewer <- readLines(system.file("contrast-evidence", "viewer.js", package = "lisaR"))
  expect_true(any(grepl("svgStatus", viewer, fixed = TRUE)))
})

test_that("the format policy of an in-run FULL tree removes unrequested formats only", {
  root <- tempfile("lisa-format-policy-")
  art <- file.path(root, "artifacts", "a")
  dir.create(art, recursive = TRUE)
  files <- c("x.png", "x.svg", "x.pdf", "x_source.tsv", "x_recipe.R", "m_kegg_base.png")
  for (f in files) writeLines("x", file.path(art, f))
  old <- options(lisaR.run_root = root)
  on.exit(options(old), add = TRUE)
  lisaR:::lisa_apply_format_policy(file.path(root, "artifacts"),
    list(png = TRUE, svg = FALSE, pdf = FALSE))
  expect_setequal(list.files(art), c("x.png", "x_source.tsv", "x_recipe.R", "m_kegg_base.png"))
})

test_that("build_LISA_report.R accepts --report-mode through its registered interface", {
  interface <- lisaR:::lisa_post_script_interfaces()[["build_LISA_report.R"]]
  expect_true("report-mode" %in% interface)
  script <- readLines(system.file("scripts", "build_LISA_report.R", package = "lisaR"))
  expect_true(any(grepl('report_mode = ""', script, fixed = TRUE)))
  expect_true(any(grepl("--report-mode conflicts", script, fixed = TRUE)))
})

test_that("run_lisa(mode = 'full') renders FULL products inside the same run", {
  bridge <- paste(deparse(body(lisaR:::run_lisa_pipeline_from_config)), collapse = "\n")
  expect_false(grepl("-full-report", bridge, fixed = TRUE))
  expect_false(grepl("render_lisa_categories(final_output_dir", bridge, fixed = TRUE))
  expect_match(bridge, ".presentation = list(", fixed = TRUE)
  post <- paste(deparse(body(lisaR:::run_lisa_pipeline_post_lisa)), collapse = "\n")
  expect_match(post, '"--report-mode", presentation_mode', fixed = TRUE)
  # The on-demand API keeps its own gallery (D2).
  on_demand <- paste(deparse(body(lisaR:::render_lisa_categories)), collapse = "\n")
  expect_match(on_demand, "lisa_extension_write_index(", fixed = TRUE)
  expect_match(on_demand, "lisa_extension_render_units(", fixed = TRUE)
})
