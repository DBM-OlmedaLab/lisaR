hommel_ui_fixture <- function() {
  n <- c(rep(100L, 6), 0L)
  d <- c(1L, 10L, 25L, 50L, 75L, 0L, NA_integer_)
  data.frame(analysis_id = "PRE", collection = "GOCC", category_id = paste0("CAT", seq_along(n)),
    category_display_name = c(rep("Duplicate label", 2), paste("Category", 3:7)),
    minimum_enriched_sets = d, n_sets_evaluable = n, minimum_enriched_pct = 100*d/n,
    n_sets_total = 120L, n_sets_missing_p = 0L, n_sets_not_eligible_or_absent = 120L-n,
    category_p_adjusted = c(rep(.01234567890123456, 5), .050000000001, NA_real_),
    confidence_level = .95, significant = !is.na(d) & d > 0,
    status = c(rep("significant", 5), "not_significant", "not_evaluable"),
    n_positive_NES = n, n_negative_NES = 0L, n_zero_NES = 0L, n_missing_NES = 0L)
}

hommel_ui_source <- function() {
  z <- hommel_ui_fixture()
  x <- data.frame(category_id = z$category_id, category_display_name = z$category_display_name,
    macrogroup_id = c(rep("M1",3),rep("M2",4)), macrogroup_name = c(rep("Group 1",3),rep("Group 2",4)), macrogroup_order = c(rep(1L,3),rep(2L,4)),
    category_order_within_macrogroup = seq_len(nrow(z)), mean_NES = c(1, -1, 2, -2, 3, NA, NA),
    n_genesets = c(2, 3, 4, 5, 6, 0, 0), color = "#447788")
  p <- tempfile(fileext = ".tsv"); on.exit(unlink(p))
  lisaR:::lisa_write_lollipop_figure_source(x, p, TRUE, "supracategory", "Fixture", figure_dpi = 40)
  lisaR:::lisa_read_figure_source_tsv(p)
}

test_that("stars use inclusive adjusted-P thresholds without rounding", {
  p <- c(0, .0001 * (1 - 1e-14), .0001, .0001 * (1 + 1e-14),
    .001, .001 * (1 + 1e-14), .01, .01 * (1 + 1e-14),
    .05, .05 * (1 + 1e-14), 1, NA_real_, NaN)
  expect_identical(lisaR:::lisa_category_significance_stars(p),
    c("****", "****", "****", "***", "***", "**", "**", "*", "*", "", "", "", ""))
  expect_error(lisaR:::lisa_category_significance_stars(c(-.1, .1)), "[0,1]", fixed = TRUE)
  expect_error(lisaR:::lisa_category_significance_stars(1.01), "[0,1]", fixed = TRUE)
  expect_error(lisaR:::lisa_category_significance_stars(Inf), "[0,1]", fixed = TRUE)
  expect_error(lisaR:::lisa_category_significance_stars(".001"), "numeric", fixed = TRUE)
  x <- hommel_ui_fixture(); x$category_p_adjusted <- p[1:7]
  # Only presentation changes: source P, N, d and significance are copied.
  x$significant <- x$category_p_adjusted <= .05
  x$minimum_enriched_sets <- ifelse(x$significant, 1L, 0L)
  x$minimum_enriched_sets[7] <- NA_integer_
  x$significant[7] <- FALSE
  x$category_p_adjusted[7] <- NA_real_
  source <- lisaR:::lisa_hommel_lollipop_source(hommel_ui_source(), x, "PRE", "GOCC")
  built <- ggplot2::ggplot_build(lisaR:::lisa_lollipop_from_source(source)$plot)
  expect_setequal(built$data[[4]]$label, c("****", "***", "**", ""))
  expect_identical(source$category_p_adjusted, x$category_p_adjusted)
})

test_that("saved category-inference figure uses P stars rather than support grades", {
  categories <- hommel_ui_fixture()
  categories$category_p_adjusted <- c(.0001, .001, .01, .05, .0500000000001, 1, NA_real_)
  categories$minimum_enriched_sets <- c(1L, 1L, 1L, 1L, 0L, 0L, NA_integer_)
  categories$significant <- c(rep(TRUE, 4), rep(FALSE, 3))
  categories$status <- c(rep("significant", 4), rep("not_significant", 2), "not_evaluable")
  members <- data.frame(analysis_id = "PRE", collection = "GOCC", category_id = categories$category_id,
    pathway = paste0("set", seq_len(nrow(categories))), NES = seq_len(nrow(categories)), pval = .01,
    hypothesis_id = paste0("set", seq_len(nrow(categories))), included_in_family = TRUE,
    raw_p_available = TRUE)
  source <- lisaR:::lisa_category_inference_source(categories, members)
  built <- ggplot2::ggplot_build(lisaR:::lisa_plot_category_inference(source))
  expect_identical(built$data[[2L]]$label, c("****", "***", "**", "*", "", "", ""))
  expect_false(any(grepl("S[1-5]|/", built$data[[2L]]$label)))
  expect_identical(source$category_p_adjusted[!duplicated(source$category_id)], categories$category_p_adjusted)
  display <- lisaR:::lisa_category_inference_display(categories)
  expect_true("Minimum enriched / evaluable sets (d/N)" %in% names(display))
})

test_that("large member NES values retain ticks in the separate descriptive panel", {
  categories <- hommel_ui_fixture()
  members <- data.frame(analysis_id = "PRE", collection = "GOCC",
    category_id = categories$category_id, pathway = paste0("set", seq_len(nrow(categories))),
    NES = c(31, 33, 35, 37, 39, 41, 43), pval = .01,
    hypothesis_id = paste0("set", seq_len(nrow(categories))),
    included_in_family = TRUE, raw_p_available = TRUE)
  source <- lisaR:::lisa_category_inference_source(categories, members)
  built <- ggplot2::ggplot_build(lisaR:::lisa_plot_category_inference(source))
  layout <- built$layout$layout
  for (i in seq_len(nrow(layout))) {
    breaks <- built$layout$panel_params[[i]]$x$get_breaks()
    if (as.character(layout$panel[[i]]) == "Category significance")
      expect_length(breaks, 0L)
    else
      expect_true(any(is.finite(breaks)))
  }
})

test_that("Hommel attachment uses three exact keys and the raw-P denominator", {
  table <- hommel_ui_fixture()
  category <- data.frame(category_id = rev(table$category_id), n_evaluable_sets = 2L,
    n_significant_sets = 1L)
  view <- structure(list(metadata = list(analysis_id = "PRE", collection = "GOCC"), categories = category),
    class = "lisa_category_evidence")
  z <- lisaR:::lisa_hommel_attach(view, table)
  expect_identical(z$categories$category_id, category$category_id)
  expect_equal(z$categories$n_sets_evaluable, rev(table$n_sets_evaluable))
  expect_equal(z$categories$n_evaluable_sets, rep(2L, 7))
  expect_equal(z$categories$n_significant_sets, rep(1L, 7))
  expect_identical(z$categories$category_p_adjusted, rev(table$category_p_adjusted))
  expect_identical(z$categories$support_grade, c("", "", "S5", "S4", "S3", "S2", "S1"))
  for (field in c("analysis_id", "collection")) {
    bad <- table; bad[[field]][1] <- "OTHER"
    expect_error(lisaR:::lisa_hommel_attach(view, bad), "scope mismatch")
  }
  bad <- table; bad$category_id[1] <- "UNKNOWN"
  expect_error(lisaR:::lisa_hommel_attach(view, bad), "no exact Hommel")
  expect_error(lisaR:::lisa_hommel_attach(view, rbind(table, table[1, ])), "Duplicate")
  view$metadata$contrast_id <- "ON-PRE"
  expect_error(lisaR:::lisa_hommel_attach(view, table), "analysis contrast")
  view$metadata$contrast_id <- NULL
  bad <- table; bad$collection <- "HALLMARKS"
  expect_error(lisaR:::lisa_hommel_attach(view, bad), "HALLMARKS")
})

test_that("new figure shows only adjusted-P stars and unchanged NES geometry", {
  old <- hommel_ui_source()
  z <- lisaR:::lisa_hommel_lollipop_source(old, hommel_ui_fixture(), "PRE", "GOCC")
  expect_identical(z$mean_NES, old$mean_NES)
  expect_identical(z$n_genesets, old$n_genesets)
  expect_identical(z$color, old$color)
  expect_identical(unique(z$annotation_variant), "hommel_support_v1")
  expect_identical(unique(z$figure_width), 13)
  rendered <- lisaR:::lisa_lollipop_from_source(z)
  b <- ggplot2::ggplot_build(rendered$plot)
  expect_equal(sum(vapply(b$layout$panel_scales_y, function(z) length(z$range$range), integer(1))), 7L)
  marks <- b$data[[4L]]
  expect_equal(nrow(marks), 7L)
  expect_identical(marks$label, c(rep("*", 5), "", ""))
  expect_equal(nrow(b$data[[7L]]), 5L)
  expect_setequal(b$data[[7L]]$x, old$mean_NES[is.finite(old$mean_NES)])
  expect_false(any(grepl("NE|ns|S[1-5]", marks$label)))
  expect_match(rendered$plot$labels$subtitle, "Adjusted category P", fixed = TRUE)
  expect_match(rendered$plot$labels$caption, "robust Hommel closed testing", fixed = TRUE)
  expect_match(rendered$plot$labels$caption, "10.1093/biomet/asz041", fixed = TRUE)
  expect_gte(rendered$plot$theme$axis.text.y$size, 11)
  expect_lt(b$layout$panel_scales_x[[1]]$range$range[[2]], 40)
  nes_range <- b$layout$panel_params[[2]]$x.range
  expect_gt(min(b$data[[7L]]$x) - nes_range[[1]], 0.06 * diff(range(b$data[[7L]]$x)))
  expect_gt(nes_range[[2]] - max(b$data[[7L]]$x), 0.06 * diff(range(b$data[[7L]]$x)))
  # All scientific means absent: still retain the full support catalogue.
  z$mean_NES <- NA_real_; z$n_genesets <- 0L
  b <- ggplot2::ggplot_build(lisaR:::lisa_lollipop_from_source(z)$plot)
  expect_equal(nrow(b$data[[4L]]), 7L)
  expect_equal(nrow(b$data[[7L]]), 0L)
})

test_that("significance facets have no numerical scale while NES retains its scale", {
  z <- lisaR:::lisa_hommel_lollipop_source(hommel_ui_source(), hommel_ui_fixture(), "PRE", "GOCC")
  # Large NES values must not be mistaken for significance layout coordinates.
  for (grouped in c(FALSE, TRUE)) for (magnitude in c(1, 30)) {
    x <- z; x$mean_NES <- x$mean_NES * magnitude
    built <- ggplot2::ggplot_build(lisaR:::lisa_plot_hommel_support(x, grouped, "supracategory"))
    layout <- built$layout$layout
    for (i in seq_len(nrow(layout))) {
      axis <- built$layout$panel_params[[i]]$x
      if (as.character(layout$panel[[i]]) == "Category significance") {
        expect_length(axis$get_breaks(), 0L)
        expect_length(axis$get_breaks_minor(), 0L)
      } else {
        expect_true(any(is.finite(axis$get_breaks())))
      }
    }
    expect_setequal(built$data[[7L]]$x, x$mean_NES[is.finite(x$mean_NES)])
  }
})

test_that("navigator keeps every row with stars and blank nonsignificant cells", {
  x <- hommel_ui_source(); x$n_genesets_significant <- x$n_genesets
  x$n_pos_genesets <- ifelse(x$mean_NES > 0, x$n_genesets, 0); x$n_pos_genesets[is.na(x$n_pos_genesets)] <- 0
  x$n_neg_genesets <- x$n_genesets - x$n_pos_genesets
  view <- lisaR:::build_lisa_category_navigation(x, "PRE", "GOCC", "core")
  view <- lisaR:::lisa_hommel_attach(view, hommel_ui_fixture())
  out <- tempfile("hommel-nav-"); on.exit(unlink(out, recursive = TRUE))
  lisaR:::render_lisa_category_navigation(view, out, "detail/index.html")
  svg <- paste(readLines(file.path(out,"category_navigation.svg")), collapse = "\n")
  expect_match(svg, 'x="1240"', fixed=TRUE)
  expect_match(svg, 'x="1012"', fixed=TRUE)
  expect_equal(length(gregexpr('class="hommel-support-cell"', svg, fixed=TRUE)[[1]]), 7L)
  cells <- regmatches(svg, gregexpr('<g class="hommel-support-cell"[^>]*>.*?</g>', svg, perl=TRUE))[[1]]
  expect_length(cells, 7L)
  expect_true(all(vapply(cells[1:5], function(cell) length(regmatches(cell,
    gregexpr('<text ',cell,fixed=TRUE))[[1]]) == 1L, logical(1))))
  expect_true(all(grepl('></g>$', cells[6:7])))
  expect_match(cells[[6]], 'data-support-status="not_significant"', fixed=TRUE)
  expect_match(cells[[7]], 'data-support-status="not_evaluable"', fixed=TRUE)
  expect_match(cells[[1]], '>*</text>', fixed=TRUE)
  expect_false(grepl('S[1-5]|1/100|≥1%', svg))
  expect_false(grepl('NA NA|NE · not evaluable|No positive minimum ·', svg))
  expect_match(svg, 'category significance: not evaluable', fixed=TRUE)
  expect_match(svg, 'category significance: not significant', fixed=TRUE)
  expect_match(svg, 'category=CAT7', fixed=TRUE)
  saved <- lisaR:::read_lisa_tsv(file.path(out,"tables","category_navigation_source.tsv"))
  expect_identical(saved$category_p_adjusted, hommel_ui_fixture()$category_p_adjusted)
  expect_identical(saved$status, hommel_ui_fixture()$status)
  expect_equal(saved$minimum_enriched_sets, hommel_ui_fixture()$minimum_enriched_sets)
  expect_equal(saved$n_sets_evaluable, hommel_ui_fixture()$n_sets_evaluable)
  env <- new.env(parent = baseenv())
  expr <- parse(file.path(out,"reproduce_category_navigation.R"))
  eval(expr[[1L]],env)
  source <- utils::read.delim(file.path(out,"tables","category_navigation_source.tsv"),colClasses="character",quote="",na.strings=NULL)
  # The portable SVG renderer is self-contained (no lisaR helper calls).
  replay <- env$lisa_category_navigation_svg(source,"Title","Subtitle")
  expect_match(replay, 'category significance: significant', fixed=TRUE)
})

test_that("portable navigator and its recipe preserve the four-star boundary", {
  x <- hommel_ui_source(); x$n_genesets_significant <- x$n_genesets
  x$n_pos_genesets <- ifelse(is.na(x$mean_NES), 0, x$n_genesets)
  x$n_neg_genesets <- 0L
  nav <- lisaR:::build_lisa_category_navigation(x, "PRE", "GOCC", "core")
  table <- hommel_ui_fixture()
  table$category_p_adjusted[1:4] <- c(.0001, .0001 * (1 + 1e-14), .001, .01)
  nav <- lisaR:::lisa_hommel_attach(nav, table)
  out <- tempfile("hommel-four-star-nav-"); on.exit(unlink(out, recursive = TRUE))
  lisaR:::render_lisa_category_navigation(nav, out, "detail/index.html")
  svg <- paste(readLines(file.path(out, "category_navigation.svg")), collapse = "\n")
  expect_match(svg, 'data-significance-stars="****"', fixed = TRUE)
  expect_match(svg, '**** ≤ 0.0001', fixed = TRUE)
  recipe <- new.env(parent = baseenv())
  eval(parse(file.path(out, "reproduce_category_navigation.R"))[[1L]], recipe)
  source <- utils::read.delim(file.path(out, "tables", "category_navigation_source.tsv"),
    colClasses = "character", quote = "", na.strings = NULL)
  replay <- recipe$lisa_category_navigation_svg(source, "Title", "Subtitle")
  expect_match(replay, 'data-significance-stars="****"', fixed = TRUE)
})

test_that("shared help explains corrected P and support without grades", {
  h <- lisaR:::lisa_hommel_help_html()
  for (text in c('id="hommel-support-help"', 'aria-labelledby="hommel-help-title"', 'data-hommel-close',
    'invented example, not a study result', 'at least 5 of the 20 sets', '5/20',
    'multiple-testing correction', 'false positives', 'at least one',
    'within one analysis', 'not exactly 5', 'more may be enriched',
    'd=0', 'N=0', 'does not mean no enrichment',
    'effect size', 'biological process is activated', 'difference between conditions',
    'robust Hommel closed testing', '10.1093/biomet/asz041'))
    expect_match(h, text, fixed=TRUE)
  for (removed in c('1/1', 'PRE', 'ON&minus;PRE', 'HALLMARKS', 'editorial', 'truncat', 'S1', 'S2', 'S3',
    'partial-conjunction', 'hommel::', 'simes=FALSE', '&mdash;', '&ndash;'))
    expect_false(grepl(removed, h, fixed=TRUE), info=removed)
  html <- '<html><head></head><body><main>Test</main></body></html>'
  once <- lisaR:::lisa_hommel_install_help(html)
  expect_identical(lisaR:::lisa_hommel_install_help(once), once)
})

test_that("refresh keeps large evidence payloads exact and updates category-only metadata", {
  root <- tempfile("hommel-large-payload-"); on.exit(unlink(root,recursive=TRUE))
  dir.create(file.path(root,"tables"),recursive=TRUE)
  table <- hommel_ui_fixture()
  payload <- list(metadata=list(analysis_id="PRE",collection="GOCC",tier="core"),
    categories=data.frame(category_id=table$category_id,n_evaluable_sets=2L),
    sentinel=strrep("evidence-data-is-preserved-",50000L))
  json <- jsonlite::toJSON(payload,auto_unbox=TRUE,dataframe="rows",digits=17)
  page <- file.path(root,"index.html")
  writeLines(paste0('<html><head></head><body><main></main><script type="application/json" id="evidence-data">',json,'</script></body></html>'),page)
  lisaR:::lisa_refresh_hommel_evidence_page(page,table)
  html <- paste(readLines(page,warn=FALSE),collapse="\n")
  marker <- '<script type="application/json" id="evidence-data">'
  start <- regexpr(marker,html,fixed=TRUE)[[1L]] + nchar(marker)
  tail <- substr(html,start,nchar(html)); end <- regexpr("</script>",tail,fixed=TRUE)[[1L]]
  z <- jsonlite::fromJSON(substr(tail,1,end-1))
  expect_identical(z$sentinel,payload$sentinel)
  expect_identical(z$categories$category_p_adjusted,table$category_p_adjusted)
  expect_equal(z$categories$n_evaluable_sets,rep(2L,nrow(table)))
  expect_identical(z$metadata$hommel_support_schema,"hommel-support-ui-v1")
  saved <- lisaR:::read_lisa_tsv(file.path(root,"tables","categories.tsv"))
  expect_identical(saved$category_p_adjusted,table$category_p_adjusted)
  expect_match(html,'id="hommel-support-help"',fixed=TRUE)
})

test_that("Hommel source reproduces PNG exactly and emits SVG and PDF", {
  skip_if_not_installed("png")
  root <- tempfile("hommel-replay-"); dir.create(root); on.exit(unlink(root,recursive=TRUE))
  z <- lisaR:::lisa_hommel_lollipop_source(hommel_ui_source(),hommel_ui_fixture(),"PRE","GOCC")
  source <- file.path(root,"source.tsv"); lisaR:::lisa_write_inference_tsv(z,source)
  saved <- lisaR:::lisa_read_figure_source_tsv(source)
  expect_identical(saved$category_p_adjusted,z$category_p_adjusted)
  for (format in c("png","svg","pdf")) lisaR:::lisa_support_save_lollipop(saved,file.path(root,paste0("original.",format)))
  e <- new.env(parent=globalenv()); e$output_png <- file.path(root,"replay.png")
  for (expr in parse(system.file("scripts","reproduce_lisa_figure.R",package="lisaR")))
    if (is.call(expr) && identical(expr[[1]],as.name("<-")) && identical(expr[[2]],as.name("render_lisa_lollipop"))) eval(expr,e)
  e$render_lisa_lollipop(saved)
  expect_equal(png::readPNG(e$output_png),png::readPNG(file.path(root,"original.png")),tolerance=0)
  expect_match(paste(readLines(file.path(root,"original.svg"),warn=FALSE),collapse="\n"),"<svg",fixed=TRUE)
  con <- file(file.path(root,"original.pdf"),"rb"); magic <- readChar(con,4); close(con)
  expect_identical(magic,"%PDF")
})

test_that("report defaults to Hommel, removes the technical section and retains compact downloads", {
  e <- new.env(parent=globalenv())
  sys.source(system.file("scripts","build_LISA_report.R",package="lisaR"),e)
  root <- tempfile("hommel-report-"); on.exit(unlink(root,recursive=TRUE))
  directory <- file.path(root,"outputs","single_de","PRE","collection_GOCC")
  dir.create(file.path(directory,"plots"),recursive=TRUE); dir.create(file.path(directory,"lisa_tables"))
  write <- lisaR:::lisa_write_inference_tsv
  write(hommel_ui_fixture(),file.path(directory,"lisa_tables","A_LISA_category_inference.tsv"))
  for (file in c("A_GSEA_lollipop_hommel_support.png","A_LISA_category_inference.png")) file.create(file.path(directory,"plots",file))
  e$report_image_formats <- function() "png"
  e$figure_card <- function(path, ...) paste0('<figure>',basename(path),'</figure>')
  e$table_html <- function(x,...) paste(unlist(x),collapse="|")
  e$short_file_copy <- function(path,...) basename(path)
  e$short_media_copy <- function(path,...) basename(path)
  e$figure_source_contract <- function(...) list(source_tsv="saved-figure.tsv.gz", recipe_r="saved-figure.R")
  page <- file.path(root,"report_pages","single_de.html")
  view <- e$report_nes_variants_section(directory,page,"Fixture")
  expect_match(view,'<option value="hommel_support" selected>',fixed=TRUE)
  expect_match(view,'data-nes-panel="hommel_support"',fixed=TRUE)
  expect_match(view,'Category significance &middot; adjusted P',fixed=TRUE)
  expect_match(view,'10.1093/biomet/asz041',fixed=TRUE)
  expect_identical(e$report_hommel_view_panel(file.path(root,"outputs","contrasts","ON-PRE","collection_GOCC"),page,"x"),"")
  html <- e$report_lisa_categories_section("Fixture",directory,character(),page)
  for (removed in c('hommel-support-technical', 'Category enrichment: minimum supported extent',
    'support-results-table', 'editorial', 'truncated', '1/1')) expect_false(grepl(removed,html,fixed=TRUE))
  expect_match(html, 'data-category-results-tsv', fixed=TRUE)
  expect_match(html, 'Complete category results (TSV)', fixed=TRUE)
  expect_match(html, 'saved-figure.tsv.gz', fixed=TRUE)
  expect_match(html, 'saved-figure.R', fixed=TRUE)
  expect_false(grepl('<table|<img[^>]+category_inference',html))
  file.create(file.path(directory, "plots", "A_GSEA_lollipop_support_source.tsv"))
  file.create(file.path(directory, "plots", "A_GSEA_lollipop_support.png"))
  html <- e$report_lisa_categories_section("Fixture",directory,character(),page)
  expect_match(html, 'class="support-historical-lollipops"', fixed=TRUE)
  expect_match(html, 'A_GSEA_lollipop_support.png', fixed=TRUE)
  expect_false(grepl('hommel-support-technical',html,fixed=TRUE))
  evidence <- file.path(root, "evidence.html"); file.create(evidence)
  e$report_inline_category_navigation <- function(...) ""
  cover <- e$report_category_evidence_section("Fixture evidence", file.path(root,"nav"),
    evidence, page, collection_dir=directory)
  expect_match(cover, 'class="panel category-evidence-cover"', fixed=TRUE)
  expect_match(cover, 'data-category-results-tsv', fixed=TRUE)
  section <- e$report_lisa_categories_section("Fixture",directory,character(),page,evidence_index=evidence)
  expect_false(grepl('data-category-results-tsv',section,fixed=TRUE))
  expect_match(section, 'class="support-historical-lollipops"', fixed=TRUE)
})
