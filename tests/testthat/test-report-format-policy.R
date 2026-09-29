format_report_env <- function(formats, tables = TRUE, recipes = TRUE) {
  env <- new.env(parent = globalenv())
  sys.source(system.file('scripts/build_LISA_report.R', package = 'lisaR'), env)
  env$report_output_policy <- lisaR:::lisa_output_policy_table(formats, tables, recipes)
  env
}

test_that('figure discovery is independent of PNG and deduplicates twins', {
  root <- tempfile('plot-discovery-'); dir.create(root)
  invisible(file.create(file.path(root, c('one.png', 'one.svg', 'two.svg', 'three.pdf', 'no.tsv'))))
  expect_setequal(basename(lisaR:::lisa_plot_files(root)), c('one.png','two.svg','three.pdf'))
  expect_setequal(basename(lisaR:::lisa_plot_files(root, formats='svg')), c('one.svg','two.svg'))
  expect_length(lisaR:::lisa_plot_files(root, formats=character()), 0)
  for (suffix in c('[.]png$', '\\.png$')) {
    expect_identical(basename(lisaR:::lisa_plot_files(root, paste0('three',suffix))), 'three.pdf')
  }
  lisaR:::write_lisa_tsv(lisaR:::lisa_output_policy_table(character()), file.path(root,'report_output_policy.tsv'))
  expect_identical(lisaR:::lisa_requested_plot_formats(root), character())
})

test_that('figure cards obey requested formats including PDF-only and no images', {
  root <- tempfile('format-cards-'); dir.create(root)
  figure <- file.path(root,'plot')
  grDevices::svg(paste0(figure,'.svg')); graphics::plot(1:3); grDevices::dev.off()
  grDevices::pdf(paste0(figure,'.pdf')); graphics::plot(1:3); grDevices::dev.off()
  page <- file.path(root,'page.html'); file.create(page)
  for (format in c('svg','pdf')) {
    env <- format_report_env(format, FALSE, FALSE)
    env$report_project_dir <- root
    html <- env$figure_card(paste0(figure,'.',format),page)
    expect_match(html,paste0('>',toupper(format),'</a>'),fixed=TRUE)
    expect_false(grepl('>PNG</a>',html,fixed=TRUE))
    expect_false(grepl('Source data|R script',html))
    if(format=='pdf') expect_match(html,'type="application/pdf"',fixed=TRUE)
  }
  env <- format_report_env(character(), FALSE, FALSE)
  expect_identical(env$figure_card(paste0(figure,'.svg'),page),'')
  expect_match(env$layer_section('Plots',character(),page,''),'not_requested',fixed=TRUE)
  env <- format_report_env(c('png','svg'), FALSE, FALSE)
  expect_error(env$figure_card(paste0(figure,'.svg'),page),'Requested figure format missing')
  env <- format_report_env('png', FALSE, FALSE)
  expect_error(env$figure_card(paste0(figure,'.svg'),page),'Requested figure format missing')
})

test_that('asset exports and required FULL cards follow the run policy', {
  for(tables in c(TRUE,FALSE)) for(recipes in c(TRUE,FALSE)) {
    env <- format_report_env('svg',tables,recipes)
    paths <- c('a.png','a.svg','a.pdf','a.tsv','a_recipe.R','a_settings.json')
    selected <- Filter(env$report_asset_requested,paths)
    expect_setequal(selected,c('a.svg',if(tables) 'a.tsv',if(recipes) c('a_recipe.R','a_settings.json')))
    formats <- c('svg',if(tables) 'tsv',if(recipes) 'r')
    attachment <- list(category_products=list(list(product='contrast_heatmap',assets=lapply(formats,function(f) list(format=f)))))
    expect_true(env$report_has_complete_category_product(attachment,'contrast_heatmap'))
  }
})

test_that('saved evidence pages inherit the containing run export policy', {
  root <- tempfile('evidence-policy-');dir.create(root)
  page <- file.path(root,'index.html')
  for(contrast in c(TRUE,FALSE)) {
    id <- if(contrast) 'contrast-evidence-data' else 'evidence-data'
    payload <- list(category_products=list(),recipes=TRUE,
      member_figure_index=list(list(category_id='A',source_tsv='a.tsv',recipe_r='a.R')))
    writeLines(paste0('<script type="application/json" id="',id,'">',
      jsonlite::toJSON(payload,auto_unbox=TRUE),'</script>',
      '<a href="reproduce_contrast_evidence.R" download>R script</a>'),page)
    env <- format_report_env('png',FALSE,FALSE)
    env$report_replace_evidence_payload(page,NULL,id)
    html <- paste(readLines(page,warn=FALSE),collapse='\n')
    expect_match(html,'"source_data":false',fixed=TRUE)
    expect_match(html,'"recipes":false',fixed=TRUE)
    expect_false(grepl('a.tsv|a.R',html,fixed=FALSE))
    expect_false(grepl('href="reproduce_contrast_evidence.R"',html,fixed=TRUE))
  }
})

test_that('FULL galleries without category sheets keep all requested assets manifest-bound', {
  root <- tempfile('full-gallery-'); dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  page <- file.path(root, 'report_pages', 'single_de.html')
  dir.create(dirname(page), recursive = TRUE)
  art <- file.path(root, 'artifacts', 'A', 'HALLMARKS', 'gene_cards')
  dir.create(art, recursive = TRUE)
  writeLines('<svg xmlns="http://www.w3.org/2000/svg"></svg>', file.path(art, 'H_SIGNAL.svg'))
  writeLines('symbol\tvalue\nA\t1', file.path(art, 'H_SIGNAL_source.tsv'))
  writeLines('plot(1)', file.path(art, 'H_SIGNAL_recipe.R'))
  for (enabled in c(TRUE, FALSE)) {
    e <- format_report_env('svg', enabled, enabled)
    title <- 'A - HALLMARKS - FULL figures'
    body <- e$report_full_gallery_section(title,
      list(GeneCards = list(path = art)), page, root)
    html <- e$page_shell('A', 'single_de', body, root, page)
    writeLines(html, page)
    expect_match(html, 'id="a-hallmarks-full-figures"', fixed = TRUE)
    expect_match(html, 'H_SIGNAL.svg', fixed = TRUE)
    expect_identical(grepl('H_SIGNAL_source.tsv', html, fixed = TRUE), enabled)
    expect_identical(grepl('H_SIGNAL_recipe.R', html, fixed = TRUE), enabled)
    refs <- e$report_category_product_reference_rows(root)
    expect_equal(nrow(refs), if (enabled) 3L else 1L)
    expect_true(all(file.exists(file.path(root, refs$source_path))))
    expect_true(all(grepl('^[0-9a-f]{64}$', refs$sha256)))
  }
})
