test_that('incremental category refresh is scoped and preserves HTML outside its element', {
  e<-new.env(parent=globalenv());sys.source(system.file('scripts/build_LISA_report.R',package='lisaR'),e)
  original <- 'BEFORE<div><div>nested</div><span>unchanged elsewhere</span></div>AFTER'
  expect_identical(e$report_replace_element(original,7L,'div','<div>linked</div>'),
    'BEFORE<div>linked</div>AFTER')
  expect_error(e$report_replace_element('x<div>unclosed',2L,'div',''), 'Unclosed')
  long <- paste0('<div>',strrep('x',1100000),'</div>TAIL')
  expect_identical(e$report_replace_element(long,1L,'div','OK'),'OKTAIL')
  root<-tempfile();dir.create(file.path(root,'config'),recursive=TRUE)
  writeLines('contrast_id\toutput_id',file.path(root,'config/contrast_index.tsv'))
  expect_message(e$report_refresh_category_links(root),'No contrasts')
  expect_false(file.exists(file.path(root,'report_pages/contrasts.html')))
  unlink(root,recursive=TRUE)
})
