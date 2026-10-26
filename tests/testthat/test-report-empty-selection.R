empty_selection_run <- function(positive = FALSE) {
  root <- tempfile("public-empty-selection-")
  scope <- file.path(root, "outputs/single_de/A/collection_HALLMARKS")
  dir.create(file.path(scope, "enrichment"), recursive = TRUE)
  dir.create(file.path(scope, "lisa_tables"))
  dir.create(file.path(root, "config"))
  write <- lisaR:::write_lisa_tsv
  summary <- data.frame(category_id = "HALLMARK_SYNTHETIC", n_genesets = as.integer(positive),
    n_genesets_mapped=1L,n_genesets_evaluable=1L,n_genesets_significant=as.integer(positive),
    mean_NES = if (positive) 1.5 else NA_real_, gsea_padj_cutoff = .05)
  gsea <- data.frame(category_id = "HALLMARK_SYNTHETIC", pathway = "HALLMARK_SYNTHETIC",
    NES = 1.5, padj = if (positive) .01 else .8)
  s <- file.path(scope, "lisa_tables/A_hallmarks_GSEA_category_summary.tsv")
  g <- file.path(scope, "enrichment/A_GSEA_hallmarks_annotated.tsv")
  write(summary, s); write(gsea, g)
  write(data.frame(key = "gsea_padj_cutoff", value = .05), file.path(root, "contract_manifest.tsv"))
  write(data.frame(analysis_id="A", species="Homo sapiens"), file.path(root, "config/de_index.tsv"))
  write(data.frame(contrast_id=character(),output_id=character()), file.path(root, "config/contrast_index.tsv"))
  list(root=root, summary=s, gsea=g)
}

test_that("empty figure hit regions keep schema, hashes and source-identity checks", {
  dir <- tempfile("public-empty-regions-"); dir.create(dir)
  x <- data.frame(category_id=character(),pathway=character(),NES=numeric(),padj=numeric())
  s <- data.frame(category_id=character(),n_genesets=integer(),mean_NES=numeric())
  v <- lisaR:::lisa_prepare_category_nes_variants(s,x,.05)
  out <- lisaR:::render_lisa_category_nes_variants(v,dir,variants="clean",formats="png",width=6,height=3,dpi=72)
  before <- lisaR:::lisa_sha256_file(out$path[[1]])
  regions <- lisaR:::lisa_category_figure_regions(out$source[[1]],out$settings[[1]])
  expect_equal(nrow(regions),0L)
  expect_identical(names(regions),c("category_id","kind","side","left","top","width","height"))
  expect_identical(lisaR:::lisa_sha256_file(out$path[[1]]),before)
  m <- jsonlite::read_json(out$settings[[1]],simplifyVector=TRUE)
  save_settings <- function(value) jsonlite::write_json(value,out$settings[[1]],auto_unbox=TRUE,digits=NA)
  invalid <- m; invalid$variant <- "invalid"; save_settings(invalid)
  expect_error(lisaR:::lisa_category_figure_regions(out$source[[1]],out$settings[[1]]),"NES variants")
  invalid <- m; invalid$distribution_layout <- "invalid"; save_settings(invalid)
  expect_error(lisaR:::lisa_category_figure_regions(out$source[[1]],out$settings[[1]]),"Unsupported NES distribution")
  invalid <- m; invalid$width <- -1; save_settings(invalid)
  expect_error(lisaR:::lisa_category_figure_regions(out$source[[1]],out$settings[[1]]),"Invalid empty NES-variant settings")
  source_bytes <- readBin(out$source[[1]],"raw",n=file.info(out$source[[1]])$size)
  lisaR:::write_lisa_tsv(data.frame(bogus=character()),out$source[[1]])
  invalid <- m; invalid$source_sha256 <- lisaR:::lisa_sha256_file(out$source[[1]]); save_settings(invalid)
  expect_error(lisaR:::lisa_category_figure_regions(out$source[[1]],out$settings[[1]]),"Incomplete empty NES-variant source")
  writeBin(source_bytes,out$source[[1]]); save_settings(m)
  cat("\n",file=out$source[[1]],append=TRUE)
  expect_error(lisaR:::lisa_category_figure_regions(out$source[[1]],out$settings[[1]]),"identity mismatch")
})

test_that("positive category geometry remains linked without changing its bitmap", {
  dir <- tempfile("public-positive-regions-"); dir.create(dir)
  x <- data.frame(category_id="SYNTHETIC",pathway="GS_SYNTHETIC",NES=1.5,padj=.01)
  s <- data.frame(category_id="SYNTHETIC",n_genesets=1L,mean_NES=1.5)
  v <- lisaR:::lisa_prepare_category_nes_variants(s,x,.05)
  out <- lisaR:::render_lisa_category_nes_variants(v,dir,variants="clean",formats="png",width=6,height=3,dpi=72)
  before <- lisaR:::lisa_sha256_file(out$path[[1]])
  regions <- lisaR:::lisa_category_figure_regions(out$source[[1]],out$settings[[1]])
  expect_gt(nrow(regions),0L)
  expect_true(all(regions$category_id=="SYNTHETIC"))
  expect_identical(lisaR:::lisa_sha256_file(out$path[[1]]),before)
})

test_that("FULL empty selections require exact saved evidence and a consistent recorded cutoff", {
  state <- lisaR:::lisa_report_selection_state
  f <- empty_selection_run()
  expect_identical(state(f$root,"A","HALLMARKS"),"completed_empty")
  expect_identical(state(empty_selection_run(TRUE)$root,"A","HALLMARKS"),"missing")
  expect_error(state(f$root,"B","HALLMARKS"),"missing or ambiguous")
  change <- function(which, column, value) {
    f <- empty_selection_run(); x <- lisaR:::read_lisa_tsv(f[[which]])
    x[[column]] <- value; lisaR:::write_lisa_tsv(x,f[[which]]); f
  }
  for (bad in list(change("summary","n_genesets",NA),
                   change("summary","n_genesets",1),
                   change("summary","gsea_padj_cutoff",.1),
                   change("summary","n_genesets_mapped",2),
                   change("gsea","padj","corrupt"),
                   change("gsea","padj",NA),
                   change("gsea","padj",1.1),
                   change("gsea","NES",NA),
                   change("gsea","analysis_id","B")))
    expect_error(state(bad$root,"A","HALLMARKS"),"LISA-REPORT-EMPTY-001")
  f <- empty_selection_run(); unlink(f$gsea)
  expect_error(state(f$root,"A","HALLMARKS"),"missing or ambiguous")
  f <- empty_selection_run(); x <- lisaR:::read_lisa_tsv(f$gsea)
  lisaR:::write_lisa_tsv(x[FALSE,,drop=FALSE],f$gsea)
  expect_error(state(f$root,"A","HALLMARKS"),"inconsistent n_genesets_mapped")
  f <- empty_selection_run(); x <- lisaR:::read_lisa_tsv(f$summary)
  lisaR:::write_lisa_tsv(x[FALSE,,drop=FALSE],f$summary)
  expect_error(state(f$root,"A","HALLMARKS"),"omits mapped members")
  f <- empty_selection_run(); file.copy(f$gsea,sub("hallmarks","other",f$gsea))
  expect_error(state(f$root,"A","HALLMARKS"),"missing or ambiguous")
  f <- empty_selection_run(); file.copy(f$summary,sub("hallmarks","other",f$summary))
  expect_error(state(f$root,"A","HALLMARKS"),"missing or ambiguous")
  f <- empty_selection_run()
  lisaR:::write_lisa_tsv(data.frame(key=c("gsea_padj_cutoff","gsea_padj_cutoff"),value=.05),file.path(f$root,"contract_manifest.tsv"))
  expect_error(state(f$root,"A","HALLMARKS"),"missing or ambiguous recorded cutoff")
  f <- empty_selection_run(); x <- lisaR:::read_lisa_tsv(f$summary); x$gsea_padj_cutoff <- NULL
  lisaR:::write_lisa_tsv(x,f$summary)
  expect_error(state(f$root,"A","HALLMARKS"),"invalid category-summary schema")
})

test_that("FULL family inventory records proven empty but rejects missing positive products", {
  f <- empty_selection_run()
  inventory <- lisaR:::lisa_full_product_inventory(f$root)
  expect_setequal(inventory$family,c("prioritized_genes","volcano","gene_heatmap"))
  expect_true(all(inventory$state=="completed_empty"))
  expect_error(lisaR:::lisa_full_product_inventory(empty_selection_run(TRUE)$root),"FULL required family absent")
  f <- empty_selection_run(); unlink(f$summary)
  expect_error(lisaR:::lisa_full_product_inventory(f$root),"missing or ambiguous")
})

test_that("the separate all-empty FULL catalog path retains its fail-closed behavior", {
  f <- empty_selection_run()
  old <- options(lisaR.run_root=f$root)
  on.exit(options(old),add=TRUE)
  report <- list(mode="full",formats=list(png=TRUE,svg=FALSE,pdf=FALSE),source_data=TRUE,recipes=TRUE)
  expect_error(lisaR:::lisa_render_run_full_products(f$root,report),"LISA-EXTENSION-009")
  expect_false(dir.exists(file.path(f$root,"artifacts")))
  options(lisaR.run_root=NULL)
  f <- empty_selection_run(); options(lisaR.run_root=f$root)
  lisaR:::write_lisa_tsv(data.frame(analysis_id=c("A","B"),species="Homo sapiens"),file.path(f$root,"config/de_index.tsv"))
  expect_error(lisaR:::lisa_render_run_full_products(f$root,report),"LISA-EXTENSION-009")
  expect_false(dir.exists(file.path(f$root,"artifacts")))
  options(lisaR.run_root=NULL)
  f <- empty_selection_run(TRUE); options(lisaR.run_root=f$root)
  expect_error(lisaR:::lisa_render_run_full_products(f$root,report),"LISA-EXTENSION-009")
})
