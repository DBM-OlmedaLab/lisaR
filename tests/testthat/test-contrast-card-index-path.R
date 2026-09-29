test_that("FULL card indexes remain portable without changing their rows or exact-card contract", {
  root <- withr::local_tempdir("card-index-")
  e <- new.env(parent = globalenv())
  script <- system.file("scripts/build_contrast_gene_level_product.R", package="lisaR")
  e$commandArgs <- function(trailingOnly=FALSE) paste0("--file=",script)
  sys.source(script,e)
  row <- data.frame(contrast_id="response_separation_by_time",
    output_id="response_separation_by_time",contrast_a="ON",contrast_b="PRE",
    label_a="ON",label_b="PRE")
  paired <- data.frame(category_id=c("CAT_A","CAT_A","CAT_B"),
    category_display_name=c("Category A","Category A","Category B"),
    symbol=c("GENE1","GENE2","GENE3"),present_A=TRUE,present_B=TRUE)
  summary <- data.frame(category_id=c("CAT_A","CAT_B"))
  e$contrast_row <- function(...) row
  e$prepare_contrast_tables <- function(...) list(paired=paired,summary=summary)
  # Test index creation, not plotting or enrichment. Marker files permit the
  # same checks on figure links without rendering scientific figures.
  e$plot_category_card <- function(df,out_png,out_pdf,cfg) {
    writeLines("test marker",out_png); writeLines("test marker",out_pdf); TRUE
  }
  e$write_card_contract <- function(...) invisible(NULL)
  e$plot_combined_heatmap <- function(...) FALSE
  parser <- e$parse_args
  for (collection in c("GOBP-C2","PATHWAYS")) {
    cfg <- parser(c("--project-dir",root,"--contrast-id",row$contrast_id,
      "--universe",collection))
    e$parse_args <- function(...) cfg
    invisible(capture.output(e$main()))
    paths <- e$paths_for(root,row,collection)
    indexes <- list.files(paths$cards_dir,pattern="index[.]tsv$",full.names=TRUE)
    expect_length(indexes,1L)
    # Include the temporary sibling, with the real Riaz root and long IDs.
    rel <- substring(indexes,nchar(root)+2L)
    physical <- file.path("D:/lisa test/riaz/out",rel)
    expect_invisible(lisaR:::lisa_assert_portable_path(c(physical,
      file.path(dirname(physical),".lisa-write-123456789abc.json")),windows=TRUE))
    index <- e$read_tsv(indexes)
    expect_identical(index$category_id,summary$category_id)
    expect_equal(index$rank,1:2)
    expect_equal(index$n_genes_total,c(2,1))
    expect_equal(index$n_genes_visible,c(2,1))
    expect_true(all(file.exists(c(index$png,index$pdf))))
    before <- tools::md5sum(indexes)
    cfg$category_id <- "CAT_A"
    e$exact_gene_card(cfg,paired,summary,paths)
    expect_true(file.exists(file.path(paths$cards_dir,"CAT_A_contrast_gene_card_index.tsv")))
    expect_identical(tools::md5sum(indexes),before)
  }
})
