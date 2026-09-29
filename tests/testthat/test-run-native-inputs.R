test_that('DESeq2 LRT objects accept explicit signed ranks but never the unsigned statistic', {
  skip_if_not_installed('DESeq2'); skip_if_not_installed('S4Vectors')
  table <- S4Vectors::DataFrame(baseMean=c(50,60,70),log2FoldChange=c(-2,1,0.5),
    lfcSE=c(.3,.2,.1),stat=c(12,10,5),pvalue=c(.001,.002,.03),padj=c(.003,.004,.03),
    signed_stat=c(-4,3,2),row.names=c('A','B','C'))
  S4Vectors::mcols(table)$description <- c('mean','effect','standard error',
    'LRT statistic: full vs reduced','p-value','adjusted p-value','signed Wald statistic')
  res <- DESeq2::DESeqResults(table)
  direction <- 'higher in treated than control'
  expect_error(lisaR:::lisa_run_collect_de(res,direction,NULL,NULL),'LISA-RUN-DE-006')
  expect_error(lisaR:::lisa_run_collect_de(res,direction,NULL,list(rank='stat')),'LISA-RUN-DE-006')
  prepared <- lisaR:::lisa_run_collect_de(res,direction,NULL,list(rank='signed_stat'))[[1]]
  expect_identical(prepared$data$rank_value,table$signed_stat)
  expect_identical(prepared$data$padj,table$padj)
})

test_that('edgeR native objects require caller-provided signed ranking and adjusted p-values', {
  table <- data.frame(logFC=c(-2,1,.5),PValue=c(.001,.002,.03),
    q=c(.003,.004,.03),signed_stat=c(-4,3,2),row.names=c('A','B','C'))
  for(cl in c('DGEExact','DGELRT','TopTags')) {
    result <- structure(list(table=table),class=cl)
    expect_error(lisaR:::lisa_run_collect_de(result,'treated-control',NULL,list(rank='signed_stat')),
      'LISA-RUN-DE-008.*topTags')
    prepared <- lisaR:::lisa_run_collect_de(result,'treated-control',NULL,list(rank='signed_stat',padj='q'))[[1]]
    expect_identical(prepared$data$rank_value,table$signed_stat)
    expect_identical(prepared$data$padj,table$q)
    expect_identical(prepared$data$log2FoldChange,table$logFC)
  }
})

test_that('Wald and limma inputs retain their published signed statistics', {
  skip_if_not_installed('DESeq2'); skip_if_not_installed('S4Vectors')
  frame <- S4Vectors::DataFrame(baseMean=c(50,60),log2FoldChange=c(-2,1),
    lfcSE=c(.3,.2),stat=c(-6,5),pvalue=c(.001,.002),padj=c(.002,.002),row.names=c('A','B'))
  S4Vectors::mcols(frame)$description <- c('mean','effect','SE','Wald statistic','p-value','adjusted p-value')
  res <- DESeq2::DESeqResults(frame)
  prepared <- lisaR:::lisa_run_collect_de(res,'treated-control',NULL,NULL)[[1]]
  expect_identical(prepared$data$rank_value,frame$stat)
  limma <- data.frame(gene=c('A','B'),logFC=c(-2,1),t=c(-6,5),P.Value=c(.001,.002),adj.P.Val=c(.002,.002))
  prepared <- lisaR:::lisa_run_collect_de(limma,'treated-control','limma',list(symbol='gene'))[[1]]
  expect_identical(prepared$data$rank_value,limma$t)
  expect_identical(prepared$data$padj,limma$adj.P.Val)
})
