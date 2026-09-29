\name{plan_lisa_outputs}
\alias{plan_lisa_outputs}
\title{Estimate lisaR report outputs before execution}
\usage{
plan_lisa_outputs(config, category_counts = NULL, bytes_per_figure = 1024^2)
}
\arguments{
  \item{config}{A YAML/JSON path or equivalent named configuration list.}
  \item{category_counts}{Optional named expected category counts by collection.}
  \item{bytes_per_figure}{Conservative byte estimate per rendered format.}
}
\value{A list of summary, collection, product, and output-policy tables.}
\description{Estimates report file counts and storage before a run is started. In standard mode this includes the GSEA overview families, one evidence matrix and one member-set chart per applicable LISA category, and a linked SVG overview per individual analysis and semantic collection. Selected NES variants, existing contrast profile families and paired evidence matrices for semantic contrasts are also included; direct HALLMARKS are not counted as LISA evidence categories. Shared gene-explorer size depends on membership and DE cardinality and is not included in the figure-size estimate. With report.category_evidence = FALSE, legacy member-gene-set figures replace the evidence matrices.}
\examples{
config <- system.file("examples", "minimal-study.json", package = "lisaR")
if (nzchar(config)) {
  estimate <- plan_lisa_outputs(config, category_counts = c(`GOBP-C2` = 25))
  estimate$summary
}
}
