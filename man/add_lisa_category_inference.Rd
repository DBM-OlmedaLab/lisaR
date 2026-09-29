\name{add_lisa_category_inference}
\alias{add_lisa_category_inference}
\title{Add category inference to a new copy of saved LISA results}
\usage{
add_lisa_category_inference(source_run, output_dir, term2gene = NULL,
  plot_formats = c("png", "svg"))
}
\arguments{
\item{source_run}{Existing complete LISA run directory.}
\item{output_dir}{New destination directory; never an existing run.}
\item{term2gene}{Exact TERM2GENE resource used by the source run. Its saved
checksum is verified. With NULL, use the saved pipeline metadata path.}
\item{plot_formats}{Formats to generate: png, svg, pdf.}
}
\value{Invisibly, the new run directory. Original fgsea provenance is preserved;
the new inference records its input hashes and method separately.}
\description{Reuse raw fgsea results without rerunning enrichment. The source
run is untouched and the destination must not exist. Rebuild the report in the
copy with the installed scripts/build_LISA_report.R. New run_lisa() runs include
this layer automatically. HALLMARKS and descriptive DE-vs-DE comparisons do not
receive category inference.}
\details{Robust closed testing uses hommel(simes=FALSE), one family per analysis
across applicable LISA category collections, and all eligible raw gene-set
P values including unclassified sets. Missing eligible P values become 1.
Identical hypotheses repeated across collections contribute their maximum
P value once. A category is significant at adjusted P <= 0.05, equivalent to
at least one member supported by the simultaneous 95 percent lower count.
The count is obtained by inversion of the same public closed partial-conjunction
tests to preserve this equivalence at machine-precision boundaries. NES
direction remains descriptive.}
\references{Goeman et al. (2019), Biometrika 106, 841--856.
\doi{10.1093/biomet/asz041}}
