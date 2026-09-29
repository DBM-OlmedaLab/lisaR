\name{run_lisa_contrast}
\alias{run_lisa_contrast}
\title{Compare two completed LISA analyses}
\usage{
run_lisa_contrast(contrast_a, contrast_b, output_dir,
  comparison_name = NULL, contrast_a_label = "Contrast A",
  contrast_b_label = "Contrast B",
  universes = lisa_default_collections(), ...)
}
\arguments{
\item{contrast_a, contrast_b}{Compatible completed single-analysis result directories.}
\item{output_dir}{Directory in which to write the contrast.}
\item{comparison_name}{Contrast identifier.}
\item{contrast_a_label, contrast_b_label}{Labels shown in tables and plots.}
\item{universes}{Gene-set collections to compare. The five options are
\code{"GOBP-C2"}, \code{"GOMF"}, \code{"GOCC"}, \code{"PATHWAYS"} and
\code{"HALLMARKS"}.}
\item{...}{Additional presentation and export arguments passed to the contrast engine.
\code{category_nes_variants} selects a nonempty unique subset of
\code{c("clean", "percentages", "direction", "dispersion")}, defaulting to all four.
Medians, directional percentages and P25--P75 are calculated separately from
the significant gene sets on each side. The existing A-minus-B delta is not
recalculated; no confidence interval or delta percentile is inferred.
Exact annotated GSEA member tables are required for distribution variants.
Legacy summary-only inputs retain their existing plots and record unavailable
variants explicitly in the supplementary NES-variant manifest.}
}
\value{A list containing \code{comparison_name}, \code{output_dir},
\code{universe}, the contrast \code{summary}, and plot/input manifests.}
\description{Compares semantic category summaries from two compatible LISA analyses.}
\examples{
# A contrast consumes two completed single-analysis directories. The
# sample project is fully offline; execution is interactive because the
# prerequisite analyses create scientific results and figures.
if (interactive()) {
  project <- lisa_init_project(tempfile("lisa-contrast-project-"))
  result <- run_lisa(file.path(project, "study.yml"))
  analyses <- utils::read.delim(
    file.path(result$output_dir, "single_de_status.tsv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  contrast <- run_lisa_contrast(
    contrast_a = analyses$output_dir[analyses$analysis_id == "response_a"],
    contrast_b = analyses$output_dir[analyses$analysis_id == "response_b"],
    output_dir = tempfile("lisa-contrast-example-"),
    comparison_name = "response_a_vs_b",
    contrast_a_label = "Response A",
    contrast_b_label = "Response B",
    universes = "GOBP-C2",
    export_formats = "tsv",
    verbose = FALSE
  )
  contrast$output_dir
  head(contrast$summary)
}
}
