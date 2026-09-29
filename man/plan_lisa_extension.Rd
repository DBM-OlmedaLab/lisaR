\name{plan_lisa_extension}
\alias{plan_lisa_extension}
\title{Plan additional selected or full report figures}
\usage{
plan_lisa_extension(source_run, selection, exact_products = NULL)
}
\arguments{
  \item{exact_products}{Optional named list of exact heatmap variants, KEGG maps or native products used by the exploration interface. Leave as \code{NULL} for the standard selected or full report.}
  \item{source_run}{A completed and checked lisaR run.}
  \item{selection}{The string \code{"full"}, a YAML/JSON selection path, or a named list.}
}
\value{A plan listing the requested figures and expected files.}
\description{Identifies selected or full category figures without repeating enrichment or changing the original run.}
\examples{
# Planning an extension requires a verified scientific run. The installed
# sample project is offline; the prerequisite run is kept interactive.
if (interactive()) {
  project <- lisa_init_project(tempfile("lisa-extension-project-"))
  result <- run_lisa(file.path(project, "study.yml"))
  selection <- list(
    selection = list(
      mode = "selected",
      categories = "SYN_SIGNAL",
      products = "member_gene_sets"
    ),
    report = list(
      mode = "selected",
      formats = list(png = TRUE, svg = FALSE, pdf = FALSE),
      source_data = TRUE,
      recipes = FALSE
    )
  )
  extension_plan <- plan_lisa_extension(result$output_dir, selection)
  extension_plan$expected_files
}
}
