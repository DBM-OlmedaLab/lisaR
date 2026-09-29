\name{render_lisa_categories}
\alias{render_lisa_categories}
\title{Create additional selected or full report figures}
\usage{
render_lisa_categories(source_run, selection, output_dir, exact_products = NULL)
}
\arguments{
  \item{exact_products}{Optional named list of exact heatmap variants, KEGG maps or native products used by the exploration interface. Leave as \code{NULL} for the standard selected or full report.}
  \item{source_run}{A completed and checked lisaR run.}
  \item{selection}{The string \code{"full"}, a YAML/JSON selection path, or a named list.}
  \item{output_dir}{A new folder for the additional report.}
}
\value{Invisibly, a summary of the created and checked report.}
\description{Creates category figures from the saved scientific tables in a new folder. It checks the new report, does not repeat enrichment and leaves the original run unchanged.}
\examples{
# Rendering uses only a verified saved run and does not repeat enrichment.
# The offline demonstration run is kept interactive because it creates figures.
if (interactive()) {
  project <- lisa_init_project(tempfile("lisa-render-project-"))
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
  rendered <- render_lisa_categories(
    result$output_dir,
    selection,
    tempfile("lisa-extension-example-")
  )
  rendered$validation_status
}
}
