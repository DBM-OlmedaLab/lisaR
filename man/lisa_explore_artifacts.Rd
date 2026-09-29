% Maintained from the corresponding R source documentation.
\name{lisa_explore_artifacts}
\alias{lisa_explore_artifacts}
\title{List validated exploration artifacts}
\usage{
lisa_explore_artifacts(ws, request = NULL)
}
\arguments{
\item{ws}{A workspace from \code{\link{lisa_explore_open}()}.}
\item{request}{Optional \code{\link{lisa_figure_request}()} to report a single figure.}
}
\value{A data frame of validated artifacts.}
\description{Reads the index and reports only artifacts that still validate. Nothing is
rendered.}
\examples{
if (interactive()) {
lisa_explore_artifacts(ws)
}
}
