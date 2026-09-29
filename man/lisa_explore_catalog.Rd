% Maintained from the corresponding R source documentation.
\name{lisa_explore_catalog}
\alias{lisa_explore_catalog}
\title{Enumerate every figure that can be requested}
\usage{
lisa_explore_catalog(ws, refresh = FALSE)
}
\arguments{
\item{ws}{A workspace from \code{\link{lisa_explore_open}()}.}
\item{refresh}{Recompute the catalog instead of reusing the cached copy.}
}
\value{A data frame with one row per requestable figure.}
\description{Reads the saved evidence of the source run and returns one row per
requestable figure. No generator is invoked and nothing is written inside the
source run. A catalog cache may be written in the separate exploration workspace.}
\examples{
if (interactive()) {
ws <- lisa_explore_open(run_dir)
head(lisa_explore_catalog(ws))
}
}
