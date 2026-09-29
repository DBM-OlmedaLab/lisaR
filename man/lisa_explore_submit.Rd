% Maintained from the corresponding R source documentation.
\name{lisa_explore_submit}
\alias{lisa_explore_submit}
\title{Submit one exact figure for generation}
\usage{
lisa_explore_submit(ws, request, background = TRUE)
}
\arguments{
\item{ws}{A workspace from \code{\link{lisa_explore_open}()}.}
\item{request}{A \code{\link{lisa_figure_request}()}.}
\item{background}{Run in a separate R process. Requires the optional \code{callr} package; set \code{FALSE} to run in the calling process.}
}
\value{A list with the job identifier, the reuse key and the resulting state.}
\description{Schedules exactly one figure. Repeating a submission that is already
validated, queued or running returns the existing result or job instead of
starting a second one.}
\examples{
if (interactive()) {
lisa_explore_submit(ws, request)
}
}
