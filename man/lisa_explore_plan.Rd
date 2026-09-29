% Maintained from the corresponding R source documentation.
\name{lisa_explore_plan}
\alias{lisa_explore_plan}
\title{Plan one exact figure without rendering it}
\usage{
lisa_explore_plan(ws, request)
}
\arguments{
\item{ws}{A workspace from \code{\link{lisa_explore_open}()}.}
\item{request}{A \code{\link{lisa_figure_request}()}.}
}
\value{A list describing the request, its state and its inputs.}
\description{Determines whether a figure is applicable, which saved inputs it needs and
what its current state is. Nothing is rendered and nothing is scheduled.}
\examples{
if (interactive()) {
lisa_explore_plan(ws, request)
}
}
