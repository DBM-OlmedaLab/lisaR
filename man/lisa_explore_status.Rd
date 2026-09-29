% Maintained from the corresponding R source documentation.
\name{lisa_explore_status}
\alias{lisa_explore_status}
\title{Report the state of every requestable figure}
\usage{
lisa_explore_status(ws)
}
\arguments{
\item{ws}{A workspace from \code{\link{lisa_explore_open}()}.}
}
\value{The catalog with a \code{state} column added.}
\description{Joins the read-only catalog with the index so an interface can render the six
product states without scheduling anything.}
\examples{
if (interactive()) {
table(lisa_explore_status(ws)$state)
}
}
