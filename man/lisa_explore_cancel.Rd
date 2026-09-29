% Maintained from the corresponding R source documentation.
\name{lisa_explore_cancel}
\alias{lisa_explore_cancel}
\title{Request cancellation of an exploration job}
\usage{
lisa_explore_cancel(ws, job_id)
}
\arguments{
\item{ws}{A workspace from \code{\link{lisa_explore_open}()}.}
\item{job_id}{The job identifier to cancel.}
}
\value{A list describing the cancellation outcome.}
\description{A queued job is cancelled at once. A running job is asked to stop and honours
the request at its next safe point; no process is ever forcibly terminated.}
\examples{
if (interactive()) {
lisa_explore_cancel(ws, job_id)
}
}
