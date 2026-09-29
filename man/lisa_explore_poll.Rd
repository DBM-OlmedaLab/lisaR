% Maintained from the corresponding R source documentation.
\name{lisa_explore_poll}
\alias{lisa_explore_poll}
\title{Report the state of exploration jobs}
\usage{
lisa_explore_poll(ws, job_id = NULL, dispatch = TRUE)
}
\arguments{
\item{ws}{A workspace from \code{\link{lisa_explore_open}()}.}
\item{job_id}{Optional job identifier; all jobs are reported when omitted.}
\item{dispatch}{Hand the free worker slot to the oldest already-submitted figure that is still waiting. Set \code{FALSE} to report state without dispatching anything.}
}
\value{A data frame with one row per job.}
\description{Reconciles each job with its owning process before reporting, so an orphaned
job is reported as failed rather than as running or available.}
\examples{
if (interactive()) {
lisa_explore_poll(ws)
}
}
