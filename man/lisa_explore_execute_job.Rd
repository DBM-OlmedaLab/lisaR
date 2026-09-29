% Maintained from the corresponding R source documentation.
\name{lisa_explore_execute_job}
\alias{lisa_explore_execute_job}
\title{Execute one queued exploration job}
\usage{
lisa_explore_execute_job(root, job_id)
}
\arguments{
\item{root}{The exploration workspace root.}
\item{job_id}{The queued job identifier.}
}
\value{Invisibly, the final job record.}
\description{Renders exactly one figure, validates the produced files and promotes the
result. This is the worker entry point; it is called in a separate process by
\code{\link{lisa_explore_submit}()} and is not normally invoked directly.}
\examples{
if (interactive()) {
lisa_explore_execute_job(ws$root, job_id)
}
}
