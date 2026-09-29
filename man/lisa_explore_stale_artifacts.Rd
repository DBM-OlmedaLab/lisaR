% Maintained from the corresponding R source documentation.
\name{lisa_explore_stale_artifacts}
\alias{lisa_explore_stale_artifacts}
\title{Report index entries that are no longer current}
\usage{
lisa_explore_stale_artifacts(ws)
}
\arguments{
\item{ws}{A workspace from \code{\link{lisa_explore_open}()}.}
}
\value{A data frame of stale entries with the reason each is stale.}
\description{Lists artifacts whose own files still validate but whose reuse key no longer
matches the key the same request would be minted with today, because the
source run, the presentation policy or the engine has changed. These entries
are deliberately hidden from \code{\link{lisa_explore_artifacts}()} and from every export;
this function exists so their exclusion is visible rather than silent.}
\examples{
if (interactive()) {
lisa_explore_stale_artifacts(ws)
}
}
