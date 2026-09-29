% Maintained from the corresponding R source documentation.
\name{lisa_explore_pump}
\alias{lisa_explore_pump}
\title{Start the next queued figure if the single worker slot is free}
\usage{
lisa_explore_pump(ws, lock_timeout_seconds = 30)
}
\arguments{
\item{ws}{A workspace from \code{\link{lisa_explore_open}()}.}
\item{lock_timeout_seconds}{Maximum time to wait for the short index-lock claim. Interactive polling can set zero so a busy or unverifiable lock is reported without blocking the Shiny event loop.}
}
\value{Invisibly, the dispatched job identifier, or \code{NULL}.}
\description{Exploration runs one costly worker at a time. A submission made while another
figure is rendering waits in the queue; this function dispatches the oldest
waiting job once the slot frees. It only ever starts work that was already
explicitly submitted, so it never turns browsing into rendering.}
\examples{
if (interactive()) {
lisa_explore_pump(ws)
}
}
