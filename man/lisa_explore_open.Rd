% Maintained from the corresponding R source documentation.
\name{lisa_explore_open}
\alias{lisa_explore_open}
\title{Open an exploration workspace over a completed run}
\usage{
lisa_explore_open(run_dir, extension_root = NULL)
}
\arguments{
\item{run_dir}{A completed lisaR run directory, used read-only.}
\item{extension_root}{Optional workspace location. It must lie outside \code{run_dir}. When omitted, a persistent user-owned directory under \code{\link[tools]{R_user_dir}()} is used, so no generated or exported material is ever placed inside the results tree that holds an accepted run. Pass an explicit path to keep the workspace somewhere else.}
}
\value{An object of class \code{lisa_explore_workspace}.}
\description{Opens a saved STANDARD run read-only and prepares a separate writable
workspace outside it. The source run is never written to.}
\examples{
# Opening a workspace requires a completed run; the sample project is
# offline, so the prerequisite run is kept interactive.
if (interactive()) {
ws <- lisa_explore_open(run_dir)
ws$root
}
}
