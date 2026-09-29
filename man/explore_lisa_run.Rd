% Generated manually for the isolated H1 UI handoff; parent integration must add
% the export declaration before package build.
\name{explore_lisa_run}
\alias{explore_lisa_run}
\title{Explore a completed lisaR report interactively}
\usage{
explore_lisa_run(run_dir, extension_root = NULL, launch.browser = TRUE)
}
\arguments{
  \item{run_dir}{A completed STANDARD run directory, opened read-only.}
  \item{extension_root}{Optional exploration workspace outside \code{run_dir}.}
  \item{launch.browser}{Whether to ask Shiny to launch a browser.}
}
\value{Invisibly, the result from the local Shiny application.}
\description{
Serves the existing report navigation and adds explicit per-category controls for
already catalogued figure products. No figure is generated while browsing. A
generation starts only after the corresponding product-specific action;
the action is delegated to \code{lisa_explore_submit}. The save action delegates
to \code{lisa_explore_export}, which copies the complete STANDARD plus available
extensions and never completes missing products.
}
\details{\code{shiny} is an optional dependency. This function reports a clear
error when it is unavailable and never installs it automatically.}
\examples{
\dontrun{
explore_lisa_run("/path/to/completed-run")
}
}
