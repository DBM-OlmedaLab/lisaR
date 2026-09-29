\name{refresh_lisa_support_grades}
\alias{refresh_lisa_support_grades}
\title{Refresh minimum-support presentation in a new copy of an inferred LISA run}
\usage{
refresh_lisa_support_grades(source_run, output_dir,
  plot_formats = c("png", "svg", "pdf"))
}
\arguments{
\item{source_run}{Complete run with saved category inference and manifest.}
\item{output_dir}{New destination outside the source run; must not exist.}
\item{plot_formats}{Any of \code{png}, \code{svg}, \code{pdf}.}
}
\value{Invisibly, the new run directory.}
\description{Copies a manifest-verified run and decorates its already calculated
category inference. No DE, fgsea, closed tests, P values, lower bounds or families
are recalculated. Historical lollipop sources and figures remain unchanged.
Adds the distinct Hommel support view and refreshes exact-key navigator and
category-evidence metadata, with reusable offline interpretation help. Repeated
presentation refreshes preserve earlier support variants and provenance archives.
Current views show the highest adjusted-category-P star threshold reached:
\code{*} at P <= 0.05, \code{**} at P <= 0.01, \code{***} at P <= 0.001,
and \code{****} at P <= 0.0001. The exact P and simultaneous 95 percent
minimum supported count d/N remain available in the category card.
Rebuild the new copy with the installed \code{scripts/build_LISA_report.R} to
create the STANDARD or FULL report and seal its new manifest.}
