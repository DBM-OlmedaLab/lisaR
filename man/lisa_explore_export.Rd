% Maintained from the corresponding R source documentation.
\name{lisa_explore_export}
\alias{lisa_explore_export}
\title{Export a static report bundle}
\usage{
lisa_explore_export(ws, destination, format = c("dir", "zip"),
                                allow_absolute_references = FALSE,
                                snapshot_artifacts = NULL,
                                progress = NULL,
                                should_cancel = NULL,
                                staging_id = NULL)
}
\arguments{
\item{ws}{A workspace from \code{\link{lisa_explore_open}()}.}
\item{destination}{Directory to create for the bundle. It must not exist.}
\item{format}{\code{"dir"} for a folder, \code{"zip"} to also produce a ZIP archive.}
\item{allow_absolute_references}{Promote the bundle even when it contains absolute HTML references or symlinks. The default \code{FALSE} refuses instead, because such a bundle does not open from another location.}
\item{snapshot_artifacts}{Optional validated artifact snapshot used by the background exporter. Leave as \code{NULL} to snapshot currently valid artifacts.}
\item{progress}{Optional callback accepting a stage and message; used to report export progress without changing the selected artifact snapshot.}
\item{should_cancel}{Optional zero-argument callback returning \code{TRUE} to request cooperative cancellation at an export checkpoint.}
\item{staging_id}{Internal per-job staging identity used by the background exporter. Interactive callers should leave this as \code{NULL}.}
}
\value{Invisibly, a list describing the exported bundle.}
\description{Copies the complete STANDARD report plus every exploration artifact that is
validated when the export begins. No figure is generated and no missing
product is completed. The result is self-contained: it opens from any
location with R and Shiny closed.}
\examples{
if (interactive()) {
lisa_explore_export(ws, file.path(ws$root, "exports", "bundle-001"))
}
}
