% Generated manually for the isolated H2 handoff.
\name{lisa_explore_configure_kegg}
\alias{lisa_explore_configure_kegg}
\alias{lisa_explore_prepare_kegg_index}
\title{Configure and prepare cached native KEGG pathway maps}
\usage{
lisa_explore_configure_kegg(
  ws,
  cache_root,
  snapshot_id,
  species,
  max_abs_log2fc = 0.5,
  color_power = 1
)

lisa_explore_prepare_kegg_index(ws, analyses = NULL, collections = NULL)
}
\arguments{
  \item{ws}{An exploration workspace returned by \code{lisa_explore_open()}.}
  \item{cache_root}{Root of an existing local KEGG snapshot cache.}
  \item{snapshot_id}{Identifier of the immutable snapshot in \code{cache_root}.}
  \item{species}{\code{"Homo sapiens"} or \code{"Mus musculus"}.}
  \item{max_abs_log2fc}{Positive finite node-fill cap passed to the native painter.}
  \item{color_power}{Positive finite colour power passed to the native painter.}
  \item{analyses}{Optional analysis identifiers to prepare.}
  \item{collections}{Optional collection identifiers to prepare.}
}
\value{
\code{lisa_explore_configure_kegg()} invisibly returns the validated recorded
declaration. \code{lisa_explore_prepare_kegg_index()} invisibly returns the
prepared pathway/category index.
}
\description{
Declares an already materialised, immutable KEGG snapshot and explicitly
prepares the pathway index used by the exploration catalogue. Both operations
are cache-only. They never download a resource, install a package, run DE, GSEA
or ORA, or paint a pathway map.
}
\details{
The prepared index is bound to the source-run manifest and the declared
snapshot. Paintable rows also record a digest of the validated pathway-list,
pathway-link, KGML and base-image bytes. If that context changes, the index is
treated as unavailable until it is prepared again.
}
\examples{
\dontrun{
ws <- lisa_explore_open("/path/to/completed-run")
lisa_explore_configure_kegg(
  ws, "/path/to/cache", "snapshot-2026", "Homo sapiens"
)
lisa_explore_prepare_kegg_index(ws)
}
}
