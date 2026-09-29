% Generated manually for the isolated H3 handoff.
\name{lisa_explore_prepare_contrast_products}
\alias{lisa_explore_prepare_contrast_products}
\title{Prepare the shared native contrast products for an exploration workspace}
\usage{
lisa_explore_prepare_contrast_products(ws, contrasts = NULL, collections = NULL)
}
\arguments{
  \item{ws}{An exploration workspace returned by \code{lisa_explore_open()}.}
  \item{contrasts}{Optional contrast identifiers to restrict preparation to.}
  \item{collections}{Optional collection identifiers to restrict preparation to.}
}
\value{
Invisibly, a data frame with one row per prepared contrast bundle, giving its
contrast, collection, prepared directory, preparation digest and the number of
ranked categories and ranked pathways it carries.
}
\description{
Builds, once per contrast and collection, the paired gene evidence and the
complete ranked category summary that the native contrast gene card, paired
gene heatmap, gene-category network and painted pathway map all read. When a
validated local KEGG snapshot is declared for the workspace it also derives the
complete ranked contrast pathway index and its category associations from that
same paired evidence.
}
\details{
This is data preparation, not analysis. It reads only manifest-bound saved
evidence of the source run and the declared local snapshot; it runs no DE, GSEA,
ORA or LISA step, writes nothing inside the source run, fetches no resource and
draws no figure. Every generator it calls is invoked in a mode that refuses to
plot, so the ranking a reader is later offered is by construction the ranking a
plot or paint would use.

Each bundle is written to
\code{<workspace>/prepared/contrast/<context>/<contrast>/<collection>/} with a
receipt recording the source-manifest hash, the exact hash of every bound
scientific input, the preparation code hash, the effective fixed settings and,
where applicable, the KEGG resource digest. One bundle is reused by the card,
heatmap, network and map requests. If any bound input, code, resource or context
changes, the bundle fails closed and preparation must be repeated; prepared rows
from different contrast endpoints, collections, species or snapshots are never
mixed.
}
\examples{
\dontrun{
ws <- lisa_explore_open("/path/to/completed-run")
lisa_explore_prepare_contrast_products(ws)
}
}
