% Maintained from the corresponding R source documentation.
\name{lisa_figure_request}
\alias{lisa_figure_request}
\title{Name one exact figure to explore}
\usage{
lisa_figure_request(unit_type, collection, product,
                                category_id = NULL,
                                analysis_id = NULL, contrast_id = NULL,
                                entity = NULL, variant = NULL)
}
\arguments{
\item{unit_type}{\code{"single_de"} or \code{"contrast"}.}
\item{collection}{The collection identifier the figure belongs to.}
\item{product}{The figure product, for example \code{"volcano"}.}
\item{category_id}{The category the figure is drawn for. Required for a category product and refused for a native collection-wide product such as \code{"contrast_paired_heatmap"}, which is one figure for the whole collection.}
\item{analysis_id}{Analysis identifier; required when \code{unit_type} is \code{"single_de"}.}
\item{contrast_id}{Contrast identifier; required when \code{unit_type} is \code{"contrast"}.}
\item{entity}{Optional concrete entity, such as a KEGG map identifier.}
\item{variant}{Optional explicit variant, such as a heatmap scale.}
}
\value{An object of class \code{lisa_figure_request}.}
\description{Builds the request object used by every exploration entry point. A request
always denotes a single figure: passing a vector to any field is an error, so
a selection can never silently widen into a Cartesian product.}
\examples{
request <- lisa_figure_request(
unit_type = "single_de",
analysis_id = "demo_analysis",
collection = "H",
category_id = "SYN_SIGNAL",
product = "volcano"
)
lisa_request_id(request)
}
