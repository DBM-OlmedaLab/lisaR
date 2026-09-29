% Maintained from the corresponding R source documentation.
\name{lisa_request_id}
\alias{lisa_request_id}
\title{Stable identifier of an exact figure request}
\usage{
lisa_request_id(request)
}
\arguments{
\item{request}{A \code{\link{lisa_figure_request}()}.}
}
\value{A short stable identifier for the requested figure.}
\description{
Derive a deterministic identifier from the canonical fields of a figure
request. The identifier labels the request; it does not render a figure
or certify that an artifact has been generated.
}
\examples{
lisa_request_id(lisa_figure_request(
unit_type = "single_de", analysis_id = "demo", collection = "H",
category_id = "SYN_SIGNAL", product = "volcano"
))
}
