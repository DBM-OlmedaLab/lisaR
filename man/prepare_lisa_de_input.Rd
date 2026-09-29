\name{prepare_lisa_de_input}
\alias{prepare_lisa_de_input}
\title{Prepare an auditable differential-expression table for LISA}
\usage{
prepare_lisa_de_input(
  data,
  method = c("DESeq2", "edgeR", "limma"),
  columns = list(),
  positive_direction,
  rank_direction = c("same_as_effect", "opposite_to_effect"),
  duplicate_policy = "error",
  na_policy = c("retain", "error"),
  matrix_scale = NULL
)
}
\arguments{
\item{data}{Data frame containing the complete tested DE results.}
\item{method}{Caller-declared table convention, not verification of producing software.}
\item{columns}{Named list mapping \code{symbol}, \code{logfc}, \code{rank},
\code{pvalue} and \code{padj} to source columns. Use \code{symbol = ".rownames"}
explicitly for non-default row names. Standard effect and P-value columns,
DESeq2 signed Wald \code{stat}, and limma \code{t} are defaults. Unsigned
DESeq2 LRT results require a separately computed signed rank. edgeR requires an explicit
signed rank column; unsigned multi-df F or LR statistics are never signed automatically.}
\item{positive_direction}{One non-empty description of what positive effects mean,
e.g. \code{"treated relative to control"}. No contrast is inferred.}
\item{rank_direction}{Whether the input signed rank is oriented like or opposite
to the effect. Opposite ranks are negated once. Comparable nonzero signs must agree.}
\item{duplicate_policy}{A LISA \code{"error"}, deterministic \code{select}, or
row-level \code{mapping_file} policy. DE statistics cannot be aggregated.
Identifier collisions are checked case-insensitively.}
\item{na_policy}{\code{"retain"} preserves missing statistics and records ranking
eligibility; \code{"error"} rejects them. Missing/empty identifiers, malformed
numbers and infinite statistics always raise errors.}
\item{matrix_scale}{Optional \code{"transformed"} or \code{"count_like"} declaration,
recorded separately. This function does not inspect or transform expression matrices.}
}
\value{
List with canonical \code{data} and unfiltered \code{full_data} columns
\code{symbol}, \code{log2FoldChange}, \code{rank_value}, \code{pvalue}, \code{padj};
input-declaration \code{receipt}; \code{column_mapping}; \code{row_audit}; and
\code{duplicate_resolution}. \code{data} differs from \code{full_data} only when
an explicitly requested duplicate policy selects or remaps rows. Both effect and
rank describe \code{positive_direction}.
}
\description{
Adapts already computed data frames using base R without loading DESeq2, edgeR
or limma. Preserves all tested rows, including non-significant and permitted
missing statistics, in \code{full_data}. Applies no adjusted-P filter, fits no
models, and infers no sample metadata. Default duplicate errors prevent implicit
row loss. Original noncanonical source columns are not copied into the output;
keep the source table as the original input artifact.
}
\examples{
results <- data.frame(gene = c("A", "B", "C"),
  log2FoldChange = c(1, -0.5, 0.1), stat = c(3, -1, 0.2),
  pvalue = c(0.001, 0.3, 0.8), padj = c(0.01, 0.6, NA))
prepared <- prepare_lisa_de_input(results, "DESeq2",
  columns = list(symbol = "gene"),
  positive_direction = "treated relative to control")
prepared$data
prepared$receipt
}
