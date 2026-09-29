\name{run_lisa}
\alias{run_lisa}
\title{Run a LISA study}
\usage{run_lisa(config = NULL, ..., de = NULL, output_dir = NULL,
  species = NULL, positive_direction = NULL, method = NULL, columns = NULL,
  contrast = NULL, contrast_label = NULL, contrast_positive_direction = NULL,
  collections = NULL, mode = NULL, png = NULL, svg = NULL, pdf = NULL,
  tables = NULL, recipes = NULL, title = NULL, workers = NULL,
  dry_run = NULL, resources = NULL, pipeline = NULL, report = NULL)}
\arguments{
\item{config}{A YAML/JSON configuration path or equivalent named R list. Leave it \code{NULL} to use the simple arguments below.}
\item{...}{Must be empty. It is present so that every simple argument is supplied by name and a second positional argument can never be taken for one of them. An unknown name is an error, never an ignored option.}
\item{de}{One already computed differential-expression result, or a named list of them whose names become the analysis identifiers. Each element may be a data frame, a path to a \code{.tsv}, \code{.txt}, \code{.tab} or \code{.csv} table, a DESeq2 \code{DESeqResults} object, or an edgeR \code{TopTags}, \code{DGEExact} or \code{DGELRT} result. Serialized \code{.rds} inputs are refused here because deserializing one is a local trust decision; make it explicitly with \code{\link{run_lisa_de}}. There is no limit of two analyses.}
\item{output_dir}{A new directory that will hold the generated \code{study.json}, the adapted inputs and receipts under \code{inputs/}, and the pipeline results under \code{results/}. An existing directory is refused, never overwritten.}
\item{species}{Species represented by the input identifiers, for example \code{"Homo sapiens"}. It is never inferred.}
\item{positive_direction}{One non-empty description of what a positive effect means, or a named list with one entry per analysis. lisaR never infers which side of a comparison is positive.}
\item{method}{Declared table convention: \code{"DESeq2"}, \code{"edgeR"} or \code{"limma"}, or a named list per analysis. This is a caller declaration, not verification of the producing software.}
\item{columns}{Named mapping of \code{symbol}, \code{logfc}, \code{rank}, \code{pvalue} and \code{padj} to source columns, or a named list of such mappings per analysis. Passed unchanged to \code{\link{prepare_lisa_de_input}}. An edgeR result requires an explicit signed \code{rank}: \code{F} and \code{LR} are unsigned and are never signed for you.}
\item{contrast}{Two analysis identifiers, minuend first: \code{c("ON", "PRE")} means \emph{ON minus PRE}. The roles may also be named, as \code{c(minuend = "ON", subtrahend = "PRE")}. A list of such pairs declares several contrasts. The orientation is never taken from the order of \code{de}.}
\item{contrast_label, contrast_positive_direction}{Optional wording for the contrast. The default restates the declared orientation and claims no statistical significance for any category.}
\item{collections}{Gene-set collections to analyse. Defaults to the registered analysis collections.}
\item{mode}{\code{"standard"} (default) or \code{"full"}. \code{"full"} writes one report in the same \code{output_dir}: the standard report plus every applicable FULL figure, rendered from the run's own saved results inside the same run (under \verb{artifacts/}) without recomputing enrichment, and sealed by the same run manifest. If a FULL product fails, the run fails. It can create thousands of files.}
\item{png, svg, pdf}{Figure formats. In this form PNG and SVG default to \code{TRUE} and PDF to \code{FALSE}. A configuration file keeps its own documented defaults, which are unchanged.}
\item{tables}{Whether to export the per-figure source-data tables. \code{TRUE} by default here.}
\item{recipes}{Whether to export the figure reproduction recipes. \code{TRUE} by default here.}
\item{title}{Optional report title.}
\item{workers}{Worker processes for the pipeline.}
\item{dry_run}{Whether to plan the run instead of computing it. A dry run writes the output plan, performs no enrichment and promises no report.}
\item{resources}{Optional resource selection. \code{"core"} or \code{"expanded"} chooses the dictionary tier; a named list with any of \code{dictionary}, \code{term2gene} and \code{category_map} names registered resources explicitly. Omitted, the registered defaults for the declared species are resolved internally and you name nothing.}
\item{pipeline, report}{Named lists of advanced overrides merged over the defaults of this form before validation. Explicit values win, including \code{report = list(formats = list(svg = FALSE))}. Supplying the same setting both directly and here is refused, even when the values agree.}
}
\value{Invisibly, a lifecycle result. A scientific execution returns its output directory and \code{gate = "PASS"}. A configuration with \code{pipeline.dry_run = TRUE} returns \code{gate = "PLAN_PASS"}, the still-free intended \code{output_dir}, and the separate \code{plan_dir} and \code{plan_path}. The simple form adds \code{study_dir}, \code{config_path}, the adapted \code{inputs}, the validated \code{analyses} and \code{contrasts}, \code{effective_config}, and \code{reports}: a table of the report pages with the path of each and whether it was actually written. Report paths are resolved by existence, never constructed.}
\description{Validates and runs a LISA study, including single analyses, declared contrasts and requested report products. Planning is stored separately and is never promoted as a completed scientific run.

This is the single entry point, and it can be called in two ways. Supply \code{config} to run an existing YAML or JSON study configuration exactly as before. Supply the simple arguments instead to run a study directly from one or more differential-expression results, without writing a configuration file and without naming any dictionary, category map or gene-set membership resource: those are resolved internally for the declared species. The simple form assembles an ordinary study, writes it where the caller can read it, and hands it to the same engine. No threshold, method, ranking orientation or scientific value differs between the two forms.

The two forms are mutually exclusive. Mixing \code{config} with any simple argument is an error rather than a silent preference for one of them.}
\details{
\code{report.category_nes_variants} is a nonempty array selecting
\code{clean}, \code{percentages}, \code{direction}, and/or \code{dispersion}. The default renders
all four. The clean view retains the mean NES, direction adds median and
directional percentages, and dispersion adds P25--P75 between-gene-set
dispersion. All views reuse the same results; no additional GSEA tests are run.
Study files must declare \code{pipeline.schema_version: "1.0.0"}.
Earlier schema versions require explicit migration in a copy of the file;
they are not silently upgraded. At least one variant must be selected,
independently of the requested output formats.

In the simple form every checkable argument is validated and every input is
adapted in memory before anything is written, and the whole study is
materialized in a staging directory beside \code{output_dir} and validated
there. A wrong flag or a missing resource therefore leaves \code{output_dir}
untouched, and correcting the argument and calling again works without
deleting anything. If the engine itself fails after the study has been
written, nothing needs rebuilding: rerun the saved configuration with
\code{run_lisa(config = file.path(output_dir, "study.json"))}.
}
\examples{
project <- lisa_init_project(tempfile("lisa-plan-project-"))
config <- file.path(project, "study.yml")
plan <- yaml::read_yaml(config)
plan$pipeline$output_dir <- "results/example"
plan$pipeline$dry_run <- TRUE
plan$pipeline$sufficiency_action <- "ignore_with_reason"
plan$pipeline$sufficiency_ignore_reason <- "Synthetic sample fixture."
yaml::write_yaml(plan, config)
result <- run_lisa(config)
result$gate

# The simple form. Running it needs the installed semantic resources, so the
# call itself stays interactive.
results <- data.frame(
  symbol = c("GENE01", "GENE02", "GENE03"),
  log2FoldChange = c(1.2, -0.8, 0.2),
  stat = c(4.1, -2.6, 0.5),
  pvalue = c(1e-4, 0.01, 0.6),
  padj = c(1e-3, 0.04, 0.8)
)
if (interactive()) {
  run <- run_lisa(
    de = results,
    output_dir = tempfile("lisa-study-"),
    species = "Homo sapiens",
    positive_direction = "higher in treated than in control"
  )
  run$reports
}
}
