\name{verify_lisa_run}
\alias{verify_lisa_run}
\title{Check a completed LISA run}
\usage{verify_lisa_run(run_dir)}
\arguments{\item{run_dir}{Path to a completed LISA results folder. A planning directory is accepted for inspection but deliberately fails the scientific verification gate.}}
\value{A list containing \code{gate}, \code{findings} and the path to the file list used for verification. Incomplete, changed or planning artifacts return \code{gate = "FAIL"} with explanatory findings. An unsafe or unreadable path can raise an error before verification.}
\description{Checks that every expected file is present and unchanged and that lisaR completed its final scientific validation. A \code{PLAN_PASS} artifact returns \code{gate = "FAIL"} with finding \code{planning_artifact_not_scientific_run}.}
\examples{
project <- lisa_init_project(tempfile("lisa-verify-project-"))
config <- file.path(project, "study.yml")
plan <- yaml::read_yaml(config)
plan$pipeline$output_dir <- "results/example"
plan$pipeline$dry_run <- TRUE
plan$pipeline$sufficiency_action <- "ignore_with_reason"
plan$pipeline$sufficiency_ignore_reason <- "Synthetic sample fixture."
yaml::write_yaml(plan, config)
planned <- run_lisa(config)
check <- verify_lisa_run(planned$plan_dir)
check$gate
check$findings
}
