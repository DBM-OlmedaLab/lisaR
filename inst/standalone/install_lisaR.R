#!/usr/bin/env Rscript
# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.


# Install lisaR and its dependencies into one R library.
#
# Usage:
#   Rscript install_lisaR.R lisaR_0.6.0.tar.gz --profile=full
#   Rscript install_lisaR.R lisaR_0.6.0.tar.gz --profile=full --library=/path/to/R-library
#   Rscript install_lisaR.R lisaR_0.6.0.tar.gz --profile=validation
#
# Profiles:
#   core        Required packages for standard lisaR analyses.
#   full        Core plus every optional analysis and reporting dependency.
#   validation  Full plus packages used only to test and check the package.

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || any(args %in% c("-h", "--help"))) {
  cat(
    "Usage: Rscript install_lisaR.R <lisaR tarball> ",
    "[--profile=core|full|validation] [--library=<directory>]\n",
    sep = ""
  )
  quit(status = if (length(args)) 0L else 2L)
}

tarball <- args[[1L]]
profile_arg <- grep("^--profile=", args, value = TRUE)
library_arg <- grep("^--library=", args, value = TRUE)
profile <- if (length(profile_arg)) sub("^--profile=", "", profile_arg[[1L]]) else "full"
library_path <- if (length(library_arg)) {
  sub("^--library=", "", library_arg[[1L]])
} else {
  Sys.getenv("R_LIBS_USER")
}

if (!profile %in% c("core", "full", "validation")) {
  stop(
    "`--profile` must be `core`, `full`, or `validation`.",
    call. = FALSE
  )
}
if (!file.exists(tarball)) {
  stop("lisaR tarball does not exist: ", tarball, call. = FALSE)
}
tarball <- normalizePath(tarball, mustWork = TRUE)
if (!nzchar(library_path)) {
  library_path <- file.path(path.expand("~"), "R", "lisaR-library")
}
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
library_path <- normalizePath(library_path, mustWork = TRUE)
# Use the requested library plus R's site and base libraries. In particular,
# do not inherit the caller's personal development library: selected profile
# packages must be installed into `library_path` and are verified there below.
.libPaths(library_path)

script_args <- commandArgs(FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)][1L])
script_dir <- dirname(normalizePath(script_file, mustWork = TRUE))
profile_path <- file.path(script_dir, "dependency-profiles.tsv")
if (!file.exists(profile_path)) {
  stop("Dependency profile is missing beside the installer: ", profile_path,
       call. = FALSE)
}
requirements <- utils::read.delim(
  profile_path,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_rows <- requirements$profile == "core"
if (profile %in% c("full", "validation")) {
  required_rows <- required_rows | requirements$profile == "full"
}
if (identical(profile, "validation")) {
  required_rows <- required_rows | requirements$profile == "validation"
}
requirements <- requirements[required_rows, , drop = FALSE]

cran_repo <- getOption("repos")[["CRAN"]]
if (is.null(cran_repo) || is.na(cran_repo) ||
    !nzchar(cran_repo) || identical(cran_repo, "@CRAN@")) {
  cran_repo <- "https://cloud.r-project.org"
}

installed_in_target <- function(package) {
  nzchar(system.file(package = package, lib.loc = library_path))
}

missing_from_target <- function(packages) {
  packages[!vapply(packages, installed_in_target, logical(1))]
}

# BiocManager provides one repository view for CRAN and Bioconductor. It is a
# bootstrap tool, not a lisaR runtime dependency, but it must itself be present
# in the requested library so the installation does not depend on a developer
# library elsewhere on the machine.
if (!installed_in_target("BiocManager")) {
  install.packages(
    "BiocManager",
    lib = library_path,
    repos = cran_repo,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}
if (!installed_in_target("BiocManager")) {
  stop("BiocManager was not available after bootstrap.", call. = FALSE)
}

old_repos <- getOption("repos")
on.exit(options(repos = old_repos), add = TRUE)
options(repos = c(CRAN = cran_repo))
repositories <- BiocManager::repositories()
available <- utils::available.packages(repos = repositories)

unknown_roots <- setdiff(requirements$package, rownames(available))
if (length(unknown_roots)) {
  stop(
    "Dependency metadata was unavailable for: ",
    paste(unknown_roots, collapse = ", "),
    call. = FALSE
  )
}

# install.packages() can regard packages in another user library as satisfying
# a dependency. For a portable standalone library that is not sufficient.
# Compute the complete recursive Depends/Imports/LinkingTo closure first, then
# explicitly install every non-base member that is absent from the target.
dependency_map <- tools::package_dependencies(
  requirements$package,
  db = available,
  which = c("Depends", "Imports", "LinkingTo"),
  recursive = TRUE
)
closure_packages <- sort(unique(c(
  requirements$package,
  unlist(dependency_map, use.names = FALSE)
)))

system_packages <- utils::installed.packages(
  lib.loc = unique(c(.Library, .Library.site))
)
system_priorities <- system_packages[, "Priority"]
base_or_recommended <- rownames(system_packages)[
  !is.na(system_priorities) &
    system_priorities %in% c("base", "recommended")
]
target_closure <- setdiff(closure_packages, base_or_recommended)

unknown_closure <- setdiff(target_closure, rownames(available))
if (length(unknown_closure)) {
  stop(
    "Repository metadata was unavailable for transitive dependencies: ",
    paste(unknown_closure, collapse = ", "),
    call. = FALSE
  )
}

closure_missing <- missing_from_target(target_closure)
if (length(closure_missing)) {
  BiocManager::install(
    closure_missing,
    lib = library_path,
    ask = FALSE,
    update = FALSE,
    force = TRUE,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

missing_closure <- missing_from_target(target_closure)
if (length(missing_closure)) {
  stop(
    "The target library is missing packages from the dependency closure: ",
    paste(missing_closure, collapse = ", "),
    call. = FALSE
  )
}

root_loadable <- vapply(requirements$package, function(package) {
  requireNamespace(package, lib.loc = library_path, quietly = TRUE)
}, logical(1))
if (any(!root_loadable)) {
  stop(
    "Installed profile packages could not be loaded: ",
    paste(requirements$package[!root_loadable], collapse = ", "),
    call. = FALSE
  )
}

install.packages(
  tarball,
  lib = library_path,
  repos = NULL,
  type = "source",
  dependencies = FALSE
)

if (!installed_in_target("lisaR")) {
  stop("lisaR was not available after installation.", call. = FALSE)
}

installed_packages <- unique(c("lisaR", "BiocManager", target_closure))
versions <- vapply(installed_packages, function(package) {
  description <- file.path(
    system.file(package = package, lib.loc = library_path),
    "DESCRIPTION"
  )
  unname(read.dcf(description, fields = "Version")[[1L]])
}, character(1))
result <- data.frame(
  package = names(versions),
  version = unname(versions),
  library = library_path,
  stringsAsFactors = FALSE
)
utils::write.table(
  result,
  file = file.path(library_path, "lisaR-installed-packages.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat("LISAR_INSTALL_GATE=PASS\n")
cat("profile=", profile, "\n", sep = "")
cat("library=", library_path, "\n", sep = "")
cat("dependency_closure=", length(target_closure), "\n", sep = "")
