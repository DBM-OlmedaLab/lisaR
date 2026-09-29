# Render the source vignettes without rebuilding the package or running
# eval=FALSE analysis examples. Usage: Rscript build_guides.R REPO OUTPUT
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2L)
root <- normalizePath(args[[1]], mustWork = TRUE)
dir.create(args[[2]], recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(args[[2]], mustWork = TRUE)
cat(R.version.string, "\n")
print(.libPaths())
stopifnot(requireNamespace("lisaR", quietly = TRUE),
          requireNamespace("rmarkdown", quietly = TRUE),
          rmarkdown::pandoc_available())
print(utils::packageVersion("lisaR"))
files <- list.files(file.path(root, "vignettes"), pattern = "[.]Rmd$",
                    full.names = TRUE)
for (path in files) {
  cat("Rendering ", basename(path), "\n", sep = "")
  rmarkdown::render(path, output_dir = out, quiet = TRUE,
                    envir = new.env(parent = globalenv()))
}
file.copy(file.path(root, "vignettes", "figures"), out, recursive = TRUE)
file.copy(files, out, overwrite = TRUE)
capture.output(sessionInfo(), file = file.path(out, "session-info.txt"))
cat("Rendered ", length(files), " guides.\n", sep = "")
