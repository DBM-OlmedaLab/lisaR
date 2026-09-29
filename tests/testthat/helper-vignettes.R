vignette_source_path <- function(name) {
  filename <- paste0(name, ".Rmd")
  candidates <- c(
    system.file("doc", filename, package = "lisaR"),
    testthat::test_path("..", "..", "vignettes", filename)
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    stop("Vignette source is unavailable: ", filename, call. = FALSE)
  }
  candidates[[1L]]
}

read_vignette_source <- function(name) {
  paste(readLines(vignette_source_path(name), warn = FALSE), collapse = "\n")
}
