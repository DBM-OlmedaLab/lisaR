# Helpers for the exploration engine tests.
#
# The heavy tests need a completed STANDARD run. Its location is supplied by the
# LISAR_EXPLORE_FIXTURE_RUN environment variable so the suite carries no machine
# specific path; without it those tests skip rather than invent a run.

# The route the shared fixture's one volcano attaches to. It is the
# category evidence sheet -- the page the collection navigator links to and the
# reader actually opens -- not the navigator itself. Defined once so every test
# asserts the same route.
explore_attach_route <- function() {
  "report_pages/evidence/response_a/GOBP-C2/index.html"
}

explore_fixture_run <- function() {
  path <- Sys.getenv("LISAR_EXPLORE_FIXTURE_RUN", "")
  if (!nzchar(path) || !dir.exists(path)) {
    testthat::skip("LISAR_EXPLORE_FIXTURE_RUN is not set to a completed run")
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

# Digest of a whole tree: relative path plus content hash for every file. Used to
# prove the immutable source run is untouched.
explore_tree_digest <- function(root) {
  files <- sort(list.files(root, recursive = TRUE, all.files = TRUE,
                           full.names = TRUE, no.. = TRUE))
  files <- files[!dir.exists(files)]
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  entries <- vapply(files, function(file) {
    paste(substring(normalizePath(file, winslash = "/"), nchar(root) + 2L),
          unname(tools::md5sum(file)))
  }, character(1))
  list(n = length(files), digest = digest::digest(unname(entries)))
}

# Instrument the native generators themselves, not the files they leave behind.
# `lisa_run_post_script` dispatches the gene-level builders (volcano, cards,
# heatmaps) and `write_lisa_gsea_category_pathway_plots` draws the gene-set and
# KEGG gene-set figures. Counting produced PNGs would not catch a generator that
# drew everything and deleted what was not requested; counting calls does.
explore_with_generator_trace <- function(code) {
  env <- new.env(parent = emptyenv())
  env$post <- list()
  env$pathway <- 0L
  assign(".LISA_EXPLORE_TRACE", env, envir = globalenv())
  on.exit({
    suppressMessages(try(untrace("lisa_run_post_script",
                                 where = asNamespace("lisaR")), silent = TRUE))
    suppressMessages(try(untrace("write_lisa_gsea_category_pathway_plots",
                                 where = asNamespace("lisaR")), silent = TRUE))
    suppressWarnings(try(rm(".LISA_EXPLORE_TRACE", envir = globalenv()),
                         silent = TRUE))
  }, add = TRUE)
  suppressMessages({
    trace("lisa_run_post_script", where = asNamespace("lisaR"), print = FALSE,
          tracer = quote({
            trace_env <- get(".LISA_EXPLORE_TRACE", envir = globalenv())
            trace_env$post <- c(trace_env$post,
                                list(list(script = script_name,
                                          args = paste(args, collapse = " "))))
          }))
    trace("write_lisa_gsea_category_pathway_plots", where = asNamespace("lisaR"),
          print = FALSE,
          tracer = quote({
            trace_env <- get(".LISA_EXPLORE_TRACE", envir = globalenv())
            trace_env$pathway <- trace_env$pathway + 1L
          }))
  })
  value <- force(code)
  list(value = value, post = env$post, pathway = env$pathway)
}

explore_generator_scripts <- function(trace) {
  vapply(trace$post, function(call) as.character(call$script), character(1))
}

explore_generator_args <- function(trace) {
  vapply(trace$post, function(call) as.character(call$args), character(1))
}

# One rendered volcano, shared by every test that needs an available artifact.
# Rendering is the expensive operation in this suite, so it happens once.
explore_shared <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    run <- explore_fixture_run()
    root <- file.path(tempdir(), "lisa-explore-shared-ws")
    if (dir.exists(root)) unlink(root, recursive = TRUE, force = TRUE)
    before <- explore_tree_digest(run)
    ws <- lisa_explore_open(run, root)
    request <- lisa_figure_request(
      unit_type = "single_de", analysis_id = "response_a",
      collection = "GOBP-C2", category_id = "SYN_SIGNAL", product = "volcano")
    started <- Sys.time()
    trace <- explore_with_generator_trace(
      lisa_explore_submit(ws, request, background = FALSE))
    cache <<- list(ws = ws, request = request, trace = trace,
                   before = before, run = run,
                   elapsed = as.numeric(difftime(Sys.time(), started,
                                                 units = "secs")))
    cache
  }
})

# One exported bundle, shared by every test that needs one.
#
# Each full export copies and hashes ~1,700 files through the guarded-copy
# layer, so exporting once per test made this suite cost roughly five of them --
# the expense the parent review flagged as finding 6. The bundle is built once
# here and reused; the assertions in each test stay completely independent and
# none of them is weakened, because they are read-only checks over a bundle that
# is identical however many times it is produced.
explore_shared_bundle <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    shared <- explore_shared()
    destination <- file.path(tempdir(), "lisa-explore-shared-bundle")
    unlink(destination, recursive = TRUE, force = TRUE)
    started <- Sys.time()
    trace <- explore_with_generator_trace(
      lisa_explore_export(shared$ws, destination))
    cache <<- list(dir = destination, bundle = trace$value, trace = trace,
                   shared = shared,
                   elapsed = as.numeric(difftime(Sys.time(), started,
                                                 units = "secs")))
    cache
  }
})

# A throwaway copy of the shared bundle, for tests that need to move or mutate
# one. Copying is far cheaper than exporting again.
explore_copy_bundle <- function(destination) {
  source_dir <- explore_shared_bundle()$dir
  unlink(destination, recursive = TRUE, force = TRUE)
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(list.files(source_dir, full.names = TRUE, all.files = TRUE,
                             no.. = TRUE),
                  destination, recursive = TRUE, copy.date = TRUE)
  stopifnot(all(ok))
  destination
}
