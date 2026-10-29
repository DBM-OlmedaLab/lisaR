# PROPOSAL ONLY. Suggested package path:
# tests/testthat/test-function-code-fingerprint.R

metadata_hash <- function(fun) lisa_function_code_sha256(list(target = fun))
metadata_source_function <- function(filename) {
  text <- "function(x = 1L, ...) { x + 2L }"
  src <- srcfilecopy(filename, text, isFile = FALSE)
  fun <- eval(parse(text = text, srcfile = src, keep.source = TRUE)[[1L]],
    envir = new.env(parent = baseenv()))
  list(fun = fun, src = src)
}

test_that("same syntax ignores distinct and later-mutated parser source environments", {
  a <- metadata_source_function("synthetic-one.R")
  b <- metadata_source_function("synthetic-two.R")
  expect_false(is.null(attr(body(a$fun), "srcref")))
  old <- function(fun) digest::digest(body(fun), algo = "sha256",
    serialize = TRUE, serializeVersion = 3L)
  # Reproduce the possible failure mechanism, not attribution of the saved runs.
  expect_false(identical(old(a$fun), old(b$fun)))
  before_old <- old(a$fun)
  before_new <- metadata_hash(a$fun)
  a$src$synthetic_cache <- list(warmed = TRUE, marker = "not code")
  expect_false(identical(before_old, old(a$fun)))
  expect_identical(metadata_hash(a$fun), before_new)
  expect_identical(metadata_hash(b$fun), before_new)
  expect_false(is.null(attr(body(a$fun), "srcref")))
  expect_match(before_new, "^[a-f0-9]{64}$")
})

test_that("formatting, execution environment and compilation are not function syntax", {
  a <- eval(parse(text = "function(x = 1L) { x + 2L }", keep.source = TRUE),
    new.env(parent = baseenv()))
  b <- eval(parse(text = "function ( x = 1L ) {\n # comment\n x + 2L\n}",
    keep.source = TRUE), new.env(parent = baseenv()))
  environment(a)$unrelated <- "cold"
  before <- metadata_hash(a)
  environment(a)$unrelated <- "warm"
  expect_identical(metadata_hash(a), before)
  expect_identical(metadata_hash(b), before)
  expect_identical(metadata_hash(compiler::cmpfun(a)), before)
  nested_a <- eval(parse(text = "function(x) { function(y = x) { y + 1L } }",
    keep.source = TRUE), new.env(parent = baseenv()))
  nested_b <- eval(parse(text = "function(x) {\n # offset\n function(y = x) {\n y + 1L\n }\n}",
    keep.source = TRUE), new.env(parent = baseenv()))
  nested_no_source <- eval(parse(text = "function(x) { function(y = x) { y + 1L } }",
    keep.source = FALSE), new.env(parent = baseenv()))
  expect_identical(metadata_hash(nested_a), metadata_hash(nested_b))
  expect_identical(metadata_hash(nested_a), metadata_hash(nested_no_source))
})

test_that("real body, defaults, argument names, ordering and literal types change hashes", {
  a <- function(x = 1L, y = 2L, ...) { x + y }
  variants <- list(
    function(x = 1L, y = 2L, ...) { x - y },
    function(x = 3L, y = 2L, ...) { x + y },
    function(z = 1L, y = 2L, ...) { z + y },
    function(y = 2L, x = 1L, ...) { x + y },
    function(x = 1, y = 2L, ...) { x + y },
    function(x = 1L, y = 2L) { x + y }
  )
  for (b in variants) expect_false(identical(metadata_hash(a), metadata_hash(b)))
  pairs <- list(
    list(function(x) x, function(x = NULL) x),
    list(function() f(a = 1L), function() f(b = 1L)),
    list(function() f(, NULL), function() f(NULL, )),
    list(function() NA_integer_, function() NA_real_),
    list(function() 1L, function() 1),
    list(function() 0, function() -0),
    list(function() "NA", function() NA_character_),
    list(function() quote(alpha), function() "alpha")
  )
  for (pair in pairs) expect_false(identical(metadata_hash(pair[[1L]]), metadata_hash(pair[[2L]])))
})

test_that("nested syntax, empty defaults and delimiter-like Unicode strings are unambiguous", {
  a <- function(x, ...) { function(y = x) list(, y, "a:2{}\u03b1") }
  expect_identical(metadata_hash(a), metadata_hash(a))
  b <- function(x, ...) { function(y = x) list(, y, "a:2{}\u03b2") }
  expect_false(identical(metadata_hash(a), metadata_hash(b)))
  # Fingerprinting must not evaluate any default or body expression.
  expect_match(metadata_hash(function(x = stop("default evaluated")) stop("body evaluated")),
    "^[a-f0-9]{64}$")
})

test_that("non-source attributes on literal AST nodes are preserved", {
  a <- function() identity(1L)
  body_a <- body(a)
  body_a[[2L]] <- structure(1L, semantic_label = "one")
  body(a) <- body_a
  b <- a
  body_b <- body(b)
  body_b[[2L]] <- structure(1L, semantic_label = "two")
  body(b) <- body_b
  expect_false(identical(metadata_hash(a), metadata_hash(b)))
})

test_that("function labels and list order are part of the explicit fingerprint contract", {
  a <- function(x) x
  b <- function(x) x + 1L
  one <- lisa_function_code_sha256(list(first = a, second = b))
  expect_false(identical(one, lisa_function_code_sha256(list(second = b, first = a))))
  expect_false(identical(one, lisa_function_code_sha256(list(renamed = a, second = b))))
})

test_that("malformed input and opaque AST objects fail closed", {
  a <- function(x) x
  bad <- list(NULL, list(), a, list(a), list("not a function"), list(target = 1L),
    list(target = sum), structure(list(target = a), class = "unexpected"),
    stats::setNames(list(a), NA_character_), stats::setNames(list(a), ""),
    stats::setNames(list(a), "bad\nname"), stats::setNames(list(a, a), c("dup", "dup")))
  for (value in bad) expect_error(lisa_function_code_sha256(value), "LISA-CODE-001")
  opaque <- function() NULL
  body(opaque) <- call("identity", new.env(parent = emptyenv()))
  expect_error(metadata_hash(opaque), "unsupported function-syntax node type: environment")
  malformed <- function() NULL
  body(malformed) <- as.call(list(as.name("function"), pairlist(x = 1L),
    as.name("x"), "not a source reference"))
  expect_error(metadata_hash(malformed), "unsupported nested-function source-reference slot")
})
