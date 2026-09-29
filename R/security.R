# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Security run-root confinement and subprocess contracts.

# Technical identifiers become path components and file-name fragments. Keep
# their grammar deliberately narrower than human-facing labels so one config
# has the same identity on Linux, Windows and macOS. In particular, do not
# transliterate or otherwise rewrite Unicode here: labels retain Unicode, while
# callers must choose an explicit ASCII technical identifier.
lisa_portable_id_pattern <- function() {
  "^[A-Za-z0-9]([A-Za-z0-9._-]*[A-Za-z0-9_-])?$"
}

lisa_windows_device_id_pattern <- function() {
  "^(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])([.]|$)"
}

lisa_safe_id <- function(id, field = "id") {
  if (length(id) != 1L) {
    stop(field, " must be one non-empty identifier.", call. = FALSE)
  }
  id <- tryCatch(as.character(id), error = function(e) NA_character_)
  if (length(id) != 1L || is.na(id) || !nzchar(id) || !nzchar(trimws(id))) {
    stop(field, " must be one non-empty identifier.", call. = FALSE)
  }
  if (grepl("[[:cntrl:]]|[/\\\\]", id) || id %in% c(".", "..")) {
    stop(
      "Unsafe ", field,
      ": identifiers cannot contain paths, traversal, or control characters.",
      call. = FALSE
    )
  }
  if (!grepl(
    lisa_portable_id_pattern(), id, perl = TRUE, useBytes = TRUE
  )) {
    stop(
      "Unsafe ", field, ": technical identifiers must use only ASCII ",
      "letters, numbers, dot, hyphen, or underscore; start with a letter or ",
      "number; and end with a letter, number, hyphen, or underscore. ",
      "Use a separate human-readable label for Unicode or spaces.",
      call. = FALSE
    )
  }
  if (grepl(
    lisa_windows_device_id_pattern(), id,
    ignore.case = TRUE, perl = TRUE, useBytes = TRUE
  )) {
    stop(
      "Unsafe ", field, ": ", id,
      " is a reserved Windows device name, including when followed by an extension.",
      call. = FALSE
    )
  }
  id
}

# External scientific keys (for example KEGG's `hsa:1`) are not technical
# identifiers and must remain unchanged in metadata and APIs. Derive only their
# path component here. Already-portable keys retain their historical cache
# names, except for the reserved encoding namespace; every other UTF-8 byte
# sequence has an injective hexadecimal representation.
lisa_external_key_path_component <- function(key, field = "external key") {
  if (length(key) != 1L) {
    stop(field, " must be one non-empty key.", call. = FALSE)
  }
  key <- tryCatch(as.character(key), error = function(e) NA_character_)
  if (length(key) != 1L || is.na(key) || !nzchar(key) ||
      !nzchar(trimws(key)) || grepl("[[:cntrl:]]", key)) {
    stop(field, " must be one non-empty key without control characters.",
         call. = FALSE)
  }

  encoding_prefix <- "lisa-key-v1-"
  portable <- tryCatch(lisa_safe_id(key, field), error = function(e) NULL)
  if (!is.null(portable) && !startsWith(portable, encoding_prefix)) {
    return(portable)
  }

  utf8_key <- enc2utf8(key)
  encoded <- paste(
    sprintf("%02x", as.integer(charToRaw(utf8_key))), collapse = ""
  )
  component <- paste0(encoding_prefix, encoded)
  lisa_safe_id(component, paste0(field, " path component"))
}

lisa_path_walk_lexically <- function(path, base = getwd()) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) ||
      grepl("[[:cntrl:]]", path)) {
    stop("Unsafe filesystem path.", call. = FALSE)
  }
  slash <- function(value) gsub("\\\\", "/", path.expand(as.character(value)))
  raw <- slash(path)
  absolute <- lisa_is_absolute_path(raw)
  if (!absolute) raw <- paste0(sub("/+$", "", slash(base)), "/", raw)

  if (grepl("^[A-Za-z]:/", raw)) {
    anchor <- substr(raw, 1L, 3L)
    remainder <- substring(raw, 4L)
  } else if (startsWith(raw, "///")) {
    anchor <- "/"
    remainder <- sub("^/+", "", raw)
  } else if (startsWith(raw, "//")) {
    unc <- strsplit(sub("^/+", "", raw), "/", fixed = TRUE)[[1]]
    if (length(unc) < 2L || any(!nzchar(unc[1:2]))) {
      stop("Unsafe UNC filesystem path.", call. = FALSE)
    }
    anchor <- paste0("//", unc[[1]], "/", unc[[2]])
    remainder <- paste(unc[-c(1L, 2L)], collapse = "/")
  } else if (startsWith(raw, "/")) {
    anchor <- "/"
    remainder <- sub("^/+", "", raw)
  } else {
    stop("Filesystem path could not be made absolute.", call. = FALSE)
  }

  pieces <- if (nzchar(remainder)) {
    strsplit(remainder, "/", fixed = TRUE)[[1]]
  } else {
    character()
  }
  # Each component's absolute spelling is its parent's spelling plus one piece.
  # Carrying that forward keeps the walk linear in the depth of the path; the
  # earlier form re-pasted the whole retained prefix once per component, which
  # profiling showed to be a measurable share of every guarded write.
  root_prefix <- if (identical(anchor, "/")) "" else sub("/+$", "", anchor)
  kept <- character()
  visited <- character()
  for (piece in pieces) {
    if (!nzchar(piece) || identical(piece, ".")) next
    if (identical(piece, "..")) {
      if (length(kept)) kept <- kept[-length(kept)]
      next
    }
    parent <- if (length(kept)) kept[[length(kept)]] else root_prefix
    current <- paste0(parent, "/", piece)
    kept <- c(kept, current)
    # Keep every lexical component, including one later followed by `..`.
    # POSIX resolves such a component before resolving `..`, so dropping it
    # before calling Sys.readlink() would recreate the original escape bug.
    visited <- c(visited, current)
  }
  components <- kept
  collapsed <- if (!length(components)) anchor else utils::tail(components, 1L)
  list(path = collapsed, visited = unique(visited), components = components)
}

lisa_path_state <- function(path) {
  state <- tryCatch(
    fs::file_info(path, follow = FALSE),
    error = function(error) {
      stop(
        "Could not inspect a filesystem entry without following links: ",
        path, ". ", conditionMessage(error), call. = FALSE
      )
    }
  )
  if (nrow(state) != 1L) {
    stop("Filesystem inspection returned an unexpected result: ", path,
         call. = FALSE)
  }
  state
}

lisa_path_is_link <- function(path) {
  state <- lisa_path_state(path)
  state_type <- as.character(state$type[[1L]])
  identical(state_type, "symlink")
}

# One no-follow inspection covering several paths at once. This is the same
# primitive and the same `follow = FALSE` contract as lisa_path_state(): fs
# reports one row per requested path, in order, and a type of NA for a path
# that does not exist, exactly as the single-path call does. Only the number of
# calls changes, and with it the per-call tibble construction that dominated
# the profile of a full run. It is deliberately NOT used where a guard depends
# on inspecting the same entry twice as two separate observations.
lisa_path_are_links <- function(paths) {
  paths <- as.character(paths)
  if (!length(paths)) return(logical())
  state <- tryCatch(
    fs::file_info(paths, follow = FALSE),
    error = function(error) {
      # Re-inspect one by one so the offending path is still named exactly as
      # a single-path inspection would have named it.
      for (path in paths) lisa_path_state(path)
      stop(
        "Could not inspect a filesystem entry without following links: ",
        paths[[1L]], ". ", conditionMessage(error), call. = FALSE
      )
    }
  )
  if (nrow(state) != length(paths)) {
    stop("Filesystem inspection returned an unexpected result: ",
         paste(paths, collapse = ", "), call. = FALSE)
  }
  observed <- as.character(state$type)
  !is.na(observed) & observed == "symlink"
}

lisa_link_delete <- function(path) {
  if (!lisa_path_is_link(path)) return(invisible(FALSE))
  tryCatch(
    fs::link_delete(path),
    error = function(error) {
      stop(
        "Could not remove a filesystem link without following it: ", path,
        ". ", conditionMessage(error), call. = FALSE
      )
    }
  )
  if (lisa_path_entry_exists(path)) {
    stop("Filesystem entry remained after no-follow link removal: ", path,
         call. = FALSE)
  }
  invisible(TRUE)
}

lisa_path_link_targets <- function(path, base = getwd()) {
  walked <- lisa_path_walk_lexically(path, base = base)
  if (!length(walked$visited)) {
    return(stats::setNames(character(), character()))
  }
  linked <- lisa_path_are_links(walked$visited)
  stats::setNames(rep("<filesystem-link>", sum(linked)),
                  walked$visited[linked])
}

lisa_path_has_symlink <- function(path, base = getwd()) {
  length(lisa_path_link_targets(path, base = base)) > 0L
}

lisa_assert_no_symlink_components <- function(path, base = getwd()) {
  links <- lisa_path_link_targets(path, base = base)
  if (length(links)) {
    stop(
      "Symbolic links are not allowed in filesystem paths: ",
      paste(names(links), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(lisa_path_walk_lexically(path, base = base)$path)
}

lisa_path_key <- function(path) {
  key <- sub("/+$", "", gsub("\\\\", "/", as.character(path)))
  if (!nzchar(key)) key <- "/"
  if (.Platform$OS.type == "windows") tolower(key) else key
}

lisa_path_within <- function(path, root) {
  path <- lisa_path_key(path)
  root <- lisa_path_key(root)
  identical(path, root) || startsWith(path, paste0(root, "/"))
}

# One directory is reachable under several absolute spellings. macOS resolves
# `/var` through a symbolic link to `/private/var`, and Windows still hands out
# 8.3 short names such as `C:/Users/RUNNER~1` beside `C:/Users/runneradmin`.
# `normalizePath()` cannot repair that on its own: with `mustWork = FALSE` it
# returns a not-yet-created path unchanged, so a planned destination keeps the
# caller spelling while its already existing parent is canonicalized. Resolve
# the longest existing ancestor and re-append the remaining lexical components
# so that both operands of a comparison name the same filesystem entry.
lisa_path_canonical <- function(path, base = getwd()) {
  walked <- lisa_path_walk_lexically(path, base = base)
  components <- walked$components
  if (!length(components)) return(walked$path)
  link <- suppressWarnings(Sys.readlink(components))
  present <- file.exists(components) | (!is.na(link) & nzchar(link))
  deepest <- if (any(present)) max(which(present)) else 0L
  if (!deepest) return(walked$path)
  prefix <- tryCatch(
    normalizePath(components[[deepest]], winslash = "/", mustWork = FALSE),
    error = function(error) components[[deepest]]
  )
  prefix <- sub("/+$", "", gsub("\\\\", "/", prefix))
  if (!nzchar(prefix)) prefix <- "/"
  if (identical(deepest, length(components))) return(prefix)
  remainder <- vapply(components[seq.int(deepest + 1L, length(components))],
                      basename, character(1), USE.NAMES = FALSE)
  paste0(prefix, if (identical(prefix, "/")) "" else "/",
         paste(remainder, collapse = "/"))
}

# Spellings that may denote an authorized root. The root itself is trusted, so
# resolving it through a platform alias is intended.
lisa_path_aliases <- function(path, base = getwd()) {
  lexical <- lisa_path_walk_lexically(path, base = base)$path
  unique(c(lexical, lisa_path_canonical(lexical, base = base)))
}

# Spellings that may denote a requested entry. Only the parent prefix is
# resolved: the final component is never followed, so a symbolic link keeps its
# own identity and stays visible to the link guards.
lisa_path_entry_aliases <- function(path, base = getwd()) {
  lexical <- lisa_path_walk_lexically(path, base = base)$path
  parent <- dirname(lexical)
  if (identical(parent, lexical)) return(lexical)
  unique(c(lexical,
           paste0(sub("/+$", "", lisa_path_canonical(parent, base = base)),
                  "/", basename(lexical))))
}

# Location of `path` relative to `root` when both name the same managed subtree
# under any platform spelling, `""` when they are the same entry, or NULL when
# `path` is genuinely outside `root`.
lisa_path_relative_within <- function(path, root, base = getwd()) {
  roots <- lisa_path_aliases(root, base = base)
  targets <- lisa_path_entry_aliases(path, base = base)
  for (candidate_root in roots) {
    trimmed <- sub("/+$", "", candidate_root)
    if (!nzchar(trimmed)) trimmed <- "/"
    offset <- if (identical(trimmed, "/")) 2L else nchar(trimmed) + 2L
    for (candidate in targets) {
      if (identical(lisa_path_key(candidate), lisa_path_key(trimmed))) {
        return("")
      }
      if (lisa_path_within(candidate, trimmed)) {
        return(substring(candidate, offset))
      }
    }
  }
  NULL
}

lisa_path_entry_exists <- function(path) {
  state <- lisa_path_state(path)
  state_type <- as.character(state$type[[1L]])
  length(state_type) == 1L && !is.na(state_type)
}

lisa_set_private_mode <- function(path, directory = FALSE, mode = NULL) {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("A private-mode target must be one non-empty path.", call. = FALSE)
  }
  path <- as.character(path)
  if (lisa_path_is_link(path)) {
    stop("Symbolic links cannot be assigned managed permissions: ", path,
         call. = FALSE)
  }
  correct_type <- if (isTRUE(directory)) {
    dir.exists(path)
  } else {
    file.exists(path) && !dir.exists(path) &&
      isTRUE(utils::file_test("-f", path))
  }
  if (!isTRUE(correct_type) || lisa_path_is_link(path)) {
    stop("Managed permission target changed or has the wrong type: ", path,
         call. = FALSE)
  }
  if (.Platform$OS.type != "windows") {
    requested_mode <- mode %||% if (isTRUE(directory)) "0700" else "0600"
    changed <- Sys.chmod(path, mode = requested_mode, use_umask = FALSE)
    if (length(changed) != 1L || !isTRUE(changed)) {
      stop("Could not assign private permissions: ", path, call. = FALSE)
    }
  }
  if (lisa_path_is_link(path)) {
    stop("Managed permission target changed while permissions were applied: ",
         path, call. = FALSE)
  }
  correct_type_after <- if (isTRUE(directory)) {
    dir.exists(path)
  } else {
    file.exists(path) && !dir.exists(path) &&
      isTRUE(utils::file_test("-f", path))
  }
  if (!isTRUE(correct_type_after)) {
    stop("Managed permission target changed while permissions were applied: ",
         path, call. = FALSE)
  }
  invisible(path)
}

lisa_prepare_external_parent <- function(path, create = FALSE) {
  walked <- lisa_path_walk_lexically(path)
  requested <- walked$path
  if (!dir.exists(requested)) {
    if (!isTRUE(create)) {
      stop("Managed destination parent does not exist: ", requested,
           call. = FALSE)
    }
    if (lisa_path_entry_exists(requested)) {
      stop("Managed destination parent is not a directory: ", requested,
           call. = FALSE)
    }
    missing <- walked$components[!vapply(walked$components, dir.exists,
                                         logical(1))]
    # `dir.create()` also reports failure when a concurrent worker won the race
    # and already created this shared parent. Two analyses of the same run
    # legitimately share `report_pages/evidence/<analysis>`, so a lost race is
    # not an error: recheck the path and only refuse when it is still absent or
    # is not a directory. Nothing else about the guard is relaxed.
    if (!dir.create(requested, recursive = TRUE, showWarnings = FALSE,
                    mode = "0700") && !dir.exists(requested)) {
      stop("Could not create the destination parent: ", requested,
           call. = FALSE)
    }
    if (length(missing)) {
      invisible(vapply(missing[dir.exists(missing)], lisa_set_private_mode,
                       character(1), directory = TRUE))
    }
  }
  # Existing system/HPC ancestors are outside lisaR's managed boundary. They
  # may legitimately be symlinks (/var -> /private/var, $SCRATCH -> a mount),
  # so canonicalize that trusted prefix instead of rejecting it.
  normalizePath(requested, winslash = "/", mustWork = TRUE)
}

lisa_managed_destination <- function(path, create_parent = FALSE) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) ||
      grepl("[[:cntrl:]]", path)) {
    stop("Unsafe managed destination.", call. = FALSE)
  }
  raw <- gsub("\\\\", "/", path.expand(as.character(path)))
  raw_parts <- strsplit(sub("^[A-Za-z]:/|^/+", "", sub("/+$", "", raw)),
                        "/", fixed = TRUE)[[1]]
  if (any(raw_parts %in% c(".", ".."))) {
    stop("Managed destinations may not contain `.` or `..` components.",
         call. = FALSE)
  }
  requested <- lisa_path_walk_lexically(path)$path
  if (identical(lisa_path_key(requested), lisa_path_key("/")) ||
      basename(requested) %in% c("", ".", "..")) {
    stop("Unsafe managed destination.", call. = FALSE)
  }
  if (lisa_path_is_link(requested)) {
    stop("Symbolic links are not allowed at the managed boundary: ",
         requested, call. = FALSE)
  }
  parent <- lisa_prepare_external_parent(dirname(requested),
                                         create = create_parent)
  target <- file.path(parent, basename(requested))
  if (lisa_path_is_link(target)) {
    stop("Symbolic links are not allowed at the managed boundary: ",
         target, call. = FALSE)
  }
  list(path = gsub("\\\\", "/", target), parent = parent,
       requested = requested)
}

lisa_run_root <- function(run_root) {
  if (length(run_root) != 1L || is.na(run_root) || !nzchar(run_root) || grepl("[[:cntrl:]]", run_root)) stop("Unsafe run_root.", call. = FALSE)
  destination <- lisa_managed_destination(run_root, create_parent = TRUE)
  root <- destination$path
  home <- normalizePath(path.expand("~"), winslash = "/", mustWork = TRUE)
  if (identical(lisa_path_key(root), lisa_path_key(home))) {
    stop("Unsafe run_root: a dedicated non-home directory is required.", call. = FALSE)
  }
  if (lisa_path_entry_exists(root) && !dir.exists(root)) {
    stop("Unsafe run_root: the managed boundary is not a directory.",
         call. = FALSE)
  }
  created <- FALSE
  if (!dir.exists(root)) {
    if (!dir.create(root, recursive = FALSE, showWarnings = FALSE,
                    mode = "0700")) {
      stop("Unsafe run_root: could not create the managed directory.",
           call. = FALSE)
    }
    created <- TRUE
  }
  if (lisa_path_is_link(root) || !dir.exists(root)) {
    stop("Unsafe run_root: symbolic links are not allowed at the managed boundary.",
         call. = FALSE)
  }
  if (created) lisa_set_private_mode(root, directory = TRUE)
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

lisa_existing_run_root <- function(run_root) {
  if (length(run_root) != 1L || is.na(run_root) || !nzchar(run_root) ||
      grepl("[[:cntrl:]]", run_root)) {
    stop("Unsafe run_root.", call. = FALSE)
  }
  destination <- lisa_managed_destination(run_root, create_parent = FALSE)
  root <- destination$path
  home <- normalizePath(path.expand("~"), winslash = "/", mustWork = TRUE)
  if (identical(lisa_path_key(root), lisa_path_key(home))) {
    stop("Unsafe run_root: a dedicated non-home directory is required.",
         call. = FALSE)
  }
  if (!dir.exists(root)) {
    stop("Run root does not exist: ", root, call. = FALSE)
  }
  if (lisa_path_is_link(root)) {
    stop("Symbolic links are not allowed at the managed boundary: ", root,
         call. = FALSE)
  }
  normalizePath(root, winslash = "/", mustWork = TRUE)
}

lisa_scan_run_tree <- function(run_root) {
  root <- lisa_existing_run_root(run_root)
  queue <- root
  rows <- list()
  while (length(queue)) {
    directory <- queue[[1]]
    queue <- queue[-1L]
    directory_state <- lisa_path_state(directory)
    directory_type <- as.character(directory_state$type[[1L]])
    if (identical(directory_type, "symlink")) {
      stop("Run directory became a symbolic link: ", directory,
           call. = FALSE)
    }
    if (!identical(directory_type, "directory") ||
        file.access(directory, 5L) != 0L) {
      stop("Run directory is absent or unreadable: ", directory, call. = FALSE)
    }
    directory_type_after <- as.character(
      lisa_path_state(directory)$type[[1L]]
    )
    if (!identical(directory_type_after, "directory")) {
      stop("Run directory became a symbolic link: ", directory,
           call. = FALSE)
    }
    children <- list.files(
      directory, recursive = FALSE, all.files = TRUE, include.dirs = TRUE,
      full.names = TRUE, no.. = TRUE
    )
    if (!length(children)) next
    for (child in children) {
      child_state <- lisa_path_state(child)
      child_type <- as.character(child_state$type[[1L]])
      if (identical(child_type, "symlink")) {
        # Do not call file.info(), dir.exists(), or recursive list.files() on a
        # link: each could follow an existing directory link outside the run.
        stop("Symbolic links are not allowed inside the run_root: ", child,
             call. = FALSE)
      }
      if (!child_type %in% c("file", "directory")) {
        stop("Run entry is non-regular or changed during inspection: ", child,
             call. = FALSE)
      }
      child_type_after <- as.character(lisa_path_state(child)$type[[1L]])
      if (!identical(child_type_after, child_type)) {
        stop("Run entry changed while the tree was inspected: ",
             child, call. = FALSE)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        path = gsub("\\\\", "/", child),
        isdir = identical(child_type, "directory"),
        bytes = if (identical(child_type, "directory")) {
          NA_real_
        } else {
          as.numeric(child_state$size[[1L]])
        },
        stringsAsFactors = FALSE
      )
      if (identical(child_type, "directory")) queue <- c(queue, child)
    }
  }
  if (!length(rows)) {
    return(data.frame(path = character(), isdir = logical(), bytes = numeric(),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  out[order(out$path), , drop = FALSE]
}

lisa_assert_run_tree_safe <- function(run_root) {
  root <- lisa_existing_run_root(run_root)
  lisa_scan_run_tree(root)
  invisible(root)
}

lisa_guarded_path <- function(path, run_root = getOption("lisaR.run_root", NULL)) {
  if (length(path) != 1L || is.na(path) || !nzchar(path) || grepl("[[:cntrl:]]", path)) stop("Unsafe filesystem path.", call. = FALSE)
  # Standalone exported writers remain usable; their immediate parent becomes the
  # explicit authorization boundary when a pipeline run root was not configured.
  if (is.null(run_root)) {
    return(lisa_managed_destination(path, create_parent = TRUE)$path)
  }
  requested_root <- lisa_path_walk_lexically(run_root)$path
  root <- lisa_run_root(run_root)
  raw <- gsub("\\\\", "/", as.character(path))
  absolute <- lisa_is_absolute_path(raw)
  if (!absolute) raw <- file.path(requested_root, raw)
  walked <- lisa_path_walk_lexically(raw)
  requested_target <- walked$path
  links <- lisa_path_link_targets(raw)
  # Compare every spelling of both roots. A link reached through a platform
  # alias of the run_root is still inside the run_root and must be rejected.
  managed_roots <- unique(c(lisa_path_aliases(requested_root),
                            lisa_path_aliases(root)))
  managed_links <- names(links)[vapply(names(links), function(link_path) {
    aliases <- lisa_path_entry_aliases(link_path)
    any(vapply(managed_roots, function(managed_root) {
      any(vapply(aliases, lisa_path_within, logical(1), root = managed_root))
    }, logical(1)))
  }, logical(1))]
  if (length(managed_links)) {
    stop("Symbolic links are not allowed inside the run_root: ",
         paste(managed_links, collapse = ", "), call. = FALSE)
  }
  relative <- lisa_path_relative_within(requested_target, requested_root)
  if (!is.null(relative)) {
    target <- if (nzchar(relative)) file.path(root, relative) else root
  } else {
    target <- requested_target
  }
  if (is.null(lisa_path_relative_within(target, root))) {
    raw_components <- strsplit(
      sub("^[A-Za-z]:/|^/+", "", gsub("\\\\", "/", raw)),
      "/", fixed = TRUE
    )[[1]]
    if (any(raw_components %in% c(".", ".."))) {
      stop("Filesystem path escapes the authorized run_root.", call. = FALSE)
    }
    # Resolve only the external parent prefix. This permits callers that used a
    # canonicalizable parent alias while still enforcing containment after it.
    destination <- lisa_managed_destination(raw, create_parent = FALSE)
    target <- destination$path
    if (is.null(lisa_path_relative_within(target, root))) {
      stop("Filesystem path escapes the authorized run_root.", call. = FALSE)
    }
  }
  target
}

lisa_guarded_dir_create <- function(path, run_root = getOption("lisaR.run_root", NULL)) {
  lisa_assert_portable_path(lisa_path_canonical(path))
  if (is.null(run_root)) {
    destination <- lisa_managed_destination(path, create_parent = TRUE)
    if (lisa_path_entry_exists(destination$path) &&
        !dir.exists(destination$path)) {
      stop("A filesystem entry blocks the required directory: ",
           destination$path, call. = FALSE)
    }
    created <- FALSE
    if (!dir.exists(destination$path)) {
      # A concurrent worker of the same run may create this shared directory
      # between the check above and the call below. Losing that race is not a
      # failure; only a still-absent path is. The link and type guards below
      # still run on whatever now occupies the path.
      created <- dir.create(destination$path, recursive = FALSE,
                            showWarnings = FALSE, mode = "0700")
      if (!created && !dir.exists(destination$path)) {
        stop("Could not create a protected directory: ", destination$path,
             call. = FALSE)
      }
    }
    if (lisa_path_is_link(destination$path)) {
      stop("Symbolic links are not allowed at the managed boundary: ",
           destination$path, call. = FALSE)
    }
    if (created) lisa_set_private_mode(destination$path, directory = TRUE)
    return(invisible(destination$path))
  }
  path <- lisa_guarded_path(path, run_root)
  root <- lisa_run_root(run_root)
  walked <- lisa_path_walk_lexically(path)
  candidates <- walked$components[vapply(walked$components, lisa_path_within,
                                         logical(1), root = root)]
  for (candidate in candidates) {
    if (lisa_path_is_link(candidate)) {
      stop("Symbolic links are not allowed inside the run_root: ", candidate,
           call. = FALSE)
    }
    if (lisa_path_entry_exists(candidate)) {
      if (!dir.exists(candidate)) {
        stop("A filesystem entry blocks the required directory: ", candidate,
             call. = FALSE)
      }
      next
    }
    if (dir.create(candidate, recursive = FALSE, showWarnings = FALSE)) {
      lisa_set_private_mode(candidate, directory = TRUE)
    } else if (!dir.exists(candidate)) {
      stop("Could not create a protected directory: ", candidate, call. = FALSE)
    } else if (lisa_path_is_link(candidate)) {
      # Another worker of this run created the shared directory first. Repeat
      # the link and type guards on the winner's entry rather than trusting it.
      stop("Symbolic links are not allowed inside the run_root: ", candidate,
           call. = FALSE)
    }
    lisa_guarded_path(candidate, root)
  }
  lisa_guarded_path(path, root)
  invisible(path)
}

lisa_atomic_temp_path <- function(parent, label = "artifact", destination = NULL,
                                  create_file = FALSE) {
  extension <- if (!is.null(destination)) tools::file_ext(destination) else ""
  fileext <- if (nzchar(extension)) paste0(".", extension) else ""
  for (attempt in seq_len(100L)) {
    candidate <- tempfile(
      paste0(".lisa-", gsub("[^A-Za-z0-9._-]", "-", label), "-"),
      tmpdir = parent, fileext = fileext
    )
    if (lisa_path_entry_exists(candidate)) next
    lisa_assert_portable_path(candidate)
    if (isTRUE(create_file)) {
      if (!file.create(candidate, showWarnings = FALSE)) next
      lisa_set_private_mode(candidate, directory = FALSE)
    }
    return(candidate)
  }
  stop("Could not reserve a protected temporary path.", call. = FALSE)
}

lisa_unlink_no_follow <- function(path) {
  if (!lisa_path_entry_exists(path)) return(invisible(FALSE))
  if (lisa_path_is_link(path)) {
    return(lisa_link_delete(path))
  }
  if (dir.exists(path)) {
    if (lisa_path_is_link(path)) {
      return(lisa_link_delete(path))
    }
    children <- list.files(
      path, recursive = FALSE, all.files = TRUE, include.dirs = TRUE,
      full.names = TRUE, no.. = TRUE
    )
    for (child in children) lisa_unlink_no_follow(child)
    # Re-check the directory itself immediately before removal. Base R cannot
    # bind this check and the removal to one file descriptor; see the TOCTOU
    # note in lisa_atomic_install().
    if (lisa_path_is_link(path)) {
      return(lisa_link_delete(path))
    } else {
      status <- unlink(path, recursive = TRUE, force = FALSE)
    }
    return(invisible(identical(as.integer(status), 0L)))
  }
  # Re-read the entry immediately before removal. A last-moment replacement
  # by a link is removed through fs::link_delete(), never through its target.
  if (lisa_path_is_link(path)) return(lisa_link_delete(path))
  status <- unlink(path, recursive = FALSE, force = FALSE)
  invisible(identical(as.integer(status), 0L))
}

lisa_remove_private_temp <- function(path) {
  lisa_unlink_no_follow(path)
  invisible(TRUE)
}

lisa_assert_regular_managed_file <- function(path, run_root = NULL,
                                             label = "Managed artifact") {
  if (lisa_path_is_link(path)) {
    stop(label, " became a symbolic link: ", path, call. = FALSE)
  }
  checked <- if (is.null(run_root)) {
    lisa_managed_destination(path, create_parent = FALSE)$path
  } else {
    lisa_guarded_path(path, lisa_existing_run_root(run_root))
  }
  regular <- file.exists(checked) && !dir.exists(checked) &&
    isTRUE(utils::file_test("-f", checked))
  if (!isTRUE(regular) || lisa_path_is_link(checked)) {
    stop(label, " must remain one regular non-link file: ", checked,
         call. = FALSE)
  }
  # Repeat the containment check after inspecting the entry. Base R cannot
  # bind these operations to one descriptor, but this avoids every known
  # direct follow of a writer-controlled link.
  if (is.null(run_root)) {
    lisa_managed_destination(checked, create_parent = FALSE)
  } else {
    lisa_guarded_path(checked, lisa_existing_run_root(run_root))
  }
  invisible(checked)
}

lisa_filesystem_recovery_abort <- function(message, source = "",
                                           destination = "", backup = "",
                                           staged = "",
                                           recovery_paths = character(),
                                           parent = NULL) {
  condition <- structure(
    list(
      message = paste0(
        "LISA-FS-RECOVERY-001 ", message,
        if (length(recovery_paths)) {
          paste0(" Recovery paths: ", paste(recovery_paths, collapse = "; "))
        } else {
          ""
        }
      ),
      call = NULL,
      code = "LISA-FS-RECOVERY-001",
      source = as.character(source %||% ""),
      destination = as.character(destination %||% ""),
      backup = as.character(backup %||% ""),
      staged = as.character(staged %||% ""),
      recovery_paths = unique(as.character(recovery_paths)),
      parent = parent
    ),
    class = c(
      "lisa_filesystem_recovery_error", "lisa_error", "error", "condition"
    )
  )
  stop(condition)
}

lisa_existing_recovery_paths <- function(...) {
  paths <- unique(as.character(unlist(list(...), use.names = FALSE)))
  paths <- paths[!is.na(paths) & nzchar(paths)]
  paths[vapply(paths, lisa_path_entry_exists, logical(1))]
}

lisa_atomic_install <- function(staged, destination, overwrite, run_root,
                                source = staged) {
  if (is.null(run_root)) {
    destination_info <- lisa_managed_destination(destination,
                                                 create_parent = FALSE)
    staged_info <- lisa_managed_destination(staged, create_parent = FALSE)
    parent <- destination_info$parent
    destination <- destination_info$path
    staged <- staged_info$path
    root <- NULL
  } else {
    root <- lisa_existing_run_root(run_root)
    parent <- lisa_guarded_path(dirname(destination), root)
    destination <- lisa_guarded_path(destination, root)
    staged <- lisa_guarded_path(staged, root)
  }
  if (!identical(lisa_path_key(dirname(staged)), lisa_path_key(parent))) {
    stop("Protected temporary paths must share the destination parent.",
         call. = FALSE)
  }
  lisa_assert_regular_managed_file(
    staged, run_root = root,
    label = "Protected temporary artifact"
  )
  if (lisa_path_entry_exists(destination) && !isTRUE(overwrite)) {
    stop("Refusing to overwrite an existing file: ", destination, call. = FALSE)
  }
  if (dir.exists(destination)) {
    stop("Atomic file installation cannot replace a directory.", call. = FALSE)
  }
  staged_receipt <- c(
    bytes = as.numeric(file.info(staged)$size),
    sha256 = lisa_sha256_file(staged)
  )
  lisa_assert_regular_managed_file(
    staged, run_root = root,
    label = "Protected temporary artifact"
  )
  if (!identical(as.numeric(file.info(staged)$size),
                 as.numeric(staged_receipt[["bytes"]]))) {
    stop("Protected temporary artifact changed while it was receipted.",
         call. = FALSE)
  }

  backup <- NULL

  if (lisa_path_entry_exists(destination)) {
    # Existing directories are deliberately not replaced by file writers or
    # copies. Guarded directory replacement is only used by guarded rename.
    backup <- lisa_atomic_temp_path(parent, "backup")
    if (!file.rename(destination, backup)) {
      stop("Could not stage the existing destination for atomic replacement.",
           call. = FALSE)
    }
  }

  install_error <- tryCatch({
    # Base R has no openat()/O_NOFOLLOW interface. These checks close the known
    # lexical and dangling-link escapes, but an attacker able to mutate the
    # same directory concurrently leaves a small residual TOCTOU window
    # between this revalidation and file.rename().
    if (is.null(root)) {
      lisa_prepare_external_parent(parent, create = FALSE)
      lisa_managed_destination(destination, create_parent = FALSE)
      lisa_managed_destination(staged, create_parent = FALSE)
    } else {
      lisa_guarded_path(parent, root)
      lisa_guarded_path(destination, root)
      lisa_guarded_path(staged, root)
    }
    lisa_assert_regular_managed_file(
      staged, run_root = root,
      label = "Protected temporary artifact"
    )
    if (!file.rename(staged, destination)) {
      stop("Atomic filesystem promotion failed: ", destination,
           call. = FALSE)
    }
    if (is.null(root)) {
      lisa_managed_destination(destination, create_parent = FALSE)
    } else {
      lisa_guarded_path(destination, root)
    }
    lisa_assert_regular_managed_file(
      destination, run_root = root, label = "Installed artifact"
    )
    installed_receipt <- c(
      bytes = as.numeric(file.info(destination)$size),
      sha256 = lisa_sha256_file(destination)
    )
    lisa_assert_regular_managed_file(
      destination, run_root = root, label = "Installed artifact"
    )
    if (!identical(as.character(installed_receipt),
                   as.character(staged_receipt))) {
      stop("Post-install verification failed.", call. = FALSE)
    }
    NULL
  }, error = identity)

  if (inherits(install_error, "error")) {
    rollback_failures <- character()
    if (lisa_path_entry_exists(destination)) {
      if (!lisa_path_entry_exists(staged)) {
        if (!isTRUE(file.rename(destination, staged))) {
          rollback_failures <- c(
            rollback_failures,
            "installed artifact could not be returned to staging"
          )
        }
      } else {
        rollback_failures <- c(
          rollback_failures,
          "both destination and staging are occupied"
        )
      }
    }
    if (!is.null(backup) && lisa_path_entry_exists(backup)) {
      if (!lisa_path_entry_exists(destination)) {
        if (!isTRUE(file.rename(backup, destination))) {
          rollback_failures <- c(
            rollback_failures,
            "the original destination backup could not be restored"
          )
        }
      } else {
        rollback_failures <- c(
          rollback_failures,
          "the destination remained occupied before backup restoration"
        )
      }
    }
    if (length(rollback_failures)) {
      lisa_filesystem_recovery_abort(
        paste0(
          conditionMessage(install_error), " Rollback was incomplete: ",
          paste(rollback_failures, collapse = "; "), "."
        ),
        source = source, destination = destination,
        backup = backup %||% "", staged = staged,
        recovery_paths = lisa_existing_recovery_paths(
          source, destination, backup, staged
        ),
        parent = install_error
      )
    }
    stop(install_error)
  }

  # The old destination remains intact under `backup` until the installed
  # bytes and digest have both been independently checked above.
  if (!is.null(backup) && lisa_path_entry_exists(backup)) {
    removed <- lisa_unlink_no_follow(backup)
    if (!isTRUE(removed) || lisa_path_entry_exists(backup)) {
      lisa_filesystem_recovery_abort(
        "The new destination was verified, but its retained backup could not be removed.",
        source = source, destination = destination, backup = backup,
        staged = staged,
        recovery_paths = lisa_existing_recovery_paths(
          source, destination, backup, staged
        )
      )
    }
  }
  invisible(destination)
}

lisa_private_tree_modes <- function(root) {
  root <- lisa_assert_run_tree_safe(root)
  if (.Platform$OS.type == "windows") return(invisible(root))
  tree <- lisa_scan_run_tree(root)
  if (nrow(tree)) {
    for (i in seq_len(nrow(tree))) {
      entry <- lisa_guarded_path(tree$path[[i]], root)
      lisa_set_private_mode(entry, directory = isTRUE(tree$isdir[[i]]))
      lisa_guarded_path(entry, root)
    }
  }
  # The managed root itself is owned by lisaR for staging/output promotion.
  lisa_set_private_mode(root, directory = TRUE)
  lisa_assert_run_tree_safe(root)
  invisible(root)
}

lisa_promote_managed_directory <- function(staging, destination,
                                           run_root = NULL) {
  if (is.null(run_root)) {
    staging <- lisa_assert_run_tree_safe(staging)
    destination_info <- lisa_managed_destination(
      destination, create_parent = TRUE
    )
    destination <- destination_info$path
  } else {
    root <- lisa_existing_run_root(run_root)
    staging <- lisa_guarded_path(staging, root)
    destination <- lisa_guarded_path(destination, root)
    lisa_assert_run_tree_safe(staging)
  }
  if (lisa_path_entry_exists(destination)) {
    stop("Refusing to overwrite an existing managed directory: ",
         destination, call. = FALSE)
  }
  lisa_private_tree_modes(staging)
  lisa_assert_run_tree_safe(staging)
  if (lisa_path_is_link(destination)) {
    stop("Symbolic links are not allowed at the managed destination.",
         call. = FALSE)
  }
  # Staging and destination are on the same parent filesystem in every lisaR
  # transaction. Promotion is one rename, with no second temporary hop.
  if (!file.rename(staging, destination)) {
    stop("Atomic managed-directory promotion failed: ", destination,
         call. = FALSE)
  }
  verification_error <- tryCatch({
    lisa_assert_run_tree_safe(destination)
    NULL
  }, error = identity)
  if (inherits(verification_error, "error")) {
    restored <- !lisa_path_entry_exists(staging) &&
      isTRUE(file.rename(destination, staging))
    if (!isTRUE(restored)) {
      lisa_filesystem_recovery_abort(
        paste0(
          "Post-promotion directory verification failed and the directory ",
          "could not be returned to staging: ",
          conditionMessage(verification_error), "."
        ),
        source = staging, destination = destination, backup = "",
        staged = staging,
        recovery_paths = lisa_existing_recovery_paths(staging, destination),
        parent = verification_error
      )
    }
    stop("Post-promotion directory verification failed: ",
         conditionMessage(verification_error), call. = FALSE)
  }
  invisible(destination)
}

lisa_guarded_write <- function(path, writer, run_root = getOption("lisaR.run_root", NULL), overwrite = TRUE) {
  path <- lisa_guarded_path(path, run_root)
  lisa_assert_portable_path(path)
  root <- if (is.null(run_root)) NULL else lisa_run_root(run_root)
  if (dir.exists(path)) stop("Refusing to replace a directory with a file: ", path, call. = FALSE)
  if (lisa_path_entry_exists(path) && !isTRUE(overwrite)) stop("Refusing to overwrite an existing file: ", path, call. = FALSE)
  parent <- if (is.null(root)) {
    lisa_prepare_external_parent(dirname(path), create = TRUE)
  } else {
    lisa_guarded_dir_create(dirname(path), root)
  }
  staged <- lisa_atomic_temp_path(parent, "write", destination = path,
                                  create_file = TRUE)
  recovery_required <- FALSE
  on.exit(if (!recovery_required && lisa_path_entry_exists(staged)) {
    lisa_remove_private_temp(staged)
  }, add = TRUE)
  writer(staged)
  lisa_assert_regular_managed_file(
    staged, run_root = root, label = "Guarded writer output"
  )
  lisa_set_private_mode(staged, directory = FALSE)
  lisa_assert_regular_managed_file(
    staged, run_root = root, label = "Guarded writer output"
  )
  if (is.null(root)) {
    lisa_managed_destination(path, create_parent = FALSE)
    lisa_managed_destination(staged, create_parent = FALSE)
  } else {
    lisa_guarded_path(parent, root)
    lisa_guarded_path(path, root)
    lisa_guarded_path(staged, root)
  }
  tryCatch(
    lisa_atomic_install(staged, path, overwrite, root, source = staged),
    lisa_filesystem_recovery_error = function(error) {
      recovery_required <<- TRUE
      stop(error)
    }
  )
  invisible(path)
}

lisa_file_copy <- function(from, to) {
  file.copy(
    from, to, overwrite = TRUE, copy.mode = FALSE, copy.date = FALSE
  )
}

lisa_guarded_copy <- function(from, to, overwrite = FALSE, run_root = getOption("lisaR.run_root", NULL)) {
  to <- lisa_guarded_path(to, run_root)
  lisa_assert_portable_path(to)
  root <- if (is.null(run_root)) NULL else lisa_run_root(run_root)
  if (dir.exists(to)) stop("Refusing to replace a directory with a file: ", to, call. = FALSE)
  if (lisa_path_entry_exists(to) && !isTRUE(overwrite)) stop("Refusing to overwrite an existing file: ", to, call. = FALSE)
  if (lisa_path_is_link(from) || !file.exists(from) ||
      dir.exists(from) || !isTRUE(utils::file_test("-f", from))) {
    stop("Guarded copy source is not one regular file: ", from, call. = FALSE)
  }
  if (lisa_path_is_link(from)) {
    stop("Guarded copy source became a symbolic link: ", from, call. = FALSE)
  }
  parent <- if (is.null(root)) {
    lisa_prepare_external_parent(dirname(to), create = TRUE)
  } else {
    lisa_guarded_dir_create(dirname(to), root)
  }
  staged <- lisa_atomic_temp_path(parent, "copy", destination = to,
                                  create_file = TRUE)
  recovery_required <- FALSE
  on.exit(if (!recovery_required && lisa_path_entry_exists(staged)) {
    lisa_remove_private_temp(staged)
  }, add = TRUE)
  lisa_assert_regular_managed_file(
    staged, run_root = root, label = "Guarded copy temporary"
  )
  copy_warnings <- character()
  copied <- withCallingHandlers(
    lisa_file_copy(from, staged),
    warning = function(warning) {
      diagnostic <- trimws(gsub(
        "[[:cntrl:]]+", " ", conditionMessage(warning)
      ))
      if (nzchar(diagnostic)) {
        copy_warnings <<- c(copy_warnings, diagnostic)
      }
      invokeRestart("muffleWarning")
    }
  )
  copy_warnings <- unique(copy_warnings)
  if (!isTRUE(copied)) {
    diagnostic <- if (length(copy_warnings)) {
      paste0(
        "; filesystem diagnostic: ", paste(copy_warnings, collapse = " | ")
      )
    } else {
      ""
    }
    stop("Guarded file copy failed: ", from, diagnostic, call. = FALSE)
  }
  if (length(copy_warnings)) {
    warning(paste(copy_warnings, collapse = " | "), call. = FALSE)
  }
  lisa_assert_regular_managed_file(
    staged, run_root = root, label = "Guarded copy output"
  )
  if (lisa_path_is_link(from) || !file.exists(from) ||
      dir.exists(from) || !isTRUE(utils::file_test("-f", from))) {
    stop("Guarded copy source changed during the copy: ", from,
         call. = FALSE)
  }
  source_mode <- as.integer(file.info(from)$mode)
  if (lisa_path_is_link(from)) {
    stop("Guarded copy source became a symbolic link: ", from, call. = FALSE)
  }
  private_mode <- as.integer(as.octmode("0600"))
  if (.Platform$OS.type != "windows" && !is.na(source_mode)) {
    lisa_set_private_mode(
      staged, directory = FALSE,
      mode = as.octmode(bitwAnd(source_mode, private_mode))
    )
  }
  lisa_assert_regular_managed_file(
    staged, run_root = root, label = "Guarded copy output"
  )
  if (is.null(root)) {
    lisa_managed_destination(to, create_parent = FALSE)
    lisa_managed_destination(staged, create_parent = FALSE)
  } else {
    lisa_guarded_path(parent, root)
    lisa_guarded_path(to, root)
    lisa_guarded_path(staged, root)
  }
  tryCatch(
    lisa_atomic_install(staged, to, overwrite, root, source = from),
    lisa_filesystem_recovery_error = function(error) {
      recovery_required <<- TRUE
      stop(error)
    }
  )
  invisible(to)
}

lisa_guarded_rename <- function(from, to, overwrite = FALSE, run_root = getOption("lisaR.run_root", NULL)) {
  from <- lisa_guarded_path(from, run_root)
  to <- lisa_guarded_path(to, run_root)
  lisa_assert_portable_path(to)
  root <- if (is.null(run_root)) NULL else lisa_existing_run_root(run_root)
  if (!lisa_path_entry_exists(from)) stop("Guarded rename source does not exist: ", from, call. = FALSE)
  if (lisa_path_entry_exists(to) && !isTRUE(overwrite)) stop("Refusing to overwrite an existing file: ", to, call. = FALSE)
  if (dir.exists(from)) {
    if (isTRUE(overwrite)) {
      stop("Guarded directory promotion does not overwrite existing destinations.",
           call. = FALSE)
    }
    lisa_promote_managed_directory(from, to, run_root = root)
    return(invisible(to))
  }
  parent <- if (is.null(root)) {
    lisa_prepare_external_parent(dirname(to), create = TRUE)
  } else {
    lisa_guarded_dir_create(dirname(to), root)
  }
  staged <- lisa_atomic_temp_path(parent, "rename", destination = to)
  moved <- FALSE
  recovery_required <- FALSE
  restore_source <- function(parent_error = NULL) {
    if (!lisa_path_entry_exists(staged)) return(invisible(TRUE))
    restored <- !lisa_path_entry_exists(from) &&
      isTRUE(file.rename(staged, from))
    if (isTRUE(restored)) return(invisible(TRUE))
    recovery_required <<- TRUE
    lisa_filesystem_recovery_abort(
      paste0(
        "Guarded rename failed and the staged source could not be ",
        "restored to its original path."
      ),
      source = from, destination = to, staged = staged,
      recovery_paths = lisa_existing_recovery_paths(from, to, staged),
      parent = parent_error
    )
  }
  on.exit({
    if (!recovery_required && !moved && lisa_path_entry_exists(staged)) {
      restore_source()
    }
    # If rollback fails, retain the staged source bytes instead of deleting
    # the only remaining copy. Successful moves leave no staged entry.
    if (!recovery_required && moved && lisa_path_entry_exists(staged)) {
      lisa_remove_private_temp(staged)
    }
  }, add = TRUE)
  if (!file.rename(from, staged)) {
    stop("Guarded file rename failed while staging source: ", from, call. = FALSE)
  }
  if (is.null(root)) {
    lisa_managed_destination(to, create_parent = FALSE)
    lisa_managed_destination(staged, create_parent = FALSE)
  } else {
    lisa_guarded_path(parent, root)
    lisa_guarded_path(to, root)
    lisa_guarded_path(staged, root)
  }
  install_error <- tryCatch({
    lisa_atomic_install(staged, to, overwrite, root, source = from)
    NULL
  },
    lisa_filesystem_recovery_error = function(error) {
      recovery_required <<- TRUE
      stop(error)
    },
    error = identity
  )
  if (inherits(install_error, "error")) {
    restore_source(install_error)
    stop(install_error)
  }
  moved <- TRUE
  invisible(to)
}

lisa_guarded_delete <- function(path, recursive = FALSE, run_root = getOption("lisaR.run_root", NULL)) {
  if (is.null(run_root)) {
    destination <- lisa_managed_destination(path, create_parent = FALSE)
    path <- destination$path
    root <- NULL
    parent <- destination$parent
  } else {
    root <- lisa_existing_run_root(run_root)
    path <- lisa_guarded_path(path, root)
    parent <- dirname(path)
  }
  if (!is.null(root) && identical(lisa_path_key(path), lisa_path_key(root))) stop("Refusing to delete the run_root.", call. = FALSE)
  if (!lisa_path_entry_exists(path)) return(invisible(path))
  if (dir.exists(path)) {
    if (!isTRUE(recursive)) stop("Refusing recursive directory deletion without `recursive = TRUE`.", call. = FALSE)
    lisa_assert_run_tree_safe(path)
  }
  if (is.null(root)) {
    lisa_prepare_external_parent(parent, create = FALSE)
    lisa_managed_destination(path, create_parent = FALSE)
  } else {
    lisa_guarded_path(parent, root)
    lisa_guarded_path(path, root)
  }
  removal_path <- path
  if (dir.exists(path) && isTRUE(recursive)) {
    quarantine <- lisa_atomic_temp_path(parent, "delete")
    if (!file.rename(path, quarantine)) {
      stop("Could not quarantine a directory before recursive deletion: ",
           path, call. = FALSE)
    }
    removal_path <- quarantine
  }
  removed <- lisa_unlink_no_follow(removal_path)
  if (!isTRUE(removed) || lisa_path_entry_exists(removal_path)) {
    stop("Guarded delete failed; retained quarantine: ", removal_path,
         call. = FALSE)
  }
  invisible(path)
}

# Run `code` with `env` exported to this process, then restore the previous
# state exactly (including variables that were previously unset).
#
# `system2(env = )` cannot be used for this. It is not portable: see ?system2,
# "On Windows, 'env' is only supported for commands such as 'R' and 'make'
# which accept environment variables on their command line". On Windows R
# appends the NAME=value tokens to the child's *command line*, and `Rscript`
# reads the first non-option argument as the script to execute -- so the child
# runs `R_LIBS=...` instead of the requested script and never starts. Exporting
# the variables here instead gives the child the same environment on every
# platform and removes shell quoting from the path entirely.
#
# The call this wraps is synchronous, and parallel workers are separate
# processes with their own environment blocks, so the window is confined to
# this process for the duration of one subprocess.
lisa_with_child_environment <- function(env, code) {
  if (!length(env)) return(code)
  requested <- names(env)
  previous <- Sys.getenv(requested, unset = NA_character_, names = TRUE)
  on.exit({
    restore <- previous[!is.na(previous)]
    if (length(restore)) do.call(Sys.setenv, as.list(restore))
    unset <- requested[is.na(previous)]
    if (length(unset)) Sys.unsetenv(unset)
  }, add = TRUE)
  do.call(Sys.setenv, as.list(env))
  code
}

lisa_run_subprocess <- function(executable, args = character(), required = TRUE, stage = "subprocess",
                                run_root = getOption("lisaR.run_root", NULL),
                                env = character()) {
  if (length(executable) != 1L || is.na(executable) || !nzchar(executable) || grepl("[[:cntrl:]]", executable)) stop("Unsafe subprocess executable.", call. = FALSE)
  args <- as.character(args)
  if (anyNA(args) || any(grepl("[[:cntrl:]]", args))) stop("Unsafe subprocess arguments.", call. = FALSE)
  if (length(env) && (is.null(names(env)) || anyNA(names(env)) ||
      any(!nzchar(names(env))) || anyDuplicated(names(env)))) {
    stop("Subprocess environment must have unique, non-empty names.", call. = FALSE)
  }
  env_names <- names(env)
  env <- as.character(env)
  names(env) <- env_names
  if (length(env) && (anyNA(env) || any(!grepl("^[A-Za-z_][A-Za-z0-9_]*$", names(env))) ||
      any(grepl("[[:cntrl:]]", env)))) {
    stop("Unsafe subprocess environment.", call. = FALSE)
  }
  if (is.null(run_root)) stop("An authorized run_root is required for subprocess diagnostics.", call. = FALSE)
  root <- lisa_run_root(run_root)
  diagnostics_root <- lisa_guarded_dir_create(file.path(root, ".lisa_subprocess"), root)
  stage_prefix <- gsub("[^A-Za-z0-9._-]+", "-", as.character(stage[[1]]))
  diagnostics_dir <- lisa_guarded_dir_create(
    tempfile(paste0(stage_prefix, "-", Sys.getpid(), "-"), tmpdir = diagnostics_root),
    root
  )
  stdout_path <- lisa_atomic_temp_path(
    diagnostics_dir, "stdout", destination = "stdout.txt", create_file = TRUE
  )
  stderr_path <- lisa_atomic_temp_path(
    diagnostics_dir, "stderr", destination = "stderr.txt", create_file = TRUE
  )
  on.exit({
    if (file.exists(stdout_path)) lisa_guarded_delete(stdout_path, run_root = root)
    if (file.exists(stderr_path)) lisa_guarded_delete(stderr_path, run_root = root)
    if (dir.exists(diagnostics_dir)) lisa_guarded_delete(diagnostics_dir, recursive = TRUE, run_root = root)
    # diagnostics_root is shared by concurrent workers.  A child may remove
    # only its own directory; removing the shared root races other subprocesses.
  }, add = TRUE)
  exit_code <- tryCatch(
    {
      lisa_assert_regular_managed_file(stdout_path, root,
                                       "Subprocess stdout")
      lisa_assert_regular_managed_file(stderr_path, root,
                                       "Subprocess stderr")
      lisa_with_child_environment(env, suppressWarnings(system2(
        executable, args = shQuote(args), stdout = stdout_path,
        stderr = stderr_path
      )))
    },
    error = function(e) structure(127L, message = conditionMessage(e))
  )
  launch_error <- attr(exit_code, "message")
  if (is.null(exit_code)) exit_code <- 0L
  lisa_assert_regular_managed_file(stdout_path, root, "Subprocess stdout")
  lisa_assert_regular_managed_file(stderr_path, root, "Subprocess stderr")
  stdout <- readLines(stdout_path, warn = FALSE)
  stderr <- readLines(stderr_path, warn = FALSE)
  if (!is.null(launch_error)) stderr <- c(stderr, launch_error)
  result <- list(executable = executable, arguments = args, exit_code = as.integer(exit_code), stdout = stdout, stderr = stderr, stage = stage)
  result$diagnostics <- paste(c(result$stdout, result$stderr), collapse = "\n")
  if (result$exit_code != 0L && isTRUE(required)) {
    stop("Required stage failed [", stage, "] with exit code ", result$exit_code, ": ", result$diagnostics, call. = FALSE)
  }
  result
}
