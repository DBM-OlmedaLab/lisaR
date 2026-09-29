# Copyright (C) 2026 Agencia Estatal Consejo Superior de Investigaciones
# Científicas (CSIC). Author: David Olmeda Casadomé.
# This file is part of lisaR, free software under the GNU General Public
# License version 3 (GPL-3). See the DESCRIPTION file and
# <https://www.gnu.org/licenses/gpl-3.0.html>.

# Optional, local, content-addressed reuse of raw fgsea results.
# This is passive typed JSON, not R serialization: no cache read can restore
# an environment, function, S4 object, or arbitrary class. Doubles are stored
# as IEEE754 bytes so JSON decimal conversion never changes scientific values.

lisa_content_cache_options <- function(cache_dir = NULL, cache_mode = "off",
                                       cache_max_bytes = 512 * 1024^2) {
  cache_mode <- match.arg(cache_mode, c("off", "readwrite", "readonly", "refresh"))
  if (!is.numeric(cache_max_bytes) || length(cache_max_bytes) != 1L ||
      is.na(cache_max_bytes) || !is.finite(cache_max_bytes) ||
      cache_max_bytes < 1 || cache_max_bytes != floor(cache_max_bytes) || cache_max_bytes > 2^53 - 1) {
    stop("LISA-CACHE-001: cache_max_bytes must be one positive whole number.", call. = FALSE)
  }
  if (!identical(cache_mode, "off") &&
      (!is.character(cache_dir) || length(cache_dir) != 1L ||
       is.na(cache_dir) || !nzchar(cache_dir) || grepl("[[:cntrl:]]", cache_dir))) {
    stop("LISA-CACHE-001: active caching requires one explicit local cache_dir.", call. = FALSE)
  }
  list(dir = cache_dir, mode = cache_mode, max_bytes = as.numeric(cache_max_bytes))
}

lisa_cache_digest <- function(object) {
  digest::digest(object, algo = "sha256", serialize = TRUE, serializeVersion = 3L)
}

lisa_fgsea_cache_versions <- function() {
  packages <- c("fgsea", "BiocParallel", "data.table", "Rcpp", "BH", "fastmatch", "parallel", "stats")
  versions <- vapply(packages, function(package) {
    tryCatch(as.character(utils::packageVersion(package)), error = function(e) "unavailable")
  }, character(1))
  defaults <- lapply(c("fgseaSimple", "fgseaMultilevel", "fgsea"), function(name) {
    fun <- get0(name, envir = asNamespace("fgsea"), inherits = FALSE)
    if (is.null(fun)) return(NULL)
    lapply(formals(fun), function(value) paste(deparse(value), collapse = "\n"))
  })
  list(R = R.version.string, platform = R.version$platform,
       packages = versions, fgsea_defaults = defaults,
       external_libraries = extSoftVersion(), locale = Sys.getlocale(),
       thread_environment = lisa_thread_environment(),
       wrapper = lisa_cache_digest(list(formals(run_fgsea_lisa), body(run_fgsea_lisa))))
}

lisa_fgsea_cache_contract <- function(pathways, ranks, min_gs_size, max_gs_size,
                                      n_threads, fgsea_nperm, fgsea_eps,
                                      random_seed, task_id, universe,
                                      versions = lisa_fgsea_cache_versions()) {
  # No intersection, sorting, or tier labels are substituted for the complete
  # selected membership structure and the actual ordered named rank vector.
  parts <- list(
    schema = "lisa_fgsea_content_v1",
    ranks = lisa_cache_digest(ranks),
    memberships = lisa_cache_digest(pathways),
    universe = lisa_cache_digest(universe),
    parameters = lisa_cache_digest(list(
      minSize = min_gs_size, maxSize = max_gs_size, nproc = n_threads,
      nperm = fgsea_nperm, eps = fgsea_eps,
      algorithm = if (is.null(fgsea_nperm)) "multilevel" else "permutation"
    )),
    rng = lisa_cache_digest(list(
      root_seed = random_seed, task_id = task_id,
      effective_seed = lisa_task_seed(random_seed, task_id),
      kind = c("L'Ecuyer-CMRG", RNGkind()[2:3])
    )),
    versions = lisa_cache_digest(versions)
  )
  list(key = lisa_cache_digest(parts), parts = parts)
}

lisa_cache_raw_hex <- function(value) {
  paste(format(value), collapse = "")
}

lisa_cache_hex_raw <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      nchar(value, type = "bytes") %% 2L || grepl("[^0-9a-f]", value)) {
    stop("Invalid cache byte encoding.")
  }
  if (!nzchar(value)) return(raw())
  starts <- seq.int(1L, nchar(value, type = "bytes"), by = 2L)
  as.raw(strtoi(substring(value, starts, starts + 1L), base = 16L))
}

lisa_cache_encode_vector <- function(value) {
  if (length(setdiff(names(attributes(value)), "names"))) stop("Unsupported cache vector attributes.")
  value_names <- names(value)
  value <- unname(value)
  type <- typeof(value)
  if (type %in% c("double", "integer", "logical")) {
    bytes <- writeBin(if (identical(type, "logical")) as.integer(value) else value,
                      raw(), size = if (identical(type, "double")) 8L else 4L,
                      endian = "little")
    data <- lisa_cache_raw_hex(bytes)
  } else if (identical(type, "character")) {
    data <- list(values = as.list(replace(value, is.na(value), "")), missing = as.list(is.na(value)))
  } else if (identical(type, "list") &&
             all(vapply(value, is.character, logical(1)))) {
    data <- lapply(value, lisa_cache_encode_vector)
  } else stop("Unsupported cache vector type.")
  list(type = type, data = data,
       names = if (is.null(value_names)) NULL else as.list(value_names))
}

lisa_cache_decode_vector <- function(encoded) {
  if (!is.list(encoded) || !identical(names(encoded), c("type", "data", "names")) ||
      !is.character(encoded$type) || length(encoded$type) != 1L) stop("Invalid cache vector schema.")
  type <- encoded$type
  if (type %in% c("double", "integer", "logical")) {
    bytes <- lisa_cache_hex_raw(encoded$data)
    size <- if (identical(type, "double")) 8L else 4L
    if (length(bytes) %% size) stop("Invalid cache numeric byte count.")
    value <- readBin(bytes, what = if (identical(type, "double")) double() else integer(),
                     n = length(bytes) / size, size = size, endian = "little")
    if (identical(type, "logical")) {
      if (any(!is.na(value) & !value %in% c(0L, 1L))) stop("Invalid cached logical values.")
      value <- as.logical(value)
    }
  } else if (identical(type, "character")) {
    data <- encoded$data
    if (!is.list(data) || !identical(names(data), c("values", "missing")) ||
        !is.list(data$values) || !is.list(data$missing) ||
        length(data$values) != length(data$missing)) stop("Invalid cache character schema.")
    value <- vapply(data$values, function(x) {
      if (!is.character(x) || length(x) != 1L || is.na(x)) stop("Invalid cached string.")
      x
    }, character(1))
    missing <- vapply(data$missing, function(x) {
      if (!is.logical(x) || length(x) != 1L || is.na(x)) stop("Invalid cached missing mask.")
      x
    }, logical(1))
    value[missing] <- NA_character_
  } else if (identical(type, "list") && is.list(encoded$data)) {
    value <- lapply(encoded$data, lisa_cache_decode_vector)
    if (!all(vapply(value, is.character, logical(1)))) stop("Only character lists may be cached.")
  } else stop("Unsupported cached vector type.")
  if (!is.null(encoded$names)) {
    vector_names <- vapply(encoded$names, function(x) {
      if (!is.character(x) || length(x) != 1L || is.na(x)) stop("Invalid cached vector names.")
      x
    }, character(1))
    if (length(vector_names) != length(value)) stop("Invalid cached vector name length.")
    names(value) <- vector_names
  }
  value
}

lisa_cache_encode_result <- function(result) {
  if (!identical(class(result), "data.frame") ||
      length(setdiff(names(attributes(result)), c("names", "row.names", "class", "lisa_notes")))) {
    stop("Only plain fgsea data frames may be cached.")
  }
  list(columns = lapply(unname(result), lisa_cache_encode_vector),
       column_names = lisa_cache_encode_vector(names(result)),
       row_names = lisa_cache_encode_vector(attr(result, "row.names")),
       notes = if (is.null(attr(result, "lisa_notes"))) NULL else
         lisa_cache_encode_vector(attr(result, "lisa_notes")))
}

lisa_cache_decode_result <- function(payload) {
  if (!is.list(payload) ||
      !identical(names(payload), c("columns", "column_names", "row_names", "notes")) ||
      !is.list(payload$columns)) stop("Invalid cached result schema.")
  columns <- lapply(payload$columns, lisa_cache_decode_vector)
  column_names <- lisa_cache_decode_vector(payload$column_names)
  row_names <- lisa_cache_decode_vector(payload$row_names)
  if (!is.character(column_names) || length(column_names) != length(columns) ||
      anyNA(column_names) || anyDuplicated(column_names) ||
      !is.null(names(row_names)) ||
      !typeof(row_names) %in% c("integer", "character") || anyNA(row_names) ||
      anyDuplicated(row_names) || any(lengths(columns) != length(row_names))) {
    stop("Invalid cached table dimensions.")
  }
  result <- structure(columns, names = column_names, class = "data.frame", row.names = row_names)
  if (!is.null(payload$notes)) {
    notes <- lisa_cache_decode_vector(payload$notes)
    if (!is.character(notes)) stop("Invalid cached notes.")
    attr(result, "lisa_notes") <- notes
  }
  result
}

lisa_cache_json <- function(object) {
  as.character(jsonlite::toJSON(object, auto_unbox = TRUE, null = "null", na = "null", digits = NA))
}

lisa_cache_slot <- function(options, key, create = FALSE) {
  # Sixteen direct-mapped slots give a hard bound on committed cache bytes
  # without a global lock or racy count-and-evict admission. Hash collisions
  # replace only a cache artifact, and can never be interpreted as a hit.
  root <- file.path(options$dir, "lisa-fgsea-v1")
  if (!dir.exists(root) && !create) return(NULL)
  if (create) tryCatch(lisa_guarded_dir_create(root, run_root = NULL), error = function(error) {
    # Another authorised caller may have created the same private root after
    # the pre-check. Still enforce the full root policy immediately below.
    if (!dir.exists(root) || lisa_path_is_link(root)) stop(error)
  })
  root <- lisa_existing_run_root(root)
  slots <- file.path(root, paste0("slot-", c(0:9, letters[1:6]), ".json"))
  for (slot in slots[lisa_path_entry_exists_vector(slots)]) {
    lisa_assert_regular_managed_file(slot, root, "Cache slot")
    if (file.info(slot)$size > floor(options$max_bytes / 32)) {
      stop("LISA-CACHE-002: cache capacity differs from existing slots; use a fresh directory or the original cache_max_bytes.", call. = FALSE)
    }
  }
  list(root = root, path = file.path(root, paste0("slot-", substr(key, 1L, 1L), ".json")),
       entry_max_bytes = floor(options$max_bytes / 32))
}

lisa_path_entry_exists_vector <- function(paths) {
  vapply(paths, lisa_path_entry_exists, logical(1))
}

lisa_cache_read <- function(slot, contract) {
  miss <- function(reason) list(hit = FALSE, reason = reason, result = NULL, artifact_sha256 = "")
  if (is.null(slot) || !file.exists(slot$path)) return(miss("absent"))
  tryCatch({
    lisa_assert_regular_managed_file(slot$path, slot$root, "Cache slot")
    info <- file.info(slot$path)
    if (info$size > slot$entry_max_bytes) return(miss("oversized"))
    # Parse one bounded byte snapshot: concurrent replacement cannot splice
    # a metadata file from one writer into a payload from another.
    bytes <- readBin(slot$path, raw(), n = info$size + 1)
    if (length(bytes) != info$size) return(miss("changed_during_read"))
    record <- jsonlite::fromJSON(rawToChar(bytes), simplifyVector = FALSE)
    if (!is.list(record) || !identical(names(record), c("schema", "key", "parts", "payload_sha256", "payload")) ||
        !identical(record$schema, "lisa_fgsea_cache_json_v1")) return(miss("malformed"))
    if (!identical(record$key, contract$key) || !identical(record$parts, contract$parts)) return(miss("content_changed_or_slot_collision"))
    if (!is.character(record$payload) || length(record$payload) != 1L ||
        !identical(record$payload_sha256, digest::digest(record$payload, algo = "sha256", serialize = FALSE))) {
      return(miss("digest_mismatch"))
    }
    result <- lisa_cache_decode_result(jsonlite::fromJSON(record$payload, simplifyVector = FALSE))
    list(hit = TRUE, reason = "exact_content_match", result = result,
         artifact_sha256 = digest::digest(bytes, algo = "sha256", serialize = FALSE))
  }, error = function(error) miss("malformed_or_unreadable"))
}

lisa_cache_write <- function(slot, contract, result) {
  tryCatch({
    if (as.numeric(utils::object.size(result)) > slot$entry_max_bytes) return("entry_exceeds_slot_capacity")
    payload <- lisa_cache_json(lisa_cache_encode_result(result))
    record <- lisa_cache_json(list(
      schema = "lisa_fgsea_cache_json_v1", key = contract$key, parts = contract$parts,
      payload_sha256 = digest::digest(payload, algo = "sha256", serialize = FALSE), payload = payload
    ))
    if (nchar(record, type = "bytes") > slot$entry_max_bytes) return("entry_exceeds_slot_capacity")
    # Reserve at most one bounded staging directory for each slot. Atomic
    # mkdir is a local writer claim, not a locking service; a busy or stale
    # claim only declines caching. A crashed writer cannot accumulate an
    # unbounded series of temporary artifacts. Half the byte budget covers
    # committed slots and half covers all possible staging records.
    staging_dir <- lisa_guarded_path(paste0(slot$path, ".pending"), slot$root)
    if (!dir.create(staging_dir, showWarnings = FALSE, mode = "0700")) {
      return("slot_writer_busy_or_stale")
    }
    on.exit(lisa_remove_private_temp(staging_dir), add = TRUE)
    staged <- file.path(staging_dir, "record.json")
    connection <- file(staged, "wb")
    tryCatch(writeBin(charToRaw(enc2utf8(record)), connection), finally = close(connection))
    lisa_set_private_mode(staged, directory = FALSE)
    lisa_assert_regular_managed_file(staged, slot$root, "Cache staging artifact")
    lisa_guarded_path(slot$path, slot$root)
    if (lisa_path_entry_exists(slot$path)) lisa_assert_regular_managed_file(slot$path, slot$root, "Cache slot")
    # A single same-filesystem rename is atomic. No check-delete-write window,
    # backup swap, or cross-worker eviction is used. Platforms unable to
    # replace an occupied file atomically simply decline this optional write.
    if (!isTRUE(suppressWarnings(file.rename(staged, slot$path)))) return("atomic_replace_unavailable")
    verified <- lisa_cache_read(slot, contract)
    if (!verified$hit) return("concurrent_replacement_or_validation_miss")
    "stored"
  }, error = function(error) "write_failed_or_unsupported_result")
}

lisa_cached_fgsea <- function(pathways, ranks, min_gs_size, max_gs_size, n_threads,
                              fgsea_nperm, fgsea_eps, random_seed, task_id,
                              universe, cache_options) {
  options <- cache_options
  qc <- data.frame(stage = "fgsea", mode = options$mode, status = "disabled",
                   reason = "cache_off", cache_key = "", artifact_sha256 = "",
                   write_status = "not_requested", elapsed_seconds = NA_real_,
                   stringsAsFactors = FALSE)
  for (component in c("ranks", "memberships", "universe", "parameters", "rng", "versions")) {
    qc[[paste0(component, "_sha256")]] <- ""
  }
  start <- proc.time()[["elapsed"]]
  contract <- NULL
  slot <- NULL
  if (!identical(options$mode, "off")) {
    contract <- lisa_fgsea_cache_contract(pathways, ranks, min_gs_size, max_gs_size,
      n_threads, fgsea_nperm, fgsea_eps, random_seed, task_id, universe)
    slot <- lisa_cache_slot(options, contract$key, create = options$mode %in% c("readwrite", "refresh"))
    cached <- if (identical(options$mode, "refresh")) {
      list(hit = FALSE, reason = "refresh_requested", artifact_sha256 = "")
    } else lisa_cache_read(slot, contract)
    qc$status <- if (cached$hit) "hit" else "miss"
    qc$reason <- cached$reason
    qc$cache_key <- contract$key
    for (component in c("ranks", "memberships", "universe", "parameters", "rng", "versions")) {
      qc[[paste0(component, "_sha256")]] <- contract$parts[[component]]
    }
    qc$artifact_sha256 <- cached$artifact_sha256
    if (cached$hit) {
      qc$elapsed_seconds <- proc.time()[["elapsed"]] - start
      return(list(result = cached$result, qc = qc))
    }
  }
  # Exactly the pre-existing seeded fgsea call: hits do not advance the caller
  # RNG and misses/off preserve the established task-seed restoration contract.
  result <- lisa_with_task_seed(random_seed, task_id,
    run_fgsea_lisa(pathways, ranks, min_gs_size, max_gs_size, n_threads, fgsea_nperm, fgsea_eps))
  if (options$mode %in% c("readwrite", "refresh")) {
    qc$write_status <- lisa_cache_write(slot, contract, result)
    if (identical(qc$write_status, "stored")) {
      qc$artifact_sha256 <- lisa_cache_read(slot, contract)$artifact_sha256
    }
  }
  qc$elapsed_seconds <- proc.time()[["elapsed"]] - start
  list(result = result, qc = qc)
}
