#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(lisaR))

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], mustWork = FALSE) else tempfile("lisaR-workers-benchmark-")
dir.create(root, recursive = TRUE, showWarnings = FALSE)
input_dir <- file.path(root, "inputs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

genes <- sprintf("GENE%03d", seq_len(200L))
write_de <- function(path, shift) {
  effect <- sin(seq_along(genes) / 7 + shift) * 2
  pvalue <- pmin(0.99, exp(-abs(effect) * 3))
  padj <- p.adjust(pvalue, method = "BH")
  write.table(data.frame(symbol = genes, log2FoldChange = effect, pvalue = pvalue, padj = padj),
              path, sep = "\t", row.names = FALSE, quote = FALSE)
}
de_a <- file.path(input_dir, "de_a.tsv")
de_b <- file.path(input_dir, "de_b.tsv")
write_de(de_a, 0)
write_de(de_b, 0.8)

collections <- c("PATHWAYS", "GOBP-C2", "GOMF", "GOCC")
prefix <- c(PATHWAYS = "REACTOME", `GOBP-C2` = "GOBP", GOMF = "GOMF", GOCC = "GOCC")
dictionary_rows <- list()
category_rows <- list()
term_rows <- list()
for (collection in collections) {
  for (set_index in seq_len(4L)) {
    gene_set <- sprintf("%s_WORKER_%02d", prefix[[collection]], set_index)
    category <- sprintf("%s_CATEGORY_%02d", gsub("-", "_", collection), set_index)
    selected <- genes[seq.int(set_index, length(genes), by = 7L)][seq_len(24L)]
    dictionary_rows[[length(dictionary_rows) + 1L]] <- data.frame(
      universe = collection, gene_set_id = gene_set, gene_set_name = gene_set,
      source_id = "synthetic", category_id = category,
      category_display_name = category, tier = "core"
    )
    category_rows[[length(category_rows) + 1L]] <- data.frame(
      category_id = category, display_name = category,
      macrogroup_id = "SYN_UNCLASSIFIED",
      macrogroup_name = collection, macrogroup_order = match(collection, collections),
      category_order_within_macrogroup = set_index, notes = "workers benchmark",
      family_id = category, pathway_id = category
    )
    term_rows[[length(term_rows) + 1L]] <- data.frame(
      gs_collection = "C5", gs_subcollection = collection, gs_name = gene_set,
      gs_exact_source = "synthetic", gene_symbol = selected
    )
  }
}
for (set_index in seq_len(4L)) {
  gene_set <- sprintf("HALLMARK_WORKER_%02d", set_index)
  selected <- genes[seq.int(set_index + 2L, length(genes), by = 7L)][seq_len(24L)]
  term_rows[[length(term_rows) + 1L]] <- data.frame(
    gs_collection = "H", gs_subcollection = "", gs_name = gene_set,
    gs_exact_source = "synthetic", gene_symbol = selected
  )
}
dictionary <- file.path(input_dir, "dictionary.tsv")
category_map <- file.path(input_dir, "category_map.tsv")
term2gene <- file.path(input_dir, "term2gene.tsv")
write.table(do.call(rbind, dictionary_rows), dictionary, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, category_rows), category_map, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(do.call(rbind, term_rows), term2gene, sep = "\t", row.names = FALSE, quote = FALSE)

de_index <- data.frame(
  analysis_id = c("A", "B"), de_path = c(de_a, de_b), species = "Homo sapiens",
  symbol_col = "symbol", logfc_col = "log2FoldChange", pvalue_col = "pvalue",
  padj_col = "padj", stringsAsFactors = FALSE
)
contrast_index <- data.frame(
  contrast_id = "A_vs_B", analysis_a = "A", analysis_b = "B",
  stringsAsFactors = FALSE
)

run_one <- function(label, workers) {
  output <- file.path(root, label)
  elapsed <- system.time(lisaR:::run_lisa_pipeline(
    de_index = de_index, contrast_index = contrast_index,
    dictionary_dir = input_dir, term2gene = term2gene,
    dictionary_path = dictionary, category_map_path = category_map,
    output_dir = output,
    collections = c("PATHWAYS", "HALLMARKS", "GOBP-C2", "GOMF", "GOCC"),
    base_dir = input_dir, lisa_dictionary = "core",
    plot_formats = "png", export_formats = "tsv",
    run_gene_level = FALSE, run_reports = FALSE, run_kegg_maps = FALSE,
    workers = workers, dry_run = FALSE
  ))[["elapsed"]]
  parallel_path <- file.path(output, "parallel_execution.tsv")
  if (!file.exists(parallel_path)) {
    stop("Benchmark run did not record parallel_execution.tsv: ", label,
         call. = FALSE)
  }
  execution <- lisaR:::read_lisa_tsv(parallel_path)
  required <- c(
    "phase", "platform", "backend_effective", "workers_requested",
    "workers_effective", "cap_reason"
  )
  if (!nrow(execution) || !all(required %in% names(execution))) {
    stop("Benchmark parallel execution receipt is incomplete: ", label,
         call. = FALSE)
  }
  execution$benchmark_label <- label
  execution$workers_argument <- as.integer(workers)
  execution$workers_requested <- as.integer(execution$workers_requested)
  execution$workers_effective <- as.integer(execution$workers_effective)
  runtime_platform <- lisaR:::lisa_runtime_platform()
  stopifnot(
    all(execution$platform == runtime_platform),
    all(execution$workers_requested == as.integer(workers)),
    all(execution$workers_effective >= 1L),
    all(execution$workers_effective <= execution$workers_requested)
  )
  if (identical(runtime_platform, "linux")) {
    expected_backend <- ifelse(
      execution$workers_effective > 1L, "multicore", "sequential"
    )
    stopifnot(identical(
      as.character(execution$backend_effective), unname(expected_backend)
    ))
  } else {
    stopifnot(
      all(execution$workers_effective == 1L),
      all(execution$backend_effective == "sequential")
    )
  }
  list(root = output, elapsed = unname(elapsed), parallel = execution)
}

runs <- lapply(c(1L, 2L, 4L), function(workers) {
  run_one(sprintf("workers-%d", workers), workers)
})
names(runs) <- c("1", "2", "4")
serial <- runs[["1"]]
serial_root <- serial$root

runtime_platform <- lisaR:::lisa_runtime_platform()
if (identical(runtime_platform, "linux")) {
  capability <- lisaR:::lisa_parallel_plan(workers = 2L, tasks = 2L)
  if (capability$workers_effective > 1L) {
    for (label in c("2", "4")) {
      execution <- runs[[label]]$parallel
      if (!any(execution$workers_effective > 1L &
               execution$backend_effective == "multicore")) {
        stop(
          "Linux exposes multicore capacity, but requested workers=", label,
          " never used the multicore backend.", call. = FALSE
        )
      }
    }
    parallel_validation <- "PASS linux_multicore_observed"
  } else {
    parallel_validation <- paste0(
      "SKIP linux_multicore_not_observable;cap_reason=", capability$cap_reason
    )
  }
} else {
  parallel_validation <- paste0(
    "PASS serial_platform_contract;platform=", runtime_platform
  )
}

parallel_contract <- do.call(rbind, lapply(runs, `[[`, "parallel"))
parallel_contract <- parallel_contract[c(
  "benchmark_label", "workers_argument", "phase", "platform",
  "backend_effective", "workers_requested", "workers_effective", "cap_reason"
)]
parallel_contract_path <- file.path(root, "benchmark_parallel_contract.tsv")
write.table(
  parallel_contract, parallel_contract_path, sep = "\t", row.names = FALSE,
  quote = FALSE
)

scientific_files <- function(path) {
  files <- list.files(file.path(path, "outputs"), recursive = TRUE, full.names = TRUE,
                      pattern = "[.]tsv$")
  files <- files[grepl("/(enrichment|lisa_tables|category_contrasts)/", files)]
  files <- files[!grepl("/qc/", files, fixed = TRUE)]
  relative <- substring(files, nchar(path) + 2L)
  stats::setNames(files, relative)
}
serial_files <- scientific_files(serial_root)
if (!length(serial_files)) {
  stop("Benchmark produced no scientific TSV files.", call. = FALSE)
}
serial_hash <- vapply(serial_files, lisaR:::lisa_sha256_file, character(1))

stable_status <- function(path, name) {
  value <- lisaR:::read_lisa_tsv(file.path(path, name))
  value[c("started_at", "ended_at", "timestamp", "output_dir", "output_path")] <- NULL
  value
}
for (workers in c("2", "4")) {
  candidate_root <- runs[[workers]]$root
  candidate_files <- scientific_files(candidate_root)
  stopifnot(identical(names(candidate_files), names(serial_files)))
  candidate_hash <- vapply(candidate_files, lisaR:::lisa_sha256_file, character(1))
  stopifnot(identical(unname(candidate_hash), unname(serial_hash)))
  stopifnot(identical(stable_status(serial_root, "single_de_status.tsv"),
                      stable_status(candidate_root, "single_de_status.tsv")))
  stopifnot(identical(stable_status(serial_root, "contrast_status.tsv"),
                      stable_status(candidate_root, "contrast_status.tsv")))
}

observed_contract <- vapply(names(runs), function(label) {
  execution <- runs[[label]]$parallel
  sprintf(
    "requested=%s,effective=%s,backend=%s",
    label,
    paste(sort(unique(execution$workers_effective)), collapse = "/"),
    paste(sort(unique(execution$backend_effective)), collapse = "/")
  )
}, character(1))

cat(sprintf(
  paste0(
    "SCIENTIFIC_EQUIVALENCE=PASS | requested_workers=1/2/4 | observed=%s | scientific_tsv=%d | ",
    "elapsed=%.3f/%.3f/%.3fs | speedup_requested_2=%.2fx | ",
    "speedup_requested_4=%.2fx | parallel_validation=%s | contract=%s | root=%s\n"
  ),
  paste(observed_contract, collapse = ";"),
  length(serial_files), runs[["1"]]$elapsed, runs[["2"]]$elapsed, runs[["4"]]$elapsed,
  runs[["1"]]$elapsed / runs[["2"]]$elapsed,
  runs[["1"]]$elapsed / runs[["4"]]$elapsed,
  parallel_validation,
  parallel_contract_path,
  root
))
