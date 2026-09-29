# Follow the three documented routes with the installed candidate on Windows.
# No upstream DE model is refitted and no historical result is overwritten.
library(lisaR)
args <- commandArgs(TRUE)
example <- args[[1L]]
phase <- if (length(args) >= 2L) args[[2L]] else "all"
stopifnot(phase %in% c("all", "standard", "full"))
stopifnot(example %in% c("sample", "riaz", "cptac"), .Platform$OS.type == "windows")
root <- "D:/lisa test"
dir.create(root, recursive=TRUE, showWarnings=FALSE)
evidence <- file.path(Sys.getenv("RUNNER_TEMP"), paste0("lisa-", example))
dir.create(evidence, recursive=TRUE, showWarnings=FALSE)
writeLines(c(R.version.string, as.character(packageVersion("lisaR")), find.package("lisaR")), file.path(evidence,"environment.txt"))
options(lisaR.dictionary_cache_root = file.path(root,"cache"))
Sys.setenv(R_LIBS=paste(.libPaths(),collapse=.Platform$path.sep),NCPUS="1")
project <- file.path(root,example)
rscript <- file.path(R.home("bin"),"Rscript")
invoke <- function(script, args=character(), label="prepare") {
  cat("TUTORIAL_STAGE_START=", label, " ", format(Sys.time(), tz="UTC"), " UTC\n", sep="")
  log <- file.path(evidence,paste0(label,".log"))
  # Stream child progress to Actions and retain an independent log.
  raw_args <- sub('^"(.*)"$', '\\1', args)
  raw_args <- sub("^'(.*)'$", "\\1", raw_args)
  capture <- function(lines, ...) {
    cat(paste(lines, collapse="\n"), "\n")
    cat(paste(lines, collapse="\n"), "\n", file=log, append=TRUE)
  }
  result <- processx::run(rscript,c("--vanilla",script,raw_args),
    error_on_status=FALSE, stdout_line_callback=capture,
    stderr_line_callback=capture)
  status <- result$status
  if (status != 0L) {
    cat(paste(tail(readLines(log, warn=FALSE), 60L), collapse="\n"), "\n")
    stop("Tutorial stage failed: ", label, "; exit status ", status)
  }
  cat("TUTORIAL_STAGE_PASS=", label, " ", format(Sys.time(), tz="UTC"), " UTC\n", sep="")
}
if (phase != "full") {
if (example == "sample") {
  project <- lisa_init_project(project)
  # The standard route uses the same saved inputs; full is a second
  # presentation product obtained from its verified results below.
} else {
  # Same fixed human resource as the reviewed examples; institutional
  # research validation, no redistribution of the membership in artifacts.
  msigdb <- prepare_lisa_msigdb_resource(accept_terms=TRUE)
  id <- if (example=="riaz") "riaz-gse91061" else "cptac-ccrcc"
  dir.create(file.path(root,"bundles"),showWarnings=FALSE)
  bundle <- install_lisa_example_bundle(id,destination=file.path(root,"bundles",example))
  if (example == "riaz") {
    Sys.setenv(LISAR_RIAZ_PROJECT_DIR=project,LISAR_RIAZ_PREPARED_BUNDLE=bundle$bundle,
      LISAR_RIAZ_RESOURCE_CACHE=file.path(root,"cache"))
    script <- system.file("examples/riaz-gse91061/scripts/11_prepare_prepared_project.R",package="lisaR")
  } else {
    Sys.setenv(LISAR_CPTAC_CCRCC_PROJECT_DIR=project,LISAR_CPTAC_CCRCC_BUNDLE=bundle$bundle,
      LISAR_CPTAC_CCRCC_RESOURCE_CACHE=file.path(root,"cache"),LISAR_CPTAC_CCRCC_PREPARE_ONLY="1")
    script <- system.file("examples/cptac-ccrcc/prepare_and_run_cptac_ccrcc.R",package="lisaR")
  }
  invoke(script)
  stopifnot(lisaR:::lisa_verify_compact_example_project(project))
  file.copy(file.path(project,"provenance","path_map.tsv"),evidence)
}
}
config <- file.path(project,"study.yml")
if (example == "sample") {
  study <- yaml::read_yaml(config)
  study$pipeline$workers <- 1L
  study$report$recipes <- TRUE
  config <- file.path(project,"run.yml")
  yaml::write_yaml(study,config)
}
validation <- validate_lisa_config(config,check_files=TRUE,strict=TRUE)
stopifnot(validation$execution_ready,
  validation$analyses==switch(example,sample=2L,riaz=2L,cptac=1L),
  validation$contrasts==switch(example,sample=1L,riaz=1L,cptac=0L))
plan <- plan_lisa_outputs(config)
write.table(plan$summary,file.path(evidence,"plan.tsv"),sep="\t",row.names=FALSE)
if (phase != "full") {
if (example == "riaz") {
  invoke(system.file("examples/riaz-gse91061/scripts/09_run_lisa_example.R",package="lisaR"),label="run")
  # Idempotent preparation must still recognise the unmodified scientific inputs.
  invoke(system.file("examples/riaz-gse91061/scripts/11_prepare_prepared_project.R",package="lisaR"),label="recheck")
} else {
  run <- run_lisa(config)
  stopifnot(identical(normalizePath(run$output_dir,winslash="/"),normalizePath(file.path(project,"out"),winslash="/")))
}
}
output <- file.path(project,"out")
check <- verify_lisa_run(output)
stopifnot(identical(check$gate,"PASS"),file.exists(file.path(output,"report_index.html")))
writeLines(c("STANDARD_VERIFIED_PASS", format(Sys.time(), tz="UTC")), file.path(evidence,"standard-verified.txt"))
if (phase == "standard") quit(status=0L)
if (example != "cptac") {
  selection <- if (example=="riaz") c("--analyses","responders_vs_pd_pre,responders_vs_pd_on",
    "--contrasts","response_separation_by_time") else c("--analyses","response_a,response_b", "--contrasts","response_a_vs_b")
  invoke(system.file("scripts/build_selected_LISA_report.R",package="lisaR"),c(
    "--source-dir",shQuote(output),"--output-dir",shQuote(file.path(project,"extra")),
    selection,"--complete-missing","true","--title",shQuote("Example results")),label="full")
  stopifnot(file.exists(file.path(project,"extra","report_index.html")))
}
if (example == "sample") {
  ws <- lisa_explore_open(output,file.path(project,"explore"))
  request <- lisa_figure_request(unit_type="single_de",analysis_id="response_a",
    collection="GOBP-C2",category_id="SYN_SIGNAL",product="volcano")
  lisa_explore_submit(ws,request,background=FALSE)
  exported <- lisa_explore_export(ws,file.path(project,"export"))
  stopifnot(exported$extensions==1L,file.exists(file.path(project,"export","report_index.html")),
    identical(verify_lisa_run(output)$gate,"PASS"))
}
# Actual output inventory, plus temporary-sibling lengths at every written
# directory. Write receipts only; do not upload data, figures or full reports.
files <- list.files(project,recursive=TRUE,full.names=TRUE,all.files=TRUE)
files <- files[!file.info(files)$isdir]
rel <- substring(files,nchar(project)+2L)
stopifnot(!any(grepl("teaching",rel,ignore.case=TRUE)))
units <- lisaR:::lisa_path_utf16_length(files)
tmp <- file.path(dirname(files),".lisa-write-123456789abc.json")
stopifnot(all(units<=240),all(lisaR:::lisa_path_utf16_length(tmp)<=240))
write.table(data.frame(path=rel,utf16=units),file.path(evidence,"paths.tsv"),sep="\t",row.names=FALSE)
if (example == "sample") {
  cfg <- yaml::read_yaml(config)
  cfg$pipeline$output_dir <- paste0("D:/",paste(rep(strrep("a",40),5),collapse="/"),"/out")
  bad <- file.path(project,"long.yml"); yaml::write_yaml(cfg,bad)
  error <- tryCatch(run_lisa(bad),error=identity)
  stopifnot(inherits(error,"error"),grepl("LISA-PATH-001",conditionMessage(error)),
    !dir.exists(cfg$pipeline$output_dir))
  writeLines(conditionMessage(error),file.path(evidence,"long-path-rejected.txt"))
}
jsonlite::write_json(list(example=example,status="PASS",platform=R.version$platform,
  package=as.character(packageVersion("lisaR")),analyses=validation$analyses,
  contrasts=validation$contrasts,files=length(files),max_utf16=max(units),
  root_has_spaces=grepl(" ",project)),file.path(evidence,"result.json"),pretty=TRUE,auto_unbox=TRUE)
cat("NATIVE_TUTORIAL_PASS=",example,"\n",sep="")
