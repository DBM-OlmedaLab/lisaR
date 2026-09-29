# A FULL acceptance contract is a list of required families, not just a count
# of whichever files happen to exist. No rendering or scientific calculation.
lisa_full_product_inventory <- function(project_dir, kegg_maps = FALSE) {
  root <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)
  de <- read_lisa_tsv(file.path(root, "config/de_index.tsv"))
  co <- read_lisa_tsv(file.path(root, "config/contrast_index.tsv"))
  result <- list()
  requested_formats <- lisa_requested_plot_formats(root, default = c("png", "pdf"))
  add <- function(owner, collection, family, dirs, pattern = "[.]png$", empty_index = character()) {
    paths <- sort(unique(unlist(lapply(dirs, function(d)
      lisa_plot_files(d, pattern=pattern, formats=requested_formats)), use.names=FALSE)))
    empty <- length(empty_index) && any(vapply(empty_index, function(p)
      file.exists(p) && nrow(read_lisa_tsv(p)) == 0L, logical(1)))
    if (!length(requested_formats)) {
      evidence <- unlist(lapply(dirs, function(d) list.files(d,
        pattern = "[.](tsv|json)$", recursive = TRUE, full.names = TRUE)), use.names = FALSE)
      if (!length(evidence) && !empty) stop("FULL required family evidence absent: ",owner,"/",collection,"/",family,call.=FALSE)
      result[[length(result)+1L]] <<- data.frame(owner=owner,collection=collection,
        family=family,state="image_formats_not_requested",path="",sha256="")
      return(invisible(NULL))
    }
    if (!length(paths) && !empty) stop("FULL required family absent: ",owner,"/",collection,"/",family,call.=FALSE)
    if (!length(paths)) {
      result[[length(result)+1L]] <<- data.frame(owner=owner,collection=collection,
        family=family,state="empty_with_evidence",path="",sha256="")
    } else for (p in paths) {
      stem <- tools::file_path_sans_ext(p)
      source <- paste0(stem,"_source.tsv")
      if (!file.exists(source)) source <- paste0(stem,".tsv")
      heatmap_stem <- sub("(_leading_edge_gene_heatmap)_(zscore|log2|none)$", "\\1", stem)
      if (!file.exists(source)) source <- paste0(heatmap_stem,"_matrix.tsv")
      recipe <- paste0(stem,"_recipe.R")
      if (!file.exists(recipe)) recipe <- paste0(stem,".R")
      if (!file.exists(recipe)) recipe <- paste0(heatmap_stem,"_recipe.R")
      # Formats follow the run's own output policy: a PDF is required only
      # when the run requested PDF (C3 / one source of formats).
      needed <- c(paste0(stem, ".", requested_formats), source, recipe)
      if (!all(file.exists(needed)))
        stop("Incomplete FULL figure contract: ",p,call.=FALSE)
      result[[length(result)+1L]] <<- data.frame(owner=owner,collection=collection,
        family=family,state="present",path=substring(p,nchar(root)+2L),sha256=lisa_sha256_file(p))
    }
  }
  for (owner in as.character(de$analysis_id)) {
    for (d in list.dirs(file.path(root,"outputs/single_de",owner),recursive=FALSE)) {
      collection <- sub("^collection_","",basename(d))
      g <- file.path(root,"outputs/gene_level/single_de",owner,basename(d))
      a <- file.path(root,"artifacts",owner,collection)
      prefix <- paste(owner,collection,"gene_level",sep="_")
      add(owner,collection,"prioritized_genes",c(file.path(g,"category_gene_cards"),file.path(a,"gene_cards")),
        empty_index=file.path(g,"category_gene_cards",paste0(prefix,"_category_gene_cards_index.tsv")))
      add(owner,collection,"volcano",c(file.path(g,"category_volcano_overlays"),file.path(a,"volcano")))
      add(owner,collection,"gene_heatmap",c(file.path(g,"leading_edge_gene_heatmaps"),file.path(a,"heatmap")))
      if (collection != "HALLMARKS") add(owner,collection,"recurrent_genes",file.path(g,"recurrent_gene_screen"))
      if (kegg_maps) add(owner,collection,"kegg_maps",file.path(g,"kegg_painter"),"_painted[.]png$",
        file.path(g,"kegg_painter",paste0(prefix,"_kegg_pathway_painter_index.tsv")))
    }
  }
  for (i in seq_len(nrow(co))) {
    owner <- paste(co$contrast_id[i],co$output_id[i],sep="_")
    for (d in list.dirs(file.path(root,"outputs/category_contrasts",owner),recursive=FALSE)) {
      collection <- sub("^collection_","",basename(d))
      g <- file.path(root,"outputs/gene_level/category_contrasts",owner,basename(d))
      add(owner,collection,"prioritized_genes",file.path(g,"contrast_category_cards"))
      add(owner,collection,"paired_gene_heatmap",g,"_paired_gene_heatmap[.]png$")
      add(owner,collection,"gene_category_network",file.path(g,"contrast_gene_category_network"))
      if(kegg_maps) add(owner,collection,"kegg_maps",file.path(g,"contrast_kegg_pathway_painter"),"_contrast_painted[.]png$")
    }
  }
  out <- if(length(result)) do.call(rbind,result) else data.frame()
  dir.create(file.path(root,"audit/full_products"),recursive=TRUE,showWarnings=FALSE)
  write_lisa_tsv(out,file.path(root,"audit/full_products/required_family_inventory.tsv"))
  out
}
