# Figure-source TSVs can carry multiline titles. Encode those fields explicitly
# rather than writing physical newlines into an unquoted TSV record.
lisa_write_figure_source_tsv <- function(x, path) {
  encoded <- names(x)[vapply(x, function(v) is.character(v) &&
    any(grepl("[\r\n\t]",v),na.rm=TRUE),logical(1))]
  for (name in encoded) x[[name]] <- vapply(x[[name]],function(v)
    as.character(jsonlite::toJSON(v,auto_unbox=TRUE,na="null")),character(1))
  if(length(encoded)) x$encoded_text_columns <- paste(encoded,collapse=";")
  if(inherits(path,"connection")) utils::write.table(x,path,sep="\t",quote=FALSE,
    row.names=FALSE,col.names=TRUE,na="") else write_lisa_tsv(x,path)
}

lisa_read_figure_source_tsv <- function(path) {
  lines <- readLines(path,warn=FALSE)
  if (!length(lines)) return(data.frame())
  header <- strsplit(lines[[1]],"\t",fixed=TRUE)[[1]]
  # Recover only the evidenced legacy KEGG format: an unquoted final subtitle
  # field containing newlines. Never guess how to repair numeric/source columns.
  if (utils::tail(header,1L) == "figure_subtitle") {
    widths <- lengths(strsplit(lines[-1],"\t",fixed=TRUE))
    if (any(widths != length(header))) {
      rebuilt <- character(); current <- NULL
      for (line in lines[-1]) {
        fields <- strsplit(line,"\t",fixed=TRUE)[[1]]
        if (length(fields) == length(header)) {
          if(!is.null(current)) rebuilt <- c(rebuilt,current)
          current <- line
        } else if (length(fields) <= 1L && !is.null(current)) {
          current <- paste0(current,"\\n",line)
        } else stop("Malformed figure TSV outside the known legacy subtitle field: ",path,call.=FALSE)
      }
      if(!is.null(current)) rebuilt <- c(rebuilt,current)
      con <- textConnection(c(lines[[1]],rebuilt));on.exit(close(con),add=TRUE)
      x <- utils::read.delim(con,sep="\t",quote="",comment.char="",check.names=FALSE)
      x$figure_subtitle <- gsub("\\n","\n",x$figure_subtitle,fixed=TRUE)
      return(x)
    }
  }
  x <- read_lisa_tsv(path)
  if ("encoded_text_columns" %in% names(x) && nrow(x)) {
    fields <- unique(as.character(x$encoded_text_columns))
    if(length(fields)!=1L || is.na(fields)) stop("Inconsistent figure text encoding.",call.=FALSE)
    fields <- strsplit(fields,";",fixed=TRUE)[[1]]
    if(!all(fields %in% names(x))) stop("Missing encoded figure fields.",call.=FALSE)
    for(name in fields) x[[name]] <- vapply(x[[name]],function(v) {
      decoded <- jsonlite::fromJSON(v)
      if(is.null(decoded)) return(NA_character_)
      if(!is.character(decoded)||length(decoded)!=1L) stop("Invalid encoded figure text.",call.=FALSE)
      decoded
    },character(1))
    x$encoded_text_columns <- NULL
  }
  x
}
