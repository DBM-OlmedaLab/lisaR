# Interactive regions for an EXISTING NES bitmap. Never save a new figure.
# Geometry is measured with its original plot function, source, settings and
# cairo device metrics on the null sink. The original PNG remains untouched.
lisa_category_figure_regions <- function(source_path, settings_path) {
  m <- jsonlite::read_json(settings_path, simplifyVector = TRUE)
  if (!identical(m$source_sha256, digest::digest(file=source_path,algo='sha256')))
    stop('Category figure source identity mismatch.',call.=FALSE)
  if (!identical(m$presentation_version,'2.0'))
    stop('Clickable regions require the recorded 2.0 category layout.',call.=FALSE)
  x <- utils::read.delim(source_path,check.names=FALSE,quote='"',comment.char='')
  prepared <- structure(list(source=x,metadata=m),class='lisa_category_nes_variants')
  # A graphics device supplies font/layout metrics; the sink emits no file.
  sink <- if (.Platform$OS.type=='windows') 'NUL' else '/dev/null'
  grDevices::png(sink,width=m$width,height=m$height,units='in',res=m$dpi,type='cairo')
  on.exit(grDevices::dev.off(),add=TRUE)
  p <- plot_lisa_category_nes_variant(prepared,m$variant)
  b <- ggplot2::ggplot_build(p); g <- ggplot2::ggplot_gtable(b)
  grid::grid.newpage();grid::grid.draw(g);grid::grid.force()
  out <- list()
  add <- function(id,kind,cx,cy,w,h,side='') {
    out[[length(out)+1L]] <<- data.frame(category_id=id,kind=kind,side=side,
      left=max(0,cx-w/2),top=max(0,cy-h/2),width=w,height=h,
      stringsAsFactors=FALSE)
  }
  for (i in seq_len(nrow(b$layout$layout))) {
    panel <- b$layout$layout[i,]; pn <- paste0('panel-',panel$COL,'-',panel$ROW)
    z <- g$layout[g$layout$name==pn,,drop=FALSE]
    if(nrow(z)!=1L) stop('Cannot locate original category panel.',call.=FALSE)
    grid::seekViewport(paste0(z$name,'.',z$t,'-',z$l,'-',z$b,'-',z$r))
    lo <- grid::deviceLoc(grid::unit(0,'npc'),grid::unit(0,'npc'),valueOnly=TRUE)
    hi <- grid::deviceLoc(grid::unit(1,'npc'),grid::unit(1,'npc'),valueOnly=TRUE)
    params <- b$layout$panel_params[[i]]
    ids <- b$layout$panel_scales_y[[panel$SCALE_Y]]$get_breaks()
    # Text labels and every visible mean/median head use the same ID. No
    # inference from display labels, sorting or colour is allowed.
    yvalues <- b$layout$panel_scales_y[[panel$SCALE_Y]]$map(ids)
    for(j in seq_along(ids)) {
      xy <- b$layout$coord$transform(data.frame(x=0,y=yvalues[j]),params)
      cy <- 1-(lo$y+xy$y*(hi$y-lo$y))/m$height
      label <- x$category_display_name[match(ids[j],x$category_id)]
      tw <- grid::convertWidth(grid::grobWidth(grid::textGrob(label,
        gp=grid::gpar(fontsize=8))), 'in',valueOnly=TRUE)
      # axis labels are right-aligned next to the original panel.
      w <- (tw+.09)/m$width
      add(ids[j],'label',lo$x/m$width-.04/m$width-w/2,cy,w,.15/m$height)
    }
    for(k in seq_along(p$layers)) {
      if(!inherits(p$layers[[k]]$geom,'GeomPoint')) next
      original <- p$layers[[k]]$data; built <- b$data[[k]]
      if(!is.data.frame(original)||!nrow(original)) next
      if(nrow(original)!=nrow(built)) stop('Point identity changed during layout.',call.=FALSE)
      hit <- which(as.character(built$PANEL)==as.character(panel$PANEL)&is.finite(built$x)&is.finite(built$y))
      for(j in hit) {
        xy <- b$layout$coord$transform(built[j,,drop=FALSE],params)
        cx <- (lo$x+xy$x*(hi$x-lo$x))/m$width
        cy <- 1-(lo$y+xy$y*(hi$y-lo$y))/m$height
        diameter <- max(2.8, built$size[j]+.8)/25.4
        add(as.character(original$category_id[j]),'head',cx,cy,
          diameter/m$width,diameter/m$height,
          if('side'%in%names(original)) as.character(original$side[j]) else '')
      }
    }
    grid::upViewport(0)
  }
  if(!length(out)) return(data.frame())
  answer <- do.call(rbind,out)
  if(anyNA(answer)||any(answer$left+answer$width>1.01)||any(answer$top+answer$height>1.01))
    stop('Original figure interaction region outside image.',call.=FALSE)
  answer
}
