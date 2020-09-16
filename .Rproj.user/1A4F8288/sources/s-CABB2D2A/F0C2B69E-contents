
### Plot heatmap column anotation
make_hcb <- function(data, var, name = NULL, fillScale = NULL, ...) {
  hcb <- ggplot(data=data, aes_string(x="index", y=1, fill=var)) +
    geom_raster() +
    #geom_tile(size=0.25, show.legend = TRUE)+
    geom_tile(show.legend = TRUE)+
    #coord_fixed() +
    scale_y_continuous(expand=c(0,0), breaks=1, labels=name) +
    scale_x_continuous(expand=c(0,0)) +
    xlab(NULL) + ylab(NULL) +
    theme(axis.title=element_blank(), 
          axis.ticks=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_text(size=8, face="plain"), 
          plot.margin=unit(c(0.1,0.5,0,0),"lines"), ...) +
    guides(fill=F)
  if(!is.null(fillScale)) hcb <- hcb + fillScale
  return(hcb)
}


###Plot heatmap of phyloseq object
plot_heatmap.2 <- function(ps, sample.label=NULL, taxa.label=NULL, ...) {
  #col_breaks = c(seq(0,0,length=30),  # for green
  #               seq(0,1,length=30),              # for yellow
  #               seq(1,3,length=60))              # for red
 # c(seq(-3,0,length=30),  # for green
  #  seq(0,1,length=30),              # for yellow
   # seq(1,3,length=60))              # for red
  
  hm <- plot_heatmap(ps, taxa.label="Species", sample.order=sample.order, taxa.order = taxa.order)
  hm <- hm + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
  
  #low = "#fee8c8"; mid = "#fdbb84"; high = "#e34a33"; trans = scales::log_trans(4); na.value = "grey90"
  #low = "#e0ecf4"; mid = "#9ebcda"; high = "#8856a7"; trans = scales::log_trans(4); na.value = "grey90"
  
  low = muted("green"); mid = "yellow"; high = muted("maroon"); trans = scales::log_trans(4); na.value = "grey90"
  # From plot_heatmap defaults
  # midpoint
  # breaks
  g_colors <- c("#ffff33", "#00ff33", "#00ffff", "#0000ff","#ff00cc","#ff0033")
  #new_gradient <- scale_fill_gradient2(low = low, mid = mid,
  # high = high, midpoint = 0.5, na.value= na.value, breaks = c(0, 0.2, 0.4, 0.8, 1), labels =c(0, 0.25, 0.5, 0.75, 1), name="Relative\nabundance") #c(0.001, 0.1, 1, 10, 100) trans = trans, 
  new_gradient <- scale_fill_gradientn(colors = g_colors, na.value= na.value, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels =c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0), name="Relative\nabundance") #c(0.001, 0.1, 1, 10, 100) trans = trans, 
  hm <- hm + theme(plot.margin=unit(c(0,0.5,0.5,0.5),"lines"), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  hm <- hm + new_gradient
  hm <- hm + geom_raster() #
  hm <- hm + ylab("Taxa")
  hm <- hm + theme(legend.position="bottom")
  
  hm$layers <- hm$layers[2] #
  return(hm)
}

###Plot annotated heatmap
heatmap2.plus <-
  function (x, Rowv = NULL, Colv = if (symm) "Rowv" else NULL, 
            distfun = dist, hclustfun = hclust, reorderfun = function(d, w) reorder(d, w), add.expr, symm = FALSE, revC = identical(Colv, "Rowv"), scale = c("row", "column", "none"), na.rm = TRUE, 
            margins = c(5, 5), ColSideColors, RowSideColors, cexRow = 0.2 + 1/log10(nr), cexCol = 0.2 + 1/log10(nc), labRow = NULL, 
            labCol = NULL, main = NULL, xlab = NULL, ylab = NULL, keep.dendro = FALSE, 
            verbose = getOption("verbose"), ...) 
  {
    scale <- if (symm && missing(scale)) 
      "none"
    else match.arg(scale)
    if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
      stop("'x' must be a numeric matrix")
    nr <- di[1]
    nc <- di[2]
    if (nr <= 1 || nc <= 1) 
      stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins) != 2) 
      stop("'margins' must be a numeric vector of length 2")
    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (is.null(Rowv)) 
      Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv)) 
      Colv <- colMeans(x, na.rm = na.rm)
    if (doRdend) {
      if (inherits(Rowv, "dendrogram")) 
        ddr <- Rowv
      else {
        hcr <- hclustfun(distfun(x))
        ddr <- as.dendrogram(hcr)
        if (!is.logical(Rowv) || Rowv) 
          ddr <- reorderfun(ddr, Rowv)
      }
      if (nr != length(rowInd <- order.dendrogram(ddr))) 
        stop("row dendrogram ordering gave index of wrong length")
    }
    else rowInd <- 1:nr
    if (doCdend) {
      if (inherits(Colv, "dendrogram")) 
        ddc <- Colv
      else if (identical(Colv, "Rowv")) {
        if (nr != nc) 
          stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
        ddc <- ddr
      }
      else {
        hcc <- hclustfun(distfun(if (symm) 
          x
          else t(x)))
        ddc <- as.dendrogram(hcc)
        if (!is.logical(Colv) || Colv) 
          ddc <- reorderfun(ddc, Colv)
      }
      if (nc != length(colInd <- order.dendrogram(ddc))) 
        stop("column dendrogram ordering gave index of wrong length")
    }
    else colInd <- 1:nc
    x <- x[rowInd, colInd]
    labRow <- if (is.null(labRow)) 
      if (is.null(rownames(x))) 
        (1:nr)[rowInd]
    else rownames(x)
    else labRow[rowInd]
    labCol <- if (is.null(labCol)) 
      if (is.null(colnames(x))) 
        (1:nc)[colInd]
    else colnames(x)
    else labCol[colInd]
    if (scale == "row") {
      x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
      sx <- apply(x, 1, sd, na.rm = na.rm)
      x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
      x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
      sx <- apply(x, 2, sd, na.rm = na.rm)
      x <- sweep(x, 2, sx, "/")
    }
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend) 1 else 0.05, 4)
    lhei <- c((if (doCdend) 1 else 0.05) + if (!is.null(main)) 0.2 else 0, 
              4)
    
    if (!missing(ColSideColors)) {
      if (!is.matrix(ColSideColors))
        stop("'ColSideColors' must be a matrix")
      if (!is.character(ColSideColors) || dim(ColSideColors)[1] != nc) 
        stop("'ColSideColors' dim()[2] must be of length ncol(x)")
      lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
      lhei <- c(lhei[1], 0.6, lhei[2]) ##control height of annotation from here
    }
    if (!missing(RowSideColors)) {
      if (!is.matrix(RowSideColors))
        stop("'RowSideColors' must be a matrix")
      if (!is.character(RowSideColors) || dim(RowSideColors)[1] != nr) 
        stop("'RowSideColors' must be a character vector of length nrow(x)")
      lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 
                                     1), lmat[, 2] + 1)
      lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
      cat("layout: widths = ", lwid, ", heights = ", lhei, 
          "; lmat=\n")
      print(lmat)
    }
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
    if (!missing(RowSideColors)) {
      par(mar = c(margins[1], 0, 0, 0.5))
      rsc=RowSideColors[rowInd,];
      rsc.colors=matrix();
      rsc.names=names(table(rsc));
      rsc.i=1;
      for(rsc.name in rsc.names){
        rsc.colors[rsc.i]=rsc.name;
        rsc[rsc==rsc.name]=rsc.i;
        rsc.i=rsc.i+1;
      }
      rsc=matrix(as.numeric(rsc), nrow=dim(rsc)[1]);
      image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
      
      if (length(colnames(RowSideColors))>0) {
        axis(1, 0:(dim(rsc)[2]-1) / (dim(rsc)[2]-1), colnames(RowSideColors), las=2, tick=FALSE);
      }
    }
    if (!missing(ColSideColors)) {
      par(mar = c(0.5, 0, 0, margins[2]))
      csc=ColSideColors[colInd,];
      csc.colors=matrix();
      csc.names=names(table(csc));
      csc.i=1;
      for(csc.name in csc.names){
        csc.colors[csc.i]=csc.name;
        csc[csc==csc.name]=csc.i;
        csc.i=csc.i+1;
      }
      csc=matrix(as.numeric(csc), nrow=dim(csc)[1]);
      image(csc, col = as.vector(csc.colors), axes = FALSE)
      
      if (length(colnames(ColSideColors))>0) {
        axis(2, 0:(dim(csc)[2]-1) / (dim(csc)[2]-1), colnames(ColSideColors), las=2, tick=FALSE);
      }
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    if (!symm || scale != "none") {
      x <- t(x)
    }
    if (revC) {
      iy <- nr:1
      ddr <- rev(ddr)
      x <- x[, iy]
    }
    else iy <- 1:nr
    image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
            c(0, nr), axes = FALSE, xlab = "", ylab = "", ...)
    axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
         cex.axis = cexCol)
    if (!is.null(xlab)) 
      mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
         cex.axis = cexRow)
    if (!is.null(ylab)) 
      mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr)) 
      eval(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 0))
    if (doRdend) 
      plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    else frame()
    par(mar = c(0, 0, if (!is.null(main)) 1 else 0, margins[2]))
    if (doCdend) 
      plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    else if (!is.null(main)) 
      frame()
    if (!is.null(main)) 
      title(main, cex.main = 1.5 * op[["cex.main"]])
    invisible(list(rowInd = rowInd, colInd = colInd, Rowv = if (keep.dendro && 
                                                                doRdend) ddr, Colv = if (keep.dendro && doCdend) ddc))
  }

###GGplot heatmap
library(ggplot2) 
library(reshape2)
library(scales)
library(plyr)

gg_heatmap <- function(mat){

  mat.m <- melt(dt.ps1.cyto.top)
  mat.m <- ddply(mat.m, .(Var2), transform, rescale = scale(value))
  #mat.m <- ddply(mat.m, .(Var2), transform,  rescale = scale(value))
  
  p <- ggplot(mat.m, aes(Var1, Var2, fill = value)) + geom_tile() + 
      scale_fill_gradient(low = "black", high = "red") + 
    #scale_x_discrete(expand = c(0, 0)) +
    #scale_y_discrete(expand = c(0, 0)) +
      theme(legend.position = "bottom",  axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(plot.margin=unit(c(0,0,0,0),"lines"),
          axis.ticks.margin = unit(0,"null")) +
      guides(fill=F)
  
  #p = ggplot(mat.m, aes(x = Var1, y = Var2, fill = rescale)) + 
  #  geom_raster()
  
  #base_size <- 9
  #p + theme_grey(base_size = base_size) + labs(x = "",  y = "") +  
    #scale_x_discrete(expand = c(0, 0)) +
    #scale_y_discrete(expand = c(0, 0)) + 
    #theme(legend.position = "none",
    #Saxis.ticks = theme_blank(), axis.text.x = theme_text(size = base_size * 0.8, angle = 330, hjust = 0, colour = "grey50"))
  return(p)
}


##Combine heatmap with column annotations
mush <- function(hmap, hcbs) {
  cbgs <- lapply(hcbs, ggplotGrob) #list of other grobs generated
  hmg <- ggplotGrob(hmap)
  # Make sure All plots have the same width in our final output
  cbWidths <- lapply(cbgs, function(x) x$widths[1:4])
  maxWidth <- do.call(unit.pmax, cbWidths)
  maxWidth <- unit.pmax(hmg$widths[1:4], maxWidth)
  
  # For visibility, set to the maximum width
  hmg$widths[1:4] <- as.list(maxWidth)
  for(i in seq_along(cbgs)) {
    cbgs[[i]]$widths[1:5] <- as.list(unit.c(maxWidth, hmg$widths[5]+hmg$widths[5]))
  }
  heights <- unit.c(unit(rep(1,length(cbgs)), "lines"), unit(1, "null"))
  rval <- do.call(arrangeGrob, args = c(cbgs, list(hmg), ncol=1, heights=list(heights
  )))
  return(rval)
}


draw_phyloseq_ggtree <- function(phyloseq) {	
  tree <- phyloseq@phy_tree
  p <- ggtree(tree, ladderize = F)
  p <- p + geom_text(subset=.(!isTip), aes(label=label), hjust=-.2, size=4)
  
  dd <- psmelt(phyloseq)
  dd <- dd[dd$Abundance > 0, ]
  data <- merge(p$data, dd, by.x="label", by.y="OTU")
  
  spacing <- 0.02
  idx <- with(data, sapply(table(node)[unique(node)], function(i) 1:i)) %>% unlist
  hjust <- spacing * idx * max(data$x)
  data$xdodge <- data$x + hjust
  
  p <- p + geom_point(data=data, aes(x=xdodge, color=SampleType,
                                     shape=Family, size=Abundance), na.rm=T) + 
    theme(legend.position="right") + scale_size_continuous(trans=log_trans(5))
  
  d2 <- unique(data[, c("x", "y", "Genus")])
  p + geom_text(data=d2, aes(label=Genus), hjust=-.3, na.rm=T, size=4)
}

ggtree.2 <- function(physeq){
  if (is.null(phy_tree(physeq, FALSE))) {
    stop("Object \"physeq\" does not have a tree slot")
  }
  else {
    tree <- phy_tree(physeq)
  }
  if (!is.rooted(tree)) {
    stop("Tree must be rooted, consider using midpoint rooting")
  }
  if (freq) {
    physeq <- transform_sample_counts(physeq, function(x) x/sum(x))
  }
  if (!is.null(group)) {
    physeq <- merge_samples(physeq, group)
  }
  if (!is.null(edge.group)) {
    x <- color_edges(physeq, edge.group, edge.method)
    edge.group <- factor(names(x$palette)[match(x$edge, x$palette)])
    tip.group <- factor(names(x$palette)[match(x$tip, x$palette)])
  }
  fattened.edges <- fattenEdges(physeq, method, width.lim, 
                                base)$edge.raw.width
  mean.fattened.edges <- colMeans(fattened.edges)
  edge.names <- colnames(fattened.edges)
  fattened.edges <- data.frame(edge.name = edge.names, t(fattened.edges))
  mean.fattened.edges <- data.frame(edge.name = edge.names, 
                                    mean.width = mean.fattened.edges)
  fattened.edges <- melt(fattened.edges, id.vars = "edge.name")
  colnames(fattened.edges)[2:3] <- c("Sample", "width")
  fattened.edges <- merge(fattened.edges, mean.fattened.edges)
  plot.phylo(phy_tree(physeq), plot = FALSE, ...)
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  ggdata <- data.frame(edge.name = edge.names, edge.x1 = lastPP$xx[lastPP$edge[, 1]], edge.x2 = lastPP$xx[lastPP$edge[, 2]], edge.y1 = lastPP$yy[lastPP$edge[, 1]], edge.y2 = lastPP$yy[lastPP$edge[, 2]])
  ggdata$edge <- lastPP$edge
  if (!is.null(edge.group)) {
    edge.group <- data.frame(edge.name = edge.names, edge.group = edge.group)
    ggdata <- merge(ggdata, edge.group)
  }
  ggdata <- merge(ggdata, fattened.edges)
  if (lastPP == "cladogram") {
    p <- ggplot(ggdata) + geom_segment(aes(x = edge.x1, y = edge.y1, 
                                           xend = edge.x2, yend = edge.y2))
  }
  if (lastPP == "phylogram") {
    p <- ggplot(ggdata) + geom_segment(aes(x = edge.x1, y = edge.y1, 
                                           xend = edge.x1, yend = edge.y2))
    p <- p + geom_segment(aes(x = edge.x1, y = edge.y2, xend = edge.x2, 
                              yend = edge.y2))
  }
  if (lastPP == "radial") {
  }
  if (lastPP == "") {
  }
}


# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
