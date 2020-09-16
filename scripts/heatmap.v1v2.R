source("Scripts/heatmap_function.R")

##Transform to relative abundances
ps2.v1v2.cst.glom.Species.ra <- transform_sample_counts(ps2.v1v2.cst.glom.Species, function(x) round(100 * x/sum(x)))

##Pick out the 25 most abundant
top25 <- names(sort(taxa_sums(ps2.v1v2.cst.glom.Species.ra), decreasing=T))[1:25]
ps2.v1v2.cst.glom.Species.ra.top25 <- prune_taxa(top25, ps2.v1v2.cst.glom.Species.ra)
taxa.order.v1v2 <- names(sort(taxa_sums(ps2.v1v2.cst.glom.Species.ra.top25)))


sample.order.v1v2 <- rownames(sample_data(ps2.v1v2.cst.glom.Species.ra.top25)[order(get_variable(ps2.v1v2.cst.glom.Species.ra.top25, "classif2"))])
hm.v1v2 <- plot_heatmap.2(ps2.v1v2.cst.glom.Species.ra.top25, taxa.label="Species", sample.order=sample.order.v1v2, taxa.order=taxa.order.v1v2)
hm.v1v2 <- hm.v1v2 + theme(axis.title.x = element_text(size=10),
                 axis.title.y = element_text(size=10),
                 axis.text.x = element_text(size=7),
                 axis.text.y = element_text(size=7),
                 plot.title = element_text(size=8),
                 legend.text = element_text(size=7),
                 legend.title = element_text(size=8),
                 # legend.margin = unit(c(0.1,0.1,0.1,0.1),"mm"),
                 # legend.key.height = unit(1, "in"),
                 legend.key.width = unit(0.15, "in"),
                 plot.margin=unit(c(0,0,0,0),"mm"))

CSTColors <- brewer.pal(6,"Paired")[c(1,2,3,4)] # Length 5 for consistency with pre-revision CT+ coloration. RColorBrewer
names(CSTColors) <- CTs
CSTColorScale <- scale_colour_manual(name = "CT", values = CSTColors[1:4])
CSTFillScale <- scale_fill_manual(name = "CT", values = CSTColors[1:4])

### CHANGING SPECIES TO TAXA ON YLABEL
labvec <- as(tax_table(ps2.v1v2.cst.glom.Species.ra.top25)[, "Species"], "character")
names(labvec) <- taxa_names(ps2.v1v2.cst.glom.Species.ra.top25)
labvec <- labvec[taxa.order.v1v2]
labvec[is.na(labvec)] <- ""

hm.v1v2 <- hm.v1v2 + scale_y_discrete("Taxa", labels = labvec)
hm.v1v2 <- hm.v1v2 + theme(axis.title = element_text(size=10))

hcbdf.v1v2 <- data.frame(sample_data(ps2.v1v2.cst.glom.Species.ra.top25))
hcbdf.v1v2 <- hcbdf[sample.order,] #resolve
hcbdf.v1v2$index <- seq(1,nsamples(ps2.v1v2.cst.glom.Species.ra.top25))
hcb.v1v2 <- make_hcb(hcbdf.v1v2, "classif2", name="CT", fillScale = CSTFillScale)
hcb.v1v2 <- hcb.v1v2 + annotate("text", x=tapply(hcbdf.v1v2$index, hcbdf.v1v2[,"classif2",drop=T], mean), y=1,
                      label=levels(hcbdf.v1v2[,"classif2",drop=T]), size=2)
bvscores.v1v2 <- make_hcb(hcbdf.v1v2, "bvscat", name="Nuggent Score", 
                     fillScale = scale_fill_manual(values=c("Negative"="white", "Intermediate"="maroon", "BV"="grey60")))
bvscores.v1v2 <- bvscores.v1v2 + theme(axis.text.y = element_text(size=8, face="bold", color="grey60"))
Fig2.v1v2 <- mush(hm.v1v2, list(bvscores, hcb)) #arrangeGrob from gridExtra
grid.newpage() #unit.pmax from grid
grid.draw(Fig2.v1v2)


####Create transition matrix###

dt.ps2.v1v2.net <- sample_data(ps2.v1v2.cst.glom.Species.ra)
dt.ps2.v1.net <- subset(dt.ps2.v1v2.net, VisitCode == 1000, select = c(SampleID, ParticipantID, VisitCode, classif2))
colnames(dt.ps2.v1.net)[4] <- "PrevCST"
dt.ps2.v2.net <- subset(dt.ps2.v1v2.net, VisitCode == 1020, select = c(SampleID, ParticipantID, VisitCode, classif2))
colnames(dt.ps2.v2.net)[4] <- "CurrCST"
dt.ps2.v1v2.netM <- cbind(dt.ps2.v2.net, dt.ps2.v1.net[, c("SampleID", "PrevCST")] )
View(dt.ps2.v1v2.netM)

ttab <- table(dt.ps2.v1v2.netM$CurrCST, dt.ps2.v1v2.netM$PrevCST)
transMat <- matrix(ttab, nrow = 4)
transMat <- transMat/rowSums(transMat)
CSTrans <- transMat
colnames(CSTrans) <- CSTs
rownames(CSTrans) <- CSTs
CST_Persist <- -1/log(diag(CSTrans))
CST_Persist

# Make Markov chain object
mcCSTTrans <- new("markovchain", states=CSTs,
              transitionMatrix = transMat, name="CSTTrans")
mcCSTTrans

# Set up igraph of the markov chain
netMC <- markovchain:::.getNet(mcCSTTrans, round = TRUE)


###Define plotting parameters and assign node colors
wts <- E(netMC)$weight/100

edgel <- get.edgelist(netMC)
elcat <- paste(edgel[,1], edgel[,2])
elrev <- paste(edgel[,2], edgel[,1])
edge.curved <- sapply(elcat, function(x) x %in% elrev)

bvscat <- table(dt.ps2.v1v2.net$classif2, dt.ps2.v1v2.net$bvscat)
rownames(bvscat) <- markovchain::states(mcCSTTrans)
colnames(bvscat) <- c("BV", "Intermediate", "Negative")
bvscat
bvscat <- bvscat/rowSums(bvscat)
vert.CSTclrs <- CSTColors


default.par <- par(no.readonly = TRUE)
# Define color scale
# Plotting function for markov chain
plotMC <- function(object, ...) {
  netMC <- markovchain:::.getNet(object, round = TRUE)
  plot.igraph(x = netMC, ...)  
}
# Color bar for the markov chain visualization, gradient in strength of preterm association
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title=NULL) {
  scale = (length(lut)-1)/(max-min)
  
  #    dev.new(width=1.75, height=5)
  
  cur.par <- par(no.readonly=T)
  par(mar=c(0,4,1,4)+0.1, oma=c(0,0,0,0)+0.1)
  par(ps = 10, cex = 0.8)
  par(tcl=-0.2, cex.axis=0.8, cex.lab = 0.8)
  plot(c(min,max), c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(1, c(0, 0.5, 1))
  for (i in 1:(length(lut)-1)) {
    x = (i-1)/scale + min
    rect(x,0,x+1/scale,10, col=lut[i], border=NA)
  }
}


pal <- colorRampPalette(c("grey50", "maroon", "magenta2"))(101)
vert.colors <- sapply(states(mcCSTTrans), function(x) pal[1+round(100*bvscat[x,"BV"])])
vert.size <- 4 + 2*sapply(states(mcCSTTrans), 
                        function(x) nrow(unique(sample_data(ps2.v1v2.cst.glom.Species)[sample_data(ps2.v1v2.cst.glom.Species)$classif2==x,"ParticipantID"])))
vert.sz <- vert.sz * 0.85
vert.font.clrs <- c("white", "white", "white", "white", "white")
# E(netMC) to see edge list, have to define loop angles individually by the # in edge list, not vertex
edge.loop.angle = c(0, 0, 0, 0, 3.14, 3.14, 0, 0, 0, 0, 3.14, 0, 0, 0, 0, 0)-0.45

layout <- matrix(c(0.6,0.95, 0.43,1, 0.3,0.66, 0.55,0.3, 0.75,0.65), nrow=5, ncol=2, byrow=T)

# Colored by association with BV
layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE), heights=c(1,10))
color.bar(pal, min=0, max=1, nticks=6, title="Fraction BV")
par(mar=c(0,1,1,1)+0.1)
edge.arrow.size=0.8    
edge.arrow.width=1.4
edge.width = (15*wts + 0.1)*0.6
edge.labels <- as.character(E(netMC)$weight/100)
edge.labels[edge.labels<0.4] <- NA  # labels only for self-loops


plotMC(mcCSTTrans, edge.arrow.size=edge.arrow.size, edge.arrow.width = edge.arrow.width,
       #       edge.label = edge.labels, edge.label.font=2, edge.label.cex=1.3, edge.label.color="black",
       # FIX EDGE LABELS FOR PUBLICATION IN POST-PROCESSING
       edge.width=edge.width, edge.curved=edge.curved, 
       vertex.color=vert.colors, vertex.size=(vert.size),
       vertex.label.font = 2, vertex.label.cex = 1,
       vertex.label.color = vert.font.clrs, vertex.frame.color = NA, 
       layout=layout, edge.loop.angle = edge.loop.angle)
par(default.par)
