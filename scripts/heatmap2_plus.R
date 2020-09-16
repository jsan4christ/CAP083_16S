
bvscat_colors <- unlist(lapply(hcbdf$bvscat , function(x){
  if(x == "Negative") 'white'
  else if(x == "Intermediate") 'pink'
  else if(x == "BV") 'purple'
}))


Chl_colors <- unlist(lapply(hcbdf$Chlamydia, function(x){
  if(x == 0) 'grey'
  else if(x == 1) 'red'
}))
Gon_colors <- unlist(lapply(hcbdf$Gonorrhoea, function(x){
  if(x == 0) 'grey'
  else if(x == 1) 'red'
}))
Tri_colors <- unlist(lapply(hcbdf$Trichomoniasis, function(x){
  if(x == 0) 'grey'
  else if(x == 1) 'red'
}))
hsv2_colors <- unlist(lapply(hcbdf$HSV.2, function(x){
  if(x == 0) 'grey'
  else if(x == 1) 'red'
}))
inf_colors <- unlist(lapply(hcbdf$Inflammation, function(x){
  if(x == 0) 'grey'
  else if(x == 1) 'red'
}))

sti_colors <- unlist(lapply(hcbdf$Any_STI, function(x){
  if(x == 0) 'grey'
  else if(x == 1) 'red'
}))

myCols <- cbind(inf_colors, Chl_colors, Gon_colors, Tri_colors, hsv2_colors,  sti_colors, bvscat_colors)
colnames(myCols) <- c("Inflammation", "Chlamydia", "Gonorrhoea", "Trichomoniasis", "HSV.2", "Any_STI", "Nuggent Score")

input <- as.matrix(t(dt.ps1.cyto.top))
heatmap.2(input, trace="none", density="none", col=bluered(20), cexRow=1, cexCol=0.2, margins = c(20,13),
          ColSideColors=myCols, na.rm = T)
pdf(file="AnnotatedHeatmap.pdf")
par(cex.main=0.8,mar=c(4,4,4,8), cin=c(1,))

heatmap2.plus(t(cytokine.mat), col=bluered(20), Rowv=NULL, Colv=NULL, cexRow=1,cexCol=0.4, margins = c(20,13),
              ColSideColors=myCols)

dev.off()



#####aheatmap###
hc <- hclust(dist(x, method = 'minkowski'), method = 'centroid')

pdf(file = "aheatmap.pdf")
aheatmap(dt.ps1.cyto.top.tr, annCol = hm_annot_all)
dev.off()

