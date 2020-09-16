###Custom functions 
## Author: San James

## Plot prevalence of bacteria/taxa accross samples
PlotPrevalence <- function(physeq, yintercept = 0.05){
  #Compute prevalence values
  prevalence.values <- apply(X = otu_table(physeq), MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
  cat("Range of prevalence is: ", range(prevalence.values))
  
  #Now add taxa table to the dataframe
  df.physeq <- data.frame(Prevalence = prevalence.values, TotalAbundance = taxa_sums(physeq), tax_table(physeq))
  
  df.physeq.phylum <- subset(df.physeq, Phylum %in% get_taxa_unique(physeq, "Phylum"))
  
  ##Plot prevalence at Phylum level
  gg.physeq.phylum.prevalence <- ggplot(df.physeq.phylum, aes(TotalAbundance, Prevalence / nsamples(physeq), color=Family)) +
    geom_hline(yintercept = yintercept, alpha = 0.5, linetype = 2) +
    geom_point(size = 3, alpha = 0.7) +
    scale_x_log10() +
    xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
    facet_wrap(~Phylum) +
    theme(legend.position="none") +
    ggtitle("Phylum Prevalence in All Samples\nColored by Family")
  
  ggplotly(gg.physeq.phylum.prevalence)
}

## Filter out low abundance taxa by Phylum prevalence
FilterByOTULevel <- function(physeq, TaxaToDrop, taxa.rank = "Pylum"){
  cat("Dropping the following features: ", TaxaToDrop, "at OTU Level: ", taxa.rank,  "\n")
  physeq.prevalence <- subset_taxa(physeq, !taxa.rank %in% TaxaToDrop)
  return(physeq.prevalence)
}

## Filter out low abundance taxa by prevalence levels
FilterByPrevalence <- function(physeq, prevalence = 0.05){
  #Thresh hold to use
  prevalenceThreshold <- prevalence *  nsamples(physeq)
  
  # Calculate feature prevalence across the data set. i.e all the samples in which an ASV is found
  df.physeq <- apply(X = otu_table(physeq), MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2), FUN = function(x){sum(x > 0)})
  
  # Add taxonomy and total read counts (all sequences from a sample) to df.prevalence
  df.physeq <- data.frame(Prevalence = df.physeq, TotalAbundance = taxa_sums(physeq), tax_table(physeq))
  
  # Define which taxa fall within the prevalence threshold
  keepTaxa <- rownames(df.physeq)[(df.physeq$Prevalence >= prevalenceThreshold)]
  
  #Drop taxa that does not meet the prevalence threshold
  filtered.physeq <- prune_taxa(keepTaxa, physeq)
  
  #return filtered physeq
  return(filtered.physeq)
}

## Perform adonis test on a physeq object and a variable from metadata 
doadonis <- function(physeq, category, strata = NULL) {
  bdist <- phyloseq::distance(physeq, "bray")
  col <- as(sample_data(physeq), "data.frame")[ ,category]
  
  # Adonis test
  adonis.bdist <- adonis(bdist ~ col, strata = strata)
  print("Adonis results:")
  print(adonis.bdist)
  
}

## Perform DiffAbb on physeq object
DiffAbundanceUsingDeseq <- function(ps, full.design, design = ~1, test = "LRT", FDR = 0.1, fitType="local") {
  require(DESeq2)
  
  # Convert phyloseq object to deseq
  ds.ps <- phyloseq_to_deseq2(ps, design = full.design)
  
  # Calculate geometric means prior to estimate size factors
  geoMeans = apply(counts(ds.ps), 1, gm_mean)
  ds.ps = estimateSizeFactors(ds.ps, geoMeans = geoMeans)
  
  # Run deseq
  dds.ps <- DESeq(ds.ps, test = test, fitType = fitType, reduced = design)
  
  # Process results
  res <- results(dds.ps)
  res <- res[order(res$padj, na.last = T), ]
  res$prank <- seq(1,nrow(res))
  res$OTU <- rownames(res)
  res <- as.data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[res$padj < FDR,]
  res
}

## Transform to relative abundance.
TransformToRelativeAbundance <- function(ps){
  require(phyloseq)
  # Call phylosecs transform sample counts function
  return(transform_sample_counts(ps, function(x) round(100 * x/sum(x))))
}


## Compute geometric mean
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


## Convert phyloseq object to edgeR
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}


## Function to cluster samples using mclust

generate_mclust_csts <- function(physeq, nClusters = 5, distMethod = "bray"){
  require(mclust)
  distObj <- phyloseq::distance(physeq, method = distMethod)
  #Mclust Model
  BIC <- mclustBIC(scale(distObj), G = nClusters)
  
  #plot(BIC) #plot bic to see how model was arrived at
  mclustObj <- Mclust(scale(distObj), x = BIC) # Model-based-clustering
  summary(mclustObj)                 # Print a summary
  
  mclustObj$modelName                # Optimal selected model ==> "VVV"
  mclustObj$G                        # Optimal number of cluster => 3
  #head(mclustObj$z, 30)              # Probality to belong to a given cluster
  #head(mclustObj$classification, 30) # Cluster assignement of each observation
  #class(mclustObj$classification)
  
  mclust.sample.cluster <- mclustObj$classification
  
  sample_data(physeq)$CSTs <- factor(mclust.sample.cluster, levels = unique(sort(mclust.sample.cluster)), labels = paste0("CST", seq(length(unique(sort(mclust.sample.cluster))))))
  
  return(physeq)
}

# Draw heatmap annotated with CSTs and STIs
cst_annotated_heatmap <- function(physeq){
  require(ggplotify)
  #physeq <- ps2.v2.cst.glom.Species.ra
  top25 <- names(sort(taxa_sums(physeq), decreasing=T))[1:16]
  physeq.top25 <- prune_taxa(top25, physeq)
  
  # Sorting vars
  taxa.order <- names(sort(taxa_sums(physeq.top25)))
  sample.order <- rownames(sample_data(physeq.top25)[order(get_variable(physeq.top25, "CSTs"))])
  
  hm <- plot_heatmap.2(physeq.top25, taxa.label="Species", sample.order=sample.order, taxa.order=taxa.order)
  hm <- hm + theme(axis.title.x = element_text(size=10),
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
  
  CSTColors <- brewer.pal(6,"Paired")[c(1,2,3,4)] # Length 5 for consistency with pre-revision CT+ coloration
  CSTs <- as.character(c("CST1", "CST2", "CST3", "CST4"))
  names(CSTColors) <- CSTs
  CSTColorScale <- scale_colour_manual(name = "CSTs", values = CSTColors[1:4])
  CSTFillScale <- scale_fill_manual(name = "CSTs", values = CSTColors[1:4])
  
  ### CHANGING SPECIES TO TAXA ON YLABEL
  df.physeq <- data.frame(sample_data(physeq.top25))
  df.physeq <- df.physeq[sample.order,] #resolve
  df.physeq$index <- seq(1, nsamples(physeq.top25))
  
  hcb <- make_hcb(df.physeq, "CSTs", name = "CSTs", fillScale = CSTFillScale)
  hcb <- hcb + annotate("text", x=tapply(df.physeq$index, df.physeq[, "CSTs", drop=T], mean), y=1,
                        label=levels(df.physeq[, "CSTs", drop=T]), size=2)
  
  Chl <- make_hcb(df.physeq, "Chlamydia", name="Chlamydia", 
                  fillScale = scale_fill_manual(values=c("0"="papayawhip", "1"="tan4")))
  Gon <- make_hcb(df.physeq, "Gonorrhoea", name="Gonorrhoea", 
                  fillScale = scale_fill_manual(values=c("0"="papayawhip", "1"="green4")))
  Tri <- make_hcb(df.physeq, "Trichomoniasis", name="Trichomoniasis", 
                  fillScale = scale_fill_manual(values=c("0"="papayawhip", "1"="yellow4")))
  Cand <- make_hcb(df.physeq, "Candidiasis", name="Candidiasis", 
                   fillScale = scale_fill_manual(values=c("0"="papayawhip", "1"="tomato4")))
  inf <- make_hcb(df.physeq, "Inflammation", name="Inflammation", 
                  fillScale = scale_fill_manual(values=c("0"="papayawhip", "1"="darkred")))
  bvscores <- make_hcb(df.physeq, "bvscat", name="Nuggent Score", 
                       fillScale = scale_fill_manual(values=c("Negative"="white", "Intermediate"="grey60", "BV"="maroon")))
  
  inf <- inf + theme(axis.text.y = element_text(size=8, face="bold", color="grey60"))
  Chl <- Chl + theme(axis.text.y = element_text(size=8, face="bold", color="tan4"))    
  Gon <- Gon + theme(axis.text.y = element_text(size=8, face="bold", color="green4"))
  Tri  <- Tri + theme(axis.text.y = element_text(size=8, face="bold", color="yellow4")) 
  Cand <- Cand + theme(axis.text.y = element_text(size=8, face="bold", color="tomato4"))
  
  Fig <- mush(hm, list(Chl, Gon, Tri, Cand, inf, bvscores, hcb))
  grid.newpage()
  grid.draw(Fig)
  #as.ggplot(Fig)
}



##DEseq
getDEsigs <- function(ps, design, FDR = 0.1, fitType="local") {
  require(DESeq2)
  pDE <- phyloseq_to_deseq2(ps, design = design)
  pDE <- DESeq(pDE, fitType=fitType)
  res <- results(pDE)
  res <- res[order(res$padj, na.last = T), ]
  res$prank <- seq(1,nrow(res))
  res$OTU <- rownames(res)
  res <- as.data.frame(res)
  res <- res[!is.na(res$padj),]
  res <- res[res$padj < FDR,]
  res
}
  
  ReorderCTs <- function(physeq, CTList){
    sample_data(physeq)$CSTs <- factor(sample_data(physeq)$CSTs, levels = CTList[[1]], labels = CTList[[2]])
    levels(sample_data(physeq)$CSTs)
    return(physeq)
  }
  
  
  cestimate_richness <- function(physeq, measures = c("Observed", "Shannon", "simpson")){
    div <- estimate_richness(physeq, measures = measures)
    div$ShannonLgT <- log10(div$Shannon + 1)
    sample_data(physeq) <- cbind(sample_data(physeq), div)
    return(physeq)
  }
  
  
  merge_low_abundance <- function(pobject, threshold=0.05){
    transformed <- transform_sample_counts(pobject, function(x) x/sum(x))
    otu.table <- as.data.frame(otu_table(transformed))
    #otu.list <- row.names(otu.table[rowMeans(otu.table) < threshold,])
    otu.list <- colnames(otu.table[colMeans(otu.table) < 0.05,])
    merged <- merge_taxa(transformed, otu.list, 1)
    for (i in 1:dim(tax_table(merged))[1]){
      if (is.na(tax_table(merged)[i,2])){
        taxa_names(merged)[i] <- "Other"
        tax_table(merged)[i,1:7] <- "Other"}
    }
    return(merged)
  }
  
  
  merge_less_than_top <- function(pobject, top=20){
    pobject <- tax_glom(pobject, taxrank = "Species")
    transformed <- transform_sample_counts(pobject, function(x) x/sum(x))
    otu.table <- as.data.frame(t(otu_table(transformed)))
    otu.sort <- otu.table[order(rowMeans(otu.table), decreasing = TRUE),]
    otu.list <- row.names(otu.sort[(top+1):nrow(otu.sort),])
    merged <- merge_taxa(transformed, otu.list, 1)
    for (i in 1:dim(tax_table(merged))[1]){
      if (is.na(tax_table(merged)[i,2])){
        taxa_names(merged)[i] <- "Other"
        tax_table(merged)[i,1:7] <- "Other"}
    }
    return(merged)
  }