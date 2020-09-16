###Differential abundance using edgeR

ps2.v1v2v3.t = transform_sample_counts(ps2.v1v2v3, function(x){x/sum(x)})
hist(log10(apply(otu_table(ps2.v1v2v3.t), 1, var)),
     xlab="log10(variance)", breaks=50,
     main="A large fraction of OTUs have very low variance")

varianceThreshold = 1e-5
keepOTUs = names(which(apply(otu_table(ps2.v1v2v3.t), 1, var) > varianceThreshold))
ps2.v1v2v3.er = filter_taxa(ps2.v1v2v3.t, function(x) mean(x) > 1e-5, TRUE)
ps2.v1v2v3.er

library()
dge = phyloseq_to_edgeR(ps2.v1v2v3.er, group="STI")
# Perform binary test
et = exactTest(dge)
# Extract values from test results
tt = topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue")
res = tt@.Data[[1]]
alpha = 0.05
sigtab = res[(res$FDR < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps2.v1v2v3.er)[rownames(sigtab), ], "matrix"))
dim(sigtab)



## Plot results

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
sigtabgen = subset(sigtab, !is.na(Genus))
# Phylum order
x = tapply(sigtabgen$logFC, sigtabgen$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels = names(x))
# Genus order
x = tapply(sigtabgen$logFC, sigtabgen$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtabgen$Genus = factor(as.character(sigtabgen$Genus), levels = names(x))
ggplot(sigtabgen, aes(x = Genus, y = logFC, color = Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5))
