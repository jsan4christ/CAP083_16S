## Author: San Emmanuel James
## Date Last Modified: 16th September, 2020
## Purpose: Convert output from PECAN (http://pecan-lab.org/speciateit) into phyloseq compartible for downstream analyses

library(Biostrings)
library(dplyr)
library(stringr)
library(DECIPHER); packageVersion("DECIPHER")
library(phangorn); packageVersion("phangorn")
library(dada2)

## Import pecan output for processing
pecan_main <- read.csv("pecan_out/andile_all_runs_dada2_abundance_table_PECAN_taxa-merged_StR_CST.csv")
pecan_asv <- as_tibble(read.csv("pecan_out/andile_all_runs_dada2_abundance_table_PECAN_asvs+taxa.csv"))

asv_taxa <- read.table("pecan_out/cmb_tx.txt")
colnames(asv_taxa) <- c("ASV", "Species")
rownames(asv_taxa) <- asv_taxa$ASV
asv_taxa$ASV <- as.character(asv_taxa$ASV)

taxa_upto_genus <- read.csv("pecan_out/andile_silva_classification.csv")
colnames(taxa_upto_genus)[1] <- "ASV"
taxa_upto_genus$ASV <- as.character(taxa_upto_genus$ASV)


seqs_asvs <- data.frame(readDNAStringSet("pecan_out/all_runs_dada2_ASV.fasta"))
colnames(seqs_asvs) <- "sequence"
seqs_asvs$ASV <- rownames(seqs_asvs)

## add sequences
taxa_upto_genus_seqs <- inner_join(seqs_asvs, taxa_upto_genus, by = "ASV")

## Add species:
taxa_table <- inner_join(taxa_upto_genus_seqs, asv_taxa, by = "ASV")
rownames(taxa_table) <- taxa_table$sequence

## Extract taxa table
taxa_table_final <- as.matrix(taxa_table[, c(3:9)])
saveRDS(taxa_table_final, "dada2-Chimera-Taxonomy/tax_table_final.RDS")

### Extract CST assignments
metadata_cst <- pecan_main[, c(1, 377:393)]
metadata_cst$visit <- sapply(str_split(metadata_cst$sampleID, "-", n = 2, simplify = FALSE), `[`, 2)
metadata_cst[is.na(metadata_cst$visit), "visit"] <- 1000
metadata_cst$sampleID <- ifelse(metadata_cst$visit == 1000, paste0(sapply(str_split(metadata_cst$sampleID, "-", n = 2, simplify = FALSE), `[`, 1), 1),
                             ifelse(metadata_cst$visit == 1020, paste0(sapply(str_split(metadata_cst$sampleID, "-", n = 2, simplify = FALSE), `[`, 1), 2),
                                    paste0(sapply(str_split(metadata_cst$sampleID, "-", n = 2, simplify = FALSE), `[`, 1), 3)))
colnames(metadata_cst)[1] <- "SampleID"
saveRDS(metadata_cst, "metadata/metadata_cst.RDS")

## Extract OTU table
dim(pecan_asv)
colnames(pecan_asv)[1] <- "SampleID"

pecan_asv_ <- pecan_asv[c(2:nrow(pecan_asv)), ]
pecan_asv_$visit <- sapply(str_split(pecan_asv_$SampleID, "-", n = 2, simplify = FALSE), `[`, 2)
pecan_asv_[is.na(pecan_asv_$visit), "visit"] <- 1000
pecan_asv_$SampleID <- ifelse(pecan_asv_$visit == 1000, paste0(sapply(str_split(pecan_asv_$SampleID, "-", n = 2, simplify = FALSE), `[`, 1), 1),
                                ifelse(pecan_asv_$visit == 1020, paste0(sapply(str_split(pecan_asv_$SampleID, "-", n = 2, simplify = FALSE), `[`, 1), 2),
                                       paste0(sapply(str_split(pecan_asv_$SampleID, "-", n = 2, simplify = FALSE), `[`, 1), 3)))

pecan_asv_$SampleID[duplicated(pecan_asv_$SampleID)]
rownames(pecan_asv_) <- pecan_asv_$SampleID

pecan_asv_$visit <- NULL 
otu_table <- as.matrix(pecan_asv_[, c(2:ncol(pecan_asv_))])
class(otu_table) <- "numeric"
colnames(otu_table) <- rownames(taxa_table)
saveRDS(otu_table, "dada2-Chimera-Taxonomy/seqtab_pecan_final.RDS")

## Construct phylogenetic tree
seqs <- getSequences(t(otu_table))
names(seqs) <- seqs # This propagates to the tip labels of the tree
##write.table(names(seqs), "seqs.txt", sep = "\n")
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,
                       processors = detectCores() - 1) ## AlignSeqs function from package Decipher
writeXStringSet(alignment, "dada2-Phangorn/aligned_seqs.fasta")

phang.align <- phyDat(as(alignment, "matrix"), type="DNA") ## phyDat function from package phangorn
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
write.tree(fit$tree, file = "dada2-Phangorn/tree.newick")

## negative edges length changed to 0!
fitGTR <- update(fit, k=4, inv=0.2) 
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
saveRDS(fitGTR, "dada2-Phangorn/phangorn.tree.RDS")
write.tree(fitGTR$tree, file = "tree.GTR.newick")
detach("package:phangorn", unload=TRUE)

## Clean up
rm(list = ls())