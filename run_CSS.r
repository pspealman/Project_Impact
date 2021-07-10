if (!requireNamespace("BiocManager", quietly = TRUE))
  
  install.packages("BiocManager")

BiocManager::install("metagenomeSeq")

library(metagenomeSeq)

OTU_read_count <- read.delim("feature-table.biom.txt")
rownames(OTU_read_count) <- OTU_read_count$OTU_ID
OTU_read_count <- OTU_read_count[-c(1)]

metaSeqObject = newMRexperiment(OTU_read_count)

metaSeqObject_CSS  = cumNorm( metaSeqObject , p=cumNormStatFast(metaSeqObject) )

OTU_read_count_CSS = data.frame(MRcounts(metaSeqObject_CSS, norm=TRUE, log=TRUE))

cor(OTU_read_count$P1, OTU_read_count_CSS$P1)
cor(OTU_read_count$M1, OTU_read_count_CSS$M1)

write.csv(OTU_read_count_CSS, "feature-table_CSS_norm.biom.txt", row.names = TRUE)