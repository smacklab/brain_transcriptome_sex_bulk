#!/usr/bin/env Rscript

# download metadata and kallisto summarized to the gene level

txi.gene = readRDS('checkpoints/kallisto_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

txi.gene$counts = txi.gene$counts[,which(colnames(txi.gene$counts) %in% meta$LID)]
txi.gene$abundance = txi.gene$abundance[,which(colnames(txi.gene$abundance) %in% meta$LID)]
txi.gene$length = txi.gene$length[,which(colnames(txi.gene$length) %in% meta$LID)]

# remove low count genes (TPM)
# genes with mean >= 10 TPM (transcripts per million) 
# for males OR females within each region

samples.by.region.and.sex = split(meta$LID, list(meta$Region,meta$sex), drop = TRUE)

tpm.cutoff = 10

keep.genes = lapply(samples.by.region.and.sex,function(x) {
  names(which(rowMeans(txi.gene$abundance[,x]) >= tpm.cutoff))
})

regions = levels(as.factor(meta$Region))
keep.new = list()
keep.genes2 = unlist(keep.genes)
for (i in 1:length(regions)){
  now = keep.genes2[grepl(regions[i], names(keep.genes2))]
  keep.new[[i]] = unique(now)
}
names(keep.new) = regions

keep = Reduce(union,keep.genes)

saveRDS(keep.new, 'checkpoints/keep_genes.rds')

# normalize

library(limma)
library(edgeR)

counts = txi.gene$counts
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
head(d0$samples)
d0 = d0[keep,]

v=voom(d0)

e = v$E[keep,rownames(meta)]

saveRDS(e, 'checkpoints/filtered_expression_matrix.rds')

