#!/usr/bin/env Rscript

############
# load data 
############

e.keep = checkpoints/readRDS('filtered_expression_matrix.rds')
keep.genes = checkpoints/readRDS('keep_genes.rds')
meta = checkpoints/readRDS('cayo_bulkbrain_combined_metadata.rds')
meta = meta[order(match(meta$LID, colnames(e.keep))),]

library(BRETIGEA)
library(biomaRt)
library(expss)
library(reshape2)
library(ggplot2)
library(edgeR)

'%!in%' <- function(x,y)!('%in%'(x,y))

## download top 1000 genes for each of the six cell types (ranked by specificity)

markers = markers_df_brain

## convert human gene names to macaque ensembl ids

hsap = useEnsembl(biomart = 'ensembl',dataset='hsapiens_gene_ensembl',mirror='www') 
gene_conv = getBM(attributes=c('ensembl_gene_id','external_gene_name','mmulatta_homolog_ensembl_gene','mmulatta_homolog_orthology_type'), mart = hsap)
gene_conv = subset(gene_conv, mmulatta_homolog_orthology_type == 'ortholog_one2one')

for (i in 1:length(markers$markers)){
  markers$id[i] = vlookup(markers$markers[i], dict = gene_conv, result_column = 3,lookup_column = 2)}
markers = markers[which(!is.na(markers$id)),]
colnames(markers) = c("name","cell","markers")
markers = markers[markers$markers %in% row.names(e.keep),]
markers = markers[markers$markers %!in% markers$markers[duplicated(markers$markers)],]
table(markers$cell)
markers = markers[,c(3,2)]

################
## estimate SPVs
################

regions = names(keep.genes)
out = list()
adj = list()

for (i in 1:length(regions)){
  m = subset(meta, Region == regions[i])
  nobatch = removeBatchEffect(e.keep[,colnames(e.keep) %in% m$LID], batch = m$Library.batch, covariates = cbind(m$RIN))
  out[[i]] = findCells(inputMat=nobatch, markers=markers, nMarker = 50, method = "SVD", scale = TRUE)
  adj[[i]] = adjustCells(inputMat = nobatch, cellSPV = out[[i]])
  }
names(out) = regions

adj_exp = do.call(cbind, adj)  
saveRDS(adj_exp, 'checkpoints/adjusted_filtered_expression_matrix.rds')

###########
# plot SPVs
###########

out2 = melt(out)
out2$Sex = vlookup(out2$Var1, meta, lookup_column = 2, result_column = 15)
out2$lib = vlookup(out2$Var1, meta, lookup_column = 2, result_column = 13)
out2$L1 = vlookup(out2$Var1, meta, lookup_column = 2, result_column = 3)

levels(out2$Var2) = c("Ast","End","Mic","Neu","Oli","OPC")
levels(out2$Sex) = c("F","M")
out2$L1 = factor(out2$L1, levels = region.levels)

ggplot(out2, aes(x=value, y=Sex, group=Sex)) + 
  geom_vline(aes(xintercept = 0), colour="grey") +
  geom_violin(aes(fill=Sex)) +
  scale_fill_manual(values = region.colors[c(3,6)]) +
  theme_classic() +
  facet_grid(L1 ~ Var2, switch = "y") +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill=NA,colour=NA),
        panel.spacing=unit(0,"cm"), axis.title.y = element_blank()) +
  coord_cartesian(clip="off") +
  theme(panel.margin = unit(.3, "lines"), legend.title = element_text(size=18),
        plot.title = element_text(size=18),
        legend.text = element_text(size=18),strip.text.x.top = element_text(size=18), 
        strip.text.y.left = element_text(angle = 0, size=18), 
        axis.text.y = element_blank(), axis.text.x= element_blank(), 
        axis.title.x=element_blank(), axis.ticks = element_blank()) +
  scale_x_continuous(breaks=c(0)) +
  ggtitle("Surrogate Proportion Variables")

# estimate sex differences in SPVs

sex_diffs = data.frame()
cell_types = levels(out2$Var2)
for (i in 1:length(cell_types)){
  datanow = subset(out2, Var2 == cell_types[i])
  mod = summary(lm(value ~ Sex, data = datanow))
  sex_diffs[i,1] = mod$coefficients[2,4]
}

rownames(sex_diffs) = cell_types
sex_diffs$padj = p.adjust(sex_diffs$V1, method='bonferroni')
sex_diffs
