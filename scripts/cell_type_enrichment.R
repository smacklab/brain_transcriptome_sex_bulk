#!/usr/bin/env Rscript

library(mashr)

#######################
# load data and results
#######################

e.keep = checkpoints/readRDS('filtered_expression_matrix.rds')
keep.genes = checkpoints/readRDS('keep_genes.rds')
meta = checkpoints/readRDS('cayo_bulkbrain_combined_metadata.rds')
meta = meta[order(match(meta$LID, colnames(e.keep))),]

# primary analyses
mash.results = checkpoints/readRDS('mashr_results.rds')

# human gtex data
mash.hsap = data/readRDS('gtex_mashr_results_sex.rds')

#################################
# extract coefficients and p vals (and select data set)
#################################

# primary analyses
mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)
mash.beta = mash.beta[,region.levels]
mash.lfsr = mash.lfsr[,region.levels]

# human gtex data
mash.beta = get_pm(mash.hsap)
mash.lfsr = get_lfsr(mash.hsap)

#################
## get gene lists (macaque analyses)
#################

ensembl.gene.names = unique(unlist(keep.genes))
mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)
all.region = numeric(length=length(mashr.genes))
names(all.region) = mashr.genes

all.region[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff) >= 1/15 * length(keep.genes)) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] > 0) >= 1/15 * length(keep.genes)
}))))] = 1
all.region[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff) >= 1/15 * length(keep.genes)) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] < 0) >= 1/15 * length(keep.genes)
}))))] = -1
table(all.region)

mb = names(all.region[which(all.region == 1)])
fb = names(all.region[which(all.region == -1)])
sb = c(mb, fb)
all = names(all.region)

#################
## get gene lists (human analyses)
#################

rownames(mash.beta) = str_sub(rownames(hsap.beta),1,15)
rownames(mash.lfsr) = str_sub(rownames(hsap.lfsr),1,15)

ensembl.gene.names = rownames(mash.beta)
mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)
all.region.fet = numeric(length=length(mashr.genes))
names(all.region.fet) = mashr.genes

all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff) >= 1/10 * length(colnames(mash.beta))) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] > 0) >= 1/10 * length(colnames(mash.beta))
}))))] = 1
all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff) >= 1/10 * length(colnames(mash.beta))) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] < 0) >= 1/10 * length(colnames(mash.beta))
}))))] = -1
table(all.region.fet)

mb = names(all.region.fet[which(all.region.fet == 1)])
fb = names(all.region.fet[which(all.region.fet == -1)])
sb = c(mb, fb)
all = names(all.region.fet)

######################################
## cell type enrichment (macaques)
######################################

library(BRETIGEA)
library(biomaRt)
library(expss)

markers = markers_df_brain

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
colnames(markers) = c('macID', 'cellName')

ct_markers_all = markers
head(ct_markers_all)
table(ct_markers_all$cellName)
cells = levels(as.factor(ct_markers_all$cellName))

celltypes = data.frame()

for (i in 1:length(cells)){
  cellnow = cells[i]
  
  sexnow = fb
  markers = unique(subset(ct_markers_all, cellName == cellnow)$macID)
  combo.cell = length(which(markers %in% sexnow))
  cell.all = markers[which(markers %in% all)]
  cell.only = length(which(cell.all %!in% sexnow))
  sex.all = length(sexnow)
  sex.only = sex.all- combo.cell
  neither = length(all) - combo.cell - cell.only - sex.only
  mat = matrix(c(combo.cell, sex.only, cell.only, neither), nrow = 2, ncol = 2)
  t1 = fisher.test(mat, alternative = 'greater')
  
  sexnow = mb
  markers = unique(subset(ct_markers_all, cellName == cellnow)$macID)
  combo.cell = length(which(markers %in% sexnow))
  cell.all = markers[which(markers %in% all)]
  cell.only = length(which(cell.all %!in% sexnow))
  sex.all = length(sexnow)
  sex.only = sex.all- combo.cell
  neither = length(all) - combo.cell - cell.only - sex.only
  mat = matrix(c(combo.cell, sex.only, cell.only, neither), nrow = 2, ncol = 2)
  t2 = fisher.test(mat, alternative = 'greater')
  
  sexnow = sb
  markers = unique(subset(ct_markers_all, cellName == cellnow)$macID)
  combo.cell = length(which(markers %in% sexnow))
  cell.all = markers[which(markers %in% all)]
  cell.only = length(which(cell.all %!in% sexnow))
  sex.all = length(sexnow)
  sex.only = sex.all- combo.cell
  neither = length(all) - combo.cell - cell.only - sex.only
  mat = matrix(c(combo.cell, sex.only, cell.only, neither), nrow = 2, ncol = 2)
  t3 = fisher.test(mat, alternative = 'greater')
  
  celltypes[i,1] = t1$p.value
  celltypes[i,2] = t1$estimate
  celltypes[i,3] = t2$p.value
  celltypes[i,4] = t2$estimate
  celltypes[i,5] = t3$p.value
  celltypes[i,6] = t3$estimate
  
}

rownames(celltypes) = cells
celltypes

celltypes$fbadj = p.adjust(celltypes$V1, method = 'BH')
celltypes$mbadj = p.adjust(celltypes$V3, method = 'BH')
celltypes$sbadj = p.adjust(celltypes$V5, method = 'BH')
celltypes = celltypes[,c(2,1,7,4,3,8,6,5,9)]
colnames(celltypes) = c('female OR', 'female p val', 'female p adj', 'male OR', 'male p val', 'male p adj', 'overall OR', 'overall p val', 'overall p adj')
celltypes

# plot and save 
plotdata = data.frame(genes = c(fb, mb))
plotdata$bias = ifelse(plotdata$genes %in% fb, "Female", "Male")
for(i in 1:length(plotdata$genes)){plotdata$`Cell Type`[i] = vlookup(plotdata$genes[i], dict = ct_markers_all, lookup_column = 1, result_column = 2)}
plotdata2 <- aggregate(.~`Cell Type`+bias, plotdata, length)
sum(subset(plotdata2, bias == 'Female')$genes)
sum(subset(plotdata2, bias == 'Male')$genes)
plotdata2$Proportion = ifelse(plotdata2$bias == "Female", plotdata2$genes/sum(subset(plotdata2, bias == "Female")$genes), plotdata2$genes/sum(subset(plotdata2, bias == "Male")$genes))
plotdata2$`Cell Type` = as.factor(plotdata2$`Cell Type`)
levels(plotdata2$`Cell Type`) = c("Astrocytes","Endothelial Cells","Microglia","Neurons","Oligodendrocytes","OPCs")

ggplot(data=plotdata2, aes(x=bias, y=Proportion, fill=`Cell Type`)) + 
  geom_bar(stat="identity",position="stack") + theme_classic() + 
  xlab("Bias") + ylab("Proportion") + scale_fill_brewer(palette='BrBG') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14))

saveRDS(plotdata2, file = 'checkpoints/macaque_sex_cell_types.rds')

######################################
## cell type enrichment (humans)
######################################

library(BRETIGEA)
library(biomaRt)
library(expss)

markers = markers_df_brain

hsap = useEnsembl(biomart = 'ensembl',dataset='hsapiens_gene_ensembl',mirror='www') 
gene_conv = getBM(attributes=c('ensembl_gene_id','external_gene_name'), mart = hsap)

for (i in 1:length(markers$markers)){
  markers$id[i] = vlookup(markers$markers[i], dict = gene_conv, result_column = 1,lookup_column = 2)}
markers = markers[which(!is.na(markers$id)),]
colnames(markers) = c("name","cell","markers")
table(markers$cell)
markers_data = markers
cells = levels(as.factor(markers$cell))

celltypes = data.frame()

for (i in 1:length(cells)){
  
  print(cells[i])
  cellnow = cells[i]
  
  sexnow = fb
  markers = unique(subset(markers_data, cell == cellnow)$markers)
  combo.cell = length(which(markers %in% sexnow))
  cell.all = markers[which(markers %in% all)]
  cell.only = length(which(cell.all %!in% sexnow))
  sex.all = length(sexnow)
  sex.only = sex.all- combo.cell
  neither = length(all) - combo.cell - cell.only - sex.only
  mat = matrix(c(combo.cell, sex.only, cell.only, neither), nrow = 2, ncol = 2)
  t1 = fisher.test(mat, alternative = 'greater')
  
  sexnow = mb
  markers = unique(subset(markers_data, cell == cellnow)$markers)
  combo.cell = length(which(markers %in% sexnow))
  cell.all = markers[which(markers %in% all)]
  cell.only = length(which(cell.all %!in% sexnow))
  sex.all = length(sexnow)
  sex.only = sex.all- combo.cell
  neither = length(all) - combo.cell - cell.only - sex.only
  mat = matrix(c(combo.cell, sex.only, cell.only, neither), nrow = 2, ncol = 2)
  t2 = fisher.test(mat, alternative = 'greater')
  
  sexnow = sb
  markers = unique(subset(markers_data, cell == cellnow)$markers)
  combo.cell = length(which(markers %in% sexnow))
  cell.all = markers[which(markers %in% all)]
  cell.only = length(which(cell.all %!in% sexnow))
  sex.all = length(sexnow)
  sex.only = sex.all- combo.cell
  neither = length(all) - combo.cell - cell.only - sex.only
  mat = matrix(c(combo.cell, sex.only, cell.only, neither), nrow = 2, ncol = 2)
  t3 = fisher.test(mat, alternative = 'greater')
  
  celltypes[i,1] = t1$p.value
  celltypes[i,2] = t1$estimate
  celltypes[i,3] = t2$p.value
  celltypes[i,4] = t2$estimate
  celltypes[i,5] = t3$p.value
  celltypes[i,6] = t3$estimate
  
}

rownames(celltypes) = cells

celltypes$fbadj = p.adjust(celltypes$V1, method = 'BH')
celltypes$mbadj = p.adjust(celltypes$V3, method = 'BH')
celltypes$sbadj = p.adjust(celltypes$V5, method = 'BH')
celltypes = celltypes[,c(2,1,7,4,3,8,6,5,9)]
colnames(celltypes) = c('female OR', 'female p val', 'female p adj', 'male OR', 'male p val', 'male p adj', 'overall OR', 'overall p val', 'overall p adj')
celltypes

# plot and save
plotdata = data.frame(genes = c(fb, mb))
plotdata$bias = ifelse(plotdata$genes %in% fb, "Female", "Male")
for(i in 1:length(plotdata$genes)){plotdata$`Cell Type`[i] = vlookup(plotdata$genes[i], dict = markers_data, lookup_column = 3, result_column = 2)}
plotdata2 <- aggregate(.~`Cell Type`+bias, plotdata, length)
sum(subset(plotdata2, bias == 'Female')$genes)
sum(subset(plotdata2, bias == 'Male')$genes)
plotdata2$Proportion = ifelse(plotdata2$bias == "Female", plotdata2$genes/sum(subset(plotdata2, bias == "Female")$genes), plotdata2$genes/sum(subset(plotdata2, bias == "Male")$genes))
plotdata2$`Cell Type` = as.factor(plotdata2$`Cell Type`)
levels(plotdata2$`Cell Type`) = c("Astrocytes","Endothelial Cells","Microglia","Neurons","Oligodendrocytes","OPCs")

ggplot(data=plotdata2, aes(x=bias, y=Proportion, fill=`Cell Type`)) + 
  geom_bar(stat="identity",position="stack") + theme_classic() + 
  xlab("Bias") + ylab("Proportion") + scale_fill_brewer(palette='BrBG') +
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14))
  theme(axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.title = element_text(size=14))

saveRDS(plotdata2, file = 'checkpoints/human_sex_cell_types.rds')
