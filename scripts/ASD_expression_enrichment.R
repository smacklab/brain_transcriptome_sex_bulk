#!/usr/bin/env Rscript

##############
# load results
#############

e.keep = checkpoints/readRDS('filtered_expression_matrix.rds')

keep.genes = checkpoints/readRDS('keep_genes.rds')
meta = checkpoints/readRDS('cayo_bulkbrain_combined_metadata.rds')
meta = meta[order(match(meta$LID, colnames(e.keep))),]

mash.results = checkpoints/readRDS('mashr_results.rds')
mash.ct.results = checkpoints/readRDS('mashr_results_cell.rds')
mash.hsap = data/readRDS('gtex_mashr_results_sex.rds')

#################################
# extract coefficients and p vals
# select data set
#################################

library(mashr)

# macaque primary analyses
mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)
mash.sbet = mash.beta / get_psd(mash.results)
mash.beta = mash.beta[,region.levels]
mash.lfsr = mash.lfsr[,region.levels]
mash.sbet = mash.sbet[,region.levels]
ensembl.gene.names = unique(unlist(keep.genes))

# macaque cell type corrected
mash.beta = get_pm(mash.ct.results)
mash.lfsr = get_lfsr(mash.ct.results)
mash.sbet = mash.beta / get_psd(mash.ct.results)
mash.beta = mash.beta[,region.levels]
mash.lfsr = mash.lfsr[,region.levels]
mash.sbet = mash.sbet[,region.levels]
ensembl.gene.names = unique(unlist(keep.genes))

# human gtex
mash.beta = get_pm(mash.hsap)
mash.lfsr = get_lfsr(mash.hsap)
hsap.berr = get_psd(mash.hsap)
mash.sbet = mash.beta / hsap.berr
library(stringr)
rownames(mash.lfsr) = str_sub(rownames(mash.lfsr),1,15)
rownames(mash.beta) = str_sub(rownames(mash.beta),1,15)
ensembl.gene.names = str_sub(rownames(mash.beta),1,15)

##############################
# load disease expression data
##############################

diseases = data/read.csv('Haney2021.csv')

# get orthologs
ortho = checkpoints/readRDS('human_macaque_one2one.rds')
diseases2 = merge(diseases, ortho, by = 'ensembl_gene_id', all.x = T)

## remove Y chrom

library(biomaRt)

# macaques
mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = ensembl.gene.names, mart = mmul)
Ychrom = subset(gene2chrom, chromosome_name == "Y")
ensembl.gene.names = ensembl.gene.names[which(ensembl.gene.names %!in% Ychrom$ensembl_gene_id)]

# humans
hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = ensembl.gene.names, mart = hsap)
Ychrom = subset(gene2chrom, chromosome_name == "Y")
ensembl.gene.names = ensembl.gene.names[which(ensembl.gene.names %!in% Ychrom$ensembl_gene_id)]

# remove chrY

mash.lfsr = mash.lfsr[which(rownames(mash.lfsr) %!in% Ychrom$ensembl_gene_id),]
mash.beta = mash.beta[which(rownames(mash.beta) %!in% Ychrom$ensembl_gene_id),]

# run if limiting to a subset of regions

reg = c("dmPFC","dlPFC","vmPFC","vlPFC","ACCg","M1","STS","V1") # cortical only
reg = c("AMY","CA3","DG","CN","Pu","LGN","VMH") # subcortical only

mash.beta = mash.beta[,which(colnames(mash.beta) %in% reg)]
mash.lfsr = mash.lfsr[,which(colnames(mash.lfsr) %in% reg)]

# gene sex biased genes

mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)

all.region.fet = numeric(length=length(mashr.genes))
names(all.region.fet) = mashr.genes

fraction.cutoff.now = 1/length(reg)

all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff.now) >= fraction.cutoff.now * length(colnames(mash.beta))) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff.now] > 0) >= fraction.cutoff.now * length(colnames(mash.beta))
}))))] = 1
all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff.now) >= fraction.cutoff.now * length(colnames(mash.beta))) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff.now] < 0) >= fraction.cutoff.now * length(colnames(mash.beta))
}))))] = -1
table(all.region.fet)

# run for macaque
all.region.join = data.frame(mmulatta_homolog_ensembl_gene = mashr.genes, direction = as.integer(all.region.fet))
all.region.do.pass = merge(all.region.join, diseases2, by='mmulatta_homolog_ensembl_gene')

# run for human
all.region.join = data.frame(ensembl_gene_id = mashr.genes, direction = as.integer(all.region.fet))
all.region.do.pass = merge(all.region.join, diseases2, by='ensembl_gene_id')

##############
# run enrichment tests
##############

# run for macaques
x = all.region.do.pass[,c(2,6)] #ASD associated
x = all.region.do.pass[,c(2,7)] #ASD up
x = all.region.do.pass[,c(2,8)] #ASD down
rownames(x) = all.region.do.pass$mmulatta_homolog_ensembl_gene

# run for humans
x = all.region.do.pass[,c(2,5)] #ASD associated
x = all.region.do.pass[,c(2,6)] #ASD up
x = all.region.do.pass[,c(2,7)] #ASD down
rownames(x) = all.region.do.pass$ensembl_gene_id

table(x)

f = table(x)[1,]
n = table(x)[2,]
m = table(x)[3,]

contingency.matrix.m = matrix(c(m[2], f[2] + n[2], m[1], f[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
contingency.matrix.f = matrix(c(f[2], m[2] + n[2], f[1], m[1] + n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))
contingency.matrix.all = matrix(c(f[2] + m[2], n[2], f[1] + m[1], n[1]), nrow=2, ncol=2, dimnames=list(c('biased','not biased'),c('disease','not disease')))

contingency.matrix.m
contingency.matrix.f
contingency.matrix.all

fisher.test(contingency.matrix.m,alternative='greater')
fisher.test(contingency.matrix.f,alternative='greater')
fisher.test(contingency.matrix.all,alternative='greater')
