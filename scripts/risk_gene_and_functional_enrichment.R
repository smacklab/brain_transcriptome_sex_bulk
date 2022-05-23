#!/usr/bin/env Rscript

# load data and results

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')

keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')
meta = meta[order(match(meta$LID, colnames(e.keep))),]

mash.results = readRDS('checkpoints/mashr_results.rds')
mash.ct.results = readRDS('checkpoints/mashr_results_cell_type.rds')

#################################
# extract coefficients and p vals 
# select data set
#################################

library(mashr)
library(stringr)

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

# human GTEx
mash.hsap = readRDS('data/gtex_mashr_results_sex.rds')
mash.beta = get_pm(mash.hsap)
mash.lfsr = get_lfsr(mash.hsap)
hsap.berr = get_psd(mash.hsap)
mash.sbet = mash.beta / hsap.berr
rownames(mash.lfsr) = str_sub(rownames(mash.lfsr),1,15)
rownames(mash.beta) = str_sub(rownames(mash.beta),1,15)
rownames(mash.sbet) = str_sub(rownames(mash.sbet),1,15)
ensembl.gene.names = rownames(mash.beta)

## remove Y chr (for macaque)

library(biomaRt)
`%!in%` = Negate(`%in%`)

mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)

mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = rownames(e.keep), mart = mmul)
Ychrom = subset(gene2chrom, chromosome_name == "Y")

mash.beta = mash.beta[which(rownames(mash.beta) %!in% Ychrom$ensembl_gene_id),]
mash.lfsr = mash.lfsr[which(rownames(mash.lfsr) %!in% Ychrom$ensembl_gene_id),]
mash.sbet = mash.sbet[rownames(sbetnow)  %!in% Ychrom$ensembl_gene_id,]
ensembl.gene.names = ensembl.gene.names[which(ensembl.gene.names %!in% Ychrom$ensembl_gene_id)]
mashr.genes = mashr.genes[which(mashr.genes %!in% Ychrom$ensembl_gene_id)]
all.region = numeric(length=length(ensembl.gene.names))
names(all.region) = ensembl.gene.names

## remove Y chr (for humans)

hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = ensembl.gene.names, mart = hsap)
Ychrom = subset(gene2chrom, chromosome_name == "Y")

mash.lfsr = mash.lfsr[which(rownames(mash.lfsr) %!in% Ychrom$ensembl_gene_id),]
mash.beta = mash.beta[which(rownames(mash.beta) %!in% Ychrom$ensembl_gene_id),]
mashr.genes = rownames(mash.beta)
names(mashr.genes) = rownames(mash.beta)

################
## gene ontology (macaques)
################

# Code male-biased genes as 1 and female-biased genes as -1 if they are significant and co-directional in at least a fraction [fraction.shared.cutoff] of regions

all.region[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(lsfrnow[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(betanow[x,][lsfrnow[x,] < fsr.cutoff] > 0) >= fraction.shared.cutoff * length(keep.genes)
}))))] = 1
all.region[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(lsfrnow[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(betanow[x,][lsfrnow[x,] < fsr.cutoff] < 0) >= fraction.shared.cutoff * length(keep.genes)
}))))] = -1
table(all.region) 

# Code male-biased genes as 1 and female-biased genes as -1 if they are significant and co-directional in at least 1 region

all.region[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(lsfrnow[x,] < fsr.cutoff) >= 1/15 * length(keep.genes)) && sum(betanow[x,][lsfrnow[x,] < fsr.cutoff] > 0) >= 1/15 * length(keep.genes)
}))))] = 1
all.region[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(lsfrnow[x,] < fsr.cutoff) >= 1/15 * length(keep.genes)) && sum(betanow[x,][lsfrnow[x,] < fsr.cutoff] < 0) >= 1/15 * length(keep.genes)
}))))] = -1
table(all.region)  

# run VISEAGO

library(ViSEAGO)

Ensembl<-ViSEAGO::Ensembl2GO("genes")

myGENE2GO<-ViSEAGO::annotate(
  "mmulatta_gene_ensembl",
  Ensembl)

BP.m<-ViSEAGO::create_topGOdata(
  geneSel=names(all.region[which(all.region == 1)]),
  allGenes=names(all.region),
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

BP.f<-ViSEAGO::create_topGOdata(
  geneSel=names(all.region[which(all.region == -1)]),
  allGenes=names(all.region),
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

test.m<-topGO::runTest(
  BP.m,
  algorithm ="parentchild",
  statistic = "fisher")

test.f<-topGO::runTest(
  BP.f,
  algorithm ="parentchild",
  statistic = "fisher")

# run merge_enrich_terms_0.05.R

BP_sResults<-merge_enrich_terms_0.05(
  Input=list(males=c("BP.m","test.m"),females=c("BP.f","test.f")))

BP_sResults2 = BP_sResults@data

myGOs<-ViSEAGO::build_GO_SS(
  gene2GO=myGENE2GO,
  enrich_GO_terms=BP_sResults)

myGOs<-ViSEAGO::compute_SS_distances(
  myGOs,
  distance="Wang")

Wang_clusters_wardD2<-ViSEAGO::GOterms_heatmap(
  myGOs,
  showIC=TRUE,
  showGOlabels=TRUE,
  GO.tree=list(
    tree=list(
      distance="Wang",
      aggreg.method="ward.D2"
    ),
    cut=list(
      dynamic=list(
        pamStage=TRUE,
        pamRespectsDendro=TRUE,
        deepSplit=2,
        minClusterSize =2))),
  samples.tree=NULL)

ViSEAGO::show_heatmap(
  Wang_clusters_wardD2,
  "GOterms")

Wang = Wang_clusters_wardD2@enrich_GOs@data
View(Wang)

clusters= unique(Wang$GO.cluster)
for(i in 1:length(clusters)){
  print(clusters[i])
  now = subset(Wang, GO.cluster == clusters[i])
  print(paste('females: ',sum(now$females.pvalue < 0.05)))
  print(paste('males: ',sum(now$males.pvalue < 0.05)))
}

#####################
## risk gene enrichment (macaques or humans)
#####################

# Import disease associations from DISEASES dataset (i.e., "Disease Ontology")

do.data = read.table('data/human_disease_associations.tsv',
                     sep='\t',
                     quote='',
                     col.names=c('protein_id','protein_name','do_id','do_name','z_score','confidence'),
                     stringsAsFactors=FALSE)

do.def = unique(subset(do.data,select=c('do_id','do_name')))

ignore.checkpoints = FALSE
if (ignore.checkpoints || !file.exists('checkpoints/disease_human_orthologs.rds')) {
  library(biomaRt)
  
  hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',
                 dataset='hsapiens_gene_ensembl')
  
  hsap.info = getBM(
    attributes=c('ensembl_gene_id','ensembl_peptide_id','external_gene_name','mmulatta_homolog_ensembl_gene','mmulatta_homolog_orthology_type'),
    mart = hsap)
  
  saveRDS(hsap.info,file='checkpoints/disease_human_orthologs.rds')
} else {
  message('Checkpoint found!\nLoading human ortholog annotations from file.')
  
  hsap.info = readRDS('checkpoints/disease_human_orthologs.rds')
}

# For matching ensembl peptides, linking genes is straightforward
do.ensembl = subset(do.data,grepl('^ENSP[0-9]{11}',protein_id))
do.ensembl = merge(do.ensembl,hsap.info,by.x='protein_id','ensembl_peptide_id',all.x=FALSE,all.y=FALSE)
do.ensembl$ensembl_peptide_id = do.ensembl$protein_id

# For everything else
do.proteinname = subset(do.data,!grepl('^ENSP[0-9]{11}',protein_id))
do.proteinname = merge(do.proteinname,hsap.info,by.x='protein_name',by.y='external_gene_name',all.x=FALSE,all.y=FALSE)
do.proteinname$external_gene_name = do.proteinname$protein_name

do.all = rbind(do.ensembl[intersect(names(do.ensembl),names(do.proteinname))],do.proteinname[intersect(names(do.ensembl),names(do.proteinname))])

do.mmul = subset(do.all,
                 nchar(mmulatta_homolog_ensembl_gene) > 0 & 
                   mmulatta_homolog_orthology_type == 'ortholog_one2one',
                 select=c('mmulatta_homolog_ensembl_gene','external_gene_name','ensembl_gene_id','ensembl_peptide_id','do_id','do_name','z_score','confidence'))

names(do.mmul) = c('ensembl_gene_id','external_gene_name','hsapiens_homolog_ensembl_gene','hsapiens_homolog_ensembl_peptide_id','do_id','do_name','z_score','confidence')

do.mmul.out = subset(do.mmul,ensembl_gene_id %in% ensembl.gene.names,select=c('ensembl_gene_id','do_id','do_name','z_score','confidence'))
rownames(do.mmul.out) = NULL
saveRDS(do.mmul.out,file='rnaseq_genes_do_d0.rds')

all.region.fet = all.region.kst = numeric(length=length(ensembl.gene.names))
names(all.region.fet) = names(all.region.kst) = ensembl.gene.names

# Code upregulated genes as 1 and downregulated genes as -1 if they are significant and co-directional in at least a fraction [fraction.shared.cutoff] of regions

all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] > 0) >= fraction.shared.cutoff * length(keep.genes)
}))))] = 1
all.region.fet[names(which(unlist(lapply(mashr.genes,function(x) {
  (sum(mash.lfsr[x,] < fsr.cutoff) >= fraction.shared.cutoff * length(keep.genes)) && sum(mash.beta[x,][mash.lfsr[x,] < fsr.cutoff] < 0) >= fraction.shared.cutoff * length(keep.genes)
}))))] = -1
table(all.region.fet)

all.region.kst[rownames(mash.sbet)] = rowMeans(mash.sbet)
all.region.join = data.frame(ensembl_gene_id = ensembl.gene.names, direction = as.integer(all.region.fet), effect = as.numeric(all.region.kst))
all.region.do = merge(all.region.join, do.mmul, by='ensembl_gene_id')

all.region.do.pass = subset(all.region.do,do_id %in% names(which(table(subset(all.region.do,confidence >= 0)$do_id) >= 10)))
all.region.do.gene.pass = unique(all.region.do.pass[c('ensembl_gene_id','direction')])

do.deg.total = as.integer(table(factor(all.region.do.pass$direction != 0,levels=c('TRUE','FALSE'))))
do.inc.total = as.integer(table(factor(all.region.do.pass$direction > 0,levels=c('TRUE','FALSE'))))
do.dec.total = as.integer(table(factor(all.region.do.pass$direction < 0,levels=c('TRUE','FALSE'))))

all.region.do.split = split(all.region.do.pass,all.region.do.pass$do_id)

library(parallel)

n.cores = detectCores() - 4

all.region.do.test = do.call(rbind,mclapply(names(all.region.do.split),function(i) {
  x = all.region.do.split[[i]]
  
  this.deg.total = as.integer(table(factor(x$direction != 0,levels=c('TRUE','FALSE'))))
  contingency.matrix.deg = matrix(rbind(this.deg.total,do.deg.total - this.deg.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))
  
  this.inc.total = as.integer(table(factor(x$direction > 0,levels=c('TRUE','FALSE'))))
  contingency.matrix.inc = matrix(rbind(this.inc.total,do.inc.total - this.inc.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))
  
  this.dec.total = as.integer(table(factor(x$direction < 0,levels=c('TRUE','FALSE'))))
  contingency.matrix.dec = matrix(rbind(this.dec.total,do.dec.total - this.dec.total),nrow=2,dimnames=list(c('in group','not in group'),c('in direction','not in direction')))
  
  deg.fet.test = fisher.test(contingency.matrix.deg,alternative='greater')
  inc.fet.test = fisher.test(contingency.matrix.inc,alternative='greater')
  dec.fet.test = fisher.test(contingency.matrix.dec,alternative='greater')
  inc.kst.test = ks.test(x$effect,subset(all.region.do.pass,do_id != i)$effect,alternative='less')
  dec.kst.test = ks.test(x$effect,subset(all.region.do.pass,do_id != i)$effect,alternative='greater')
  
  data.frame(
    do_id = unique(x$do_id),
    do.size = sum(this.deg.total),
    deg.n = this.deg.total[1],
    inc.n = this.inc.total[1],
    dec.n = this.dec.total[1],
    deg.fet.score = deg.fet.test$estimate,
    inc.fet.score = inc.fet.test$estimate,
    dec.fet.score = dec.fet.test$estimate,
    inc.kst.score = inc.kst.test$statistic,
    dec.kst.score = dec.kst.test$statistic,
    deg.fet.pval = deg.fet.test$p.value,
    inc.fet.pval = inc.fet.test$p.value,
    dec.fet.pval = dec.fet.test$p.value,
    inc.kst.pval = inc.kst.test$p.value,
    dec.kst.pval = dec.kst.test$p.value
  )
},mc.cores=n.cores))

all.region.do.test = within(all.region.do.test,{
  dec.kst.qval = p.adjust(dec.kst.pval,'fdr')
  inc.kst.qval = p.adjust(inc.kst.pval,'fdr')
  dec.fet.qval = p.adjust(dec.fet.pval,'fdr')
  inc.fet.qval = p.adjust(inc.fet.pval,'fdr')
  deg.fet.qval = p.adjust(deg.fet.pval,'fdr')
})

all.region.do.results = merge(all.region.do.test,do.def,by='do_id')
all.region.do.results$dataset = 'DISEASES'
all.region.do.results$set = 'union'
all.region.do.results$region = 'all'
View(all.region.do.results)
subset(all.region.do.results, inc.kst.qval < 0.05)[,c('do_name','inc.kst.score','inc.kst.qval')]
subset(all.region.do.results, dec.kst.qval < 0.05)[,c('do_name','inc.kst.score','inc.kst.qval')]

#write.csv(all.region.do.results, file = "checkpoints/manual_do_results_macaque.csv")
#write.csv(all.region.do.results, file = "checkpoints/manual_do_results_human.csv")


