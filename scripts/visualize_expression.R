#!/usr/bin/env Rscript

## load normalized expression data

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

## Remove batch effects

library(limma)

meta = meta[order(match(meta$LID, colnames(e.keep))),]
nobatch = removeBatchEffect(e.keep, batch = meta$Library.batch, covariates = cbind(meta$RIN))

#######
## UMAP 
#######

library(umap)
library(ggplot2)
library(RColorBrewer)

a = umap(t(nobatch), n_neighbors = 200, min_dist = 0.5, metric = 'manhattan')

a.umap = data.frame(as.data.frame(a$layout), Region = meta$Region, Age=meta$exact_age_years, RIN=meta$RIN, Batch = meta$Sequencing.batch, LibBatch = meta$Library.batch, Sex = meta$sex)

a.umap$Region = factor(a.umap$Region, levels = region.levels, labels = region.levels)
mycolors <- region.colors
myshapes <- region.shapes

# plot by region 
p = ggplot(a.umap,aes(x=V1,y=V2,color=Region,shape=Region)) + 
  scale_color_manual(values = mycolors) + geom_point(size=2.5) + 
  theme_classic() + xlab('UMAP 1') + ylab('UMAP 2') + 
  theme(axis.ticks=element_blank(),axis.text=element_blank()) + 
  scale_shape_manual(values = myshapes) +
  theme(axis.title=element_text(size=18), legend.text = element_text(size=18), legend.title = element_text(size=18))
plot(p)

# plot by sex
levels(a.umap$Sex) = c("F","M")
p = ggplot(a.umap,aes(x=V1,y=V2,color=Sex)) + 
  scale_color_manual(values = mycolors[c(3,6)]) + geom_point(size=2.5) + 
  theme_classic() + xlab('UMAP 1') + ylab('UMAP 2') + 
  theme(axis.ticks=element_blank(),axis.text=element_blank()) + 
  scale_shape_manual(values = myshapes) +
  theme(axis.title=element_text(size=18), legend.text = element_text(size=18), legend.title = element_text(size=18))
plot(p)

###########
## t-SNE
###########

library(Rtsne)

a = Rtsne(t(nobatch), dims = 3, perplexity=30, verbose=TRUE, max_iter = 1000)
plot.tsne = as.data.frame(a$Y)
plot.tsne$Region = factor(meta$Region, levels = region.levels, labels = region.levels)
plot.tsne$Sex = meta$sex
levels(plot.tsne$Sex) = c("F","M")

#plot by region
p = ggplot(plot.tsne,aes(x=V1,y=V2,color=Region,shape=Region)) + 
  scale_color_manual(values = mycolors) + geom_point(size=2.5) + 
  scale_shape_manual(values = myshapes) +
  theme_classic() + xlab('Dim 1') + ylab('Dim 2') + 
  theme(axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_text(size=18), legend.text = element_text(size=18), legend.title = element_text(size=18))
plot(p)

#plot by sex
p = ggplot(plot.tsne,aes(x=V1,y=V2,color=Sex)) + 
  scale_color_manual(values = mycolors[c(3,6)]) + geom_point(size=2.5) + 
  theme_classic() + xlab('Dim 1') + ylab('Dim 2') + 
  theme(axis.ticks=element_blank(),axis.text=element_blank(),axis.title=element_text(size=18), legend.text = element_text(size=18), legend.title = element_text(size=18))
plot(p)

###########
## PCA
###########

pca = prcomp(cor(nobatch))
plot.pca = as.data.frame(pca$x)
plot.pca$Sex = factor(meta$sex)
levels(plot.pca$Sex) = c("F","M")
plot.pca$Region = factor(meta$Region, levels = region.levels, labels = region.levels)

# plot by region
p = ggplot(plot.pca,aes(x=PC1,y=PC2,color=Region,shape=Region)) + 
  scale_shape_manual(values = myshapes) + scale_color_manual(values = mycolors) + 
  geom_point(size=2.5) + theme_classic() + xlab('PC 1') + ylab('PC 2') + 
  theme(axis.ticks=element_blank(),axis.text=element_blank(),axis.title = element_text(size=18), legend.text = element_text(size=18), legend.title = element_text(size=18))
plot(p)

# plot by sex
p = ggplot(plot.pca,aes(x=PC1,y=PC2,color=Sex)) + 
  scale_color_manual(values = mycolors[c(3,6)]) + geom_point(size=2.5) + 
  theme_classic() + xlab('PC 1') + ylab('PC 2') + 
  theme(axis.ticks=element_blank(),axis.text=element_blank(), axis.title = element_text(size=18), legend.text = element_text(size=18), legend.title = element_text(size=18))
plot(p)

#########################
## Variance partitioning
#########################

library(variancePartition)
library(expss)
library(BiocParallel)

n.cores = detectCores() - 4
param <- SnowParam(n.cores, "SOCK", progressbar=TRUE)

design = as.formula(paste('~',paste(c('exact_age_years','RIN','(1|Region)','(1|sex)','(1|ordinal.rank)','(1|Library.batch)','(1|Individual)'),collapse=' + ')))
varPart <- fitExtractVarPartModel(e.keep, design, meta, BPPARAM=param)

vp = sortCols(varPart)
vp = vp[,c(1,3,5,6,7,8)]

mean(vp$Region)
mean(vp$exact_age_years)
mean(vp$ordinal.rank)
mean(vp$sex)
mean(vp$Residuals)

library(biomaRt)

mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = rownames(e.keep), mart = mmul)

library(expss)
library(dplyr)
library(reshape2)

vp.chrom = vp
for(i in 1:length(rownames(vp))){
  vp.chrom$chrom[i] = vlookup(rownames(vp.chrom)[i], dict = gene2chrom, lookup_column = 1, result_column = 3)}
vp.chrom$chrom2 = ifelse(vp.chrom$chrom == 'X' | vp.chrom$chrom == 'Y', vp.chrom$chrom, 'Autosomal')
data.frame(vp.chrom %>% group_by(chrom2) %>% summarise(mean_sex = mean(sex)))

obj = data.frame(vp)
col=c(region.colors[c(1:ncol(obj)-1)], "grey85")
obj$gene <- rownames(obj)
data.plot <- melt(obj, id="gene")
data.plot$value <- data.plot$value * 100
ylim = c(0, max(data.plot$value))

g <- ggplot(data.plot, aes(x=variable, y=value)) + geom_boxplot()
split <- split(data.plot, data.plot$variable)
ld <- layer_data(g)
outliers <- lapply(seq_along(split), function(i) {
  box <- ld[ld$group == i, ]
  data <- split[[i]]
  data <- data[data$value > box$ymax , ]
  data
})
outliers <- do.call(rbind, outliers)
for(i in 1:length(outliers$gene)){
  outliers$chrom[i] = vlookup(outliers$gene[i], dict = gene2chrom, lookup_column = 1, result_column = 3)}
outliers$chrom2 = ifelse(outliers$chrom == 'X' | outliers$chrom == 'Y', outliers$chrom, 'Autosomal')

ggplot(data=data.plot, aes(x=variable, y=value)) + 
  geom_violin( scale="width", aes(alpha = 0.8, fill = factor(variable))) + 
  ylab('') + xlab('') + ylim(ylim) + theme_classic() + 
  geom_boxplot(width=0.07, fill="grey", outlier.colour=NA) + 
  scale_fill_manual(values=col) +
  theme(legend.position="none") +
  theme(plot.title=element_text(hjust=0.5)) +
  theme(axis.text.x = element_text(size  = 18,angle = 20,hjust = 1, vjust = 1)) +
  theme(axis.text = element_text(size=18)) +
  geom_point(data = outliers, aes(shape = chrom2, size = chrom2)) +
  scale_shape_manual(values = c(16, 4, 8)) +
  scale_size_manual(values = c(2, 3.5, 3.5)) +
  scale_x_discrete(labels=c("Region" = "Region", "Individual" = "Individual", "exact_age_years" = "Age", "ordinal.rank" = "Rank", "sex" = "Sex", "Residuals" = "Residuals"))

############
## Hierarchical clustering
############

library(pvclust)

# cluster regions

regions = names(keep.genes)

lr = meta[,c('LID','Region')]
lr.split = split(lr,lr$Region)

e.mean = do.call(cbind,lapply(regions,function(r) {
  matrix(apply(nobatch[,lr.split[[r]]$LID],1,mean),ncol=1,dimnames=list(rownames(nobatch),r))}))

all.clust.boot = pvclust(e.mean, method.dist='correlation',method.hclust='average',nboot=1000)
plot(all.clust.boot)
