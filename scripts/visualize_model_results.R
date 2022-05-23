#!/usr/bin/env Rscript

#######################
# load data and results
#######################

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')
meta = meta[order(match(meta$LID, colnames(e.keep))),]

emma.results = readRDS('checkpoints/emma_results.rds')
emma.ct.results = readRDS('checkpoints/emma_results_cell_type.rds')
mash.results = readRDS('checkpoints/mashr_results.rds')
mash.ct.results = readRDS('checkpoints/mashr_results_cell_type.rds')

library(biomaRt)
mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = rownames(e.keep), mart = mmul)
Ychrom = subset(gene2chrom, chromosome_name == "Y")$ensembl_gene_id
Xchrom = subset(gene2chrom, chromosome_name == "X")$ensembl_gene_id

#################################
# extract coefficients and p vals
#################################

emma.beta = emma.results[,paste('beta',predictor,sep='.'),]
emma.pval = emma.results[,paste('pval',predictor,sep='.'),]
emma.qval = apply(emma.pval,2,function(x) p.adjust(x,'fdr'))
emma.beta = emma.beta[,region.levels]
emma.pval = emma.pval[,region.levels]
emma.qval = emma.qval[,region.levels]

emma.beta.ct = emma.ct.results[,paste('beta',predictor,sep='.'),]
emma.pval.ct = emma.ct.results[,paste('pval',predictor,sep='.'),]
emma.qval.ct = apply(emma.pval.ct,2,function(x) p.adjust(x,'fdr'))
emma.beta.ct = emma.beta.ct[,region.levels]
emma.pval.ct = emma.pval.ct[,region.levels]
emma.qval.ct = emma.qval.ct[,region.levels]
get_pm=function(x){x$result$PosteriorMean}
get_lfsr=function(x){x$result$lfsr}

mash.beta = get_pm(mash.results)
mash.lfsr = get_lfsr(mash.results)
mash.beta = mash.beta[,region.levels]
mash.lfsr = mash.lfsr[,region.levels]

mash.beta.ct = get_pm(mash.ct.results)
mash.lfsr.ct = get_lfsr(mash.ct.results)
mash.beta.ct = mash.beta.ct[,region.levels]
mash.lfsr.ct = mash.lfsr.ct[,region.levels]

########################
## inspect and visualize
########################

# Significant genes per region
apply(emma.qval,2,function(x) sum(x < fdr.cutoff,na.rm=TRUE))
apply(emma.qval.ct,2,function(x) sum(x < fdr.cutoff,na.rm=TRUE))
apply(mash.lfsr,2,function(x) sum(x < fsr.cutoff))
apply(mash.lfsr.ct,2,function(x) sum(x < fsr.cutoff))

# Genes that are significant in just one region
apply(mash.lfsr[apply(mash.lfsr,1,function(x) sum(x < fsr.cutoff)) == 1,],2,function(x) sum(x < fsr.cutoff))
apply(mash.lfsr.ct[apply(mash.lfsr.ct,1,function(x) sum(x < fsr.cutoff)) == 1,],2,function(x) sum(x < fsr.cutoff))
apply(mash.lfsr.t[apply(mash.lfsr.t,1,function(x) sum(x < fsr.cutoff)) == 1,],2,function(x) sum(x < fsr.cutoff))

##########################
## plot counts of sex-biased genes
##########################

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(reshape2)

# select data to plot

# mashr vs mashr cell-type

b.mash1 = data.frame(expand.grid(rownames(mash.beta),colnames(mash.beta)), qval=as.numeric(mash.lfsr), beta=as.numeric(mash.beta))
b.mash2 = data.frame(expand.grid(rownames(mash.beta.ct),colnames(mash.beta.ct)), qval=as.numeric(mash.lfsr.ct), beta=as.numeric(mash.beta.ct))
b.mash1$qval.signed = with(b.mash1,ifelse(beta>0,1,-1) * qval)
b.mash2$qval.signed = with(b.mash2,ifelse(beta>0,1,-1) * qval)

model.counts.combined = rbind(
  within(melt(tapply(b.mash1$qval.signed,b.mash1$Var2,function(x) -sum(abs(x) < fsr.cutoff & x < 0,na.rm=TRUE))), {method='Sex-biased genes'; direction='down'} ),
  within(melt(tapply(b.mash1$qval.signed,b.mash1$Var2,function(x) sum(abs(x) < fsr.cutoff & x >= 0,na.rm=TRUE))), {method='Sex-biased genes'; direction='up'} ),
  within(melt(tapply(b.mash2$qval.signed,b.mash2$Var2,function(x) -sum(abs(x) < fsr.cutoff & x < 0,na.rm=TRUE))), {method='Cell-type corrected'; direction='down'} ),
  within(melt(tapply(b.mash2$qval.signed,b.mash2$Var2,function(x) sum(abs(x) < fsr.cutoff & x >= 0,na.rm=TRUE))), {method='Cell-type corrected'; direction='up'} ))
model.counts.combined$method = factor(model.counts.combined$method, levels=c('Sex-biased genes','Cell-type corrected'))

# emma vs mashr

b.emma = data.frame(expand.grid(rownames(emma.beta),colnames(emma.beta)), qval=as.numeric(emma.qval), beta=as.numeric(emma.beta))
b.mash = data.frame(expand.grid(rownames(mash.beta),colnames(mash.beta)), qval=as.numeric(mash.lfsr), beta=as.numeric(mash.beta))
b.emma$qval.signed = with(b.emma,ifelse(beta>0,1,-1) * qval)
b.mash$qval.signed = with(b.mash,ifelse(beta>0,1,-1) * qval)

model.counts.combined = rbind(
  within(melt(tapply(b.emma$qval.signed,b.emma$Var2,function(x) -sum(abs(x) < fdr.cutoff & x < 0,na.rm=TRUE))), {method='EMMA'; direction='down'} ),
  within(melt(tapply(b.emma$qval.signed,b.emma$Var2,function(x) sum(abs(x) < fdr.cutoff & x >= 0,na.rm=TRUE))), {method='EMMA'; direction='up'} ),
  within(melt(tapply(b.mash$qval.signed,b.mash$Var2,function(x) -sum(abs(x) < fsr.cutoff & x < 0,na.rm=TRUE))), {method='MASHR'; direction='down'} ),
  within(melt(tapply(b.mash$qval.signed,b.mash$Var2,function(x) sum(abs(x) < fsr.cutoff & x >= 0,na.rm=TRUE))), {method='MASHR'; direction='up'} ))

# emma cell-type vs mashr cell-type

b.emma = data.frame(expand.grid(rownames(emma.beta.ct),colnames(emma.beta.ct)), qval=as.numeric(emma.qval.ct), beta=as.numeric(emma.beta.ct))
b.mash = data.frame(expand.grid(rownames(mash.beta.ct),colnames(mash.beta.ct)), qval=as.numeric(mash.lfsr.ct), beta=as.numeric(mash.beta.ct))
b.emma$qval.signed = with(b.emma,ifelse(beta>0,1,-1) * qval)
b.mash$qval.signed = with(b.mash,ifelse(beta>0,1,-1) * qval)

model.counts.combined = rbind(
  within(melt(tapply(b.emma$qval.signed,b.emma$Var2,function(x) -sum(abs(x) < fdr.cutoff & x < 0,na.rm=TRUE))), {method='EMMA (cell-type corrected)'; direction='down'} ),
  within(melt(tapply(b.emma$qval.signed,b.emma$Var2,function(x) sum(abs(x) < fdr.cutoff & x >= 0,na.rm=TRUE))), {method='EMMA (cell-type corrected)'; direction='up'} ),
  within(melt(tapply(b.mash$qval.signed,b.mash$Var2,function(x) -sum(abs(x) < fsr.cutoff & x < 0,na.rm=TRUE))), {method='MASHR (cell-type corrected)'; direction='down'} ),
  within(melt(tapply(b.mash$qval.signed,b.mash$Var2,function(x) sum(abs(x) < fsr.cutoff & x >= 0,na.rm=TRUE))), {method='MASHR (cell-type corrected)'; direction='up'} ))

# plot

ylimit = ceiling(with(model.counts.combined,max(abs(value)))/100) * 100

ggplot(
  model.counts.combined,aes(Var1,value,fill=Var1,alpha=direction)) +
  geom_bar(stat='identity') +
  scale_fill_manual(name='Region',values=region.colors) +
  scale_alpha_manual(values=c(0.75,1)) +
  scale_y_continuous(
    limits = c(-ylimit,ylimit),
    breaks = c(-ylimit,-ylimit*0.5,0,ylimit*0.5,ylimit),
    labels = c(formatC(ylimit,width=5,flag=' '),'F',formatC(0,width=5,flag=' '),'M',formatC(ylimit,width=5,flag=' '))
  ) +
  theme_classic(base_size=12) +
  theme(
    axis.ticks.y = element_line(linetype=c(1,0,1,0,1)),
    axis.text.x = element_text(
      angle = -45, hjust = 0, vjust = 1,
    ),
    axis.text.y=element_text(
      face = c('plain','bold','plain','bold','plain'),
      #			size = axis.text.size * c(1,2,1,2,1),
      angle = c(0,90,0,90,0), hjust=0.5
    )
  ) +
  facet_wrap(~method,nrow=2) + 
  theme(strip.text.x = element_text(size = 18), axis.text=element_text(size=18), axis.title=element_text(size=18), strip.background = element_blank()) + 
  theme(legend.position='none') +
  xlab('Regions') +
  ylab('Number of genes') 

#############################
## plot scatterplots of betas
#############################

# select data to plot

# emma vs emma cell-type

b.emma1 = data.frame(expand.grid(rownames(emma.beta),colnames(emma.beta)), qval=as.numeric(emma.qval), beta=as.numeric(emma.beta))
b.emma2 = data.frame(expand.grid(rownames(emma.beta.ct),colnames(emma.beta.ct)), qval=as.numeric(emma.qval.ct), beta=as.numeric(emma.beta.ct))
b.emma1$qval.signed = with(b.emma1,ifelse(beta>0,1,-1) * qval)
b.emma2$qval.signed = with(b.emma2,ifelse(beta>0,1,-1) * qval)

betas = cbind(b.emma1, b.emma2)
betas = betas[,-c(6,7)]
p = ggplot(aes(x=beta, y=beta.1), data = betas) + 
  facet_wrap(~Var2) + geom_point(aes(color=Var2)) + 
  scale_color_manual(values = region.colors) + 
  theme_classic() + xlab("EMMA") + ylab("EMMA (cell-type corrected)") + 
  theme(legend.position = "none", strip.text = element_text(size=18), axis.title = element_text(size=18), axis.text = element_text(size=14))
plot(p)

# emma versus mashr

b.emma = data.frame(expand.grid(rownames(emma.beta),colnames(emma.beta)), qval=as.numeric(emma.qval), beta=as.numeric(emma.beta))
b.mash = data.frame(expand.grid(rownames(mash.beta),colnames(mash.beta)), qval=as.numeric(mash.lfsr), beta=as.numeric(mash.beta))

betas = cbind(b.emma, b.mash)
betas = betas[,-c(6,7)]
p = ggplot(aes(x=beta, y=beta.1), data = betas) + 
  facet_wrap(~Var2) + geom_point(aes(color=Var2)) + 
  scale_color_manual(values = region.colors) + 
  theme_classic() + xlab("EMMA") + ylab("MASHR") + 
  theme(legend.position = "none", strip.text = element_text(size=18), axis.title = element_text(size=18), axis.text = element_text(size=14))
plot(p)

# emma cell-type versus mashr cell-type

b.emma = data.frame(expand.grid(rownames(emma.beta.ct),colnames(emma.beta.ct)), qval=as.numeric(emma.qval.ct), beta=as.numeric(emma.beta.ct))
b.mash = data.frame(expand.grid(rownames(mash.beta.ct),colnames(mash.beta.ct)), qval=as.numeric(mash.lfsr.ct), beta=as.numeric(mash.beta.ct))

betas = cbind(b.emma, b.mash)
betas = betas[,-c(6,7)]
p = ggplot(aes(x=beta, y=beta.1), data = betas) + 
  facet_wrap(~Var2) + geom_point(aes(color=Var2)) + 
  scale_color_manual(values = region.colors) + 
  theme_classic() + xlab("EMMA (cell-type corrected)") + ylab("MASHR (cell-type corrected)") + 
  theme(legend.position = "none", strip.text = element_text(size=18), axis.title = element_text(size=18), axis.text = element_text(size=14))
plot(p)

# mashr versus mashr cell type

b.mash1 = data.frame(expand.grid(rownames(mash.beta),colnames(mash.beta)), qval=as.numeric(mash.lfsr), beta=as.numeric(mash.beta))
b.mash2 = data.frame(expand.grid(rownames(mash.beta.ct),colnames(mash.beta.ct)), qval=as.numeric(mash.lfsr.ct), beta=as.numeric(mash.beta.ct))
b.emma1$qval.signed = with(b.emma1,ifelse(beta>0,1,-1) * qval)
b.emma2$qval.signed = with(b.emma2,ifelse(beta>0,1,-1) * qval)

betas = cbind(b.mash1, b.mash2)
betas = betas[,-c(6,7)]
p = ggplot(aes(x=beta, y=beta.1), data = betas) + 
  facet_wrap(~Var2) + geom_point(aes(color=Var2)) + 
  scale_color_manual(values = region.colors) + 
  theme_classic() + xlab("MASHR") + ylab("MASHR (cell-type corrected)") + theme(legend.position = "none", strip.text = element_text(size=14), axis.title = element_text(size=14), axis.text = element_text(size=14))
plot(p)

cor.test(betas$beta, betas$beta.1, method = 'spearman')

###################################################
## plot counts of significant effects per N regions
###################################################

library(plyr)

# select data set

b.mash = data.frame(expand.grid(rownames(mash.beta),colnames(mash.beta)), qval=as.numeric(mash.lfsr), beta=as.numeric(mash.beta))
b.mash = data.frame(expand.grid(rownames(mash.beta.ct),colnames(mash.beta.ct)), qval=as.numeric(mash.lfsr.ct), beta=as.numeric(mash.beta.ct))

# prep plot

b.mash$qval.signed = with(b.mash,ifelse(beta>0,1,-1) * qval)
b.mash.split.genes = split(b.mash,b.mash$Var1)

region.combinations = table(unlist(lapply(b.mash.split.genes,function(x) {
  paste(subset(x,qval < fsr.cutoff & beta != 0)$Var2,collapse='-')
})))
region.combinations.up = table(unlist(lapply(b.mash.split.genes,function(x) {
  paste(subset(x,qval < fsr.cutoff & beta > 0)$Var2,collapse='-')
})))
region.combinations.down = table(unlist(lapply(b.mash.split.genes,function(x) {
  paste(subset(x,qval < fsr.cutoff & beta < 0)$Var2,collapse='-')
})))
region.combinations.all = rbind(data.frame(region.combinations.up), data.frame(region.combinations.down))
table(region.combinations.all$Var1 %in% data.frame(region.combinations)$Var1)
all.combo = unique(region.combinations.all$Var1)
all.combo = as.character(all.combo[which(all.combo != "")])
region.combinations.all = ddply(region.combinations.all, .(Var1), numcolwise(sum))
region.combinations = region.combinations.all$Freq
names(region.combinations) = region.combinations.all$Var1
region.combinations = region.combinations[as.logical(nchar(names(region.combinations)))]
region.combinations = region.combinations[order(region.combinations,decreasing=TRUE)]
region.combinations.inc = region.combinations.dec = integer(length(region.combinations))
names(region.combinations.inc) = names(region.combinations.dec) = names(region.combinations)

region.combinations.inc[names(region.combinations)] = table(unlist(lapply(b.mash.split.genes,function(x) {
  paste(subset(x,qval < fsr.cutoff & beta > 0)$Var2,collapse='-')
})))[names(region.combinations)]
region.combinations.dec[names(region.combinations)] = table(unlist(lapply(b.mash.split.genes,function(x) {
  paste(subset(x,qval < fsr.cutoff & beta < 0)$Var2,collapse='-')
})))[names(region.combinations)]

region.combinations.inc[is.na(region.combinations.inc)] = 0
region.combinations.dec[is.na(region.combinations.dec)] = 0

region.combinations.sum = region.combinations.inc + region.combinations.dec
region.combinations.sum = sort(region.combinations.sum,decreasing=TRUE)
region.combinations = region.combinations[names(region.combinations.sum)]
region.combinations.inc = region.combinations.inc[names(region.combinations.sum)]
region.combinations.dec = region.combinations.dec[names(region.combinations.sum)]

region.combinations.inc[is.na(region.combinations.inc)] = 0
region.combinations.dec[is.na(region.combinations.dec)] = 0

width.of.bars = 0.8
region.combinations.results = do.call(rbind,lapply(1:length(region.combinations),function(i) {
  x = names(region.combinations)[i]
  n = unlist(lapply(strsplit(names(region.combinations),'-'),length))[i]
  count.all = as.integer(region.combinations[x])
  count.inc = as.integer(region.combinations.inc[x])
  count.dec = as.integer(region.combinations.dec[x])
  out = integer(length(region.levels))
  names(out) = region.levels
  out[unlist(strsplit(x,split='-'))] = 1
  out = rbind(
    data.frame(
      combination=i,
      n_regions = n,
      share_region = n >= fraction.shared.cutoff * length(region.levels),
      region=factor(region.levels,levels=region.levels),
      region_sig = factor(ifelse(as.logical(out),names(out),NA),levels=region.levels),
      value = 1,
      xmin = seq(1,length(region.levels)) - (width.of.bars)/2,
      xmax = seq(1,length(region.levels)) + (width.of.bars)/2,
      ymin = i - (width.of.bars)/2,
      ymax = i + (width.of.bars)/2,
      chart = 'meta'
    ),
    data.frame(
      combination=i,
      n_regions = n,
      share_region = n >= fraction.shared.cutoff * length(region.levels),
      region=NA,
      region_sig=NA,
      value=count.all,
      xmin = NA,
      xmax = NA,
      ymin = NA,
      ymax = NA,
      chart='count_all'
    ),
    data.frame(
      combination=i,
      n_regions = n,
      share_region = n >= fraction.shared.cutoff * length(region.levels),
      region=NA,
      region_sig=NA,
      value=count.inc,
      xmin = NA,
      xmax = NA,
      ymin = NA,
      ymax = NA,
      chart='count_increase'
    ),
    data.frame(
      combination=i,
      n_regions = n,
      share_region = n >= fraction.shared.cutoff * length(region.levels),
      region=NA,
      region_sig=NA,
      value=count.dec,
      xmin = NA,
      xmax = NA,
      ymin = NA,
      ymax = NA,
      chart='count_decrease'
    )
  )
  out
}))

mb = subset(region.combinations.results, chart == "count_increase")[,c(2,6)]
mb = ddply(mb, .(n_regions), numcolwise(sum))
mb$Bias = 'M'
fb = subset(region.combinations.results, chart == "count_decrease")[,c(2,6)]
fb = ddply(fb, .(n_regions), numcolwise(sum))
fb$Bias = 'F'

pie = rbind(fb, mb)
pie = pie[order(pie$n_regions),]
pie$n_regions = as.factor(pie$n_regions)
pie$label = paste(pie$n_regions, paste("(",pie$Bias,")",sep=""), sep = " ")
pie$prop = pie$value/sum(pie$value)
pie$alpha = ifelse(pie$Bias == "F",0.8,1)
pie

library(moonBook)
library(webr)
library(grid)
library(tidyr)
library(RColorBrewer)

# stacked bar plot

pie$Bias = revalue(pie$Bias, c("F"="Female", "M"="Male"))

ggplot(pie, aes(x=n_regions, y = value, fill = Bias)) + geom_bar(position="stack", stat="identity") +
  theme_classic() + scale_fill_manual(values=region.colors[c(3,6)]) +
  xlab("# of regions") + ylab("# of sex-biased genes") +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18), legend.text = element_text(size=18), legend.title = element_text(size=18))

####################
## correlation plots
####################

library(corrplot)
library(corrgram)
library(RColorBrewer)

colnew <- colorRampPalette(region.colors[c(5,6)])

# select data to plot

cornow = cor(emma.beta, use = "complete.obs", method = 'spearman')
cornow = cor(emma.beta.ct, use = "complete.obs", method = 'spearman')
cornow = cor(mash.beta, use = "complete.obs", method = 'spearman')
cornow = cor(mash.beta.ct, use = "complete.obs", method = 'spearman')

corrplot(cornow, col=COL2('BrBG', 10), 
         #type='lower', 
         #method='ellipse', 
         method = 'square',
         tl.col = "black", 
         tl.cex = 1.3, 
         cl.cex = 1.3)

#################################################
## chromosome distribution, enrichment, and plots
#################################################

library(biomaRt)

mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = rownames(e.keep), mart = mmul)

b.emma1 = data.frame(expand.grid(rownames(emma.beta),colnames(emma.beta)), qval=as.numeric(emma.qval), beta=as.numeric(emma.beta))
b.emma2 = data.frame(expand.grid(rownames(emma.beta.ct),colnames(emma.beta.ct)), qval=as.numeric(emma.qval.ct), beta=as.numeric(emma.beta.ct))
b.mash1 = data.frame(expand.grid(rownames(mash.beta),colnames(mash.beta)), qval=as.numeric(mash.lfsr), beta=as.numeric(mash.beta))
b.mash2 = data.frame(expand.grid(rownames(mash.beta.ct),colnames(mash.beta.ct)), qval=as.numeric(mash.lfsr.ct), beta=as.numeric(mash.beta.ct))

library(expss)

b.emma1$chrom = vlookup(b.emma1$Var1, dict = gene2chrom, lookup_column = 1, result_column = 3)
b.emma2$chrom = vlookup(b.emma2$Var1, dict = gene2chrom, lookup_column = 1, result_column = 3)
b.mash1$chrom = vlookup(b.mash1$Var1, dict = gene2chrom, lookup_column = 1, result_column = 3)
b.mash2$chrom = vlookup(b.mash2$Var1, dict = gene2chrom, lookup_column = 1, result_column = 3)

for(i in 1:length(b.emma1$Var1)){if(is.na(b.emma1$beta[i]) == TRUE) {b.emma1$bias[i]=NA} else if (b.emma1$beta[i]<0) {b.emma1$bias[i]="Female"} else {b.emma1$bias[i]="Male"}}
for(i in 1:length(b.emma2$Var1)){if(is.na(b.emma2$beta[i]) == TRUE) {b.emma2$bias[i]=NA} else if (b.emma2$beta[i]<0) {b.emma2$bias[i]="Female"} else {b.emma2$bias[i]="Male"}}
for(i in 1:length(b.mash1$Var1)){if(is.na(b.mash1$beta[i]) == TRUE) {b.mash1$bias[i]=NA} else if (b.mash1$beta[i]<0) {b.mash1$bias[i]="Female"} else {b.mash1$bias[i]="Male"}}
for(i in 1:length(b.mash2$Var1)){if(is.na(b.mash2$beta[i]) == TRUE) {b.mash2$bias[i]=NA} else if (b.mash2$beta[i]<0) {b.mash2$bias[i]="Female"} else {b.mash2$bias[i]="Male"}}

chrom.b.emma1 = unique(subset(b.emma1, qval < fdr.cutoff)[,c(1,5,6)])
chrom.b.emma2 = unique(subset(b.emma2, qval < fdr.cutoff)[,c(1,5,6)])
chrom.b.mash1 = unique(subset(b.mash1, qval < fsr.cutoff)[,c(1,5,6)])
chrom.b.mash2 = unique(subset(b.mash2, qval < fsr.cutoff)[,c(1,5,6)])

chroms = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y")
chrom.b.emma1b = chrom.b.emma1[chrom.b.emma1$chrom %in% chroms,]
chrom.b.emma2b = chrom.b.emma2[chrom.b.emma2$chrom %in% chroms,]
chrom.b.mash1b = chrom.b.mash1[chrom.b.mash1$chrom %in% chroms,]
chrom.b.mash2b = chrom.b.mash2[chrom.b.mash2$chrom %in% chroms,]

chrom.b.emma1b$chrom = factor(chrom.b.emma1b$chrom, levels = chroms)
chrom.b.emma2b$chrom = factor(chrom.b.emma2b$chrom, levels = chroms)
chrom.b.mash1b$chrom = factor(chrom.b.mash1b$chrom, levels = chroms)
chrom.b.mash2b$chrom = factor(chrom.b.mash2b$chrom, levels = chroms)

d.emma1 <- aggregate(.~chrom+bias, chrom.b.emma1b, length)
d.emma2 <- aggregate(.~chrom+bias, chrom.b.emma2b, length)
d.mash1 <- aggregate(.~chrom+bias, chrom.b.mash1b, length)
d.mash2 <- aggregate(.~chrom+bias, chrom.b.mash2b, length)

# select data to plot

now=d.mash1
now=d.mash2

# calculate proportions

totals = data.frame(table(subset(gene2chrom, chromosome_name %in% chroms)$chromosome_name))
colnames(totals) = c('chrom','Var1')
totals$non = 0
for(i in 1:length(totals$chrom)){
  fcount = subset(now, bias == "Female" & chrom == totals$chrom[i])$Var1
  mcount = subset(now, bias == "Male" & chrom == totals$chrom[i])$Var1
  totals$non[i] = totals$Var1[i] - (fcount + mcount)
} # ignore error
totals[22,3] = 0
totals.format = data.frame(chrom = totals$chrom, bias = "None", Var1 = totals$non)
d.mash1b = rbind(now, totals.format)
d.mash1c = d.mash1b
for (i in 1:length(d.mash1c$chrom)){
  d.mash1c$total[i] = subset(totals, chrom == d.mash1c$chrom[i])$Var1}
d.mash1c$prop = d.mash1c$Var1 / d.mash1c$total

ggplot(data=d.mash1c, aes(x=chrom, y=prop, fill=bias)) + 
  geom_bar(stat="identity",position="stack") + 
  theme_classic() + xlab("Chromosome") + ylab("Proportion") + 
  scale_fill_manual(values=c(region.colors[3],region.colors[6],'#999999')) +
  labs(fill = "Bias") +
  theme(axis.text=element_text(size=14), axis.title=element_text(size=18), legend.text = element_text(size=18), legend.title = element_text(size=18))

# describe distribution 

length(chrom.b.mash1$Var1) # 561
length(chrom.b.mash2$Var1) # 427 (5 genes different directions)

length(chrom.b.mash1$Var1) / length(rownames(e.keep))
length(chrom.b.mash2$Var1) / length(rownames(e.keep))

length(subset(chrom.b.mash1, bias == 'Male')$Var1) / length(chrom.b.mash1$Var1)
length(subset(chrom.b.mash1, bias == 'Female')$Var1) / length(chrom.b.mash1$Var1)
5 / (length(chrom.b.mash2$Var1)-5)
(length(subset(chrom.b.mash2, bias == 'Male')$Var1)-5) / (length(chrom.b.mash2$Var1)-5)
(length(subset(chrom.b.mash2, bias == 'Female')$Var1)-5) / (length(chrom.b.mash2$Var1)-5)

length(subset(chrom.b.mash1, chrom == "X")$Var1) / length(chrom.b.mash1$Var1)
length(subset(chrom.b.mash1, chrom == "Y")$Var1) / length(chrom.b.mash1$Var1)
length(subset(chrom.b.mash1, chrom != "X" & chrom != "Y")$Var1) / length(chrom.b.mash1$Var1)
length(subset(chrom.b.mash2, chrom == "X")$Var1) / length(chrom.b.mash2$Var1)
length(subset(chrom.b.mash2, chrom == "Y")$Var1) / length(chrom.b.mash2$Var1)
length(subset(chrom.b.mash2, chrom != "X" & chrom != "Y")$Var1) / length(chrom.b.mash2$Var1)

## chromosome enrichment

now = b.mash1
now = b.mash2

chrom_en = data.frame()

for (j in 1:length(chroms)){
  df = now
  df = df[which(!is.na(df$beta)),]
  all <- df$Var1
  genes_in_chr <- length(unique(df[df$chrom == chroms[j],]$Var1))
  genes_total <- length(unique(all))
  
  hits <- unique(df[df$qval < fsr.cutoff,]$Var1)
  hits_total <- length(hits)
  hits_in_chr <- length(unique(df[df$chrom == chroms[j] & df$qval < fsr.cutoff,]$Var1))
  mat <- matrix(c(hits_in_chr,genes_in_chr-hits_in_chr,hits_total-hits_in_chr,genes_total-hits_total-genes_in_chr+hits_in_chr),nrow=2,ncol=2)
  fr <- fisher.test(mat, alternative="greater")
  all = data.frame(chrom=chroms[j], or=fr$estimate[["odds ratio"]], p=fr$p.value, group='all')
  
  hits <- unique(df[df$qval < fsr.cutoff & df$beta>0,]$Var1)
  hits_total <- length(hits)
  hits_in_chr <- length(unique(df[df$chrom == chroms[j] & df$qval < fsr.cutoff & df$beta>0,]$Var1))
  mat <- matrix(c(hits_in_chr,genes_in_chr-hits_in_chr,hits_total-hits_in_chr,genes_total-hits_total-genes_in_chr+hits_in_chr),nrow=2,ncol=2)
  fr <- fisher.test(mat, alternative="greater")
  male = data.frame(chrom=chroms[j], or=fr$estimate[["odds ratio"]], p=fr$p.value, group='male')
  
  hits <- unique(df[df$qval < fsr.cutoff & df$beta<0,]$Var1)
  hits_total <- length(hits)
  hits_in_chr <- length(unique(df[df$chrom == chroms[j] & df$qval < fsr.cutoff & df$beta<0,]$Var1))
  mat <- matrix(c(hits_in_chr,genes_in_chr-hits_in_chr,hits_total-hits_in_chr,genes_total-hits_total-genes_in_chr+hits_in_chr),nrow=2,ncol=2)
  fr <- fisher.test(mat, alternative="greater")
  female = data.frame(chrom=chroms[j], or=fr$estimate[["odds ratio"]], p=fr$p.value, group='female')
  
  chrom_en = rbind(chrom_en, all, male, female)
  
}
chrom_en$padj = p.adjust(chrom_en$p, method = 'bonferroni')
chrom_en$padj = round(chrom_en$padj, digits = 3)
chrom_en

##################################
## compare betas for sex chromosome vs autosomal genes
##################################

sex_ef = subset(b.mash1, qval < fsr.cutoff & chrom != "Y")
sex_ef = subset(b.mash2, qval < fsr.cutoff & chrom != "Y")

for (i in 1:length(sex_ef$Var1)) {if(sex_ef$chrom[i] == "X") {sex_ef$chrom2[i] = "X"} else {sex_ef$chrom2[i] = "Autosome"}}

p = ggplot(data = sex_ef, aes(x = chrom2, y = beta)) + 
  geom_violin(adjust = 0.5) + theme_classic() + 
  geom_hline(aes(yintercept = 0), color = "grey", linetype="dashed") + 
  ylab("Beta") + xlab("Chromosome") +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 14))
p_build = ggplot2::ggplot_build(p)$data[[1]]
p_build = transform(p_build, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
p_build <- rbind(plyr::arrange(transform(p_build, x = xminv), y), plyr::arrange(transform(p_build, x = xmaxv), -y))
p_build$Bias <- ifelse(p_build$y >= 0,'Male','Female')
p_build$group1 <- with(p_build,interaction(factor(group),factor(Bias)))
p_fill = ggplot() + geom_polygon(data = p_build, aes(x = x,y = y,group = group1,fill = Bias)) + 
  scale_fill_manual(values = c(region.colors[c(3,6)])) + theme_classic() +
  xlab("Chromosome") + ylab("Sex effect (Beta)") +
  geom_hline(aes(yintercept = 0), color = "grey", linetype="dashed") + ylim(c(-1,1)) +
  scale_x_continuous(breaks=seq(from=1,to=2, by = 1), labels=c("Autosome","X")) +
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 18), legend.text = element_text(size=18), legend.title = element_text(size=18))
p_fill

for(i in 1:length(region.levels)){
  print(paste("female # is ",length(unique(subset(sex_ef, Var2 == region.levels[i] & bias == "Female")$Var1))))
  print(paste("male # is ",length(unique(subset(sex_ef, Var2 == region.levels[i] & bias == "Male")$Var1))))
}

length(unique(sex_ef$Var1)) # 9 Y chrom genes removed

mean(abs(subset(sex_ef, beta >0)$beta))
mean(abs(subset(sex_ef, beta <0)$beta))
t.test(x = abs(subset(sex_ef, beta >0)$beta), y = abs(subset(sex_ef, beta <0)$beta))

mean(subset(sex_ef, chrom2 == "X" & beta < 0)$beta) # X female
mean(subset(sex_ef, chrom2 == "X" & beta > 0)$beta) # X male
mean(subset(sex_ef, chrom2 == "Autosome" & beta < 0)$beta) # autosomal female
mean(subset(sex_ef, chrom2 == "Autosome" & beta > 0)$beta) # autosomal male

length(unique(subset(sex_ef, chrom2 == "X" & beta < 0)$Var1))
length(unique(subset(sex_ef, chrom2 == "X" & beta > 0)$Var1))
length(unique(subset(sex_ef, chrom2 == "Autosome" & beta < 0)$Var1))
length(unique(subset(sex_ef, chrom2 == "Autosome" & beta > 0)$Var1))

sex_ef2 = sex_ef
sex_ef2$key = paste(sex_ef$bias, sex_ef$chrom2)
mod = aov(abs(beta) ~ key, data = sex_ef2)
summary(mod)
tuk = TukeyHSD(mod)
tuk$key[,'p adj']
