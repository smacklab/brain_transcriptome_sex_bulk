#!/usr/bin/env Rscript

# load data and results

library(plyr)

resid_exp = readRDS(file='checkpoints/residual_expression.rds)

keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

for(i in 1:length(resid_exp)){
  resid_exp[[i]] = data.frame(resid_exp[[i]])
  resid_exp[[i]]$ROWNAMES  <- rownames(resid_exp[[i]])}
resid_exp_all <- join_all( resid_exp, by="ROWNAMES", type="full" )
rownames(resid_exp_all) <- resid_exp_all$ROWNAMES; resid_exp_all$ROWNAMES <- NULL

meta = meta[order(match(meta$LID, colnames(resid_exp_all))),]

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')

txi.gene = readRDS('checkpoints/kallisto_genes.RDS')
tpm.gene = txi.gene$abundance
tpm.gene = tpm.gene[rownames(tpm.gene) %in% rownames(e.keep),]

emma.results = readRDS('checkpoints/emma_results.rds')

#############
## pleiotropy
#############

## calculate tissue specificity per gene

regions = region.levels

pleioT = data.frame()
for (i in 1:length(tpm.gene[,1])){
  geneexp = tpm.gene[i,]
  geneexp[which(geneexp < 1)] = 1
  regionexp = data.frame()
  for (j in 1:length(regions)){
    m = subset(meta, Region == regions[j])
    regionexp[j,1] = mean(geneexp[which(names(geneexp) %in% m$LID)])}
  regionexp[,2] = regionexp[,1]
  maxE = max(regionexp[,2], na.rm = TRUE)
  if (maxE == -Inf) {maxE = NA} else (maxE = maxE)
  for (k in 1:length(regionexp[,1])){
    regionexp[k,3] = 1 - log(regionexp[k,2])/log(maxE)}
  lengthnow = length(regionexp[,3][!is.na(regionexp[,3])])-1
  if(lengthnow == -1) {lengthnow = NA} else (lengthnow=lengthnow)
  pleioT[i,1] = sum(regionexp[,3], na.rm = TRUE) / lengthnow
}
row.names(pleioT) = row.names(tpm.gene)

min(pleioT[,1]) ##lower=evenly distributed
max(pleioT[,1]) ##higher=greater tissue specificity
mean(pleioT[,1])
sd(pleioT[,1])
hist(pleioT[,1])

saveRDS(pleioT, file = 'checkpoints/tissue_specificity_per_gene.rds')

## plot tissue specificity ~ sex differences in expression

pleioT = readRDS(file = 'checkpoints/tissue_specificity_per_gene.rds')

mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl') 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = rownames(e.keep), mart = mmul)
Ychrom = subset(gene2chrom, chromosome_name == "Y")$ensembl_gene_id
Xchrom = subset(gene2chrom, chromosome_name == "X")$ensembl_gene_id

library(ggplot2)
library(expss)
library(RColorBrewer)

pleioTlim2 = data.frame(pleioT)
for (i in 1:length(pleioTlim2[,1])){
  m = subset(meta, sex == "m")$LID
  f = subset(meta, sex == "f")$LID
  mexp = t(resid_exp_all[rownames(pleioTlim2)[i],colnames(resid_exp_all)%in%m])[,1]
  fexp = t(resid_exp_all[rownames(pleioTlim2)[i],colnames(resid_exp_all)%in%f])[,1]
  pleioTlim2$m_to_f[i] = abs(mean(mexp[!is.na(mexp)]) - mean(fexp[!is.na(fexp)]))
}
hist(pleioTlim2$m_to_f, breaks = 1000)

pleioTlim2 = subset(pleioTlim2, rownames(pleioTlim2) %!in% Ychrom)
pleioTlim2$bias = ifelse(pleioTlim2$m_to_f<0, 'Female', 'Male')

cor.test(x=pleioTlim2$m_to_f, y=log(pleioTlim2$V1), method = 'spearman') 

pleio_mov = data.frame(seq = seq(from = 0, to = round(max(pleioTlim2$m_to_f)+0.1,1), by = 0.05))
for(i in 1:length(pleio_mov$seq)){
  pleio_mov$avg[i] = mean(subset(pleioTlim2, m_to_f >= pleio_mov$seq[i] & m_to_f < pleio_mov$seq[i+1])$V1)}
colnames(pleio_mov) = c("m_to_f","V1")
summary(lm((V1) ~ m_to_f, data = pleio_mov[!is.na(pleio_mov$V1),]))

pleio_combo = rbind(cbind(pleio_mov,type="avg") , cbind(V1 = pleioTlim2$V1, m_to_f = pleioTlim2$m_to_f, type = "obs"))
pleio_combo$V1 = as.numeric(pleio_combo$V1)
pleio_combo$m_to_f = as.numeric(pleio_combo$m_to_f)

ggplot(pleio_combo) +
  geom_point(aes(x = m_to_f, y = V1, color = type)) +
  geom_point(data = subset(pleio_combo, type == 'avg'), aes(x = m_to_f, y = V1, color = type)) +
  scale_color_manual(values = c("black","grey"), labels = c("Moving Average", "Observed"), name = "") +
  geom_smooth(data=subset(pleio_combo, type == 'avg'), aes(x = m_to_f, y = V1), formula=y~x, method='lm') +
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0), labels = scales::number_format(accuracy = 0.01), breaks = seq(from = 0, to = 3, by = 0.5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0.01,1.1)) +
  xlab("Absolute sex difference in mean expression") + ylab("Tissue specificity") +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18), legend.text= element_text(size=18), legend.position = 'bottom', legend.direction = "horizontal")

################
## LOF
################

## based on the ratio of LoF to synonymous mutations
## lower  LoFtool percentiles = genes more intolerant to functional variation

lof = read.csv('data/LoFtools.csv')

library(biomaRt)
library(expss)

hsap = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl') 
gene_conv2 = getBM(attributes=c('ensembl_gene_id','external_gene_name','mmulatta_homolog_ensembl_gene','mmulatta_homolog_orthology_type'), mart = hsap)
gene_conv2 = subset(gene_conv2, mmulatta_homolog_orthology_type == 'ortholog_one2one')
for(i in 1:length(lof$Gene)){lof$mac_gene[i] = vlookup(lof$Gene[i], dict = gene_conv2, lookup_column = 2, result_column = 3)}

lof_dat = data.frame(rownames = rownames(resid_exp_all))
for (i in 1:length(lof_dat[,1])){
  m = subset(meta, sex == "m")$LID
  f = subset(meta, sex == "f")$LID
  mexp = t(resid_exp_all[lof_dat$rownames[i],colnames(resid_exp_all)%in%m])[,1]
  fexp = t(resid_exp_all[lof_dat$rownames[i],colnames(resid_exp_all)%in%f])[,1]
  lof_dat$m_to_f[i] = abs(mean(mexp[!is.na(mexp)]) - mean(fexp[!is.na(fexp)]))}
lof_dat = subset(lof_dat, m_to_f != "Inf")
hist(lof_dat$m_to_f, breaks = 1000)

for(i in 1:length(lof_dat$rownames)){lof_dat$lof[i] = vlookup(lof_dat$rownames[i], dict = lof, lookup_column = 7, result_column = 4)}
lof_dat = lof_dat[which(!is.na(lof_dat$lof)),]
rownames(lof_dat) = lof_dat$rownames
lof_dat = lof_dat[,c(-1)]
lof_dat = subset(lof_dat, rownames(lof_dat) %!in% Ychrom)

library(ggplot2)

cor.test(x=lof_dat$m_to_f, y=(lof_dat$lof), method = 'spearman') 
lof_mov = data.frame(seq = seq(from = 0, to = round(max(lof_dat$m_to_f)+0.1,1), by = 0.05))
for(i in 1:length(lof_mov$seq)){
  lof_mov$avg[i] = mean(subset(lof_dat, m_to_f >= lof_mov$seq[i] & m_to_f < lof_mov$seq[i+1])$lof)}
colnames(lof_mov) = c("m_to_f","lof")
summary(lm((lof) ~ m_to_f, data = lof_mov[!is.na(lof_mov$lof),])) 

lof_combo = rbind(cbind(lof_mov,type="avg") , cbind(lof = lof_dat$lof, m_to_f = lof_dat$m_to_f, type = "obs"))
lof_combo$lof = as.numeric(lof_combo$lof)
lof_combo$m_to_f = as.numeric(lof_combo$m_to_f)

ggplot(lof_combo) +
  geom_point(aes(x = m_to_f, y = lof, color = type)) +
  geom_point(data = subset(lof_combo, type == 'avg'), aes(x = m_to_f, y = lof, color = type)) +
  scale_color_manual(values = c("black","grey"), labels = c("Moving Average", "Observed"), name = "") +
  geom_smooth(data = subset(lof_combo, type == 'avg'), aes(x = m_to_f, y = lof), method='lm') +
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0), limits = c(0,1.6), breaks = c(0,0.5,1,1.5), labels = scales::number_format(accuracy = 0.01)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0.01,1.1)) +
  xlab("Absolute sex difference in mean expression") + ylab("LOF percentile")  +
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18), legend.text = element_text(size=18))

########################
## genetic variance
########################

gv_cont = emma.results[,'vu',]
gv_cont = melt(gv_cont)
gv_cont$m_to_f = 0
m = subset(meta, sex == "m")$LID
f = subset(meta, sex == "f")$LID
for(i in 1:length(gv_cont$Var1)){
  now = resid_exp_all[rownames(resid_exp_all) == gv_cont$Var1[i],colnames(resid_exp_all) %in% subset(meta, Region == gv_cont$Var2[i])$LID]
  dif = abs(mean(now[names(now)%in%m]) - mean(now[names(now)%in%f]))
  gv_cont$m_to_f[i] = dif}
gv_cont = gv_cont[which(!is.na(gv_cont$value)),]
gv_cont = gv_cont[which(!is.na(gv_cont$m_to_f)),]
gv_cont$log_vu = log(gv_cont$value)
hist(gv_cont$log_vu, breaks = 1000)
gv_cont = subset(gv_cont, Var1 %!in% Ychrom)

gv_mov1 = data.frame(seq = seq(from = 0, to = round(max(gv_cont$m_to_f)+0.1,1), by = 0.1))
for(i in 1:length(gv_mov1$seq)){
  gv_mov1$avg_log[i] = mean(log(subset(gv_cont, m_to_f >= gv_mov1$seq[i] & m_to_f < gv_mov1$seq[i+1] & log_vu > -15)$value))}
colnames(gv_mov1) = c("m_to_f","log_vu")
summary(lm(log_vu ~ m_to_f, data = gv_mov1[!is.na(gv_mov1$log_vu),]))
corsub1 = subset(gv_cont, log_vu > -15)
cor.test(x=corsub1$m_to_f, y=corsub1$log_vu, method = 'spearman') 

gv_mov2 = data.frame(seq = seq(from = 0, to = round(max(gv_cont$m_to_f)+0.1,1), by = 0.1))
for(i in 1:length(gv_mov2$seq)){
  gv_mov2$avg_log[i] = mean(log(subset(gv_cont, m_to_f >= gv_mov2$seq[i] & m_to_f < gv_mov2$seq[i+1] & log_vu < -15)$value))}
colnames(gv_mov2) = c("m_to_f","log_vu")
summary(lm(log_vu ~ m_to_f, data = gv_mov2[!is.na(gv_mov2$log_vu),]))
corsub2 = subset(gv_cont, log_vu < -15)
cor.test(x=corsub2$m_to_f, y=corsub2$log_vu, method = 'spearman') # 0.3024174 < 2.2e-16

gv_combo = rbind(cbind(gv_mov1,type="avg1") , cbind(gv_mov2,type="avg2"), cbind(log_vu = log(gv_cont$value), m_to_f = gv_cont$m_to_f, type = "obs"))
gv_combo$log_vu = as.numeric(gv_combo$log_vu)
gv_combo$m_to_f = as.numeric(gv_combo$m_to_f)

ggplot(gv_combo) +
  geom_point(aes(x = m_to_f, y = log_vu, color = type)) +
  geom_point(data = subset(gv_combo, type == 'avg1'), aes(x = m_to_f, y = log_vu, color = type)) +
  geom_point(data = subset(gv_combo, type == 'avg2'), aes(x = m_to_f, y = log_vu, color = type)) +
  scale_color_manual(values = c("black","red","grey"), labels = c("Upper (moving avg)", "Lower (moving avg)", "Observed"), name = "") +
  geom_smooth(data = subset(gv_combo, type == 'avg1'), aes(x = m_to_f, y = log_vu), method='lm') +
  geom_smooth(data = subset(gv_combo, type == 'avg2'), aes(x = m_to_f, y = log_vu), method='lm') +
  theme_classic() + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0), limits=c(-25,5), breaks = c(-25,-20,-15,-10,-5,0)) +
  xlab("Absolute sex difference in mean expression") + ylab("Log genetic variance") +
  theme(legend.position = 'bottom', legend.direction = "horizontal") +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18), legend.text = element_text(size=18))
