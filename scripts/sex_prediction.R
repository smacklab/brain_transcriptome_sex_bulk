#!/usr/bin/env Rscript

#######################
# load data and results
#######################

library(plyr)

resid_exp = readRDS(file='checkpoints/residual_expression.rds')

keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')

for(i in 1:length(resid_exp)){
  resid_exp[[i]] = data.frame(resid_exp[[i]])
  resid_exp[[i]]$ROWNAMES  <- rownames(resid_exp[[i]])}
resid_exp_all <- join_all( resid_exp, by="ROWNAMES", type="full" )
rownames(resid_exp_all) <- resid_exp_all$ROWNAMES; resid_exp_all$ROWNAMES <- NULL

meta = meta[order(match(meta$LID, colnames(resid_exp_all))),]

`%!in%` = Negate(`%in%`)

#################
## structure data
#################

datanow = data.frame(t(resid_exp_all))
table(rownames(datanow) == meta$LID)
datanow2 = datanow
datanow2$sex = meta$sex
outcomeName <- 'sex'

##################
## select gene set
##################

library(biomaRt)

mmul = useEnsembl(biomart = 'ensembl',dataset='mmulatta_gene_ensembl',mirror = "www") 
gene2chrom = getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name'), filters = 'ensembl_gene_id', values = rownames(resid_exp_all), mart = mmul)
Ychrom = subset(gene2chrom, chromosome_name == "Y")
Xchrom = subset(gene2chrom, chromosome_name == "X")

# X + autosomes
datanow3 = datanow2[,which(colnames(datanow2) %!in% Ychrom$ensembl_gene_id)]
predictorsNames <- names(datanow3)[names(datanow3) != outcomeName]
genelist = 'X_and_autosomes'

# autosomes only
datanow3 = datanow2[,which(colnames(datanow2) %!in% Ychrom$ensembl_gene_id)]
datanow3 = datanow3[,which(colnames(datanow3) %!in% Xchrom$ensembl_gene_id)]
predictorsNames <- names(datanow3)[names(datanow3) != outcomeName]
genelist = 'autosomes_only'

# X only
datanow3 = datanow2[,which(colnames(datanow2) %in% Xchrom$ensembl_gene_id)]
datanow3$sex = datanow2$sex
predictorsNames <- names(datanow3)[names(datanow3) != outcomeName]
genelist = 'X_only'

###################
## LOOCV per region
###################

library(caret)
library(pROC)
library(gbm)

objControl <- trainControl(method='LOOCV', classProbs = TRUE, summaryFunction = twoClassSummary)

gbmGrid <-  expand.grid(interaction.depth = c(1, 3, 5, 9), 
                        n.trees = (1:5)*50, 
                        shrinkage = 0.1,
                        n.minobsinnode = 5)

mods_out = list()
predictions_out = list()

for (i in 1:length(region.levels)){
  print(region.levels[i])
  regiondata = datanow3[subset(meta, Region == region.levels[i])$LID,]
  regiondata = regiondata[,colSums(is.na(regiondata))<nrow(regiondata)]
  predictorsNames = names(regiondata)[names(regiondata) != outcomeName]
  
  objModel <- train(sex ~ ., 
                    data = regiondata, 
                    method='gbm', 
                    trControl=objControl,
                    tuneGrid = gbmGrid,
                    metric = "ROC")
  
  mods_out[[i]] = objModel
}

names(mods_out) = region.levels

saveRDS(mods_out, file = paste('checkpoints/gbm_loocv_',genelist,'.rds',sep=""))

############################################
## investigate outputs for individual models
############################################

# select gene list

genelist = 'X_and_autosomes'
genelist = 'X_only'
genelist = 'autosomes_only'

mods_out = readRDS(paste('checkpoints/gbm_loocv_',genelist,'.rds',sep=""))

library(tibble)
library(caret)

gbmImp = list()
preds = list()
mats = list()
stats = data.frame(V1 = vector(mode = 'numeric', length = 19))
p1 = list()
p2 = list()

for (i in 1:length(region.levels)){
  
  print(region.levels[i])
  print(mods_out[[i]])
  
  prednow = data.frame(mods_out[[i]]$pred)
  depthnow = as.numeric(mods_out[[i]]$bestTune[2])
  treesnow = as.numeric(mods_out[[i]]$bestTune[1])
  prednow = subset(prednow, interaction.depth == depthnow & n.trees == treesnow)
  prednow$LID = subset(meta, Region == region.levels[i])$LID 
  prednow$Individual = subset(meta, Region == region.levels[i])$Individual
  preds[[i]] = prednow
  
  impnow = varImp(mods_out[[i]], scale = TRUE, useModel = TRUE)
  gbmImp[[i]] = data.frame(impnow$importance)
  gbmImp[[i]]$gene = rownames(gbmImp[[i]])
  rownames(gbmImp[[i]]) = NULL
  p1[[i]] = plot(impnow, top = 40)
  
  con.mat = confusionMatrix(preds[[i]][,1], preds[[i]][,2])
  p2[[i]] = fourfoldplot(con.mat$table, conf.level = 0, margin = 1, main = "Confusion Matrix")
  mats[[i]] = con.mat$table
  ROCnow =  subset(mods_out[[i]]$results, interaction.depth == depthnow & n.trees == treesnow)$ROC
  names(ROCnow) = 'ROC'
  stats[,i] = data.frame(c(con.mat$overall, con.mat$byClass, ROCnow))[,1]
}

names(gbmImp) = names(preds) = names(mats) = colnames(stats) = region.levels
rownames(stats) = names(c(con.mat$overall,con.mat$byClass, ROCnow))

## model fit statistics

mean(as.numeric(stats['ROC',]))
min(stats['ROC',])
max(stats['ROC',])
mean(as.numeric(stats['Accuracy',]))
min(stats['Accuracy',])
max(stats['Accuracy',])
mean(as.numeric(stats['Specificity',]))
min(stats['Specificity',])
max(stats['Specificity',])
mean(as.numeric(stats['Sensitivity',]))
min(stats['Sensitivity',])
max(stats['Sensitivity',])

# number of influential genes

inf_genes = data.frame(count=c(1:15))
for (i in 1:length(region.levels)){
  inf_genes$count[i] = length(subset(gbmImp[[region.levels[i]]], Overall > 0)$gene)
}
inf_genes$count
mean(inf_genes$count)

# prediction probabilities (compare individuals, regions)

library(reshape2)
library(tidyr)
library(dplyr)
library(data.table)

preds.all = rbindlist(preds, use.names=TRUE, fill=TRUE, idcol="region")
preds.all$Individual = with(preds.all, reorder(Individual,f,median))
preds.all$abs_pred = ifelse(preds.all$obs == 'f', preds.all$f, preds.all$m)
preds.all$match = ifelse(preds.all$pred == preds.all$obs, 'black', 'red')
preds.all$region = factor(preds.all$region, levels = region.levels)
ggplot(preds.all, aes(x=Individual, y=abs_pred, fill = obs)) +
  scale_fill_manual(values=region.colors[c(3,6)], breaks = c("f","m"), labels = c("F","M")) +
  geom_boxplot() + theme_classic() + ylab('Prediction Probability') +
  geom_point(size = 2, color=preds.all$match) + theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.text.x = element_blank(), axis.title = element_text(size=18), axis.text = element_text(size=18), legend.text = element_text(size=18), legend.title = element_blank())
ggplot(preds.all, aes(x=region, y=abs_pred)) +
  scale_fill_manual(values=c('#666666')) +
  geom_boxplot() + theme_classic() + ylab('Prediction Probability') +
  geom_point(size = 2, color=preds.all$match) + theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18), axis.title.x = element_blank())

# investigate variability

metanow = unique(meta[,c('Individual','exact_age_years','rank.scaled','ordinal.rank','sex')])

m = preds.all
#m = data.frame(m %>% group_by(Individual) %>% summarise(mean_pred = mean(abs_pred)))
#m = data.frame(m %>% group_by(Individual) %>% summarise(mean_pred = sd(abs_pred)))
m = merge(m, metanow, by = 'Individual')

# age effects
mod = lm(mean_pred ~ exact_age_years, data = m)
summary(mod)
mod = lm(mean_pred ~ exact_age_years, data = subset(m, sex=='f'))
summary(mod)
mod = lm(mean_pred ~ exact_age_years, data = subset(m, sex=='m'))
summary(mod)
ggplot(m, aes(x=exact_age_years, y=mean_pred, color=sex)) + geom_point() + 
  theme_classic() + geom_smooth(method='lm') +
  xlab('Age (years)') + 
  #ylab('Mean prediction probability') +
  ylab('Prediction probability variation (sd)') +
  scale_color_manual(values=region.colors[c(3,6)], breaks = c('f','m'), labels = c('F','M')) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18), legend.text = element_text(size=18), legend.title = element_blank())

# rank effects
mod = lm(mean_pred ~ rank.scaled, data = m2)
summary(mod)
mod = lm(mean_pred ~ rank.scaled, data = subset(m2, sex=='f'))
summary(mod)
mod = lm(mean_pred ~ rank.scaled, data = subset(m2, sex=='m'))
summary(mod)
ggplot(m2, aes(x=rank.scaled, y=mean_pred, color=sex)) + geom_point() + 
  theme_classic() + geom_smooth(method='lm') +
  xlab('Age (years)') + 
  ylab('Mean prediction probability') +
  scale_color_manual(values=region.colors[c(3,6)], breaks = c('f','m'), labels = c('F','M')) +
  theme(axis.text = element_text(size=18), axis.title = element_text(size=18), legend.text = element_text(size=18), legend.title = element_blank())

# gene importance (compare genes, regions)

library(expss)
library(reshape2)

gbmImp.all = melt(gbmImp)
gbmImp.all$L1 = factor(gbmImp.all$L1, levels = region.levels)
table(subset(gbmImp.all, value>0)$L1)

library(expss)

gbmImp.all2 = gbmImp.all %>% spread(L1, value)
gbmImp.all2[gbmImp.all2 == 0] <- NA
gbmImp.all2$Sum = rowSums(gbmImp.all2[,3:17], na.rm = TRUE)
gbmImp.all2$Mean = rowMeans(gbmImp.all2[,3:17], na.rm = TRUE)
for (i in 1:length(gbmImp.all2$gene)) {gbmImp.all2$chrom[i] = vlookup(gbmImp.all2$gene[i], 
  dict = gene2chrom, lookup_column = 1, result_column = 3)}

nonzero = subset(gbmImp.all2, Sum>0)
dim(nonzero)[1]
chroms = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","X","Y")
nonzero2 = subset(nonzero, chrom %in% chroms)
chrom_means = data.frame(nonzero %>% group_by(chrom) %>% summarise(mean_sum = mean(Sum)))
chrom_means
nonzero$chrom2 = ifelse(nonzero$chrom == "X","X","Autosome")
table(nonzero$chrom2)
chrom_means2 = data.frame(nonzero %>% group_by(chrom2) %>% summarise(mean_sum = mean(Sum)))
chrom_means2

nonzero$count = rowSums(nonzero[3:17] > 0, na.rm = TRUE)
nonzero = nonzero[order(nonzero$count),]
table(nonzero$count)/length(nonzero$count)
table(nonzero$chrom2,nonzero$count)

unique.infl = subset(nonzero, count == 1)$gene
unique.exp = names(table(melt(keep.genes)$value)[which(table(melt(keep.genes)$value)>=13)])
unique.combo = unique.infl[which(unique.infl %in% unique.exp)]
length(unique.combo)
length(unique.combo)/length(unique.infl)

# compare X chromosome model gene importance to humans (from Oliva et al. 2020)

oliva = read.csv('data/Oliva_data.csv')
orths = readRDS('checkpoints/human_macaque_one2one.rds')

for (i in 1:length(nonzero$gene)){
  nonzero$ENSEMBL_gene_id[i] = vlookup(nonzero$gene[i], dict = orths, lookup_column = 3, result_column = 1)}

comp = merge(oliva, nonzero, by = 'ENSEMBL_gene_id')
comp = comp[which(!is.na(comp$Sum.of.relative.sex.predictivity)),]

cor.test(x = log10(comp$Sum.of.relative.sex.predictivity), y = log10(comp$Sum), method = 'spearman')
ggplot(data = comp, aes(x = log10(Sum), y = log10(Sum.of.relative.sex.predictivity))) + geom_point() + 
  geom_smooth(method = 'lm') + theme_classic() +
  ylab('Relative Importance (log10) in Oliva et al.') + xlab('Relative Importance (log10) in Current Study') +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18))

cor.test(x = log10(comp$Avg.of.relative.sex.predictivity...100.), y = log10(comp$Mean), method = 'spearman')
ggplot(data = comp, aes(x = log10(Mean), y = log10(Avg.of.relative.sex.predictivity...100.))) + geom_point() + 
  geom_smooth(method = 'lm') + theme_classic() +
  ylab('Relative Importance (log10) in Oliva et al.') + xlab('Relative Importance (log10) in Current Study') +
  theme(axis.title = element_text(size=18), axis.text = element_text(size=18))

