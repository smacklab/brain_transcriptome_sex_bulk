#!/usr/bin/env Rscript

###############################
## TO BE RUN AFTER emma_model.R
###############################

##load data

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')

keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')
meta = meta[order(match(meta$LID, colnames(e.keep))),]

################
## get residuals
################

# Put together k and z matrices

animals = sort(unique(meta$Individual))

# Calculate k matrix (pairwise relatedness)

k = matrix(0,nrow=length(animals),ncol=length(animals),dimnames=list(animals,animals))

kinship = read.delim('cayo_brain_lcmlkin_results.txt',stringsAsFactors=FALSE)
kinship = kinship[kinship$Ind1 %in% animals & kinship$Ind2 %in% animals,]
i = as.matrix(kinship[,c('Ind1','Ind2')])
k[i[,1:2]] = kinship$pi_HAT
k[i[,2:1]] = kinship$pi_HAT
diag(k) = 1

# Calculate z matrix (matching libraries to genotype)
z = matrix(0,nrow=nrow(meta),ncol=nrow(k))
rownames(z) = meta$LID
colnames(z) = rownames(k)
i = as.matrix(meta[,c('LID','Individual')])
z[i] = 1

library(parallel)
library(doParallel)

# Within region

regions = names(keep.genes)

# Initialize output
EMMA_RNA_random_effects = vector('list',length(regions))
names(EMMA_RNA_random_effects) = regions
n.cores = detectCores() - 4

## generate uhat matrix of random effects

for (i in 1:length(regions)) {
  message('Now analyzing region ',regions[i])
  m = droplevels(subset(meta,Region %in% regions[i]))
  z.this = z[rownames(z) %in% m$LID,colnames(z) %in% m$Individual]
  e.this = e.keep[keep.genes[[regions[i]]],colnames(e.keep) %in% m$LID]
  k.this = k[rownames(k) %in% m$Individual, colnames(k) %in% m$Individual]
  
  c.this = model.covariates[apply(m[,model.covariates],2,function(x) length(unique(x))) > 1]
  design = model.matrix(as.formula(paste('~',paste(c.this,collapse=' + '))), data=m)
  clus = parallel::makeCluster(n.cores, setup_strategy = "sequential")
  registerDoParallel(cores=n.cores)  
  clusterExport(clus,varlist=c('e.this','k.this','z.this','design'),envir=environment())
  
  EMMA_RNA_random_effects[[regions[i]]] = t(parApply(clus,e.this,1,function(y) {
    require(EMMREML)
    
    emma=emmreml(y = y,X = design,Z = z.this,K = k.this,varbetahat = T,varuhat = T,PEVuhat = T,test = T)
    
    return(t(emma$uhat))
  }))
  stopCluster(clus)
  colnames(EMMA_RNA_random_effects[[i]])=colnames(k.this)
}

## just remove everything BUT the sex effects

reduced_expression_matrix = list()
design_mat = list()

for (i in 1:length(regions)){
  m = droplevels(subset(meta,Region %in% regions[i]))
  z.this = z[rownames(z) %in% m$LID,]
  e.this = e.keep[keep.genes[[regions[i]]],colnames(e.keep) %in% m$LID]
  c.this = model.covariates[apply(m[,model.covariates],2,function(x) length(unique(x))) > 1]
  
  design = model.matrix(as.formula(paste('~',paste(c.this,collapse=' + '))), data=m)
  
  fullmodel = data.frame(out[[i]])
  reduced_expression_matrix[[i]] = e.this -
    fullmodel$beta.exact_age_years%*%t(design[,"exact_age_years"]) -
    fullmodel$beta.RIN%*%t(design[,"RIN"]) -
    fullmodel$beta.ordinal.rank.L%*%t(design[,"ordinal.rank.L"]) -
    fullmodel$beta.ordinal.rank.Q%*%t(design[,"ordinal.rank.Q"]) -
    if("Library.batch2018-06-09" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2018.06.09%*%t(design[,"Library.batch2018-06-09"])} else {0} -
    if("Library.batch2018-06-12" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2018.06.12%*%t(design[,"Library.batch2018-06-12"])} else {0} -
    if("Library.batch2018-06-20" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2018.06.20%*%t(design[,"Library.batch2018-06-20"])} else {0} -
    if("Library.batch2018-07-09" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2018.07.09%*%t(design[,"Library.batch2018-07-09"])} else {0} -
    if("Library.batch2018-07-11" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2018.07.11%*%t(design[,"Library.batch2018-07-11"])} else {0} -
    if("Library.batch2018-07-21" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2018.07.21%*%t(design[,"Library.batch2018-07-21"])} else {0} -
    if("Library.batch2018-07-25" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2018.07.25%*%t(design[,"Library.batch2018-07-25"])} else {0} -
    if("Library.batch2019-09-30" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2019.09.30%*%t(design[,"Library.batch2019-09-30"])} else {0} -
    if("Library.batch2019-10-01" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2019.10.01%*%t(design[,"Library.batch2019-10-01"])} else {0} -
    if("Library.batch2019-10-03" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2019.10.03%*%t(design[,"Library.batch2019-10-03"])} else {0} -
    if("Library.batch2019-10-04" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2019.10.04%*%t(design[,"Library.batch2019-10-04"])} else {0} -
    if("Library.batch2019-10-17" %in% colnames(design) == TRUE){fullmodel$beta.Library.batch2019.10.17%*%t(design[,"Library.batch2019-10-17"])} else {0}
  design_mat[[i]] = model.matrix(~0+Individual,data=m)
  colnames(design_mat[[i]])=gsub("Individual","",colnames(design_mat[[i]]))
  design_mat[[i]]=as.matrix(design_mat[[i]][,colnames(EMMA_RNA_random_effects[[i]])])
  design_mat[[i]]=design_mat[[i]][rownames(m),]
  
}

## remove random (individual) effects

for (i in 1:length(regions)){
  for (x in 1:ncol(design_mat[[i]])) {
    reduced_expression_matrix[[i]]=reduced_expression_matrix[[i]]-EMMA_RNA_random_effects[[i]][,x]%*%t(design_mat[[i]][,x])}}

names(reduced_expression_matrix) = regions

saveRDS(reduced_expression_matrix, file='checkpoints/residual_expression.rds')

