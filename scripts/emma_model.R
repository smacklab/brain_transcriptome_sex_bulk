#!/usr/bin/env Rscript

########################
## choose data to upload and analyze
########################

# gene level data (primary analyses)

e.keep = readRDS('checkpoints/filtered_expression_matrix.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')
meta = meta[order(match(meta$LID, colnames(e.keep))),]

# cell type corrected gene level data

e.keep = readRDS('checkpoints/adjusted_filtered_expression_matrix.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')
meta = readRDS('checkpoints/cayo_bulkbrain_combined_metadata.rds')
meta = meta[order(match(meta$LID, colnames(e.keep))),]
# change model covariates in _input_options.R file
# model.covariates = c('exact_age_years','sex','ordinal.rank')

#############
## run models
#############

# Put together k and z matrices

animals = sort(unique(meta$Individual))

# Calculate k matrix (pairwise relatedness)

k = matrix(0,nrow=length(animals),ncol=length(animals),dimnames=list(animals,animals))

kinship = read.delim('cayo_brain_lcmlkin_results.txt',stringsAsFactors=FALSE)
i = as.matrix(kinship[,c('Ind1','Ind2')])
k[i[,1:2]] = kinship$pi_HAT
k[i[,2:1]] = kinship$pi_HAT
diag(k) = 1

k2 = matrix(0,nrow=length(animals),ncol=length(animals),dimnames=list(animals,animals))
diag(k2) = 1

z = matrix(0,nrow=nrow(meta),ncol=nrow(k))
rownames(z) = meta$LID
colnames(z) = rownames(k)
i = as.matrix(meta[,c('LID','Individual')])
z[i] = 1

library(parallel)
library(doParallel)

regions = names(keep.genes)

# Initialize output
out = vector('list',length(regions))
names(out) = regions

n.cores = detectCores() - 4

# Design model covariates

for (i in 1:length(regions)) {
  message('Now analyzing region ',regions[i])
  m = droplevels(subset(meta,Region %in% regions[i]))
  z.this = z[rownames(z) %in% m$LID,]
  e.this = e.keep[keep.genes[[regions[i]]],colnames(e.keep) %in% m$LID]

  # Drop covariates if they are uniform across dataset
  c.this = model.covariates[apply(m[,model.covariates],2,function(x) length(unique(x))) > 1]

  design = model.matrix(as.formula(paste('~',paste(c.this,collapse=' + '))), data=m)
  
  clus = parallel::makeCluster(n.cores, setup_timeout = 0.5)
  registerDoParallel(cores=n.cores)  
  clusterExport(clus,varlist=c('e.this','k','k2','z.this','design'),envir=environment())
  
  out[[regions[i]]] = t(parApply(clus,e.this,1,function(y) {
    require(EMMREML)
    
    emma=emmreml(y = y,X = design,Z = z.this,K = k,varbetahat = T,varuhat = T,PEVuhat = T,test = T)
    
    p = emma$pvalbeta
    varb = emma$varbetahat
    b = emma$betahat
    vu = emma$Vu
    ve = emma$Ve
    
    c(b,varb,p[,"none"],vu,ve)
  }))
  
  stopCluster(clus)
  
  colnames(out[[regions[i]]])[(ncol(design) * 0 + 1):(ncol(design) * 1)] = paste('beta',colnames(design),sep='.')
  colnames(out[[regions[i]]])[(ncol(design) * 1 + 1):(ncol(design) * 2)] = paste('bvar',colnames(design),sep='.')
  colnames(out[[regions[i]]])[(ncol(design) * 2 + 1):(ncol(design) * 3)] = paste('pval',colnames(design),sep='.')
  colnames(out[[regions[i]]])[(ncol(design) * 3 + 1)] = 'vu'
  colnames(out[[regions[i]]])[(ncol(design) * 3 + 2)] = 've'
}

regions.dimnames = list(genes = Reduce(union,lapply(out,rownames)), outputs = Reduce(union,lapply(out,colnames)), regions = regions)
regions.dim = unlist(lapply(regions.dimnames,length))
regions.numeric = numeric(Reduce(`*`,regions.dim))
regions.numeric[!regions.numeric] = NA
regions.array = array(unname(regions.numeric),dim=unname(regions.dim),dimnames=unname(regions.dimnames))

for (i in 1:length(out)) {
  foo = reshape2::melt(out[[i]])
  foo$Var3 = names(out)[i]
  j = as.matrix(foo[,paste0('Var',1:3)])
  regions.array[j] = foo$value
  rm(foo)
}

# saveRDS(regions.array,file='checkpoints/emma_results.rds')
# saveRDS(regions.array,file='checkpoints/emma_results_cell_type.rds')
