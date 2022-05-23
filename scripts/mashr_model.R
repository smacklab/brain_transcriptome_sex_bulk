#!/usr/bin/env Rscript

predictor='sexm'

########################
## choose data to upload
########################

# gene level
emma.results = readRDS('checkpoints/emma_results.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')

# cell type corrected gene level
emma.results = readRDS('checkpoints/emma_results_cell_type.rds')
keep.genes = readRDS('checkpoints/keep_genes.rds')

#############
## run models
#############

Bhat = emma.results[,paste('beta',predictor,sep='.'),]
Shat = sqrt(emma.results[,paste('bvar',predictor,sep='.'),])

# For missing data, set beta to 0 and standard error to 100
# See https://github.com/stephenslab/mashr/issues/17
Bhat[is.na(Bhat)] = 0
Shat[is.na(Shat)] = 100

library(mashr)

# Create the mashr data object
mash.data = mash_set_data(Bhat,Shat)

# select strong signals

m.1by1 = mash_1by1(mash.data)
strong.subset = get_significant_results(m.1by1, thresh = 0.05)

# # The code below is the preferred method (computes covariance matrix on a null dataset to learn correlation structure among null tests)

# Get random subset (random choose half of all genes)
random.subset = sample(1:nrow(Bhat),ceiling(nrow(Bhat)/2))

# Set temporary objects in order to estimate null correlation structure
temp = mash_set_data(Bhat[random.subset,],Shat[random.subset,])
temp.U.c = cov_canonical(temp)
Vhat = estimate_null_correlation(temp,temp.U.c)

mash.random = mash_set_data(Bhat[random.subset,],Shat[random.subset,],V=Vhat)
mash.strong = mash_set_data(Bhat[strong.subset,],Shat[strong.subset,], V=Vhat)

# Perform PCA and extreme deconvolution to obtain data-driven covariances
U.pca = cov_pca(mash.strong,5)
U.ed = cov_ed(mash.strong, U.pca)

# fit mash to the random tests using both data-driven and canonical covariances
# estimates mixture proportions
U.c = cov_canonical(mash.random)
m.r = mash(mash.random, Ulist = c(U.ed,U.c), outputlevel = 1)

# fit mash to all data using the above mash fit
m = mash(mash.data, g=get_fitted_g(m.r), fixg=TRUE)

#saveRDS(m,file='checkpoints/mashr_results.rds')
#saveRDS(m,file='checkpoints/mashr_results_cell_type.rds')

print("done")

