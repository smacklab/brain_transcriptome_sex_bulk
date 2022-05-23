#!/usr/bin/env Rscript

# Predictor of interest (set to column name in metadata file)
predictor = 'sexm'
sex.variable = 'sex'

# Predictor of interest (set to preferred label)
predictor.label = 'Sex'

# Transcripts-per-million cutoff
tpm.cutoff = 10

# Cutoffs 
fdr.cutoff = 0.1
fsr.cutoff = 0.05

# Covariates to include in the EMMA models
model.covariates = c('exact_age_years','sex','RIN','ordinal.rank','Library.batch')

# Covariates to include in the EMMA models (cell type corrected)
#model.covariates = c('exact_age_years','sex','ordinal.rank')

# Fraction of regions considered shared
fraction.shared.cutoff = 13/15 

# Fraction of regions considered unique
fraction.unique.cutoff = 1/15 

# Set factor
region.levels = c('dmPFC', 'dlPFC', 'vmPFC', 'vlPFC', 'ACCg', 'M1', 'STS', 'V1', 'AMY', 'CA3', 'DG', 'CN', 'Pu', 'LGN', 'VMH')

# Set colors (extension of palette Dark2)
region.colors = c('#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02', '#a6761d', '#666666', '#ed1c24', '#00aeef', '#86328c', '#00aaad', '#ac7eaf', '#2e3192', '#6a3e14')

# set colors for sex
sex.colors = c("#E69F00", "#0072B2")

# Set shapes
region.shapes = c(21,24,22,16,17,15,21,24,22,16,17,15,21,24,22)

# function
`%!in%` = Negate(`%in%`)
