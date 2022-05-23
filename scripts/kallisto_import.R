#!/usr/bin/env Rscript

library(tximport)
library(rhdf5)
library(biomaRt)
library(stringr)

samplesF=scan('female_ids.txt', what = "")
samplesM=scan('male_ids.txt', what = "")
samples=c(samplesF, samplesM)

files = file.path('kallisto',samples,'abundance.h5')
names(files) = samples

filesF = files[names(files) %in% samplesF]
filesM = files[names(files) %in% samplesM]

## this imports the hd5 kallisto-mapped files
txi.kallistoF = tximport(filesF, type = 'kallisto', txOut = TRUE)
txi.kallistoM = tximport(filesM, type = 'kallisto', txOut = TRUE)

txi.kallisto = list(abundance=as.matrix(merge(txi.kallistoF$abundance,txi.kallistoM$abundance,by="row.names",all=TRUE)),counts=as.matrix(merge(txi.kallistoF$counts,txi.kallistoM$counts,by="row.names",all=TRUE)),length=as.matrix(merge(txi.kallistoF$length,txi.kallistoM$length,by="row.names",all=TRUE)),countsFromAbundance=c("no"))
rownames(txi.kallisto$abundance) = txi.kallisto$abundance[,1]
rownames(txi.kallisto$counts) = txi.kallisto$counts[,1]
rownames(txi.kallisto$length) = txi.kallisto$length[,1]
txi.kallisto$abundance = txi.kallisto$abundance[,-1]
txi.kallisto$counts = txi.kallisto$counts[,-1]
txi.kallisto$length = txi.kallisto$length[,-1]
txi.kallisto$abundance[is.na(txi.kallisto$abundance)] = 0
txi.kallisto$counts[is.na(txi.kallisto$counts)] = 0
txi.kallisto$length[is.na(txi.kallisto$length)] = 0
class(txi.kallisto$abundance) = 'numeric'
class(txi.kallisto$counts) = 'numeric'
class(txi.kallisto$length) = 'numeric'

## convert to counts per gene
mmul = useMart(biomart = 'ENSEMBL_MART_ENSEMBL',dataset='mmulatta_gene_ensembl') 

## get a file that matches ensembl transcripts (with version appended) to ensembl genes
tx2gene.chrom = getBM(attributes=c('ensembl_transcript_id_version', 'ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'), filters = 'ensembl_transcript_id_version', values = rownames(txi.kallisto$abundance), mart = mmul)

## summarize kallisto mapped data to the gene level
txi.gene = summarizeToGene(txi.kallisto, tx2gene.chrom)

saveRDS(txi.gene, file = 'checkpoints/kallisto_genes.rds')
