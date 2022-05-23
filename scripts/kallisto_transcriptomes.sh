#!/bin/bash

mkdir -p genomes

# Fetch fasta from Ensembl
wget -O genomes/Mmul_10.cdna.fa.gz \
	ftp://ftp.ensembl.org/pub/release-99/fasta/macaca_mulatta/cdna/Macaca_mulatta.Mmul_10.cdna.all.fa.gz
  
# Unzip fasta
gunzip genomes/Mmul_10.cdna.fa.gz 

# Index reference sequence 
samtools faidx genomes/Mmul_10.cdna.fa

# Get all non-Y transcripts (for females)
grep ENSMMUT genomes/Mmul_10.cdna.fa | grep -v "primary_assembly:Mmul_10:Y" | cut -f1 -d " " | sed -e "s/>//g" > genomes/transcripts_noY

# Generate new fasta file for females
samtools faidx -r genomes/transcripts_noY genomes/Mmul_10.cdna.fa > genomes/Mmul_10.cdna_noY.fa

# Index female transcriptome for kallisto
kallisto index -i genomes/Mmul_10.cdna_noY.idx genomes/Mmul_10.cdna_noY.fa

# Get all non-PAR transcripts (for males)
grep ENSMMUT genomes/Mmul_10.cdna.fa | grep -v "CD99 antigen" | cut -f1 -d " " | sed -e "s/>//g" > genomes/transcripts_noPAR

# Generate new fasta file for males
samtools faidx -r genomes/transcripts_noPAR genomes/Mmul_10.cdna.fa > genomes/Mmul_10.cdna_noPAR.fa

# Index male transcriptome for kallisto
kallisto index -i genomes/Mmul_10.cdna_noPAR.idx genomes/Mmul_10.cdna_noPAR.fa
