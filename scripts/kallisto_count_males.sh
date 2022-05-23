#!/bin/bash

index=genomes/Mmul_10.cdna_noPAR.idx

mkdir -p kallisto

sample=`sed -n ${SLURM_ARRAY_TASK_ID}p male_ids.txt`

kallisto quant -i $index -t $SLURM_CPUS_ON_NODE -o kallisto/${sample} fastq/${sample}.R1.fastq.gz fastq/${sample}.R2.fastq.gz
