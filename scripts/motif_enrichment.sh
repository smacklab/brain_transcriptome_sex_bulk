#!/bin/bash

module load contrib/homer/4.10

findMotifs.pl list-sex-biased-genes.txt mmul10 output_directory -start -1000 -end 300 -bg list-all-expressed-genes.txt -mset vertebrates -nomotif -nogo

findMotifs.pl list-sex-biased-genes-cell-type-corrected.txt mmul10 output_directory -start -1000 -end 300 -bg list-all-expressed-genes.txt -mset vertebrates -nomotif -nogo
