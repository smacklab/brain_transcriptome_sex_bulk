# Brain bulk-tissue RNA-seq analysis pipeline

**Repository for the brain transcriptomic sex differences project from Cayo Santiago macaques (bulk-tissue RNA-seq)**

This repository contains scripts used in the analysis of sex effects in bulk-tissue RNA-seq data for the Cayo Santiago rhesus macaque population. 

Note that we ran most steps on the University of Washington ([Mox](https://wiki.cac.washington.edu/display/hyakusers/Hyak+mox+Overview)) high-performance computing clusters. We have aimed to generalize the code here by removing system-specific references to installed software and modules. Instead, we document required software and version numbers below (excluding standard Unix programs and R). For HPC systems, the required scripts and binaries must be in the PATH. The easiest way to do this is to use an existing module or to install your own. In these cases, the modules should be loaded prior to running the appropriate code below.

As Mox uses the [slurm](https://slurm.schedmd.com/documentation.html) scheduler, most code below should run on slurm systems with little or no modification. For non-slurm HPC systems, slurm scripts and environmental variables will need to be adjusted, though hopefully without too much hassle.

We ran most analysis steps using [R](https://cran.r-project.org) (v4.1). We recommend the following utility or visualization packages to extend base R's functionality.

| Package                                                              | Description                                        |
| -----------                                                          | -----------                                        |
| [tidyverse](https://www.tidyverse.org/)                              | utilities for data manipulation and visualization  |
| [reshape2](https://cran.r-project.org/web/packages/reshape2)         | data manipulation                                  |
| [abind](https://cran.r-project.org/web/packages/abind)               | combining multidimensional arrays                  |
| [XML](https://cran.r-project.org/web/packages/XML)                   | parsing XML files                                  |
| [jsonlite](https://cran.r-project.org/web/packages/jsonlite)         | parsing JSON files                                 |
| [ggrastr](https://cran.r-project.org/web/packages/ggrastr)           | rasterizing big data visualizations                |
| [ggtext](https://cran.r-project.org/web/packages/ggtext)             | rendering text                                     |
| [ggrepel](https://cran.r-project.org/web/packages/ggrepel)           | avoiding overplotting of labels                    |
| [ggbeeswarm](https://cran.r-project.org/web/packages/ggbeeswarm)     | beeswarm plot support                              |
| [egg](https://cran.r-project.org/web/packages/egg)                   | advanced layouts for visualization                 |
| [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer) | color advice for plots                             |
| [viridis](https://cran.r-project.org/web/packages/viridis)           | perceptually uniform color gradients               |
| [doParallel](https://cran.r-project.org/web/packages/doParallel)     | support for parallel computing                     |
| [future](https://cran.r-project.org/web/packages/future)             | support for parallel computing                     |

More specialized R packages are listed with their specific scripts below.

# Inputs

The following files are expected:

* Demultiplexed pair-end fastq files should be compressed with gzip and placed in the `fastq/` folder with the naming convention `<library ID>.R1.fastq.gz` (read 1) and `<library ID>.R2.fastq.gz` (read 2).

* An animal metadata file should be placed in `data/cayo_brain_bulk_metadata_animals.tsv`

* A library metadata file should be placed in `data/cayo_brain_bulk_metadata_technical.tsv`

* An animal social metrics file should be placed in `data/social_metrics.csv`

* Results for an analysis of sex-biased gene expression in humans (GTEx v8) should be placed in ```data/gtex_mashr_results_sex.rds```

# Pipeline

### Map reads with splice-aware aligner

* ***Required software***: [STAR](https://github.com/alexdobin/STAR) (v2.5), [SAMtools](http://www.htslib.org/) (v1.9), [GATK](https://gatk.broadinstitute.org/) (v4.1.2.0)

```
# Download and index genome
sbatch scripts/star_index.sh

# Map each library using STAR
sbatch --array=1-$(tail -n+2 data/cayo_brain_bulk_metadata_technical.tsv | wc -l | xargs) scripts/star_map.sh
```

### Merge alignments per genotype

* ***Required software***: [SAMtools](http://www.htslib.org/) (v1.9)

```
# Merge bam files for each genotype
sbatch --array=1-$(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs) scripts/samtools_merge.sh
```

### Call and filter genotypes

* ***Required software***: [SAMtools](http://www.htslib.org/) (v1.9), [GATK](https://gatk.broadinstitute.org/) (v4.1.2.0), [VCFtools](https://vcftools.github.io/) (v0.1.16)

```
# Split variants into chromosomes (per genotype)
sbatch --array=1-$(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs) scripts/samtools_split.sh

# Clean variants
sbatch --array=1-$(($(tail -n+2 data/cayo_brain_bulk_metadata_animals.tsv | wc -l | xargs)*20)) scripts/gatk_clean_reads.sh

# Call variants
sbatch --array=1-20 scripts/gatk_call_variants.sh

# Filter variants
sbatch --array=1-20 scripts/gatk_filter_variants.sh

# Concatenate variants across chromosomes
scripts/vcftools_concat.sh
```

### Compute relatedness with lcMLkin

* ***Required software***: [VCFtools](https://vcftools.github.io/) (v0.1.16), [lcMLkin](https://github.com/COMBINE-lab/maximum-likelihood-relatedness-estimation) (v20190218)

```
# Thin variants
scripts/vcftools_thin_variants.sh

# Call kinship
scripts/lcmlkin_call_kinship.sh
```

### Calculate sequencing stats

* ***Required software***: [SAMtools](http://www.htslib.org/) (v1.9), [ea-utils](https://expressionanalysis.github.io/ea-utils/) (v1.04.807), [GNU parallel](https://www.gnu.org/software/parallel/) (v20171122)

```
# Calculate sequencing stats across libraries
sbatch scripts/sequencing_stats_parallel.sh

# Summarize sequencing stats
scripts/sequencing_stats_summarize.sh
```

### Import and clean sample and library metadata

```
# Read, format, and clean animal and library metadata
# Create male and female sample lists
scripts/clean_metadata.R
```

### Quantify transcripts with kallisto

* ***Required software***: [kallisto](https://pachterlab.github.io/kallisto) (v0.43.1)

```
# Create sex specific transcriptomes and index
scripts/kallisto_transcriptomes.sh

# Count transcripts
sbatch --array=1-$(wc -l checkpoints/male_ids.txt | cut -d ' ' -f 1) scripts/kallisto_count_males.sh
sbatch --array=1-$(wc -l checkpoints/female_ids.txt | cut -d ' ' -f 1) scripts/kallisto_count_females.sh
```

### Import expression data

* ***Key libraries***: [tximport](https://doi.org/doi:10.18129/B9.bioc.tximport), [rdf5](https://doi.org/doi:10.18129/B9.bioc.rhdf5), [biomaRt](https://doi.org/doi:10.18129/B9.bioc.biomaRt)

```
# Import kallisto results into R
scripts/kallisto_import.R
```

### Filter expression matrix

* ***Key libraries***: [biomaRt](https://doi.org/doi:10.18129/B9.bioc.biomaRt), [limma](https://doi.org/doi:10.18129/B9.bioc.limma)

```
# Apply filters to gene expression dataset
scripts/filter_expression.R
```

### Visualize data (pre-modeling)

* ***Key libraries***: [variancePartition](https://doi.org/doi:10.18129/B9.bioc.variancePartition), [umap](https://cran.r-project.org/web/packages/umap), [dendextend](https://cran.r-project.org/web/packages/dendextend), [ape](https://cran.r-project.org/web/packages/ape), [phangorn](https://cran.r-project.org/web/packages/phangorn)

```
# Visualize data
scripts/visualize_expression.R
```

### Fit linear mixed model(s)

* ***Key libraries***: [EMMREML](https://cran.r-project.org/web/packages/EMMREML)

```
# Fit linear mixed effect model
scripts/emma_model.R

# Get residual expression values (for later analyses)
scripts/residual_expression
```

### Apply adaptive shrinkage

* ***Key libraries***: [mashr](https://cran.r-project.org/web/packages/mashr)

```
# Refine effects and significance values with adaptive shrinkage
scripts/mashr_model.R
```

### Cell type enrichment

* **Key libraries:** BRETIGEA, biomaRt

```
# Calculate cell type enrichment for sex-biased genes in macaques (and human GTEx data)
scripts/cell_type_enrichment.R
```

### Adjust macaque expression data for cell type proportions 

* **Key libraries:** BRETIGEA, biomaRt

```
# Fit linear mixed effect model
scripts/adjust_expression.R
```

### Fit linear mixed model(s) to cell type corrected data

* **Key libraries:** EMMREML

```
# Fit linear mixed effect model
scripts/emma_model.R
```

### Apply adaptive shrinkage to cell type corrected data

* **Key libraries:** mashr

```
# Refine effects and significance values with adaptive shrinkage
scripts/mashr_model.R
```

### Visualize and describe model results

* **Key libraries:** biomaRt, ggplot2

```
# Visualize and describe sex-biased gene distributions
# Includes chromosome enrichment analysis
scripts/visualize_model_results.R
```

### Disease, motif, and functional enrichment analyses

* **Required software**: Homer (v4.10)
* **Key libraries:** ViSEAGO, mashr, biomaRt

```
# Perform disease, motif, and functional enrichment analyses for sex-biased genes in macaques and humans
scripts/risk_gene_and_functional_enrichment.R
scripts/ASD_expression_enrichment.R
scripts/motif_enrichment.sh
```

### Human vs. rhesus macaque comparisons

* **Key libraries:** biomaRt

```
# Compare sex effects in humans and macaques
# Compare cell type enrichment and disease enrichment results in humans and macaques
scripts/gtex_comparison.R
```

### Sex prediction

* **Key libraries:** caret

```
# Fit gradient boosted models per region in macaques
scripts/sex_prediction.R
```

### Identify potential drivers of sex-biased gene expression

```
# Correlate sex-biased expression in macaques with tissue specificity, loss of function, and genetic variance for expression
scripts/evolutionary_mechanisms.R
```
