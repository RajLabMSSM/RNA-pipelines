

differential-pipeline
---------------------

Jack Humphrey 2019

A collection of workflows to perform differential analyses on RNA-seq data, wrappped into a Snakemake pipeline.

## Supported workflows

*Differential gene expression (DGE)*

* RSEM -> Limma Voom

* Kallisto -> Sleuth

*Differential transcript expression (DTE)*

* RSEM -> Limma Voom

* Kallisto -> Sleuth

*Differential transcript usage*

* RSEM -> Limma diffSplice


## Dependencies

snakemake - see below to install through conda

R packages:

* optparse
* readr
* tximport
* limma
* edgeR
* dplyr
* ggplot2
* ggrepel
* UpSetR
* janitor
* stringr
* sleuth
* DT
* knitr

Most of these packages should be installed on R/3.6.0 already.

## Installation

```
git clone git@github.com:RajLabMSSM/differential-pipeline

conda create -n differential-pipeline snakemake
```

## Running on example data

It is recommended you run this first to verify that the pipeline works.

```
cd differential-pipeline # enter directory
conda activate differential-pipeline
snakemake -s Snakefile --configfile example/config.yaml
```

## Running on your own data

Two input files are needed:

### *samples.tsv*

This table specifies the samples to be used in the analysis, with each sample as a row.
This should be tab-delimited and have at least three mandatory columns:

1. sample - the unique sample ID

2. rapid_path - the absolute path on Chimera to the sample's RAPiD output folder. example: `/sc/orga/projects/als-omics/TDP-43_knockdown_data/Hammell/RAPiD-nf-0.3.0/tdp_1` 

3. condition - which of the two conditions a sample belongs to.

### *config.yaml*

This specifies the options used by the Snakemake pipeline. Important options:

* dataCode = "my_experiment" - the experiment name to prepend to all results files
* outFolder = "results" - the subfolder to put the results in
* refCondition = "control" - the name of the reference condition
* altCondition = "case" - the name of the alternate condition

## HTML Report

An optional module takes the output of all differential analyses and creates an HTML report using Rmarkdown

An example HTML report can be seen [here](https://rajlabmssm.github.io/differential-pipeline/report.html)


