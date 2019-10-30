# qSVA-pipeline

Jack Humphrey

Snakemake pipeline for extracting degradation matrices from a set of BAM files in order to calculate quality surrogate variables (qSVs)

Based on [Jaffe et al 2017, PNAS](https://www.pnas.org/content/pnas/114/27/7130.full.pdf)  and https://github.com/nellore/region_matrix

## Installation and conda setup

1. clone the repo

2. create new conda environment with dependencies

`conda create -n qSVA-pipeline bioconda::wiggletools bioconda::snakemake python=3.6`

## Usage

### Input files

samples.tsv - a list of bam names (without the .bam suffix)

config.yaml - configuration file, sets output folder etc

lib_sizes.tsv - a column file with first column "samples" and second column "Aligned" the number of aligned reads per sample

metadata - path to samples.tsv
bamSuffix - default is ".bam", may want more complicated suffix like ".Aligned.Sorted.Clean.bam" etc


### Usage

* serial execution

`snakemake -s Snakefile --configfile test/config.yaml`

* parallel execution

 `./snakejob test/config.yaml`


