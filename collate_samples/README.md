# Collate samples

Create matrices of gene and isoform expression values following a RAPiD run.

Written in snakemake and R (3.6.0)

R dependencies:

* tximport
* limma
* edgeR
* DESeq2

## Usage

1. clone repo

2. edit the snakefile to specify `inFolder` - full path to top directory of RAPiD run, and `outFolder` - full path to where you want tables to be saved.

3. load snakemake (within a conda environment)

4a. local execution: snakemake -s Snakefile

4b. parallel execution: ./snakejob (experimental)
