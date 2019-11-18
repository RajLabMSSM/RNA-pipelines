# kallisto-pipeline

Jack Humphrey 2019

A snakemake workflow to:

* generate multiple Kallisto indexes from different transcript lists
*  run Kallisto to estimate transcript abundance with bootstraps
* run Sleuth to call differential isoform usage between groups

In addition this will contain projects using the results of this pipeline:

1. Estimate novel isoforms from AD brain isoseq reference

2. Estimate intron retention in human brains using KeepMeAround (Pimentel, 2019, unpublished)


## dependencies

snakemake
kallisto v0.45.0 - already installed on chimera

## testing

1. clone repo
2. load dependencies
3. run `snakemake -s Snakefile --configfile example/config.yaml -pr`

## usage

samples.tsv - list of sample IDs

config.yaml - see `example/config.yaml`

### assumed use case
you have a set of FASTQ files with the forms <sample>.R1.fastq.gz and <sample>.R2.fastq.gz for each <sample>. `samples.tsv` is a single column text file with a column named `sample` with a set of sample IDs.

you also have one or more transcript sets (in FASTA format) that you want to quantify for your samples. Specify these in the `config.yaml`.




