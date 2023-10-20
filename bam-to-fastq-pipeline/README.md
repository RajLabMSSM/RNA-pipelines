# bam-to-fastq-pipeline

Jack Humphrey, Ben Muller, Kailash BP

Towfique Raj Lab

2020-2023

## Rationale

Convert BAMs back to FASTQs for large number of samples, maintaining file IDs, using Picard.

This pipeline formerly used Bazam, but has been replaced with Picard as Bazam does not respect stranded RNA-seq libraries.

## Testing

conda activate snakemake

snakemake -s Snakefile --configfile test/test_config.yaml

## Config options

inFolder: where the BAMs are, path must end in "/" - these can be symlinks

outFolder: where to output FASTQs, path must end in "/"

bamSuffix: string to strip off the end of sampleIDs when naming FASTQs
	eg sample1.coord.sorted.aligned.unique.bam, the bamSuffix is ".coord.sorted.aligned.unique.bam". The resulting FASTQ will then be sample1.fastq.gz

paired: `True` if paired-end sequencing, `False` if single-end. Default is `True`.

