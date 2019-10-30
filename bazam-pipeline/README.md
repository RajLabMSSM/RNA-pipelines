# bazam-pipeline

Convert BAMs to FASTQs for large number of samples, maintaining file IDs

## conda recipe

conda create -n bazam-pipeline -c bioconda python=3.6 bazam=1.0.1-0 snakemake=5.7.4


## Testing

conda activate bazam-pipeline

snakemake -s Snakefile --configfile test/test_config.yaml

## Config options

inFolder: where the BAMs are
outFolder: where to output FASTQs

bamSuffix: string to strip off the end of sampleIDs when naming FASTQs
	eg sample1.coord.sorted.aligned.unique.bam, the bamSuffix is ".coord.sorted.aligned.unique.bam". The resulting FASTQ will then be sample1.fastq.gz


