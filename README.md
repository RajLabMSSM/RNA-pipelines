# RNA-pipelines

A collection of snakemake pipelines for performing common analyses with RNA-seq data.

Each pipeline has:
	* a README
	* a conda recipe for dependencies
	* a test dataset
	* a snakemake script
	* parallel execution via a snakejob script

## Contents

*bazam-pipeline* - convert BAMs to FASTQs as quickly as possible for a set of BAMs

*RAPiD-tips* - documentation on how to run the [RAPiD pipeline](https://github.com/genomely/RAPiD-nf)

*collate_samples* - group together RSEM outputs by gene and transcript for a set of files processed by RAPiD, group together QC outputs into a single table

*qSVA-pipeline* - run qSVA to find RNA quality-associated surrogate variables (qSVs)

*kallisto-pipeline* - experimental workflows for generating custom kallisto indexes and running kallisto on lots of samples
