# Downloading data from SRA

Jack Humphrey 2023

Fastest way to download RNA-seq data in FASTQ format from the Sequence Read Archive (SRA)

## Recipe

1. Fetches the SRA format files from the web

2. Uses `fasterq-dump` from SRAtoolkit to rapidly convert them to FASTQ

3. Then runs parallel gzip `pigz` to compress to FASTQ.GZ

4. Removes all temporary files along the way to save space - this is a key bottleneck when downloading files


Default behaviour is to download the technical reads as well - for single-cell RNA-seq these are the I1 sample index reads made by Cellranger mkfastq.

In some cases the I1 reads are uploaded to SRA as the R1s, and the barcodes or cDNA reads are labelled as R3 (technical reads).

Example:

`sh download_files.sh example_accessions.txt ` 

The first positional argument is a text file containing a set of SRA accessions (IDs starting with SRR).


