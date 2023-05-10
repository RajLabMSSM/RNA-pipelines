
---
# Salmon 
---
This is a pipeline to run Salmon to quantify transcript abundances from RNA-seq data [Salmon git](https://combine-lab.github.io/salmon/about/)

## Input files

---

You will need to provide the following parameters in your config.yaml file:
1. inFolder: _Path to fastq data_ 
2. outFolder: _Folder for where to output the abundance estimates_
3. IndexOutFolder: _Folder for where to put the Salmon index_
4. metadata: _A dataframe containing 3 columns: list of sample IDs to use, in format of columns sample, R1, R2_
	* Sample_1, sample1_R1.fastq, sample1_R2.fastq
5. Annotation file: _GTF file_ 
6. stranded: _The way the library was prepared, use "A" for salmon to infer_
	* Unstranded: -l IU, Reverse stranded, -l ISR, Forward stranded, -l ISF , Unknown, salmon will infer: -l A 
Threads: _Specify threads to run, recommended is 10_
6. refType: "gencode"
	* reference type - only specified if both genome and fasta are GENCODE 
	* if other than GENCODE, provide ""
7. fastaType: _Name for fasta type_
	* ex: "isoseq"
	* **Specify which FASTA files to use as indexes**
8. genome_fastaFiles: _Path to genome file_
	* **Must be zipped!**
9. transcriptome_fastaFiles: _Path to transcriptome file_
	* fasta.gz or fa.gz format

## Output files

---

Salmon will output a quant.sf file containing the quantification of transcript abundances for each sample provided. The pipeline will produce a Rdata file with transcript counts and transcript tpm dataframes. Additionally a tpms folder will be created containing tsv files at the transcript and gene levels. 
