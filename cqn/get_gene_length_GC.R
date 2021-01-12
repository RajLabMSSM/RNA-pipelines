#!/usr/bin/env Rscript
library("biomaRt")

ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)

df <- getBM(attributes = c("ensembl_gene_id", "start_position", "end_position", "percentage_gene_gc_content"), mart = ensembl)

df$length <- as.numeric(df$end_position) - as.numeric(df$start_position)

library("readr")

write_tsv(df, path = "data/homo_sapiens_ensembl_gc_length.tsv.gz")

