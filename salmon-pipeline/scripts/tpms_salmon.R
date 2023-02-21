
#Erica Brophy
#02-21-23
#make tpm files from Rdata

library(dplyr)
library(tibble)
library(stringr)
library(tidyverse)
library(readr)
library(readxl)


args <- commandArgs(trailingOnly = TRUE)
#print(args)
outFolder <- args[1]
print(outFolder)
annotations <- args[2]
print(annotations)
matrixData <- args[3]
print(matrixData)

cohort <- stringr::str_split_fixed(outFolder, "/", 9)[,9]
print(cohort)


#find Rdata file 
#matrixData <- list.files(path = outFolder, pattern = ".Rdata",full.names = TRUE, recursive = TRUE)
load(matrixData)

#annotaions file
annotations <- read_table2(annotations)
annotations <- as.data.frame(annotations)
names(annotations)[4] <- "transcript_id"
names(annotations)[5] <- "gene_id"

gene_id <- annotations$gene_id 
transcript_id <- annotations$transcript_id

df <- as.data.frame(gene_id) 
df$transcript_id <- transcript_id


tx_tpm <- tibble::rownames_to_column(tx_tpm, "transcript_id")
tx_counts <- tibble::rownames_to_column(tx_counts, "transcript_id")

tx_tpm_outFile <- paste0(outFolder, "/tpms/", cohort, "/", cohort, "_transcript_tpm.tsv")
write_tsv(tx_tpm, path = tx_tpm_outFile )
#readr::write_tsv(tx_tpm, "outFolder/tpms/{cohort}/{cohort}_transcript_tpm.tsv")

tx_counts_outFile <- paste0(outFolder, "/tpms/", cohort, "/", cohort, "_transcript_counts.tsv")
write_tsv(tx_counts, path = tx_counts_outFile )
#readr::write_tsv(tx_counts, "outFolder/tpms/{cohort}/{cohort}_transcript_counts.tsv")


df_tpm <- left_join(tx_tpm, df, by = "transcript_id") %>% dplyr::select(gene_id, everything()) %>% dplyr::select(-transcript_id)
gene_tpm <- df_tpm %>% group_by(gene_id) %>% summarise_all(sum)

df_count <- left_join(tx_counts, df, by = "transcript_id") %>% dplyr::select(gene_id, everything()) %>% dplyr::select(-transcript_id)
gene_count <- df_count %>% group_by(gene_id) %>% summarise_all(sum)


gene_tpm_outFile <- paste0(outFolder, "/tpms/", cohort, "/", cohort, "_gene_tpm.tsv")
write_tsv(gene_tpm, path = gene_tpm_outFile )
#readr::write_tsv(gene_tpm, "outFolder/tpms/{cohort}/{cohort}_gene_tpm.tsv")


gene_counts_outFile <- paste0(outFolder, "/tpms/", cohort, "/", cohort, "_gene_counts.tsv")
write_tsv(gene_count, path = gene_counts_outFile )
#readr::write_tsv(gene_count, "outFolder/tpms/{cohort}/{cohort}_gene_counts.tsv")

