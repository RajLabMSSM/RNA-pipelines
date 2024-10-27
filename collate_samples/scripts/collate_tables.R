# find all gene and transcript files in a RAPiD directory
# collate them into a table with tximport
# save as RData
# execution
# Rscript collate_tables.R <inFolder> <outFolder>

args <- commandArgs(trailingOnly=TRUE)
inFolder <- args[1]
outFolder <- args[2]

library(tximport)

gene_files <- list.files(path = inFolder, pattern = ".genes.results", full.names = TRUE, recursive = TRUE)

tx_files <- list.files(path = inFolder, pattern = ".isoforms.results", full.names = TRUE, recursive = TRUE)

# remove any files that are from the "work/" directory - all should have the string "rsem" in the folder name
gene_files <- gene_files[ grepl("rsem", gene_files) ]
tx_files <- tx_files[ grepl("rsem", tx_files) ]

if(length(gene_files) > 0){
  message(paste0(" * collating genes from ", length(gene_files), " samples" ))
  gene_matrix <- tximport(files = gene_files, type = "rsem", countsFromAbundance="scaledTPM"  )
  
  genes_counts <- data.frame(gene_matrix$counts, check.names = F)
  genes_tpm <- data.frame(gene_matrix$abundance, check.names = F)
  
  names(genes_counts) <- gsub(".genes.results", "", basename(gene_files))
  names(genes_tpm) <- gsub(".genes.results", "", basename(gene_files))
  
  
  save(genes_counts, genes_tpm, file = paste0(outFolder, "gene_matrix.RData"))
  
}else{
    stop( paste0("No gene files found in ", inFolder) )
}

if(length(tx_files) > 0){
  message(paste0(" * collating transcripts from ", length(tx_files), " samples" ))
  tx_matrix <- tximport(tx_files, type = "rsem", txIn = TRUE, txOut = TRUE) 
  
  tx_counts <- data.frame(tx_matrix$counts, check.names = F)
  tx_tpm <- data.frame(tx_matrix$abundance, check.names = F)
  
  names(tx_counts) <- gsub(".isoforms.results", "", basename(tx_files))
  names(tx_tpm) <- gsub(".isoforms.results", "", basename(tx_files))
  
  
  save(tx_counts, tx_tpm, file = paste0(outFolder, "tx_matrix.RData"))
  
}else{
    stop( paste0("No transcript files found in ", inFolder) )
}


