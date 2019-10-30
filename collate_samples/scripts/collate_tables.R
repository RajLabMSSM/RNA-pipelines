# find all gene and transcript files in a RAPiD directory
# collate them into a table with tximport
# save as RData

args <- commandArgs(trailingOnly=TRUE)

outFolder <- args[1]

library(tximport)

gene_files <- list.files(pattern = ".genes.results", full.names = TRUE, recursive = TRUE)

tx_files <- list.files(pattern = ".isoforms.results", full.names = TRUE, recursive = TRUE)


if(length(gene_files) > 0){
  gene_matrix <- tximport(files = gene_files, type = "rsem", countsFromAbundance="scaledTPM"  )
  
  genes_counts <- data.frame(gene_matrix$counts, check.names = F)
  genes_tpm <- data.frame(gene_matrix$abundance, check.names = F)
  
  names(genes_counts) <- gsub(".genes.results", "", basename(gene_files))
  names(genes_tpm) <- gsub(".genes.results", "", basename(gene_files))
  
  
  save(genes_counts, genes_tpm, file = paste0(outFolder, "gene_matrix.RData"))
  
}

if(length(tx_files) > 0){
  tx_matrix <- tximport(tx_files, type = "rsem", txIn = TRUE, txOut = TRUE) 
  
  tx_counts <- data.frame(tx_matrix$counts, check.names = F)
  tx_tpm <- data.frame(tx_matrix$abundance, check.names = F)
  
  names(tx_counts) <- gsub(".tx.results", "", basename(tx_files))
  names(tx_tpm) <- gsub(".tx.results", "", basename(tx_files))
  
  
  save(tx_counts, tx_tpm, file = paste0(outFolder, "tx_matrix.RData"))
  
}




