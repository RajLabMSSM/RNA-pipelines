# find all salmon abundance estimates for a given reference
# assume in subdirectories of <outFolder>/<sample>/<index>/quant.sf
# get sample ID from the structure above
# collate them into a table with tximport
# save as RData
# execution

args <- commandArgs(trailingOnly=TRUE)
outFolder <- args[1]
index <- args[2]
cohort <- args[3]

#cohort <- stringr::str_split_fixed(outFolder, "/", 9)[,9]
print(cohort)

library(tximport)

message(" * Collating Salmon outputs" )
message( paste0("* Using files quantified with index: ", index) )

# first find all quant.sf files
tx_files <- list.files(path = outFolder, pattern = "quant.sf", full.names = TRUE, recursive = TRUE)

# take just those within the provided index folder
tx_files <- tx_files[ grepl(index, tx_files) ]
print(tx_files[1:20])
sample_ids <- stringr::str_split_fixed(tx_files, "/", 12)[,10]
print(sample_ids[1:20])


stopifnot(length(sample_ids) == length(tx_files) )
names(tx_files) <- sample_ids

print(sample_ids[1:20])

if(length(tx_files) > 0){
  message(paste0(" * collating transcripts from ", length(tx_files), " samples" ))
  tx_matrix <- tximport(tx_files, type = "salmon", txOut = TRUE) 
  
  tx_counts <- data.frame(tx_matrix$counts, check.names = F)
  tx_tpm <- data.frame(tx_matrix$abundance, check.names = F)
  
  names(tx_counts) <- sample_ids
  names(tx_tpm) <- sample_ids
 
 # names(tx_matrix$infReps) <- sample_ids
 
  save(tx_counts, tx_tpm, sample_ids, file = paste0(outFolder,"/", cohort , "_", index, "_salmon_counts.RData"))

  
  #save(tx_matrix, sample_ids, file = paste0(outFolder, index, "_kallisto_matrix.RData"))  

}else{
    stop( paste0("No transcript files found in ", outFolder) )
}


