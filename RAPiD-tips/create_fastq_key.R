
# read in folder full of FASTQ files and produce a table 
# split out sample name and find R1 and R2
fq <- list.files("fastq/", full.names =  TRUE)
sample <- basename(stringr::str_split_fixed(fq, "_", 4)[,1])

pairs <- split(fq, sample)

# iterate through pairs and assign f1 and f2
assignPair <- function(pair){
    f1 = pair[ grepl("R1", pair) ]
    f2 = pair[ grepl("R2", pair) ]
    
    data.frame(f1 = paste(f1, collapse = ","), f2 = paste(f2, collapse = ","), stringsAsFactors = FALSE )
}

pair_df <- purrr::map_df(pairs, assignPair, .id = "sample")

readr::write_tsv(pair_df, "sample_fastq_table.tsv")
