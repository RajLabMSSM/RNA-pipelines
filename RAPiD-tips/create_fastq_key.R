library(optparse)
option_list <- list(
    make_option(c('-i', '--inFolder' ), help='The full path to the folder that contains the fastQ files', default = "fastq"),
    make_option(c('-e', '--endName1' ), help="the identifier in the file name which refers to the name of the end 1 (R1, END1, etc.)", default = "R1"),
    make_option(c('-f', '--endName2' ), help="the identifier in the file name which refers to the name of the end 2 (R2, END2, etc.)", default = "R2")
)
option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)
inFolder <- opt$inFolder
endName1 <- opt$endName1
endName2 <- opt$endName2

# read in folder full of FASTQ files and produce a table 
# split out sample name and find R1 and R2
fq <- list.files(inFolder, full.names =  TRUE)
# ORIG CODE
#sample <- basename(stringr::str_split_fixed(fq, "_", 4)[,1])

# TEST1 remove "fastq.gz" and take basename (problem -- keeps R1 and R2 and doesn't split)
#del_suffix <- stringr::str_remove(fq, ".fastq.gz")
#sample <- basename(del_suffix)

# TEST2 str_split on the first "." and keep everything before it (problem -- not robust enough -- only works when "end name" is after the first ".")
#sample <- basename(stringr::str_split(fq, "\\.", simplify=TRUE)[,1])

# TEST3 str_split on the "end name" option (R1, END1, etc.) and keep the first thing -- test on multiple datasets before pushing to Master
temp <- basename(stringr::str_split(fq, endName1, simplify=TRUE)[,1])
sample <- basename(stringr::str_split(temp, endName2, simplify=TRUE)[,1])

pairs <- split(fq, sample)

# iterate through pairs and assign f1 and f2
assignPair <- function(pair){
    f1 = pair[ grepl("R1", pair) ]
    f2 = pair[ grepl("R2", pair) ]
    
    data.frame(f1 = paste(f1, collapse = ","), f2 = paste(f2, collapse = ","), stringsAsFactors = FALSE )
}

pair_df <- purrr::map_df(pairs, assignPair, .id = "sample")

readr::write_tsv(pair_df, "sample_fastq_table.tsv")
