
# Collate all STAR-Fusion output files together

library(tidyverse)
library(optparse)


option_list <- list(
        make_option(c('--outFolder'), help='', default = "example")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


output <- opt$output



message(" * Collating STAR-FUSION outputs" )

# first find all abundance.h5 files
files <- list.files(path = outFolder, pattern = "star-fusion.fusion_predictions.tsv", full.names = TRUE, recursive = TRUE)

# take just those within the provided index folder
print(files[1:20])

sample_ids <- basename(dirname(files) )

stopifnot(length(sample_ids) == length(files) )

names(files) <- sample_ids

print(sample_ids[1:20])

# read in all files
# make unique ID for fusions
data <- map( files, ~{
    d <- read_tsv(.x)
    d$ID <- paste0( d$`#FusionName`, d$LeftBreakpoint, d$RightBreakpoint, sep = ",")
    return(d)
})

# get out unique set of fusions

ffpm <- map2(data, names(data), ~{
    d <- select(.x, ID, FFPM) 
    names(d)[2] <- .y
}) %>% reduce( full_join, by = "ID")

meta <- map_df(data, ~{
    select(.x, ID, annots)
}) %>% distinct()



write_tsv(ffpm, path = count_file)

write_tsv(meta, path




