
# Collate all STAR-Fusion output files together

library(tidyverse)
library(optparse)


option_list <- list(
        make_option(c('--outFolder'), help='', default = "example"),
        make_option(c('--dataCode'), help='', default = "test")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


outFolder <- opt$outFolder
dataCode <- opt$dataCode

count_file <- paste0(outFolder, "/", dataCode, "_fusion_counts.tsv.gz")
meta_file <- paste0(outFolder, "/", dataCode,  "_fusion_meta.tsv.gz")

message(" * Collating STAR-FUSION outputs" )

# first find all abundance.h5 files
files <- list.files(path = outFolder, pattern = "star-fusion.fusion_predictions.tsv", full.names = TRUE, recursive = TRUE)

# take just those within the provided index folder
print(head(files) )
#print(files[1:20])

sample_ids <- basename(dirname(files) )

stopifnot(length(sample_ids) == length(files) )

names(files) <- sample_ids

#print(sample_ids[1:20])

# read in all files
# make unique ID for fusions
message(" * reading in files")
data <- map( files, ~{
    d <- read_tsv(.x)
    d$ID <- paste( d$`#FusionName`, d$LeftBreakpoint, d$RightBreakpoint, sep = "|")
    return(d)
})

#########save.image('test.RData')
# get out unique set of fusions

ffpm <- map2(data, names(data), ~{
    d <- select(.x, ID, FFPM) 
    names(d)[2] <- .y
    return(d)
}) %>% reduce( full_join, by = "ID")

meta <- map_df(data, ~{
    select(.x, ID, annots)
}) %>% distinct()


message("* combining counts!")
head(ffpm)

message("* combining metadata!")
head(meta)

##write out 
message(" * writing ", count_file)
write_tsv(ffpm, file = count_file)
message(" * writing ", meta_file)
write_tsv(meta, file = meta_file)




