library(rjson)
library(argparse)
library(purrr)
library(stringr)
library(tibble)
library(dplyr)
library(readr)

args <- commandArgs(trailingOnly=TRUE)
outFolder <- (args[1])
dataCode <- args[2]

message(" * collating all run_info.json files")

# find all run_info.json files
all_json <- list.files(pattern = "run_info.json", recursive = TRUE, full.names = TRUE)

# read in all run_info files into a dataframe with purrr
json_df <- purrr::map_df( all_json,~{fromJSON(file = .x) } )

split_json <- str_split_fixed(all_json, pattern = "/", n =5 )

# make sure splitting is done correctly
stopifnot( all(split_json[,5] ==  "run_info.json") )

samples <- split_json[,3]

references <- split_json[,4]

# add sample and reference to table - use for downstream analysis
df  <- dplyr::bind_cols( tibble(sample = samples, ref = references), json_df)

# write out
outFile <- paste0(outFolder, "/", dataCode, "_run_info_collated.tsv")
write_tsv(df, path = outFile )
