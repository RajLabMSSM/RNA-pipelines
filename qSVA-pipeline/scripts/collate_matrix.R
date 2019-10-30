library(optparse)
options(echo = TRUE)
# parse args
option_parser=OptionParser(
  usage="%prog [options] ",
  option_list=list(
    make_option( "--code", help = "the analysis code" ),
    make_option( c("-o","--outFolder"), 
                 default="output/", 
                 help="The name and path of where to output the degradation matrix [%default]"
	),
   make_option( c("-l", "--lib_size_file"),
		help = "a tab-separated column file with sample and aligned, correponding to the same sample names as in samples.tsv"
	),
  make_option( c("-s", "--sample_file"),
	help = "a list of sample names"
	),
  make_option( c("-r", "--read_length"),
	default = 100,
	help = "read length for your BAM files"
	)
))

opt <- parse_args(option_parser)


code <- opt$code
outFolder <- opt$outFolder
lib_sizes_file <- opt$lib_size_file
samples_file <- opt$sample_file
read_length <- opt$read_length
library(sva)
library(readr)
# samples.tsv

# libSizes

# outFolder

# read in samples
samples <- read_tsv(samples_file)

names(samples)[1] <- "sample"

# read in lib sizes
lib_sizes <- read_tsv(lib_sizes_file)

names(lib_sizes)[1] <- "sample"
names(lib_sizes)[2] <- "Aligned"


# match lib size
samples$lib_size <- lib_sizes$Aligned[ match(samples$sample, lib_sizes$sample) ]

#print(head(samples))

samples$file <- paste0(outFolder, samples$sample, ".degradation_matrix.tsv")

#print(samples)

deg_matrix <- 
	read.degradation.matrix(
		covFiles = samples$file, 
		sampleNames = samples$sample, 
		totalMapped = samples$lib_size, 
		readLength = read_length,
		normFactor = 8e+07, 
		type = "region_matrix_single", BPPARAM = bpparam())


outFile <- paste0(outFolder, code ,".degradation_matrix.RData")

save(deg_matrix, file = outFile)

# calculate qSVs?
