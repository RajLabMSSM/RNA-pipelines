# filter a FASTA file to only include particular transcripts

# currently used to only retain sequences of polyadenylated transcripts
inFile <- "/sc/orga/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v32.transcripts.fa.gz"

refFile <- "/sc/orga/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v32.metadata.PolyA_feature.gz"
# which column contains transcript IDs?
ref_column <- 1

outFile <- "gencode.v32.polyadenylated_transcripts.fa.gz"


library(Biostrings)
library(stringr)

f <- readDNAStringSet(inFile)


# keep only polyadenylated transcripts?
headers <- as.data.frame(str_split_fixed(names(f), pattern = "\\|", n = 8), stringsAsFactors = FALSE )

transcripts <- headers$V1
genes <- headers$V2

message(paste0("* read ", length(f), " transcript sequences\n" ) )
message(paste0("* ", length(unique(transcripts)), " transcripts from ", length(unique(genes)), " genes") )

# read in reference table
ref <- readr::read_tsv(refFile, col_names = FALSE)

ref_transcripts <- unique(ref[, ref_column, drop = TRUE])

filtered_fasta <- f[ transcripts %in% ref_transcripts,] 

message(paste0("* kept ", length(filtered_fasta), " transcript sequences\n" ) )

writeXStringSet(filtered_fasta, outFile, compress = TRUE)
