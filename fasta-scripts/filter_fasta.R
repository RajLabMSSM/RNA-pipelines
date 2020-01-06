inFile <- "/sc/orga/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v32.pc_transcripts.fa.gz"

outFile <- "gencode.v32.pc_transcripts_min_100_aa_transcripts.fa.gz"

min_aa_length <- 100

library(Biostrings)
library(stringr)

f <- readDNAStringSet(inFile)


# keep only polyadenylated transcripts?
headers <- as.data.frame(str_split_fixed(names(f), pattern = "\\|", n = 8), stringsAsFactors = FALSE )

transcripts <- headers$V1
genes <- headers$V2

message(paste0("* read ", length(f), " transcript sequences\n" ) )
message(paste0("* ", length(unique(transcripts)), " transcripts from ", length(unique(genes)), " genes") )

# extract CDS coords

# header stored in names(f)
# parse headers - extract CDS coords
# each header has a UTR5, CDS and UTR3 section delimited by "|"
# each section has the numerical coordinates x-y
# not every transcript has a UTR5 so the CDS entry isn't always in the same place
# instead cut header string at "CDS" - the numbers x-y immediately following are the coordinates
headers <- as.data.frame(str_split_fixed(names(f), pattern = "CDS:", n = 2), stringsAsFactors = FALSE)[,2]
# split again on "|"
headers <- as.data.frame(str_split_fixed(headers, pattern = "\\|", n = 2), stringsAsFactors = FALSE )[,1]
# split x-y coordinates on "-"
coords <- as.data.frame(str_split_fixed(headers, pattern = "\\-", n = 2), stringsAsFactors = FALSE)
# any non numeric value will be set to NA
coords$cds_start <- as.numeric(coords$V1)
coords$cds_end <- as.numeric(coords$V2)

# any missing values set to 0
coords[ is.na(coords) ] <- 0 

coords$cds_length <- coords$cds_end - coords$cds_start
# counting starts at 0 so add 1
coords$aa_length <- (coords$cds_length + 1) / 3

# append amino acid length to fasta headers
names(f) <- paste0(names(f), "AA:", coords$aa_length, "|" )

# remove any transcripts with an amino acid length less than 100.

to_remove <- coords$aa_length < min_aa_length

filtered_fasta <- f[ !to_remove ]

message( paste0("* removing ", sum(to_remove), " sequences for having CDS less than ", min_aa_length) )

writeXStringSet(filtered_fasta, outFile, compress = TRUE)
