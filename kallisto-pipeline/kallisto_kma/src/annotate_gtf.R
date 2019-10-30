library('data.table')
library('dplyr')
library('tidyr')

# pre-condition: generate_introns.py has been run and 'intron_to_transcripts.txt' has
# been generated
INTRON_TO_TRANSCRIPTS <- '../kma_out/intron_to_transcripts.txt'
#
# post-condition: generates gtf file with introns
OUT_FILE <- '../index/intron.gtf'

itt <- fread(INTRON_TO_TRANSCRIPTS, header = TRUE, data.table = FALSE)
itt <- separate(itt, intron_extension, c('chr', 'coordinates'), sep = ':', remove = FALSE)
itt <- separate(itt, coordinates, c('start', 'stop'), sep = '-')
itt <- select(itt, -target_id)
itt <- unique(itt)


# example line:
# 1       protein_coding  exon    881782  881925  .       -       .       gene_id "ENSG00000188976"; transcript_id "ENST00000327044"; exon_number "15"; gene_name "NOC2L"; gene_source "ensembl_havana";
format_gtf <- function(chr, start, stop, strand, gene, target_id, intron_name) {
  start <- as.integer(start) + 1
  stop <- as.integer(stop)
  # the '-' strand reverses everything so, since we are defining everything relative to
  # the start of the reference, just call everything positive strand
  strand <- '+'
  a <- sprintf(paste0('%s\tprotein_coding\ttranscript\t%d\t%d\t.\t%s\t.\tgene_id "%s";',
      ' transcript_id "%s"; exon_number "1"; gene_name "%s"; gene_source "ensembl_havana";'),
    chr, start, stop, strand, gene, target_id, intron_name
    )
  b <- sprintf(paste0('%s\tprotein_coding\texon\t%d\t%d\t.\t%s\t.\tgene_id "%s";',
      ' transcript_id "%s"; exon_number "1"; gene_name "%s"; gene_source "ensembl_havana";'),
    chr, start, stop, strand, gene, target_id, intron_name
    )
  paste0(a, '\n', b)
}

# TODO: make sure intron coordinates are correct with index file
# XXX: target name DOES NOT match FASTA
intron_gtf <- with(itt, Map(format_gtf, chr, start, stop, strand, gene, intron_extension,
    intron))
intron_gtf <- as.character(intron_gtf)
write(as.character(intron_gtf), file = OUT_FILE, sep = '\n')
