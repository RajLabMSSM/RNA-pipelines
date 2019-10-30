# import gene matrix
# remove low count genes
# apply voom normalisation
# save normalised counts
# apply DESEq normalisation
# save normalised counts
library(edgeR)
library(limma)
library(DESeq2)

args <- commandArgs(trailingOnly=TRUE)

inFile <- args[1]

outFolder <- args[2]

message("loading counts")

# gene_matrix.Rdata
load(inFile)

# remove low count genes

cpm <- cpm(genes_counts)
# CPM >= 1 in at least 10% of samples
keep.exp <- rowSums(cpm > 1) >= (0.1 * ncol(genes_counts) )

genes_counts <- genes_counts[ keep.exp, ]

message( paste( nrow(genes_counts), "pass CPM filter" ) )

message("performing voom normalisation")

counts_voom <- limma::voom(genes_counts)
genes_counts_voom <- counts_voom$E

save(genes_counts, genes_counts_voom, file = paste0(outFolder, "genes_counts_voom.RData") )

genes_counts_rounded <- round(genes_counts)

save.image("~/test.RData")

message("DESeq normalisation")

support <- data.frame(sample = names(genes_counts_rounded), condition = 1, stringsAsFactors = FALSE)

dds <- DESeqDataSetFromMatrix(countData = genes_counts_rounded,
                               colData = support,
                               design = ~ 1)

dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)

save(dds, vsd, file = paste0(outFolder, "genes_counts_deseq.RData") )

