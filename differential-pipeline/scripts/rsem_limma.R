suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(optparse))

library(readr)

t2g_file <- "~/GENCODE/gencode.v30.tx2gene.tsv"

t2g <- read_tsv(t2g_file)

option_list <- list(
    make_option(c('--dataCode'), help='', default = "example"),
	make_option(c('--outFolder'), help='', default = "results/example/"),
	make_option(c('--metadata'), help='', default = "samples.tsv"),
	make_option(c('--refCondition'), help = '', default = "control"),
	make_option(c('--altCondition'), help = '', default = "case" )
)

parser <- OptionParser(usage = "%prog [options] genes.results1 genes.results2 .. genes.resultsN", option_list=option_list)


arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

outFolder <- opt$outFolder
dataCode <- opt$dataCode
metadata <- opt$metadata
refCondition <- opt$refCondition
altCondition <- opt$altCondition

files <- arguments$args

debug_file <- paste0(outFolder, dataCode, "_debug.RData")

print(debug_file)
save.image(file = debug_file)


# parse arguments

# figure out model design from support
support <- readr::read_tsv(metadata)


# remove any samples which are neither ref or alt condition
support <- dplyr::filter(support, condition %in% c(refCondition, altCondition) )

support$condition <- factor(support$condition, levels = c(refCondition, altCondition) )


design <- model.matrix( ~ condition , data = support )

# add covariates if present
# if any columns in metadata are not sample, rapid_path or condition, assume they are covariates
covariates <- names(support)[ !names(support) %in% c("sample", "condition", "rapid_path")]

if( length(covariates) > 0 ){

  full_mod_string <- paste("~", paste(c("condition", covariates), collapse = " + "))

  design <- model.matrix( formula(full_mod_string ), data = support )

  message("covariates present in metadata")
  message(paste0("fitting model ", full_mod_string))

}


## DIFFERENTIAL GENE EXPRESSION
message("1. differential gene expression")
genes.files <- paste0(support$rapid_path, "/Processed/RAPiD/rsem/", support$sample, ".genes.results")

genes.rsem <- tximport(genes.files, type = "rsem", countsFromAbundance="scaledTPM")

# get count matrix

x <- data.frame(genes.rsem$counts, check.names = F)

names(x) <- gsub(".genes.results", "", basename(genes.files))

cpm = cpm(x)

head(x)

# harmonise counts with samples
if( all(support$sample %in% names(x) ) ){
	x <- x[,support$sample ]
	cpm <- cpm[,support$sample ]
}else{
	message("samples in metadata not found in count matrix!")
	quit()
}
# filter out lowly expressed genes
keep.exp <- rowSums(cpm > 1) >= ceiling(0.3*ncol(x))
x <- x[keep.exp,]

# save count matrix
gene_count_matrix_file <- paste0(outFolder, "/", dataCode, "_gene_matrix.tsv")
gene_count_matrix <- tibble::rownames_to_column(x, var = "EnsemblID")

write_tsv(gene_count_matrix, path = gene_count_matrix_file)

# calculate normalisation factors
norm <- calcNormFactors(x, method = "TMM") #normalized counts with TMM

# run differential gene expression
dge <- DGEList(counts=x, samples=support, norm.factors = norm)
v <- voom(dge, design)
vfit <- lmFit(v, design)
efit <- eBayes(vfit)
topTable(efit, coef=ncol(design))
summary(efit, coef=ncol(design))
summary(decideTests(efit, p.value  = "0.05"))

coefficient <- paste0("condition", altCondition)
res_genes <- data.frame(topTable(efit, coef = coefficient, number = nrow(x), sort.by="p"), check.names = FALSE)
res_genes <- tibble::rownames_to_column(res_genes, var = "EnsemblID")

res_genes$gene <- t2g$GENENAME[match(res_genes$EnsemblID, t2g$GENEID) ]

res_genes <- dplyr::select(res_genes, EnsemblID, gene, everything() )

# save results and RData
res_genes_file <- paste0(outFolder, "/", dataCode, "_genes_results.txt")
readr::write_tsv(res_genes, path = res_genes_file)

## DIFFERENTIAL TRANSCRIPT EXPRESSION
message("2. differential transcript expression")
tx.files <- paste0(support$rapid_path, "/Processed/RAPiD/rsem/", support$sample, ".isoforms.results")

tx.rsem <- tximport(tx.files, type = "rsem", txIn = TRUE, txOut = TRUE)

# get count matrix

x <- data.frame(tx.rsem$counts, check.names = F)

names(x) <- gsub(".isoforms.results", "", basename(tx.files))

cpm = cpm(x)

head(x)

# harmonise counts with samples
if( all(support$sample %in% names(x) ) ){
	x <- x[,support$sample ]
	cpm <- cpm[,support$sample ]
}else{
	message("samples in metadata not found in count matrix!")
	quit()
}

# save count matrix
tx_count_matrix <- tibble::rownames_to_column(x, var = "transcriptID")
tx_count_matrix_file <- paste0(outFolder, "/", dataCode, "_tx_matrix.tsv")

write_tsv(tx_count_matrix, path = tx_count_matrix_file)

# filter out lowly expressed transcripts
keep.exp <- rowSums(cpm > 1) >= ceiling(0.3*ncol(x))
x <- x[keep.exp,]

# calculate normalisation factors
norm <- calcNormFactors(x, method = "TMM") #normalized counts with TMM

# run differential transcript expression
dge <- DGEList(counts=x, samples=support, norm.factors = norm)
v <- voom(dge, design)
vfit <- lmFit(v, design)
efit <- eBayes(vfit)
topTable(efit, coef=ncol(design))
summary(efit, coef=ncol(design))
summary(decideTests(efit, p.value  = "0.05"))

coefficient <- paste0("condition", altCondition)
res_tx <- data.frame(topTable(efit, coef = coefficient, number = nrow(x), sort.by="p"), check.names = FALSE)
res_tx <- tibble::rownames_to_column(res_tx, var = "EnsemblID")

res_tx$gene <- t2g$GENENAME[match(res_tx$EnsemblID, t2g$TXNAME) ]

res_tx <- dplyr::select(res_tx, EnsemblID, gene, everything() )


# save results and RData
res_tx_file <- paste0(outFolder, "/", dataCode, "_transcripts_results.txt")

readr::write_tsv(res_tx, path = res_tx_file)

## DIFFERENTIAL TRANSCRIPT USAGE
message("3. differential transcript usage")

# use the limma diffSplice method
# must match genes to transcripts
tx_order <- rownames(dge$counts)
gene_order <- t2g$GENENAME[ match(tx_order, t2g$TXNAME)]

# set gene name in dge object
dge$genes <- data.frame(GeneID = gene_order, txid = tx_order)

#redo fit
v <- voom(dge, design, plot=TRUE)
fit <- lmFit(v, design)

# test for differential transcript usage
ex <- diffSplice(fit, geneid="GeneID")

# get out per-gene test values
dtu_per_gene <- topSplice(ex, coef=2, test="simes", number = Inf)

# get out per-transcript tests (with fold changes!)
dtu_per_transcript <- topSplice(ex, coef=2, test="t", number = Inf)

# save results
res_dtu_gene_file <- paste0(outFolder, "/", dataCode, "_dtu_gene_results.txt")
res_dtu_tx_file <- paste0(outFolder, "/", dataCode, "_dtu_transcript_results.txt")

write_tsv(dtu_per_gene, path = res_dtu_gene_file)
write_tsv(dtu_per_transcript, path = res_dtu_tx_file)

## END

res_rdata <-  paste0(outFolder, "/", dataCode, "_results.RData")
save.image(file = res_rdata)
