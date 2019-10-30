options(echo=TRUE)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sleuth))
suppressPackageStartupMessages(library(optparse))

option_list <- list(
   	make_option(c('--t2g'), help ='', default = "a file"),
	make_option(c('--rapidVersion'), help ='', default = "RAPiD"),
	make_option(c('--dataCode'), help='', default = "example"),
	make_option(c('--outFolder'), help='', default = "results/example/"),
	make_option(c('--metadata'), help='', default = "samples.tsv"),
	make_option(c('--refCondition'), help = '', default = "control"),
	make_option(c('--altCondition'), help = '', default = "case" )
)

parser <- OptionParser(usage = "%prog [options] genes.results1 genes.results2 .. genes.resultsN", option_list=option_list)


arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

t2g <- opt$t2g
rapidVersion <- opt$rapidVersion
outFolder <- opt$outFolder
dataCode <- opt$dataCode
metadata <- opt$metadata
refCondition <- opt$refCondition
altCondition <- opt$altCondition

#files <- arguments$args


debug_file <- paste0(outFolder, dataCode, "_debug.RData")

#print(debug_file)
save.image(file = debug_file)

#quit()


# figure out model design from support
support <- readr::read_tsv(metadata)

#head(support)
#print(refCondition)
#print(altCondition)

# remove any samples which are neither ref or alt condition
support <- dplyr::filter(support, condition %in% c(refCondition, altCondition) )

support$condition <- factor(support$condition, levels = c(refCondition, altCondition) )

# set path to kallisto abundance files
support$path <- paste0(support$rapid_path, "/Processed/", rapidVersion, "/kallisto/abundance.h5")
support$rapid_path <- NULL

print(support$path)


if( !all(file.exists(support$path) ) ){
message("some files are missing")
quit()
}

# map transcripts to genes

# take this directly from one of the abundance files
example_kallisto <- read_kallisto_h5(support$path[1], read_bootstrap = FALSE, max_bootstrap = NULL )
abundances <- example_kallisto$abundance
ttg <- as.data.frame(stringr::str_split_fixed(abundances$target_id, pattern = "\\|", n = 3)[,-3])
names(ttg) <- c("ext_tx", "ext_gene")
ttg$target_id <- abundances$target_id


# build sleuth object
so <- sleuth_prep(support, target_mapping = ttg, aggregation_column = 'ext_gene', extra_bootstrap_summary = TRUE)


# create full and reduced models
reduced_mod <- formula("~ 1")
full_mod <- formula("~ condition")

# if covariates are present

covariates <- names(support)[ !names(support) %in% c("sample", "condition", "path")]

if( length(covariates) > 0 ){
        reduced_mod <- model.matrix( formula(paste("~", c(covariates), collapse = " + ") ), data = support )
  full_mod_string <- paste("~", paste(c("condition", covariates), collapse = " + "))
	full_mod <- model.matrix( formula(full_mod_string ), data = support )
  message("covariates present in metadata")
  message(paste0("fitting model ", full_mod_string))
}


# fit models
so <- sleuth_fit(so, reduced_mod, 'reduced')
so <- sleuth_fit(so, full_mod, 'full')

# likelihood ratio test
so <- sleuth_lrt(so, 'reduced', 'full')

# get per-gene results
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
#sleuth_table_gene <- dplyr::filter(sleuth_table_gene, qval <= 0.05)

# get per-transcript results
sleuth_table_tx <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE, pval_aggregate = FALSE)
#sleuth_table_tx <- dplyr::filter(sleuth_table_tx, qval <= 0.05)

# write results and save RData
res_file_gene <- paste0(outFolder, "/", dataCode, "_results_gene.txt")
res_file_tx <- paste0(outFolder, "/", dataCode, "_results_tx.txt")
res_rdata <-  paste0(outFolder, "/", dataCode, "_results.RData")

readr::write_tsv(sleuth_table_gene, path = res_file_gene)
readr::write_tsv(sleuth_table_tx, path = res_file_tx)

save.image(file = res_rdata)
