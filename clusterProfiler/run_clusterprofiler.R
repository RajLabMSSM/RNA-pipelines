#BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "GOSemSim"))
#install.packages("tidyverse")

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOSemSim)
# create GO term semamtic similarity databases for each type
hsGO_BP <- godata('org.Hs.eg.db', ont="BP")
hsGO_MF <- godata('org.Hs.eg.db', ont="MF")
hsGO_CC <- godata('org.Hs.eg.db', ont="CC")



# df - a limma differential expression data.frame with log_fc, geneid (EnsemblID), and adj_p_val
# min.p - to subset by adjusted P value 
# n - to take the top N genes as the gene set
# direction - one of "up" or "down" to take upregulated or downregulated genes respectively
# simplify - whether to reduce the number of GO terms using semantic similarity
# similarity_cutoff - how to decide similar GO terms
runCPGO <- function(df, n = NULL, min.p = 0.05,direction = "up", simplify = TRUE, similarity_cutoff = 0.7){
  # universe is all gene IDs present in df
  uni <- df$geneid
  
  # pvalue dplyr::filtering
  df <- dplyr::filter(df, adj_p_val < min.p)
  
  # split by direction or not
  stopifnot( direction %in% c("up", "down"))
  if( direction == "up"){
    df <- dplyr::filter(df, log_fc > 0)
  }else{
    df <- dplyr::filter(df, log_fc < 0)
  }
  
  if( !is.null(n) ){
    stopifnot( is.numeric(n) & n > 0)
    df <- df[ order(df$p_value),]
    df <- head(df, n) 
  }
  
  if(nrow(df) == 0){ stop("0 genes remain for testing")}
  
  set <- df$geneid
  #return(list(set = set, uni = uni))
  print(paste0(" * ", length(set), " genes for testing; ", length(uni), " genes in background" ))
  
  # run GOenrichment separately for each type
  # for some reason simplify doesn't work with "ALL" types together
  print(" * clusterProfiler GO enrichment - biological process")
  bp <- enrichGO(gene          = set,
                universe      = uni,
                OrgDb         = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
  print(" * clusterProfiler GO enrichment - molecular function")
  mf <-  enrichGO(gene          = set,
                universe      = uni,
                OrgDb         = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont           = "MF",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
  print(" * clusterProfiler GO enrichment - cell component")
  cc <-  enrichGO(gene          = set,
                universe      = uni,
                OrgDb         = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
        readable      = TRUE)
  results <- list(BP = as.data.frame(bp), MF = as.data.frame(mf), CC = as.data.frame(cc) )
  
  # simplify
  if( simplify == TRUE){
    print( " * collapsing GO terms by semantic similarity")
    mf_simple <- simplify(mf,  cutoff = similarity_cutoff, measure = "Wang", semData = hsGO_MF )
    bp_simple <- simplify(bp,  cutoff = similarity_cutoff, measure = "Wang", semData = hsGO_BP )
    cc_simple <- simplify(cc,  cutoff = similarity_cutoff, measure = "Wang", semData = hsGO_CC )
  
    simplify_res <- 
      list(BP_simple = as.data.frame(bp_simple), CC_simple = as.data.frame(cc_simple), MF_simple = as.data.frame(mf_simple) ) 
    results <- c(results, simplify_res)
    }
  
    return(results)
}
