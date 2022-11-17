### Nov.20, 2020
# plot user defined genes by loding RData

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(BiocParallel))

visualize_dexseq_results <- function(dxr, contrast_name, geneName, anno){
  # get geneID from geneName with exon.df
  geneID <- anno$Gene[which(anno$Name == geneName)] %>% unique()
  
  # normalized_counts.pdf
  name_tmp <-
    paste(contrast_name, geneName, "normalized_counts.pdf", sep = ".")
  name_tmp <- gsub("\\+", "_", name_tmp)
  name_tmp <- gsub("\\.", "_", name_tmp)
  name_tmp <- gsub("_pdf$", ".pdf", name_tmp)
  print(paste("saving: ", outDir, name_tmp, sep = ""))
  try(plotDEXSeq_norCounts(dxr, geneID, outDir, name_tmp))
  
  ## relative_exon_usage.pdf
  name_tmp <-
    paste(contrast_name, geneName, "relative_exon_usage.pdf", sep = ".")
  name_tmp <- gsub("\\+", "_", name_tmp)
  name_tmp <- gsub("\\.", "_", name_tmp)
  name_tmp <- gsub("_pdf$", ".pdf", name_tmp)
  print(paste("saving: ", outDir, name_tmp, sep = ""))
  try(plotDEXSeq_relative_exon_usage(dxr, geneID, outDir, name_tmp))
}

geneName <- 'Runx1'
RData.path <- "../DEXSeq/contrast1/contrast1.RData"

load(RData.path)
outDir <- paste(dirname(RData.path), "/", sep = "")
visualize_dexseq_results(dxr, contrast_name = name, geneName = geneName, anno=anno)

