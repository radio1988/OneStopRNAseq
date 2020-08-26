options(show.error.locations = TRUE)
options(error=traceback)

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(BiocParallel))
sessionInfo()

BPPARAM = MulticoreParam(workers=8) # test

readExcel <- function(fname){
  df <- readxl::read_xlsx( fname, na='NA', sheet=1)  # na term important
  df <- data.frame(df)  #important
  return (df)
}

plotDEXSeq_norCounts <- function(dxr, gene, outDir, name_tmp) {
	pdf(paste(outDir, name_tmp, sep=""))
	plotDEXSeq(dxr, gene, expression=FALSE, norCounts=TRUE, displayTranscripts=TRUE,
	  legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
	dev.off()
}

plotDEXSeq_relative_exon_usage <- function(dxr, gene, outDir, name_tmp) {
	pdf(paste(outDir, name_tmp, sep=""))
	plotDEXSeq(dxr, gene, expression=FALSE, splicing=TRUE, displayTranscripts=TRUE,
	  legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
	dev.off()
}

plotDispEstsWrapper <- function(dxd, outDir, name){
	pdf(paste(outDir, name, ".disp.pdf", sep=""))
	plotDispEsts( dxd )
	dev.off()
}

plotMAWrapper <- function(dxd, ourDir, name){
	pdf(paste(outDir, name, ".MA.pdf", sep=""))
	plotMA( dxr, cex=0.8 ,  alpha = maxFDR) # contain NA error message
	dev.off()
}

outDir <- "./DEXSeq/" # test
dir.create(outDir, showWarnings = FALSE)
args = commandArgs(trailingOnly=TRUE)

# Read ARGS
# todo: remove countFile from args
if (length(args) < 2){
    # development
    metaFile <- 'meta/meta.xlsx'  
    contrastFile <- './meta/contrast.de.xlsx'
    gffFile <- '/project/umw_mccb/OneStopRNAseq/kai/DEXSeq/gff/gencode.v34.primary_assembly.annotation.DEXSeq.gff'
    gffFile <- 'gencode.v34.primary_assembly.annotation.DEXSeq.gff'
    #countFile <- 'DEXSeq_count/N052611_Alb_Dex_count.txt DEXSeq_count/N052611_Alb_count.txt DEXSeq_count/N052611_Dex_count.txt DEXSeq_count/N052611_untreated_count.txt DEXSeq_count/N061011_Alb_Dex_count.txt DEXSeq_count/N061011_Alb_count.txt DEXSeq_count/N061011_Dex_count.txt DEXSeq_count/N061011_untreated_count.txt DEXSeq_count/N080611_Alb_Dex_count.txt DEXSeq_count/N080611_Alb_count.txt DEXSeq_count/N080611_Dex_count.txt DEXSeq_count/N080611_untreated_count.txt DEXSeq_count/N61311_Alb_Dex_count.txt DEXSeq_count/N61311_Alb_count.txt DEXSeq_count/N61311_Dex_count.txt DEXSeq_count/N61311_untreated_count.txt'
  }else{
    # production
    metaFile <- args[1]
    contrastFile <- args[2]
    gffFile <- args[3]
    #countFile <-paste( unlist(args[4:length(args)]), collapse=' ')  
  }



min_count_per_exon <- 0  # test, maybe not valid for DEXSeq if larger than 1
maxFDR <- 0.05  # todo: read config
#minLFC <- 0.585

#countFile <- gsub(' +',' ',countFile) 
print(">>> Parameters: ")
print(paste("path:", getwd()))
paste("metaFile:", metaFile)
paste("contrastFile:", contrastFile)
paste("gffFile:", gffFile)


cat("\n\nreading metaFile:\n")
meta.df <- readExcel(metaFile)
print(meta.df)
cat("\n\nreading contrastFile:\n")
contrast.df <- readExcel(contrastFile)
print(contrast.df)

countFile <- paste('DEXSeq_count/', meta.df$SAMPLE_LABEL, "_count.txt", sep='')
print("countFile:")
print(countFile)
if (!all(file.exists(countFile))){
	stop("Some countFiles Can't be found!!!")
}

## name1 should be treated group
for (i in 1:dim(contrast.df)[2]) { # test
  # parse names
  name1 <- contrast.df[1, i]
  name2 <- contrast.df[2, i]
  name1 <- gsub(" ", "", name1)
  name2 <- gsub(" ", "", name2)
  name1 <- gsub(";$", "", name1)
  name2 <- gsub(";$", "", name2)
  name1s <- strsplit(name1, ";") [[1]]
  name2s <- strsplit(name2, ";") [[1]]
  name1 <- gsub(";", ".", name1)
  name2 <- gsub(";", ".", name2)
  name <- paste(name1, name2, sep = "_vs_")
  cat(paste("\n\n>>> for contrast", i, ":", name, "\n"))

  sampleTable <- data.frame(row.names = meta.df$SAMPLE_LABEL,
                            condition = meta.df$GROUP_LABEL,
                            batch     = meta.df$BATCH)

  idx <- sampleTable$condition %in% c(name1, name2)
  sampleTableSubset <- sampleTable[idx, ]
  countFilesSubset <- countFile[idx]
  
  # Print sampleTable:
  print('countFilesSubset:')
  print(countFilesSubset)
  cat("\nsampleTableSubset: \n")
  print(sampleTableSubset)
  cat("\n")
    
  # Read data
  print("reading data..")
  dxd = DEXSeqDataSetFromHTSeq(
    countFilesSubset,
    sampleData=sampleTableSubset,
    design= ~ sample + exon + condition:exon,
    flattenedfile=gffFile )
  dxd$condition <- relevel(dxd$condition, ref = name2)

  # Filter data (skipped to keep figures in results correct)
  if (min_count_per_exon > 0){
    print(paste("Removing bins/exons with less than ", min_count_per_exon, "reads"))
    dxd <- dxd[rowSums(featureCounts(dxd)) >= min_count_per_exon, ]
  }

  print(dxd)
  #print(head(geneIDs(dxd), 3))
  print(head(featureCounts(dxd), 5))
  
  print("Estimate size factors..")
  dxd = estimateSizeFactors(dxd)
  
  ## if only one batch, don't apply model, otherwise DEXSeq error:
  print("Estimate Dispersion and performing statistical test..")
  if (length(unique(sampleTableSubset$BATCH)) == 1) {
  	print("No Batch Effect")
    dxd = estimateDispersions(dxd, BPPARAM=BPPARAM)
    dxd = testForDEU(dxd, BPPARAM=BPPARAM)
    dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
  } else {
  	print("With Batch Effect")
    formulaFullModel    =  ~ sample + exon + batch:exon + condition:exon
    formulaReducedModel =  ~ sample + exon + batch:exon 
    dxd = estimateDispersions(dxd, formula = formulaFullModel, BPPARAM=BPPARAM)
    dxd = testForDEU(dxd, 
                     reducedModel = formulaReducedModel, 
                     fullModel = formulaFullModel, 
                     BPPARAM=BPPARAM)  
    dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
  }
  
  # Summarize result:
  cat("\nSummarize result:\n")
  dxr = DEXSeqResults(dxd)
  save.image(file=paste(outDir, "contrast", i, ".RData", sep=''))
  
  # Some statistics:
  cat(paste("\n\nnum of DE exons with FDR <", maxFDR, "\n"))
  print(table(dxr$padj < maxFDR))
  cat(paste("\n\nnum of DE genes with FDR <", maxFDR, "\n"))
  print(table(tapply(dxr$padj < maxFDR, dxr$groupID, any)))
  
  # Save result to xlsx:
  fname <- paste(outDir, "DEXSeq_", name, ".xlsx", sep = "")
  ## filter out rows with no padj value, otherwise, the excel might be too huge.
  df_dxr <- data.frame(dxr)
  df_dxr <- df_dxr[!is.na(df_dxr$padj), ]
  colnames(df_dxr) <- gsub('countData.', '', colnames(df_dxr))
  print(paste("Saving results to:", fname))
  try(WriteXLS(df_dxr, row.names = F, fname)) # if only nrow > 0
  
  # QC plots:
  try(plotDispEstsWrapper(dxd, outDir, name))
  try(plotMAWrapper(dxd, ourDir, name))

  # DEXSeqHTML(dxr, path=outDir, file=paste(name, ".report.html", sep=""), 
  #        fitExpToVar="condition", FDR=maxFDR) # not working properlly
  
  # Visualization:
  ## plot for top 5 genes ranked by pvalue
  gene_temp <- df_dxr[order(df_dxr$pvalue),]
  geneList  <- gene_temp$groupID %>% unique()
  
  for (k in 1:5) {
    gene <- geneList[k]

    name_tmp <- paste(name, gene, "top", k, "normalized_counts.pdf", sep=".")
    name_tmp <- gsub("\\+", "_", name_tmp)
    name_tmp <- gsub("\\.", "_", name_tmp)
    name_tmp <- gsub("_pdf$", ".pdf", name_tmp)
    print(paste("saving: ", outDir, name_tmp, sep=""))
    try(plotDEXSeq_norCounts(dxr, gene, outDir, name_tmp))

    ## relative_exon_usage.pdf
    name_tmp <- paste(name, gene, "top", k, "relative_exon_usage.pdf", sep=".")
    name_tmp <- gsub("\\+", "_", name_tmp)
    name_tmp <- gsub("\\.", "_", name_tmp)
    name_tmp <- gsub("_pdf$", ".pdf", name_tmp)
    print(paste("saving: ", outDir, name_tmp, sep=""))
    try(plotDEXSeq_relative_exon_usage(dxr, gene, outDir, name_tmp))
  }
}
