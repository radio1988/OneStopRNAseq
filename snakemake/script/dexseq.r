options(show.error.locations = TRUE)
options(error=traceback)

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(BiocParallel))
sessionInfo()

BPPARAM = MulticoreParam(workers=4)

readExcel <- function(fname){
  df <- readxl::read_xlsx( fname, na='NA', sheet=1)  # na term important
  df <- data.frame(df)  #important
  return (df)
}

outDir <- "./DEXSeq/"
args = commandArgs(trailingOnly=TRUE)


if (length(args) < 2){
    # development
    metaFile <- 'meta/meta.xlsx'  
    contrastFile <- './meta/contrast.de.xlsx'
    gffFile <- '/project/umw_mccb/OneStopRNAseq/kai/DEXSeq/gff/gencode.v34.primary_assembly.annotation.DEXSeq.gff'
    countFile <- 'DEXSeq_count/N052611_Alb_Dex_count.txt DEXSeq_count/N052611_Alb_count.txt DEXSeq_count/N052611_Dex_count.txt DEXSeq_count/N052611_untreated_count.txt DEXSeq_count/N061011_Alb_Dex_count.txt DEXSeq_count/N061011_Alb_count.txt DEXSeq_count/N061011_Dex_count.txt DEXSeq_count/N061011_untreated_count.txt DEXSeq_count/N080611_Alb_Dex_count.txt DEXSeq_count/N080611_Alb_count.txt DEXSeq_count/N080611_Dex_count.txt DEXSeq_count/N080611_untreated_count.txt DEXSeq_count/N61311_Alb_Dex_count.txt DEXSeq_count/N61311_Alb_count.txt DEXSeq_count/N61311_Dex_count.txt DEXSeq_count/N61311_untreated_count.txt'
  }else{
    # production
    metaFile <- args[1]
    contrastFile <- args[2]
    gffFile <- args[3]
    countFile <-paste( unlist(args[4:length(args)]), collapse=' ')  
  }


# min_count_per_exon <- 10000  # test, maybe not valid for DEXSeq if larger than 1
maxFDR <- 0.1
#minLFC <- 0.585

countFile <- gsub(' +',' ',countFile) 
print(">>> Parameters: ")
print(paste("path:", getwd()))
paste("metaFile:", metaFile)
paste("contrastFile:", contrastFile)
paste("gffFile:", gffFile)
print("countFile:")
print(countFile)

cat("\n\nreading metaFile:\n")
meta.df <- readExcel(metaFile)
print(meta.df)
cat("\n\nreading contrastFile:\n")
contrast.df <- readExcel(contrastFile)
print(contrast.df)

## name1 should be treated group
tem <- dim(contrast.df)[2]
for (i in 1:tem) {
  name1 <- contrast.df[1,i]
  name2 <- contrast.df[2,i]
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
  
  sample <- factor(meta.df$SAMPLE_LABEL)
  batch <- factor(meta.df$BATCH)
  
  row_names     <- c()
  condition_col <- c()
  batch_col     <- c()
  countFiles    <- unlist(strsplit(countFile, " "))
  #countFiles    <- paste("./", unlist(strsplit(countFile, " ")), sep="")
  countFilesSubset <- c()
  
  for (i2 in 1:length(sample)) {
    for (j in 1:length(name1s)) {
      if (grepl(name1s[j], sample[i2], fixed=TRUE)) {
        row_names <- c(row_names, as.character(sample[i2]))
        condition_col <- c(condition_col, name1) # GROUP1
        batch_col <- c(batch_col, as.character(batch[i2]))
        countFilesSubset <- c(countFilesSubset, countFiles[i2])
      }
    }
    for (l in 1:length(name2s)) {
      if (grepl(name2s[l], sample[i2], fixed=TRUE)) {
        row_names <- c(row_names, as.character(sample[i2]))
        condition_col <- c(condition_col, name2) # GROUP2, check when ; in contrast GROUP labels
        batch_col <- c(batch_col, as.character(batch[i2]))
        countFilesSubset <- c(countFilesSubset, countFiles[i2])
      }
    }
  }
  sampleTable <- data.frame(row.names = row_names,
                            condition = condition_col,
                            batch     = batch_col)
  
  # Print sampleTable:
  cat("\nsample table: \n")
  print(sampleTable)
  cat("\n")
  
  # Run DEXseq analysis below for each sampleTable:
  

  
  # Read in data for DEXSeq analysis:
  print(paste("reading:", countFilesSubset))
  dxd = DEXSeqDataSetFromHTSeq(
    countFilesSubset,
    sampleData=sampleTable,
    design= ~ sample + exon + condition:exon,
    flattenedfile=gffFile )

  print(dxd)

  # print(paste("Removing bins/exons with less than ", min_count_per_exon, "reads"))
  # dxd <- dxd[rowSums(featureCounts(dxd)) >= min_count_per_exon, ]  # test: correct ?
  # print(dxd)
  #print(head(geneIDs(dxd), 3))
  print(head(featureCounts(dxd), 5))
  
  print("Estimate size factors..")
  dxd = estimateSizeFactors(dxd)
  
  ## if only one batch, don't apply model, otherwise DEXSeq error:
  print("Estimate Dispersion and performing statistical test..")
  if (unique(batch) == 1) {
    dxd = estimateDispersions(dxd, BPPARAM=BPPARAM)
    dxd = testForDEU(dxd, BPPARAM=BPPARAM)
    dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
  } else {
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
  
  # Some statistics:
  ## How many DEXs with a FDR of 10%?
  cat(paste("\n\nnum of DE exons with FDR <", maxFDR, "\n"))
  print(table(dxr$padj < maxFDR))
  ## How many genes?
  cat(paste("\n\nnum of DE genes with FDR <", maxFDR, "\n"))
  print(table(tapply(dxr$padj < maxFDR, dxr$groupID, any)))
  
  # Save result to xlsx:
  fname <- paste(outDir, "DEXSeq_", name, ".xlsx", sep = "")
  ## filter out rows with no padj value, otherwise, the excel might be too huge.
  df_dxr <- data.frame(dxr)
  df_dxr <- df_dxr[!is.na(df_dxr$padj), ]
  print(paste("Saving results to:", fname))
  WriteXLS(df_dxr, row.names = F, fname)
  
  # PlotMA plot:
  # plotMA(dxr) # comment out becase for small testset might pop error
  pdf(paste(outDir, name, ".disp.pdf", sep=""))
  plotDispEsts( dxd )
  dev.off()

  pdf(paste(outDir, name, ".MA.pdf", sep=""))
  plotMA( dxr, cex=0.8 ,  alpha = maxFDR) # contain NA error message
  dev.off()
  
  # Visualization:
  ## plot for top 5 genes ranked by pvalue
  gene_temp <- df_dxr[order(df_dxr$pvalue),]
  geneList  <- gene_temp$groupID %>% unique()
  
  for (k in 1:5) {
    gene <- geneList[k]
    #print(paste("Below is for contrast: ", name, sep=""))
    
    ## plotDEXseq
    # name_tmp <- paste(name, gene, "top", k, "standard.pdf", sep=".")
    # print(paste("saving: ", outDir, name_tmp, sep=""))
    # pdf(paste(outDir, name_tmp, sep=""))
    # plotDEXSeq( dxr, gene, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
    # dev.off()
    # ## Transcripts (most useful)
    # name_tmp <- paste(name, gene, "top", k, "transcript.pdf", sep=".")
    # print(paste("saving: ", outDir, name_tmp, sep=""))
    # pdf(paste(outDir, name_tmp, sep=""))
    # plotDEXSeq(dxr, gene, displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
    # dev.off()
    ## Normalized counts
    name_tmp <- paste(name, gene, "top", k, "normalized_counts.pdf", sep=".")
    print(paste("saving: ", outDir, name_tmp, sep=""))
    pdf(paste(outDir, name_tmp, sep=""))
    plotDEXSeq(dxr, gene, expression=FALSE, norCounts=TRUE,displayTranscripts=TRUE,
               legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
    dev.off()
    # ## Fitted splicing
    # name_tmp <- paste(name, gene, "top", k, "fitted_splicing.pdf", sep=".")
    # print(paste("saving: ", outDir, name_tmp, sep=""))
    # pdf(paste(outDir, name_tmp, sep=""))
    # plotDEXSeq(dxr, gene, expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2)
    # dev.off()
  }
}
