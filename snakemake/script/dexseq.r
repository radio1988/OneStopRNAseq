options(show.error.locations = TRUE)
options(error=traceback)

suppressPackageStartupMessages(library(readxl))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(WriteXLS))
suppressPackageStartupMessages(library(BiocParallel))
sessionInfo()


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



args = commandArgs(trailingOnly=TRUE)

# Read ARGS
# todo: remove countFile from args
if (length(args) < 2){
    # development
    metaFile <- './meta/meta.xlsx'  
    contrastFile <- './meta/contrast.as.xlsx'
    gffFile <- "/project/umw_mccb/genome/Mus_musculus_UCSC_mm10/gencode.vM25.primary_assembly.annotation.gtf.dexseq.gff"  # test
    annoFile <- "https://raw.githubusercontent.com/hukai916/Collections/master/gencode.vM21.annotation.txt"
    maxFDR <- 0.1
    minLFC <- 0.585
    threads <- 4
    MIN_GENE_COUNT <- 20
    idx_contrast <- 1
    #countFile <- 'DEXSeq_count/N052611_Alb_Dex_count.txt DEXSeq_count/N052611_Alb_count.txt DEXSeq_count/N052611_Dex_count.txt DEXSeq_count/N052611_untreated_count.txt DEXSeq_count/N061011_Alb_Dex_count.txt DEXSeq_count/N061011_Alb_count.txt DEXSeq_count/N061011_Dex_count.txt DEXSeq_count/N061011_untreated_count.txt DEXSeq_count/N080611_Alb_Dex_count.txt DEXSeq_count/N080611_Alb_count.txt DEXSeq_count/N080611_Dex_count.txt DEXSeq_count/N080611_untreated_count.txt DEXSeq_count/N61311_Alb_Dex_count.txt DEXSeq_count/N61311_Alb_count.txt DEXSeq_count/N61311_Dex_count.txt DEXSeq_count/N61311_untreated_count.txt'
  }else{
    # production
    metaFile <- args[1]
    contrastFile <- args[2]
    gffFile <- args[3]
    annoFile <- args[4]
    maxFDR <- as.numeric(args[5])
    minLFC <- as.numeric(args[6])
    threads <- as.numeric(args[7])
    MIN_GENE_COUNT <- as.numeric(args[8])  # todo: from args
    idx_contrast <- as.numeric(args[9])  # index of contrast, e.g. 1,2,3
    #countFile <-paste( unlist(args[4:length(args)]), collapse=' ')  
  }


print(">>> Parameters: ")
print(paste("path:", getwd()))
paste("metaFile:", metaFile)
paste("contrastFile:", contrastFile)
paste("gffFile:", gffFile)
paste("annoFile:", annoFile)
paste("maxFDR:", maxFDR)
paste("minLFC:", minLFC)
paste("threads:", threads)
paste("MIN_GENE_COUNT:", MIN_GENE_COUNT)
paste("index of contrast:", idx_contrast)






paste("idx_contrast:", idx_contrast )


BPPARAM = MulticoreParam(workers=threads) 
outDir <- paste("./DEXSeq/contrast",idx_contrast,"/", sep="")
dir.create( "./DEXSeq/", showWarnings = FALSE)
dir.create( outDir, showWarnings = FALSE)

# Importing annotation from GitHub (must be raw, not zipped)
getAnnotation <- function(urlpath) {
  tmp <- tempfile()
  download.file(urlpath, destfile = tmp, method = 'auto')
  return(read.table(tmp, sep="\t", header = TRUE))
}

paste("Getting data from", annoFile)

if (grepl('https://', annoFile)){
  print("downloading remote annoFile")
  anno <- getAnnotation(annoFile)
}else{
  print("reading local annoFile")
  anno <- read.table(annoFile, sep = "\t", header = T)
}
print("Dimention of annotation table: ")
dim(anno)
head(anno)



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




# Main code: for each contrast
## name1 should be treated group

# parse names
name1 <- contrast.df[1, idx_contrast]
name2 <- contrast.df[2, idx_contrast]
cat(paste("name1:", name1, "\nname2:", name2))
name1 <- gsub(" ", "", name1)
name2 <- gsub(" ", "", name2)
name1 <- gsub(";$", "", name1)
name2 <- gsub(";$", "", name2)
name1s <- strsplit(name1, ";") [[1]]
name2s <- strsplit(name2, ";") [[1]]
name1 <- gsub(";", "_", name1)
name2 <- gsub(";", "_", name2)
name <- paste(name1, name2, sep = "_vs_")
cat(paste("\n\n>>> for contrast", idx_contrast, ":", name, "\n"))
cat(paste("name1:", name1, "\nname2", name2, "\n"))

# create sampleTable
sampleTable <- data.frame(
  row.names = meta.df$SAMPLE_LABEL,
  condition = meta.df$GROUP_LABEL,
  batch     = meta.df$BATCH
)
sampleTable

# filter sampleTable
idx <- sampleTable$condition %in% c(name1s, name2s)
sampleTableSubset <- sampleTable[idx,]
countFilesSubset <- countFile[idx]
sampleTableSubset
countFilesSubset

# change condition names for complex condition contrast
if (length(name1s) > 1) {
  levels(sampleTableSubset$condition)[levels(sampleTableSubset$condition) %in% name1s] <-
    name1
  sampleTableSubset$condition[sampleTableSubset$condition %in% name1s] <-
    name1
  # todo: fix batch
}
if (length(name2s) > 1) {
  levels(sampleTableSubset$condition)[levels(sampleTableSubset$condition) %in% name2s] <-
    name2
  sampleTableSubset$condition[sampleTableSubset$condition %in% name2s] <-
    name2
}
sampleTableSubset


# Print sampleTable:
print('countFilesSubset:')
print(countFilesSubset)
cat("\nsampleTableSubset: \n")
print(sampleTableSubset)
cat("\n")

# Read data
print("reading data..")
dxd <- DEXSeqDataSetFromHTSeq(
  countFilesSubset,
  sampleData = sampleTableSubset,
  design = ~ sample + exon + condition:exon,
  flattenedfile = gffFile
)
dxd$condition <- relevel(dxd$condition, ref = name2)
print("Data summary:")
print(dxd)

# Filter lowly expressed genes
cat(paste("removing genes with less than", MIN_GENE_COUNT, "counts per sample\n"))
if (MIN_GENE_COUNT > 0) {
  jump <- round(dim(dxd)[2] / 2)
  gene_count <- rowSums(counts(dxd)[, c(1, 1 + jump)])
  large_gene_idx <-
    gene_count >= MIN_GENE_COUNT * dim(sampleTableSubset)[1]  # per sample
  #large_gene_idx <- gene_count >= MIN_GENE_COUNT
  cat(paste("there are", sum(large_gene_idx), "genes will survive the filter\n"))
  if (sum(large_gene_idx) > 30){
    dxd <- dxd[large_gene_idx,]
    print("filtering performed")
  }else{
    warning("not enough genes has enough reads, filtering stopped")
    warning("please try setting a smaller MIN_GENE_COUNT parameter, or perform deeper sequencing")
    quit(status=1)
  }
}
print("Data summary after fitering:")
print(dxd)
print(head(featureCounts(dxd), 5))

print("Estimate size factors..")
dxd = estimateSizeFactors(dxd)

## if only one batch, don't apply model, otherwise DEXSeq error:
print("Estimate Dispersion and performing statistical test..")
if (length(unique(sampleTableSubset$batch)) == 1) {
  print("No Batch Effect")
  dxd = estimateDispersions(dxd, BPPARAM = BPPARAM)
  dxd = testForDEU(dxd, BPPARAM = BPPARAM)
  dxd = estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM =
                                  BPPARAM)
} else {
  print("With Batch Effect")
  formulaFullModel    =  ~ sample + exon + batch:exon + condition:exon
  formulaReducedModel =  ~ sample + exon + batch:exon
  dxd = estimateDispersions(dxd, formula = formulaFullModel, BPPARAM =
                              BPPARAM)
  dxd = testForDEU(
    dxd,
    reducedModel = formulaReducedModel,
    fullModel = formulaFullModel,
    BPPARAM = BPPARAM
  )
  dxd = estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM =
                                  BPPARAM)
}

# Summarize result:
cat("\nSummarize result:\n")
dxr = DEXSeqResults(dxd)
save.image(file = paste(outDir, "contrast", idx_contrast, ".RData", sep =
                          ''))

# Some statistics:
cat(paste("\n\nnum of DE exons with FDR <", maxFDR, "\n"))
print(table(dxr$padj < maxFDR))
cat(paste("\n\nnum of DE genes with FDR <", maxFDR, "\n"))
print(table(tapply(dxr$padj < maxFDR, dxr$groupID, any)))

# Save result to xlsx:
fname <- paste(outDir, "DEXSeq_", name, ".xlsx", sep = "")
## filter out rows with no padj value, otherwise, the excel might be too huge.
df_dxr <- data.frame(dxr)
## Annotate results
if (sum(anno[, 1] %in% dxr$groupID) < 1) {
  warning("!!! Annotation file and count file have no ID in common")
  warning("The results will be unannotated")
  warning("Please Double check Annotation file")
  print("count table ID:")
  print(row.names(cts)[1:2])
  print("anno table ID:")
  print(anno[1:2, 1])
}
df_dxr <- merge(
  df_dxr,
  anno,
  by.x = 'groupID',
  by.y = 1,
  sort = F,
  all.1 = T
)

# save results
df_dxr <- df_dxr[!is.na(df_dxr$padj),]
colnames(df_dxr) <- gsub('countData.', '', colnames(df_dxr))
print(paste("Saving results to:", fname))
try(WriteXLS(df_dxr, row.names = F, fname))
# if only nrow > 0

# save sig results
lfc_idx <-
  abs(df_dxr[, grepl("log2fold_", colnames(df_dxr))]) > minLFC
fdr_idx <- df_dxr$padj < maxFDR
df_dxr.sig <- df_dxr[lfc_idx & fdr_idx,]
fname.sig <- paste(outDir, "DEXSeq_", name, ".sig.xlsx", sep = "")
print(paste("Saving significant results to:", fname.sig))
try(WriteXLS(df_dxr.sig, row.names = F, fname.sig))
# if only nrow > 0

# QC plots:
try(plotDispEstsWrapper(dxd, outDir, name))
try(plotMAWrapper(dxd, ourDir, name))

# DEXSeqHTML(dxr, path=outDir, file=paste(name, ".report.html", sep=""),
#        fitExpToVar="condition", FDR=maxFDR) # not working properlly

# Visualization:
## plot for top 5 genes ranked by pvalue
if (dim(df_dxr.sig)[1] > 0) {
  gene_temp <- df_dxr.sig[order(df_dxr.sig$pvalue), ]
} else{
  gene_temp <- df_dxr[order(df_dxr$pvalue), ]
}
geneList  <- gene_temp$groupID %>% unique()

for (k in 1:5) {
  gene <- geneList[k]
  gene_symbol <-  anno$Name[anno$Gene == gene]
  if (is.na(gene_symbol)) {
    gene_symbol <- gene
  }
  gene_symbol <- gsub(" ", "", gene_symbol)
  
  name_tmp <-
    paste(name, gene_symbol, "top", k, "normalized_counts.pdf", sep = ".")
  name_tmp <- gsub("\\+", "_", name_tmp)
  name_tmp <- gsub("\\.", "_", name_tmp)
  name_tmp <- gsub("_pdf$", ".pdf", name_tmp)
  print(paste("saving: ", outDir, name_tmp, sep = ""))
  try(plotDEXSeq_norCounts(dxr, gene, outDir, name_tmp))
  
  ## relative_exon_usage.pdf
  name_tmp <-
    paste(name, gene_symbol, "top", k, "relative_exon_usage.pdf", sep = ".")
  name_tmp <- gsub("\\+", "_", name_tmp)
  name_tmp <- gsub("\\.", "_", name_tmp)
  name_tmp <- gsub("_pdf$", ".pdf", name_tmp)
  print(paste("saving: ", outDir, name_tmp, sep = ""))
  try(plotDEXSeq_relative_exon_usage(dxr, gene, outDir, name_tmp))
}


