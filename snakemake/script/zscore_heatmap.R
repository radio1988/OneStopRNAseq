library(WriteXLS)
library(pheatmap)

readExcel <- function(fname){
  df <- readxl::read_xlsx( fname, na='-', sheet=1)  # na term important
  df <- data.frame(df)  #important
  return (df)
}

zscore <- function(matrix){
  return( t(scale(t(matrix))))
}

Heatmap <- function(df, nclass=4, fname="heatmap", main="title"){
  # Heatmap Pre-plot to get hclust
  nclass <- 4
  p <- pheatmap(df, 
                main = "P300 Down, RNA Sig, H3K27me3Sig (156 Genes) ",
                cluster_cols = F, 
                cutree_rows = nclass, 
                show_rownames = F, 
                width = 4, height = 6)
  
  # Print out gene classification (https://www.biostars.org/p/287512/)
  # cutree manually
  gene_classes <- sort(cutree(p$tree_row, k=nclass))
  gene_classes <- data.frame(gene_classes) 
  classification <- merge(gene_classes, df, by = 0)
  row.names(classification) <- classification$Row.names
  #Re-order original data (genes) to match ordering in heatmap (top-to-bottom)
  idx <- rownames(df[p$tree_row[["order"]],])
  classification <- classification[idx,] 
  WriteXLS(classification, row.names = F,
           paste(fname,"gene_class.xlsx", sep = "."))
  
  # Get classification labled heatmap
  annotation_row = data.frame(
    class = classification$gene_classes
  )
  row.names(annotation_row) <- row.names(classification)
  annotation_row$class <- as.character(annotation_row$class)
  
  p2 <- pheatmap(df, 
                 main = main,
                 cluster_cols = F, 
                 cutree_rows = nclass, 
                 show_rownames = F, 
                 annotation_row = annotation_row,
                 filename = paste(fname, "pdf", sep = "."), 
                 width = 4, height = 6)
}


# read data and zscore calculation
df <- readExcel("k4.sig.86.xlsx")

df1 <- df[,c(2:5)]
row.names(df1) <- df$Name
head(df1)
tail(df1)
df1.zscore <- zscore(df1)

df2 <- df[, c(6:9)]
row.names(df2) <- df$Name
head(df2)
tail(df2)
df2.zscore <- zscore(df2)
head(df2.zscore)

df.zscore <- cbind(df1.zscore, df2.zscore)
head(df.zscore)


Heatmap(df.zscore, nclass = 4, 
        fname = "k4.86.zheatmap", main = "P300 Down, RNA Sig, H3K4me3 Sig (86)")
