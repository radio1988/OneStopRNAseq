library(readr)
library(ggplot2)
library(tidyverse)

setClass('gsea_edb', slots=list(gsea_db='character',  # m1.all.v2022
                                edb_path = 'character', # "gsea/TN_AKT2KO_vs_TN_WT/m1.all.v2022.1.Mm.symbols.gmt.GseaPreranked/edb/results.edb"
                                rnk_file='character', # "../DESeq2_II/gsea/TN_AKT2KO_vs_TN_WT.tar.gz"
                                df='data.frame')) #  [1] "RANKED_LIST"      "TEMPLATE"         "GENESET"  ...

ReadGseaEdb <- function(tar_file, edb_path){
  # tar_file = 'contrast.tar.gz', 
  # edb_path = 'gsea/CSC_1um_12h_vs_CSC_DMSO_24h/df.all.v2022.1.Hs.symbols.gmt.GseaPreranked/edb/results.edb'
  # return gsea_data (df) # sorted by decreasing ABS_NES 
  
  # decompress
  untar(tarfile = tar_file, files = edb_path, exdir = tempdir())
  edb_file_path <- file.path(tempdir(), edb_path)
  
  # read
  lines <- readLines(edb_file_path)
  filtered_lines <- grep("^  ", lines, value = TRUE)
  gsea_data <- read.table(text=filtered_lines, header = F, sep = " ")
  gsea_data <- gsea_data[, c(-1,-2,-3)]
  colnames(gsea_data) <- c('RANKED_LIST', 'TEMPLATE', 'GENESET', 'ES', 'NES', 'NP', 'FDR', 'FWER', 
                           'RND_ES', 'HIT_INDICES','ES_PROFILE', 'RANK_AT_ES', 'RANK_SCORE_AT_ES')
  
  # fix value strings
  for (name in colnames(gsea_data)) {
    redundant <- paste0(name,'=')
    gsea_data[,name] <- gsub(redundant, '', gsea_data[,name])
  }
  
  gsea_data$GENESET <- gsub('gene_sets.gmt#', '', gsea_data$GENESET )
  gsea_data$RANKED_LIST <- gsub('.rnk', '', gsea_data$RANKED_LIST)
  
  # set numeric data types
  numeric_cols <- c('ES', 'NES', 'NP', 'FDR', 'FWER')
  gsea_data[, numeric_cols] <- lapply(gsea_data[, numeric_cols], as.numeric)
  gsea_data[, 'RANK_AT_ES'] <- as.integer(gsea_data[, 'RANK_AT_ES'] )
  
  # Sort
  gsea_data$ABS_NES <- abs(gsea_data$NES)
  gsea_data <- gsea_data[order(gsea_data$ABS_NES, decreasing = T), ]
  
  return (gsea_data) # data frame, see colnames
}

# PARAMS
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  print("Options provided:")
  print(args)
  data_files <- args # a list of '../top10.rnk.txt.tar.gz'..
}else {
  stop("No gsea.tar.gz files were provided.")
}

max.q <- 0.05
e <- 0.001
max.n <- 100 # max num of gene sets in plots, if 0, skips filtering

# PREP
dir.create('multiBubble', showWarnings = FALSE)
dir.create('fdr_nes', showWarnings = FALSE)
dir.create('dataTables', showWarnings = FALSE)
dir.create('mergedTables', showWarnings = FALSE)
n <- length(data_files)

edbpaths1 <- grep("edb/results.edb$", untar(data_files[1], list = T), value = TRUE)
# assume edbpaths are the same among files

for (i in 1:length(edbpaths1)){ # 19 msigdbs
  EDBS <- list()
  EDBS.sig <- list()
  
  for (j in 1:length(data_files)){ # 3 tar.gz files
    rnk_file <- data_files[j]
    # "../../../../DIT/gsea/WT_HFD_vs_WT_MCD.tar.gz"
    
    edbpaths <- grep("edb/results.edb$", untar(rnk_file, list = T), value = TRUE)
    edbpaths <- sort(edbpaths) # keep order for different tar.gz
    edb_path <-  edbpaths[i]
    # "gsea/WT_HFD_vs_WT_MCD/m1.all.v2022.1.Mm.symbols.gmt.GseaPreranked/edb/results.edb"
    
    name <- basename(paste0(unlist(strsplit(x = edb_path, split = '\\.'))[1:3], collapse = '.'))
    # "m1.all.v2022"
    
    print(paste0('gsea_db', i, ': ', name))
    print(paste0('rnk_file', j, ": ", rnk_file))
    print(paste0('edb_path: ', edb_path))
    
    df <- ReadGseaEdb(rnk_file, edb_path) # sorted by FDR
    write.csv(df, paste0('dataTables/', j, '.', basename(rnk_file),'.', name, '.csv'), row.names = F)
    
    df.sig <- df[df$FDR < max.q, ]
    
    edb <- new('gsea_edb', gsea_db = name, edb_path = edb_path, rnk_file = rnk_file, df = df)
    edb.sig <- new('gsea_edb', gsea_db = name, edb_path = edb_path, rnk_file = rnk_file, df = df.sig)
    
    EDBS[[j]] <- edb
    EDBS.sig[[j]] <- edb.sig
  }
  
  # check gsea_db the same for different tar.gz files
  GSEA_DBS <- unlist(lapply(EDBS, function(edb) edb@gsea_db))
  if (length(unique(GSEA_DBS)) != 1){ 
    warning(paste0(GSEA_DBS, sep = '; '))
    stop('gsea_db not the same for all rnk_files')
  }else{
    print(unique(GSEA_DBS))
  }
  
  # PREP
  DFS <-  lapply(EDBS, function(edb) edb@df)
  DFS.SIG <- lapply(EDBS.sig, function(edb) edb@df)
  
  genesets.anysig <- unique(unlist(lapply(EDBS.sig, function(edb) edb@df$GENESET)))
  
  df.merged <- DFS %>% reduce(inner_join, by='GENESET')
  
  if ( length(grep('^FDR', colnames(df.merged))) > 1){
    df.merged$MIN_FDR <- apply(df.merged[, grep('^FDR', colnames(df.merged))], 1, min)
  }else{
    df.merged$MIN_FDR <- df.merged[, grep('^FDR', colnames(df.merged))]
  }
  
  if (length(grep('^ABS_NES', colnames(df.merged))) > 1){
    df.merged$MAX_ABS_NES <- apply(df.merged[, grep('^ABS_NES', colnames(df.merged))], 1, max)
  }else{
    df.merged$MAX_ABS_NES <-df.merged[, grep('^ABS_NES', colnames(df.merged))]
  }
  
  df.merged <- df.merged[order(df.merged$MAX_ABS_NES, decreasing = T),]
  df.union <- dplyr::filter(df.merged, GENESET %in% genesets.anysig)
  df.union[, grep('ABS_NES', colnames(df.union))]
  write.csv(df.union, paste0('mergedTables/union.', name, '.csv'), row.names = F)
  
  print(paste("There are", dim(df.merged)[1], 'genesets sig in ≥1 GSEA'))
  print(paste("There are", dim(df.union)[1], 'such genesets has data in all GSEA'))
  
  if (max.n > 0 & dim(df.union)[1] > max.n){
    df.union <- df.union[1:max.n, ]
  }
  print(paste("There are", dim(df.union)[1], 'genesets plotted'))
  
  # Multiple Bubble Plot Prep
  FDR <-  unname(unlist(df.union[, grep("FDR", colnames(df.union))[-(n+1)]]))
  logFDR <- unlist(lapply(unname(unlist(df.union[, grep("FDR", colnames(df.union))[-(n+1)]])), function(x) -log10(x+e)))
  logFDR_ <- apply(data.frame(logFDR, rep(0, length(logFDR))), 1, max)
  # data.frame(FDR, logFDR, logFDR_)
  data <- data.frame(
    GeneSet = rep(df.union$GENESET, n),
    Sig = logFDR_,
    EnrichmentScore =  unname(unlist(df.union[, grep("^NES", colnames(df.union))])),
    Comparison =  unname(unlist(df.union[, grep("^RANKED_LIST", colnames(df.union))]))
  )
  
  # Plotting
  mplot <- ggplot(
    data, aes(y = GeneSet, x = EnrichmentScore, size = Sig, color = Comparison)) +
    geom_point(alpha = 0.7) +
    scale_radius(limits = c(0, 3), range = c(1,3)) +
    labs(title = "GSEA Bubble Plot",
         y = "Gene Set",
         x = "Normalized Enrichment Score",
         size = paste0('Significance: -log10(FDR+',e,')'),
         color = "Comparison") +
    theme_minimal() 
  
  ggsave(paste0('multiBubble/', name, '.MultiBubblePlot.pdf'), 
         mplot,    
         width = max(nchar(data$GeneSet))/10 + 5,
         height = nrow(data)/22 + 1, 
         limitsize = FALSE)
  
  stat.p <- ggplot(data, aes(x=EnrichmentScore, y=Sig)) + 
    geom_point() + geom_hline(yintercept=1.3, linetype="dashed", 
                              color = "red", size=2)
  ggsave(paste0('fdr_nes/', name, '.FDR_NES.pdf'), stat.p, width = 5, height = 5)
}






