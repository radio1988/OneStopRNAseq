if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = 'http://cran.us.r-project.org')
}

if (!requireNamespace("KernSmooth", quietly = TRUE)) {
  install.packages("KernSmooth", repos = 'http://cran.us.r-project.org')
}


if (!requireNamespace("CleanUpRNAseq", quietly = TRUE)) {
  BiocManager::install(pkgs = "haibol2016/CleanUpRNAseq" , type = 'source', aks = FALSE)
}

library(CleanUpRNAseq)

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = c("output", "message"))

print(paste("Working Directory:", getwd()))  # main root

meta <- read.csv(snakemake@input[['meta']])
print(meta)

print("PAIR_END:")
print(snakemake@config[["PAIR_END"]])

x <- create_diagnostic_plot(
  isPairedEnd = snakemake@config[["PAIR_END"]],
  gtf = snakemake@input[['gtf']],
  metadata = meta,
  out_dir = 'CleanUpRNAseqQC/',
  ensdb_sqlite = snakemake@input[['ensdb']],
  threads = snakemake@threads,
  verbose = F
)

sink()
