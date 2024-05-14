if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos = 'http://cran.us.r-project.org')
}

if (!requireNamespace("CleanUpRNAseq", quietly = TRUE)) {
  install.packages('./workflow/envs/CleanUpRNAseq/',
                   repos = NULL,
                   type = "source")  # works
  # BiocManager::install(
  #   "CleanUpRNAseq",
  #   update = TRUE,
  #   ask = T,
  #   checkBuilt = FALSE,
  #   force = FALSE,
  #   version = BiocManager::version()
  # )
}

library(CleanUpRNAseq)

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = c("output", "message"))

message("Working Directory:", getwd())

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
