log <- file(snakemake@log[[1]], open="wt")
sinkall(log)

message("Working Directory:", getwd())

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos='http://cran.us.r-project.org')
}

if (!requireNamespace("CleanUpRNAseq", quietly = TRUE)) {
install.packages('./workflow/envs/CleanUpRNAseq/', repos = NULL, type="source")  # works
}

library(CleanUpRNAseq)

meta <- read.csv(snakemake@input[['meta']])
print(meta)

x <- create_diagnostic_plot(
  gtf=snakemake@input[['gtf']],
  metadata = meta,
  out_dir='CleanUpRNAseqQC/',
  ensdb_sqlite = snakemake@input[['ensdb']],
  threads = snakemake@threads,
  verbose = F
)

sink()
