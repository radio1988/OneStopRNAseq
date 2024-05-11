message("Working Directory:", getwd())

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos='http://cran.us.r-project.org')
}

if (!requireNamespace("CleanUpRNAseq", quietly = TRUE)) {
install.packages('./CleanUpRNAseq/', repos = NULL, type="source")  # works
}

library(CleanUpRNAseq)
library(readr)

meta <- read_table("./meta.txt")
print(meta)

create_diagnostic_plot(
  gtf=snakemake@input[['gtf']],
  metadata = snakemake@input[['meta']],,
  out_dir='CleanUpRNASeqQC/plots',
  ensdb_sqlite = snakemake@input[['sqlite']],
  threads = snakemake@threads,
  verbose = F
)
