if (!requireNamespace("pehatmap", quietly = TRUE)) {
  install.packages("pehatmap", repos='http://cran.us.r-project.org')
}

working_dir <- getwd()
print(working_dir)

message("Working Directory:", getwd())


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
