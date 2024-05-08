library(CleanUpRNAseq)
library(readr)

meta <- read_table("./meta.txt")
print(meta)

create_diagnostic_plot(
  gtf=snakemake@input[['gtf']],
  metadata = snakemake@input[['meta']],,
  out_dir='CleanUpRNASeqQC',
  ensdb_sqlite = snakemake@input[['sqlite']],
  threads = snakemake@threads,
  verbose = F
)
