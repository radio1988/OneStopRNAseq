message("Working Directory:", getwd())

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos='http://cran.us.r-project.org')
}

if (!requireNamespace("CleanUpRNAseq", quietly = TRUE)) {
install.packages('./CleanUpRNAseq/', repos = NULL, type="source")  # works
}

log <- file(snakemake@log[[1]], open="wt")
sink(log)

library(CleanUpRNAseq)

make_ensdb(
  gtf=snakemake@input[['gtf']],
  outfile = snakemake@output[['ensdb']]
  organism_latin_name = "NA",
  genome_version = "NA",
  Ensembl_release_version = "NA"
)

#todo: correct version numbers
