log <- file(snakemake@log[[1]], open="wt")
sink(log)

message("Working Directory:", getwd())  # main workdir

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap", repos='http://cran.us.r-project.org')
}

if (!requireNamespace("CleanUpRNAseq", quietly = TRUE)) {
  BiocManager::install(pkgs = "haibol2016/CleanUpRNAseq" , type = 'source', aks = FALSE)
}

library(CleanUpRNAseq)

make_ensdb(
  gtf=snakemake@input[['gtf']],
  outfile = snakemake@output[['ensdb']],
  organism_latin_name = "NA",
  genome_version = "NA",
  Ensembl_release_version = "NA"
)

#todo: correct version numbers
