if (!requireNamespace("pehatmap", quietly = TRUE)) {
  install.packages("pehatmap", repos='http://cran.us.r-project.org')
}

working_dir <- getwd()
print(working_dir)

message("Working Directory:", getwd())


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
