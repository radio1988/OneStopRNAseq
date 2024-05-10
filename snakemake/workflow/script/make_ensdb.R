print('test')


library(CleanUpRNAseq)



make_ensdb(
  gtf=snakemake@input[['gtf']],
  outfile = snakemake@output[['ensdb']]
  organism_latin_name = "NA",
  genome_version = "NA",
  Ensembl_release_version = "NA"
)

#todo: correct version numbers
