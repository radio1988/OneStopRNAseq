log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = c("output", "message"))

library(CleanUpRNAseq)

print(paste("Working Directory:", getwd()))

meta <- readRDS(snakemake@input[['meta_rds']])
qc <- readRDS(snakemake@input[['qc_rds']])

meta.df <- read.csv(snakemake@input[['meta_tab']])
print(meta.df)

print(paste("PAIR_END:", snakemake@config[["PAIR_END"]]))
if (snakemake@config[["PAIR_END"]]){
    corrected <- correct_for_contamination(
        is_stranded_library = F,
        unstranded_metadata = meta.df,
        salmon_summary = qc$salmon_summary,
        correction_method = 'Global',
        featurecounts_summary = qc$featurecounts_summary,
        saf_list = qc$saf_list,
        ensdb_sqlite = snakemake@input[['ensdb']])
} else{
    corrected <- correct_for_contamination(
        is_stranded_library = T,
        stranded_metadata = meta.df,
        salmon_summary = qc$salmon_summary,
        correction_method = 'Global',
        featurecounts_summary = qc$featurecounts_summary,
        saf_list = qc$saf_list,
        ensdb_sqlite = snakemake@input[['ensdb']])
}

write.csv(corrected, snakemake@output[[1]], quote=F)


close(log)
