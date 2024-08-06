log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = c("output", "message"))

library(CleanUpRNAseq)

print(paste("Working Directory:", getwd()))

meta <- readRDS(snakemake@input[['meta_rds']])
qc <- readRDS(snakemake@input[['qc_rds']])

meta.df <- read.csv(snakemake@input[['meta_tab']])
print(meta.df)

strand_string <- scan(snakemake@input[['strandness']], what = character(), nlines = 1)
print(paste('strand_string:', strand_string))
if (strand_string == "feature_count/counts.s0.strict.txt.summary" | strand_string == "feature_count/counts.s0.liberal.txt.summary") {
    stranded <- FALSE
} else if (strand_string == "feature_count/counts.s1.strict.txt.summary" | strand_string == "feature_count/counts.s1.liberal.txt.summary"){
    stranded <- TRUE
} else if (strand_string == "feature_count/counts.s2.strict.txt.summary" | strand_string == "feature_count/counts.s2.liberal.txt.summary"){
    stranded <- TRUE
} else {
    stop("strand_string not recognized")
}

if (stranded){
    corrected <- correct_for_contamination(
        is_stranded_library = T,
        stranded_metadata = meta.df,
        salmon_summary = qc$salmon_summary,
        correction_method = 'Global',
        featurecounts_summary = qc$featurecounts_summary,
        saf_list = qc$saf_list,
        ensdb_sqlite = snakemake@input[['ensdb']])
} else{
    corrected <- correct_for_contamination(
        is_stranded_library = F,
        unstranded_metadata = meta.df,
        salmon_summary = qc$salmon_summary,
        correction_method = 'Global',
        featurecounts_summary = qc$featurecounts_summary,
        saf_list = qc$saf_list,
        ensdb_sqlite = snakemake@input[['ensdb']])
}

write.csv(corrected, snakemake@output[[1]], quote=F)


close(log)
