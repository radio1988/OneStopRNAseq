import sys

output = open(snakemake.output[0], "w")
log = open(snakemake.log[0], "w")


# INTRON and gDNA correction, only can choose one at most
if snakemake.config['START'] == 'FASTQ' and snakemake.config["INTRON"] and snakemake.config["CleanUpRNAseqCorrection"]:
    message = "INTRON mode and CleanUpRNAseqCorrection incompatible. \
    gDNA correction is only possible for exon level rnaseq quantification"
    print(message, file=log)
    sys.exit(message)

# SPECIES and Analysis options
# Not now


# meta.csv
# todo: sample name can't contain count, abundance, length
