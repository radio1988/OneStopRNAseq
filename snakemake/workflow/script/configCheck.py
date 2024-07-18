import sys

output = open(snakemake.output[0], "w")
log = open(snakemake.log[0], "w")


# INTRON and gDNA correction, only can choose one at most
if snakemake.config['START'] == 'FASTQ' and snakemake.config["INTRON"] and snakemake.config["CleanUpRNAseqCorrection"]:
    message = "INTRON mode and CleanUpRNAseqCorrection incompatible. \
    gDNA correction is only possible for exon level rnaseq quantification"
    print(message, file=log)
    sys.exit(message)
else:
    print("INTRON {} and CleanUpRNAseqCorrection {} compatible".format(
        snakemake.config['INTRON'], snakemake.config['CleanUpRNAseqCorrection'],
        file = output
    ))


# SPECIES and Analysis options
# Not now
