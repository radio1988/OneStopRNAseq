

# CHECK CONFIG.YAML
if config["MODE"] == 'liberal' and config["CleanUpRNAseqCorrection"]:
    pass


if config["INTRON"] and config["CleanUpRNAseqCorrection"]:
    pass


def CheckTrimmedFiles_Input():
    L = []
    if config['STRAND'] == 0:
        L = ["trimmed/{sample}.R1.fastq.gz".format(sample) for sample in SAMPLES]
    else:
        L = ["trimmed/{sample}.R1.fastq.gz".format(sample) for sample in SAMPLES]
        L.extend(["trimmed/{sample}.R2.fastq.gz".format(sample) for sample in SAMPLES])

# CHECK FILES
rule CheckTrimmedFiles:
    input:
        CheckTrimmedFiles_Input
    output:
        'meta/CheckFiles.txt'
    script:
        "../script/CheckFiles.py"
