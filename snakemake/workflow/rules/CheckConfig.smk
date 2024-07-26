import sys

# CHECK CONFIG.YAML
if config["MODE"] == 'liberal' and config["CleanUpRNAseqCorrection"]:
    pass

if config["INTRON"] and config["CleanUpRNAseqCorrection"]:
    pass


def CheckTrimmedFiles_Input(config, SAMPLES):
    L = []
    if config['PAIR_END']:
        L = ["trimmed/{}.R1.fastq.gz".format(sample) for sample in SAMPLES]
        L.extend(["trimmed/{}.R2.fastq.gz".format(sample) for sample in SAMPLES])
    else:
        L = ["trimmed/{}.fastq.gz".format(sample) for sample in SAMPLES]
    return (L)


# CHECK FILES
rule CheckTrimmedReadFiles:
    """check trimmed.fastq.gz for all samples, pass if all non empty (>100bytes)"""
    input:
        CheckTrimmedFiles_Input(config, SAMPLES)
    output:
        'fastqc/CheckFile/all.txt'
    log:
        "fastqc/CheckFile/all.log"
    script:
        "../script/CheckFileSizes.py"


def CheckTrimmedReadFile_Input(wildcards):
    if config['START'] == 'FASTQ':
        if config['PAIR_END']:
            return ["trimmed/{wildcards.sample}.R1.fastq.gz", "trimmed/{wildcards.sample}.R2.fastq.gz"]
        else:
            return "trimmed/{wildcards.sample}.fastq.gz"
    elif config['START'] == 'BAM':
        return "mapped_reads/{wildcards.sample}.bam"
    else:
        sys.exit("error")



rule CheckTrimmedReadFile:
    """check trimmed.fastq.gz for each sample, pass if non empty (>100bytes)"""
    input:
        CheckTrimmedReadFile_Input
    output:
        'fastqc/CheckFile/CheckFile.{sample}.txt'
    log:
        'fastqc/CheckFile/CheckFile.{sample}.log'
    script:
        "../script/CheckFileSizes.py"
