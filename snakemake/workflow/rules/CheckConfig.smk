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
rule CheckTrimmedFiles:
    """check trimmed.fastq.gz for all samples, pass if all non empty (>100bytes)"""
    input:
        CheckTrimmedFiles_Input(config, SAMPLES)
    output:
        'meta/CheckFile/all.txt'
    log:
        "meta/CheckFile/all.log"
    script:
        "../script/CheckFileSizes.py"



rule CheckTrimmedFile:
    """check trimmed.fastq.gz for each sample, pass if non empty (>100bytes)"""
    input:
        ["trimmed/{sample}.R1.fastq.gz", "trimmed/{sample}.R2.fastq.gz"] \
            if config['PAIR_END'] \
            else "trimmed/{sample}.fastq.gz"
    output:
        'meta/CheckFile/CheckFile.{sample}.txt'
    log:
        'meta/CheckFile/CheckFile.{sample}.log'
    script:
        "../script/CheckFileSizes.py"
