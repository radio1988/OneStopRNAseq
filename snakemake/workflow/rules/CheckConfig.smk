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
    input:
        CheckTrimmedFiles_Input(config,SAMPLES)
    output:
        'meta/CheckFiles.txt'
    script:
        "../script/CheckFileSizes.py"



rule CheckTrimmedFile:
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
