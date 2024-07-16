

# CHECK CONFIG.YAML
if config["MODE"] == 'liberal' and config["CleanUpRNAseqCorrection"]:
    pass


if config["INTRON"] and config["CleanUpRNAseqCorrection"]:
    pass


def CheckTrimmedFiles_Input(config, SAMPLES):
    L = []
    if config['STRAND'] == 0:
        L = ["trimmed/{}.R1.fastq.gz".format(sample) for sample in SAMPLES]
    else:
        L = ["trimmed/{}.R1.fastq.gz".format(sample) for sample in SAMPLES]
        L.extend(["trimmed/{}.R2.fastq.gz".format(sample) for sample in SAMPLES])

# CHECK FILES
rule CheckTrimmedFiles:
    input:
        CheckTrimmedFiles_Input(config, SAMPLES)
    output:
        'meta/CheckFiles.txt'
    script:
        "../script/CheckFiles.py"
