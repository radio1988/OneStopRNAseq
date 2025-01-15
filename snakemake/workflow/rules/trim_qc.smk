import sys

if config['PAIR_END']:
    rule fastp_pe:
        """
        Trim paired-end reads using fastp
        default min-len 15
        -5, -3, window 4 , min 20Q, max 5N
        adaptors will be trimmed
        output is Phred33
        """
        input:
            r1="fastq/{sample,[A-Za-z0-9_-]+}.R1.fastq.gz",
            r2="fastq/{sample,[A-Za-z0-9_-]+}.R2.fastq.gz",
            r1_r2_check="fastqc/details_raw/{sample}.r1r2_checked"  # if r1 != r2, abort
        output:
            r1=temp("trimmed/{sample}.R1.fastq.gz"),
            r2=temp("trimmed/{sample}.R2.fastq.gz"),
            r1_unpaired=temp("trimmed/unpaired/{sample}.R1.fastq.gz"),
            r2_unpaired=temp("trimmed/unpaired/{sample}.R2.fastq.gz"),
            html="trimmed/report/{sample}.html",
            json="trimmed/report/{sample}.json"
        conda:
            "../envs/fastp.yaml"
        params:
            extra="--adapter_fasta " + config['ADAPTORS'] + " --detect_adapter_for_pe",
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1000
        threads:
            4
        log:
            "trimmed/log/{sample}.trim.log"
        benchmark:
            "trimmed/log/{sample}.trim.log.benchmark"
        shell:
            "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            --unpaired1 {output.r1_unpaired} --unpaired2 {output.r2_unpaired} \
            -j {output.json} -h {output.html} -R {wildcards.sample} \
            -w {threads} {params.extra} &> {log}"
    # rule Trimmomatic_PE:
    #         input:
    #             r1="fastq/{sample,[A-Za-z0-9_-]+}.R1.fastq.gz",
    #             r2="fastq/{sample,[A-Za-z0-9_-]+}.R2.fastq.gz"
    #         output:
    #             r1=temp("trimmed/{sample}.R1.fastq.gz"),
    #             r2=temp("trimmed/{sample}.R2.fastq.gz"),
    #             r1_unpaired=temp("trimmed/unpaired/{sample}.R1.fastq.gz"),
    #             r2_unpaired=temp("trimmed/unpaired/{sample}.R2.fastq.gz"),
    #         params:
    #             trimmer=["ILLUMINACLIP:" + config[
    #                 'ADAPTORS'] + ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16 TOPHRED33"],
    #         resources:
    #             mem_mb=lambda wildcards, attempt: attempt * 1000
    #         threads:
    #             4
    #         log:
    #             "trimmed/log/{sample}.trim.log"
    #         benchmark:
    #             "trimmed/log/{sample}.trim.log.benchmark"
    #         wrapper:
    #             "v3.13.1/bio/trimmomatic/pe"
else:
    rule fastp_se:
        input:
            sample=["fastq/{sample,[A-Za-z0-9_-]+}.fastq.gz"]
        output:
            trimmed=temp("trimmed/{sample}.fastq.gz"),
            failed=temp("trimmed/failed/{sample}.fastq.gz"),
            html="trimmed/{sample}.html",
            json="trimmed/{sample}.json"
        log:
            "trimmed/log/{sample}.trim.log"
        benchmark:
            "trimmed/log/{sample}.trim.benchmark"
        params:
            extra="--adapter_fasta " + config['ADAPTORS'],
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2000
        threads:
            4
        conda:
            "../envs/fastp.yaml"
        shell:
            "fastp -i {input.sample} -o {output.trimmed} \
            --failed_out {output.failed} \
            -j {output.json} -h {output.html} -R {wildcards.sample} \
            -w {threads} {params.extra} &> {log}"
    # rule Trimmomatic_SE:
    #     input:
    #         "fastq/{sample,[A-Za-z0-9_-]+}.fastq.gz",
    #     output:
    #         temp("trimmed/{sample}.fastq.gz")
    #     params:
    #         trimmer=["ILLUMINACLIP:" + config[
    #             'ADAPTORS'] + ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16 TOPHRED33"],
    #     resources:
    #         mem_mb=lambda wildcards, attempt: attempt * 1000
    #     threads:
    #         4
    #     log:
    #         "trimmed/log/{sample}.trim.log"
    #     benchmark:
    #         "trimmed/log/{sample}.trim.log.benchmark"
    #     wrapper:
    #         "v3.13.1/bio/trimmomatic/se"


rule FastQC_Raw:
    input:
        ["fastq/{sample}.R1.fastq.gz", "fastq/{sample}.R2.fastq.gz"] if config['PAIR_END'] \
            else "fastq/{sample}.fastq.gz"
    output:
        ["fastqc/details_raw/{sample}.R1_fastqc.html", "fastqc/details_raw/{sample}.R2_fastqc.html"] \
            if config['PAIR_END'] else \
            "fastqc/details_raw/{sample}_fastqc.html",
        ["fastqc/details_raw/{sample}.R1_fastqc.zip", "fastqc/details_raw/{sample}.R2_fastqc.zip"] \
            if config['PAIR_END'] else \
            "fastqc/details_raw/{sample}_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        4
    log:
        "fastqc/details_raw/log/{sample}.log"
    benchmark:
        "fastqc/details_raw/log/{sample}.benchmark"
    shell:
        "fastqc -t {threads} {input} -o fastqc/details_raw &> {log};"


if config['PAIR_END']:
    rule Check_R1_R2:
        input:
            r1="fastqc/details_raw/{sample}.R1_fastqc.zip",
            r2="fastqc/details_raw/{sample}.R2_fastqc.zip"
        output:
            "fastqc/details_raw/{sample}.r1r2_checked",
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1000
        threads:
            1
        log:
            "fastqc/details_raw/log/{sample}.check.log"
        shell:
            "python ../script/check_r1_r2.py {input.r1} {input.r2} &> {log} && touch {output}"


rule MultiQC_Raw:
    input:
        [expand("fastqc/details_raw/{sample}.R1_fastqc.zip",sample=SAMPLES),
         expand("fastqc/details_raw/{sample}.R2_fastqc.zip",sample=SAMPLES)] \
            if config['PAIR_END'] else \
            expand("fastqc/details_raw/{sample}_fastqc.zip",sample=SAMPLES)
    output:
        "fastqc/raw_reads_multiqc_report.html"
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    threads:
        1
    log:
        "fastqc/raw_reads_multiqc_report.log"
    benchmark:
        "fastqc/raw_reads_multiqc_report.benchmark"
    shell:
        "multiqc {input} -f -n {output} &> {log};"

rule FastQC_Trimmed:
    input:
        ["trimmed/{sample}.R1.fastq.gz", "trimmed/{sample}.R2.fastq.gz"] \
            if config['PAIR_END'] else \
            "trimmed/{sample}.fastq.gz"
    output:
        ["fastqc/details_trimmed/{sample}.R1_fastqc.html",
         "fastqc/details_trimmed/{sample}.R2_fastqc.html"] \
            if config['PAIR_END'] else \
            "fastqc/details_trimmed/{sample}_fastqc.html",
        ["fastqc/details_trimmed/{sample}.R1_fastqc.zip",
         "fastqc/details_trimmed/{sample}.R2_fastqc.zip"] \
            if config['PAIR_END'] else \
            "fastqc/details_trimmed/{sample}_fastqc.zip"
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        4
    log:
        "fastqc/details_trimmed/log/{sample}.log"
    benchmark:
        "fastqc/details_trimmed/log/{sample}.benchmark"
    shell:
        "fastqc -t {threads} {input} -o fastqc/details_trimmed &> {log};"


rule MultiQC_Trimmed:
    input:
        [expand("fastqc/details_trimmed/{sample}.R1_fastqc.zip",sample=SAMPLES),
         expand("fastqc/details_trimmed/{sample}.R2_fastqc.zip",sample=SAMPLES)] \
            if config['PAIR_END'] else \
            expand("fastqc/details_trimmed/{sample}_fastqc.zip",sample=SAMPLES)
    output:
        "fastqc/trimmed_reads_multiqc_report.html"
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        "fastqc/trimmed_reads_multiqc_report.log"
    benchmark:
        "fastqc/trimmed_reads_multiqc_report.benchmark"
    shell:
        """
        multiqc {input} -f -n {output} &> {log};
        """

rule FastQC_Zip:
    input:
        "fastqc/trimmed_reads_multiqc_report.html",
        "fastqc/raw_reads_multiqc_report.html"
    output:
        "fastqc/fastqc.zip"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    threads:
        1
    shell:
        "rm -f fastqc/fastqc.zip && [ -d fastqc ] && zip -rq  fastqc/fastqc.zip fastqc/"


def CheckTrimmedFiles_Input(config, SAMPLES):
    L = []
    if config['PAIR_END']:
        L = ["trimmed/{}.R1.fastq.gz".format(sample) for sample in SAMPLES]
        L.extend(["trimmed/{}.R2.fastq.gz".format(sample) for sample in SAMPLES])
    else:
        L = ["trimmed/{}.fastq.gz".format(sample) for sample in SAMPLES]
    return L


# CHECK FILES
rule CheckTrimmedReadFiles:
    """check trimmed.fastq.gz for all samples, pass if all non empty (>100bytes)"""
    input:
        CheckTrimmedFiles_Input(config, SAMPLES)
    output:
        'fastqc/CheckFile/all.txt'
    log:
        "fastqc/CheckFile/all.log"
    benchmark:
        "fastqc/CheckFile/all.benchmark"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    threads:
        1
    script:
        "../script/CheckFileSizes.py"


def CheckTrimmedReadFile_Input(wildcards):
    if config['START'] == 'FASTQ':
        if config['PAIR_END']:
            return ["trimmed/{sample}.R1.fastq.gz", "trimmed/{sample}.R2.fastq.gz"]
        else:
            return "trimmed/{sample}.fastq.gz"
    elif config['START'] == 'BAM':
        return "mapped_reads/{sample}.bam"
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
