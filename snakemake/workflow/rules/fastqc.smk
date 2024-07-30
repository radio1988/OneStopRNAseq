rule FastQC1:
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


rule FastQC_MultiQC1:
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

rule FastQC2:
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


rule FastQC_MultiQC2:
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

rule fastqc_compress:
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
