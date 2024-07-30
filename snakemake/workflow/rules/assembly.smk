"""
HISAT2 -> StringTie
"""


rule HISAT2_BUILD:
    """
    guided by GTF

    with only genome.fa, only takes ~6GB RAM, 20min
    $ hisat2-build -p 16 genome.fa genome

    with transcripts, takes 160GB RAM, 1h (human t2t, 2023)
    $ hisat2-build -p 8 --exon genome.exon --ss genome.ss genome.fa genome_tran
    http://daehwankimlab.github.io/hisat2/howto/

    """
    input:
        genome=config['GENOME'],
        gtf=config['GTF']
    output:
        ss=config['GTF'] + "ss",
        exon=config['GTF'] + "exon",
        index=config['GENOME'] + ".1.ht2",
        flag=touch(config['GENOME'] + '.hisat2_build.finished')
    conda:
        "../envs/hisat2.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 10000
    threads:
        16
    log:
        config['GENOME'] + '.hisat2_build.log'
    benchmark:
        config['GENOME'] + '.hisat2_build.benchmark'
    shell:
        """
        hisat2_extract_splice_sites.py -v {input.gtf} > {output.ss} 2> {log}
        hisat2_extract_exons.py -v {input.gtf} > {output.exon} 2>> {log}
        hisat2-build -p {threads} \
        {input.genome} {input.genome} --ss {output.ss} --exon {output.exon} 2>> {log}
        """

rule HISAT2:
    input:
        genome=config['GENOME'],
        index=config['GENOME'] + '.hisat2_build.finished',
        reads=["trimmed/{sample}.R1.fastq.gz", "trimmed/{sample}.R2.fastq.gz"] \
            if config['PAIR_END'] else \
            "trimmed/{sample}.fastq.gz",
        strandness="meta/strandness.detected.txt"
    output:
        bam="hisat2/{sample}.bam",
        bai="hisat2/{sample}.bam.bai"
    conda:
        "../envs/hisat2.yaml"
    params:
        strand=get_strandness_for_hisat2_PE("meta/strandness.detected.txt"),# todo: SE?
        pe="-1 trimmed/{sample}.R1.fastq.gz -2 trimmed/{sample}.R2.fastq.gz" \
            if config['PAIR_END'] else \
            "-U trimmed/{sample}.fastq.gz"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        16
    log:
        "hisat2/log/{sample}.log"
    benchmark:
        "hisat2/log/{sample}.benchmark"
    shell:
        """
        hisat2 -x {input.genome} -p {threads} --dta-cufflinks \
        {params.pe} {params.strand} 2> {log} | \
        samtools sort -@ 2 -o {output.bam} 2>> {log}

        samtools index {output.bam}
        """


rule StringTie:  # todo: group level
    input:
        bam="hisat2/{sample}.bam",
        gtf=config['GTF'],
        strandness='meta/strandness.detected.txt'
    output:
        "stringtie/{sample}.stringtie.gtf"
    conda:
        "../envs/stringtie.yaml"
    params:
        strand=get_strandness_for_stringtie('meta/strandness.detected.txt')
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        16
    log:
        "stringtie/log/{sample}.stringtie.log"
    benchmark:
        "stringtie/log/{sample}.stringtie.benchmark"
    shell:
        """
        stringtie -G {input.gtf} {params.strand} -o {output} {input.bam} 2> {log}
        """

rule StringTie_Merge:
    input:
        stringties=expand("stringtie/{sample}.stringtie.gtf",sample=SAMPLES),
        gtf=config['GTF']
    output:
        "stringtie/stringtie.merged.gtf"
    conda:
        "../envs/stringtie.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        16
    log:
        "stringtie/log/stringtie.merged.gtf.log"
    benchmark:
        "stringtie/log/stringtie.merged.gtf.benchmark"
    shell:
        """
        stringtie --merge -G {input.gtf} -o {output} {input.stringties} 2> {log}
        """
