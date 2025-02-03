"""
DEXSeq: slow, not well annotated
"""

DEXSeq_GFF = config['GTF'] + ".dexseq.gff"

rule DEXSeq_GFF_Prep:
    input:
        config['GTF']
    output:
        DEXSeq_GFF
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    threads:
        1
    log:
        DEXSeq_GFF + ".log"
    benchmark:
        DEXSeq_GFF + ".benchmark"
    conda:
        "../envs/bioconda_misc.yaml"
    shell:
        "python workflow/script/dexseq_prepare_annotation.py -r no {input} {output} &> {log}"


rule DEXSeq_Count:
    """
    Count reads for DEXSeq
    """
    input:
        bam="sorted_reads/{sample}.bam", # todo: maybe -r name (sort by name) is faster, low priority
        strandFile="meta/strandness.detected.txt",
        gff=DEXSeq_GFF
    output:
        "DEXSeq_count/{sample}_count.txt"
    conda:
        "../envs/bioconda_misc.yaml"
    params:
        strand=get_strandness_for_dexseq('meta/strandness.detected.txt'),
        readType=('-p yes -r pos' if config['PAIR_END'] else ' '),
        MAPQ=10
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    log:
        "DEXSeq_count/log/{sample}_count.txt.log"
    benchmark:
        "DEXSeq_count/log/{sample}_count.txt.benchmark"
    shell:
        "python workflow/script/dexseq_count.py -f bam -a {params.MAPQ} \
        {params.readType} {params.strand} {DEXSeq_GFF} {input.bam} {output} &> {log}"


if config['DEXSEQ_ANALYSIS']:
    rule DEXSeq:
        input:
            count_files=expand("DEXSeq_count/{sample}_count.txt",sample=SAMPLES),
            meta=config['META'],
            contrast=config['CONTRAST_AS'],
            gffFile=config['GTF'] + ".dexseq.gff",
            strandness="meta/strandness.detected.txt"
        output:
            "DEXSeq/contrast{ascn}/contrast{ascn}.RData"
        conda:
            "../envs/dexseq.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 3000
        threads:
            12
        params:
            rmd="'workflow/script/dexseq.r'",
            annoFile=config['ANNO_TAB'],
            max_fdr=config['MAX_FDR'],
            min_lfc=config['MIN_LFC'],
            min_gene_count=config['MIN_GENE_COUNT'] if 'MIN_GENE_COUNT' in config else 100
        log:
            "DEXSeq/contrast{ascn}/DEXSeq.log"
        benchmark:
            "DEXSeq/contrast{ascn}/DEXSeq.benchmark"
        shell:
            """
            mkdir -p DEXSeq 
            cp workflow/script/dexseq.r DEXSeq/ &> {log}
            Rscript workflow/script/dexseq.r {input.meta} {input.contrast} {input.gffFile} {params.annoFile} \
            {params.max_fdr} {params.min_lfc} {threads} {MIN_GENE_COUNT} {wildcards.ascn} &>> {log}
            D=DEXSeq/contrast{wildcards.ascn}
            rm -f $D.zip && [ -d $D ] && zip -rq $D.zip $D/ &>> {log}
            """
