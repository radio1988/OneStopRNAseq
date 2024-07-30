rule featureCounts_EXON:
    '''
    exon read count
    recognizes config.yaml: 
        - MODE: strict, liberal
        - config['PAIR_END']: True, False
    '''
    input:
        bams=mapped_bam_inputs(config,SAMPLES),
        gtf=config["GTF"]
    output:
        COUNT="feature_count/counts.s{strand}.{mode}.txt",
        summary="feature_count/counts.s{strand}.{mode}.txt.summary"
    params:
        pe='-p' if config['PAIR_END'] else ' ',
        mode='-Q 20 --minOverlap 1 --fracOverlap 0 -B -C' \
            if config['MODE'] == 'strict' else \
            '-M -Q 0 --primary --minOverlap 1 --fracOverlap 0'
    conda:
        "../envs/subread.yaml"
    priority:
        100
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
    threads:
        4
    log:
        "feature_count/log/counts.s{strand}.{mode}.txt.log"
    benchmark:
        "feature_count/log/counts.s{strand}.{mode}.txt.benchmark"
    shell:
        """
        featureCounts -a {input.gtf} -o {output.COUNT} \
        -T {threads} -g gene_id -s {wildcards.strand} \
        {params.pe}  {params.mode} \
        {input.bams} > {log} 2>&1
        """

rule featureCounts_EXON_multiqc:
    input:
        "feature_count/counts.s{strand}.{mode}.txt.summary"
    output:
        "feature_count/counts.s{strand}.{mode}.txt.summary.multiqc_report.html"
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        "feature_count/log/counts.s{strand}.{mode}.txt.summary.multiqc_report.html.log"
    benchmark:
        "feature_count/log/counts.s{strand}.{mode}.txt.summary.multiqc_report.html.benchmark"
    shell:
        """
        multiqc {input} -f -n {output} > {log} 2>&1;
        """


rule featureCounts_GENE:
    '''
    gene read count, including exon and intron reads
    the larger one of gene-count/exon-count for each gene, is used as gene-count
    recognizes config.yaml: 
        - config['MODE']: strict, liberal
        - config['PAIR_END']: True, False
    '''
    input:
        bams=mapped_bam_inputs(config,SAMPLES),
        exon_counts="feature_count/counts.s{strand}.{mode}.txt",
        gtf=config["GTF"]
    output:
        ct='feature_count_gene_level/counts.s{strand}.{mode}.txt',
        summary='feature_count_gene_level/counts.s{strand}.{mode}.txt.summary'
    params:
        pe='-p -B -C ' if config['PAIR_END'] else ' ',
        mode='-Q 20 ' if config['MODE'] == 'strict' else '-M --primary -Q 0'
    conda:
        "../envs/subread.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
    threads:
        4
    log:
        'feature_count_gene_level/log/counts.s{strand}.{mode}.log'
    benchmark:
        'feature_count_gene_level/log/counts.s{strand}.{mode}.benchmark'
    shell:
        """
        # gene level count
        featureCounts -a {input.gtf} -o {output.ct} -T {threads} \
        -t gene -g gene_id -s {wildcards.strand} \
        --minOverlap 1 --fracOverlap 0 \
        {params.pe} {params.mode} {input.bams} \
        > {log} 2>&1

        # get larger count (exon vs gene) for each gene
        python workflow/script/get_max_of_featureCountsTable.py \
        {input.exon_counts} {output.ct}  \
        >> {log} 2>&1
        """


rule featureCounts_GENE_multiqc:
    input:
        'feature_count_gene_level/counts.s{strand}.{mode}.txt.summary'
    output:
        "feature_count_gene_level/counts.s{strand}.{mode}.txt.summary.multiqc_report.html"
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        "feature_count_gene_level/log/counts.s{strand}.{mode}.txt.summary.multiqc_report.html.log"
    benchmark:
        "feature_count_gene_level/log/counts.s{strand}.{mode}.txt.summary.multiqc_report.html.benchmark"
    shell:
        """
        multiqc {input} -f -n {output} > {log} 2>&1;
        """


rule Strand_Detection:
    input:
        strand_detection_input(config)
    output:
        "feature_count/counts.strict.txt" if config['MODE'] == 'strict' else "feature_count/counts.liberal.txt",
        "meta/strandness.detected.txt"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    log:
        "meta/log/strandness.detected.txt.log"
    script:
        "../script/strandness_detection.py"


rule Strand_Detection_GeneLevel:
    input:
        ["feature_count_gene_level/counts.s{}.{}.txt.summary". \
             format(strand,config['MODE']) for strand in config['STRAND']]
    output:
        "feature_count_gene_level/counts.{}.txt".format(config['MODE']),
        "feature_count_gene_level/strandness.detected.txt"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    log:
        "feature_count_gene_level/log/strandness.detected.txt.log"
    script:
        "../script/strandness_detection.ge.py"


# rule SalmonTE_prep
if config['PAIR_END']:
    rule SalmonTE_Prep:
        input:
            r1=expand("trimmed/{sample}.R1.fastq.gz",sample=SAMPLES),
            r2=expand("trimmed/{sample}.R2.fastq.gz",sample=SAMPLES)
        output:
            r1=temp(expand("fastq_salmon/{sample}_1.fastq.gz",sample=SAMPLES)),
            r2=temp(expand("fastq_salmon/{sample}_2.fastq.gz",sample=SAMPLES))
        log:
            "fastq_salmon/SalmonTE_prep.log"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1000
        threads:
            1
        run:
            for i, f in enumerate(input):
                os.system("ln -s " + "../" + input[i] + " " + output[i])

rule SalmonTE:
    # only start from Fastq
    input:
        reads=expand("fastq_salmon/{sample}_{n}.fastq.gz",n=[1, 2],sample=SAMPLES) \
            if config['PAIR_END'] else \
            expand("trimmed/{sample}.fastq.gz",sample=SAMPLES),
        raw_reads=expand("trimmed/{sample}.R1.fastq.gz",sample=SAMPLES) \
            if config['PAIR_END'] else \
            expand("trimmed/{sample}.fastq.gz",sample=SAMPLES),
        raw_reads2=expand("trimmed/{sample}.R2.fastq.gz",sample=SAMPLES) \
            if config['PAIR_END'] else \
            expand("trimmed/{sample}.fastq.gz",sample=SAMPLES)
    output:
        "SalmonTE_output/EXPR.csv"
    conda:
        "../envs/salmonte.yaml"  # test
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    params:
        ref=config['TE_REFERENCE'],
    threads:
        16
    log:
        "log/SalmonTE.log"
    benchmark:
        "log/SalmonTE.benchmark"
    shell:
        """
        rm -rf SalmonTE_output/
        python workflow/envs/SalmonTE/SalmonTE.py --version >> {log}
        python workflow/envs/SalmonTE/SalmonTE.py quant \
        --reference={params.ref} --exprtype=count \
        --num_threads={threads} \
        {input.reads} > {log} 2>&1
        """

rule Merge_TE_and_GE:
    input:
        gene=Merge_TE_and_Gene_input(config),
        te="SalmonTE_output/EXPR.csv"
    output:
        "feature_count_gene_level/TE_included.txt" \
            if config['INTRON'] \
            else "feature_count/TE_included.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000,
    threads:
        1
    log:
        "feature_count_gene_level/log/TE_included.txt.log" \
            if config['INTRON'] \
            else "feature_count/log/TE_included.txt.log"
    benchmark:
        "feature_count_gene_level/log/TE_included.txt.benchmark" \
            if config['INTRON'] \
            else "feature_count/log/TE_included.txt.benchmark"
    shell:
        """
        python workflow/script/merge_featureCount_and_SalmonTE.py \
        {input.gene} {input.te} {output} > {log} 2>&1;
        """
