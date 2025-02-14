"""
Align (STAR/HISAT2), sort, qc, bigWig
"""

if config['ALIGNER'] == 'STAR' and config['START'] == 'FASTQ':  # only when FASTQ, to skip the align DAG, for BAM start
    rule STAR_Index:
        input:
            fa=config['GENOME'],
            gtf=config["GTF"],
        output:
            config['INDEX'] + "/SAindex"
        params:
            index=config['INDEX']
        conda:
            "../envs/star.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2400
        threads:
            16
        log:
            config['INDEX'] + "/SAindex.log"
        benchmark:
            "log/star_idx/star_idx.benchmark"
        shell:
            """
            whoami > {log};

            STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {params.index} \
            --genomeFastaFiles {input.fa} \
            --sjdbGTFfile {input.gtf} &>> {log}
            """


    rule STAR:
        input:
            index=config['INDEX'] + "/SAindex",
            gtf=config["GTF"],
            reads=["trimmed/{sample}.R1.fastq.gz", "trimmed/{sample}.R2.fastq.gz"] \
                if config['PAIR_END'] else \
                "trimmed/{sample}.fastq.gz"
        output:
            log="mapped_reads/{sample}.Log.final.out",
            bam=temp("mapped_reads/{sample}.bam"),
        params:
            index=config['INDEX']
        conda:
            "../envs/star.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 2400
        threads:
            16
        log:
            "mapped_reads/log/{sample}.bash.log"
        benchmark:
            "mapped_reads/log/{sample}.benchmark"
        shell:
            """
            STAR --runThreadN {threads} --genomeDir {params.index} --sjdbGTFfile {input.gtf} \
            --readFilesCommand zcat --readFilesIn {input.reads} \
            --outFileNamePrefix mapped_reads/{wildcards.sample}. \
            --outFilterType BySJout \
            --outMultimapperOrder Random \
            --outFilterMultimapNmax 200 \
            --alignSJoverhangMin 8 \
            --alignSJDBoverhangMin 3 \
            --outFilterMismatchNmax 999 \
            --outFilterMismatchNoverReadLmax 0.05 \
            --alignIntronMin 20 \
            --alignIntronMax 1000000 \
            --alignMatesGapMax 1000000 \
            --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
            --outSAMstrandField None \
            --outSAMtype BAM Unsorted \
            --quantMode GeneCounts \
            --outReadsUnmapped Fastx \
            > {log} 2>&1

            mv mapped_reads/{wildcards.sample}.Aligned.out.bam mapped_reads/{wildcards.sample}.bam &>> {log}
            pigz -p {threads} -f mapped_reads/{wildcards.sample}.Unmapped.out.mate*  &>> {log}
            """

    rule STAR_MultiQC:
        input:
            expand("mapped_reads/{sample}.Log.final.out",sample=SAMPLES)
        output:
            log=expand("bam_qc/STAR_Align_summary/{sample}.Log.final.out",sample=SAMPLES),
            report="bam_qc/STAR_Align_summary_multiqc_report.html",
        conda:
            "../envs/fastqc.yaml"
        resources:
            mem_mb=2000
        threads:
            1
        log:
            "bam_qc/STAR_Align_summary_multiqc_report.html.log",
        shell:
            """
            mkdir -p bam_qc/STAR_Align_summary/
            cp {input} bam_qc/STAR_Align_summary/
            multiqc {input} -f -n {output.report} &> {log}
            """

if config['ALIGNER'] == 'HISAT2':
    rule HISAT2_UNSTRANDED:
        input:
            genome=config['GENOME'],
            index=config['GENOME'] + '.hisat2_build.finished',
            reads=["trimmed/{sample}.R1.fastq.gz", "trimmed/{sample}.R2.fastq.gz"] \
                if config['PAIR_END'] else \
                "trimmed/{sample}.fastq.gz",
        output:
            bam=temp("mapped_reads/{sample}.bam")
        conda:
            "../envs/hisat2.yaml"
        params:
            strand=" ",# all set to unstranded for now(test)
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
            {params.pe} {params.strand} > {output.bam} 2> {log} 
            """


rule SAMtools_sort:
    input:
        mapped_bam_single_input(config)
    output:
        protected("sorted_reads/{sample}.bam")
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    threads:
        4
    log:
        "sorted_reads/log/{sample}.bam.log"
    benchmark:
        "sorted_reads/log/{sample}.bam.benchmark"
    shell:
        """
        rm -f {output}.tmp.*.bam;
        samtools sort -@ {threads} -m 3G {input} -o {output} &> {log}
        """

rule SAMtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 500
    threads:
        1
    log:
        "sorted_reads/log/{sample}.bam.bai.log"
    benchmark:
        "sorted_reads/log/{sample}.bam.bai.benchmark"
    shell:
        "samtools index -@ {threads} {input} &> {log}"


rule SAMtools_QC:
    input:
        bam="sorted_reads/{sample}.bam",
        bai="sorted_reads/{sample}.bam.bai"
    output:
        stats="bam_qc/stats/{sample}.stats.txt",
        idxstats="bam_qc/idxstats/{sample}.idxstats.txt",
        flagstat="bam_qc/flagstat/{sample}.flagstat.txt"
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    threads:
        2
    log:
        idxstats="bam_qc/stats/log/{sample}.stats.log",
        stats="bam_qc/idxstats/log/{sample}.idxstats.log",
        flagstat="bam_qc/flagstat/log/{sample}.flagstat.log"
    shell:
        """
        samtools idxstats {input.bam} > {output.idxstats} 2> {log.idxstats} &
        samtools stats {input.bam} > {output.stats} 2> {log.stats} &
        samtools flagstat {input.bam} > {output.flagstat} 2> {log.flagstat} &
        wait
        """

rule SAMtools_QC_MultiQc:
    input:
        stats=expand("bam_qc/stats/{sample}.stats.txt",sample=SAMPLES),
        idxstats=expand("bam_qc/idxstats/{sample}.idxstats.txt",sample=SAMPLES),
        flagstat=expand("bam_qc/flagstat/{sample}.flagstat.txt",sample=SAMPLES)
    output:
        stats="bam_qc/samtools_stats_multiqc_report.html",
        idxstats="bam_qc/samtools_idxstats_multiqc_report.html",
        flagstat="bam_qc/samtools_flagstat_multiqc_report.html",
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=2000
    threads:
        1
    log:
        "bam_qc/log/SAMtools_MultiQc.log"
    shell:
        """
        multiqc {input.stats} -f -n {output.stats} &> {log}
        multiqc {input.idxstats} -f -n {output.idxstats} &>> {log}
        multiqc {input.flagstat} -f -n {output.flagstat} &>> {log}
        """


rule QoRTs:
    input:
        bam="sorted_reads/{sample}.bam",
        gtf=config['GTF'],
        length="meta/read_length.max.txt",
        check=checkFileInput
    output:
        temp(directory("bam_qc/QoRTs/{sample}"))
    params:
        pe=' ' if config['PAIR_END'] else '--singleEnded',
        aligner="--minMAPQ 60" if config['ALIGNER'] == 'HISAT2' else "",# only support STAR and HISAT2
        length=get_read_length("meta/read_length.max.txt")
    conda:
        "../envs/java11.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16500
    threads:
        1
    log:
        "bam_qc/QoRTs/log/{sample}.log"
    benchmark:
        "bam_qc/QoRTs/log/{sample}.benchmark"
    shell:
        """
        java -Xmx16G -jar workflow/envs/hartleys-QoRTs-099881f/QoRTs.jar QC {params.pe} \
        {params.aligner} \
        --maxReadLength {params.length} \
        --maxPhredScore 45 \
        {input.bam} {input.gtf} bam_qc/QoRTs/{wildcards.sample}/ &> {log}
        """


rule QoRTs_MultiPlot:
    input:
        expand("bam_qc/QoRTs/{sample}",sample=SAMPLES)
    output:
        "bam_qc/QoRTs_MultiPlot/plot-basic.pdf"
    conda:
        "../envs/qorts.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000,
    threads:
        2
    log:
        "bam_qc/QoRTs_MultiPlot/log/QoRTs_MultiPlot.log"
    benchmark:
        "bam_qc/QoRTs_MultiPlot/log/QoRTs_MultiPlot.benchmark"
    shell:
        """
        python workflow/script/meta_to_decoder.py {config[META]} > {log} 2>&1;
        Rscript workflow/script/QoRT.R > {log} 2>&1;
        mkdir -p bam_qc/QoRTs_MultiPlot/details/;
        mv bam_qc/QoRTs_MultiPlot/plot-sample* bam_qc/QoRTs_MultiPlot/details/;
        """

rule compress_bam_qc
    input:
        "bam_qc/QoRTs_MultiPlot/plot-basic.pdf",
        "bam_qc/STAR_Align_summary_multiqc_report.html" if config['ALIGNER'] == 'STAR' and config[
            'START'] == 'FASTQ' else "bam_qc/samtools_stats_multiqc_report.html",# todo hisat2
        "bam_qc/samtools_stats_multiqc_report.html",
        "bam_qc/samtools_idxstats_multiqc_report.html",
        "bam_qc/samtools_flagstat_multiqc_report.html"
    output:
        "bam_qc/bam_qc.zip"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2000
    threads:
        1
    shell:
        """
        D=bam_qc 
        rm -f $D/$D.zip && [ -d $D ] && zip -rq  $D/$D.zip $D/
        """

rule bamCoverage:
    input:
        bam="sorted_reads/{sample}.bam",
        bai="sorted_reads/{sample}.bam.bai",
        check=checkFileInput
    output:
        "bigWig/{sample}.{mode}.cpm.bw"
    params:
        "--minMappingQuality 20" if config['MODE'] == "strict" else " "
    conda:
        "../envs/deeptools.yaml"
    threads:
        4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,# 16 GB for tricky cases
    log:
        "bigWig/log/{sample}.{mode}.cpm.bw.log"
    benchmark:
        "bigWig/log/{sample}.{mode}.cpm.bw.benchmark"
    shell:
        """
        bamCoverage --bam {input.bam} -o  {output} --numberOfProcessors {threads} \
        --outFileFormat bigwig --normalizeUsing CPM --binSize 50 \
        {params} > {log} 2>&1
        """
