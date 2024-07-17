if config['PAIR_END']:
    rule Trimmomatic_PE:
        input:
            r1="fastq/{sample,[A-Za-z0-9_-]+}.R1.fastq.gz",
            r2="fastq/{sample,[A-Za-z0-9_-]+}.R2.fastq.gz"
        output:
            r1=temp("trimmed/{sample}.R1.fastq.gz"),
            r2=temp("trimmed/{sample}.R2.fastq.gz"),
            r1_unpaired=temp("trimmed/unpaired/{sample}.R1.fastq.gz"),
            r2_unpaired=temp("trimmed/unpaired/{sample}.R2.fastq.gz"),
        params:
            trimmer=["ILLUMINACLIP:" + config[
                'ADAPTORS'] + ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16 TOPHRED33"],
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1000
        threads:
            4
        log:
            "trimmed/log/{sample}.trim.log"
        benchmark:
            "trimmed/log/{sample}.trim.log.benchmark"
        wrapper:
            "v3.13.1/bio/trimmomatic/pe"
else:
    rule Trimmomatic_SE:
        input:
            "fastq/{sample,[A-Za-z0-9_-]+}.fastq.gz",
        output:
            temp("trimmed/{sample}.fastq.gz")
        params:
            trimmer=["ILLUMINACLIP:" + config[
                'ADAPTORS'] + ":2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:16 TOPHRED33"],
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1000
        threads:
            4
        log:
            "trimmed/log/{sample}.trim.log"
        benchmark:
            "trimmed/log/{sample}.trim.log.benchmark"
        wrapper:
            "v3.13.1/bio/trimmomatic/se"


rule FastQC1:
    input:
        ["fastq/{sample}.R1.fastq.gz", "fastq/{sample}.R2.fastq.gz"] if config['PAIR_END'] else \
            "fastq/{sample}.fastq.gz"
    output:
        ["fastqc/details_raw/{sample}.R1_fastqc.html", "fastqc/details_raw/{sample}.R2_fastqc.html"] if config['PAIR_END'] else \
            "fastqc/details_raw/{sample}_fastqc.html",
        ["fastqc/details_raw/{sample}.R1_fastqc.zip", "fastqc/details_raw/{sample}.R2_fastqc.zip"] if config['PAIR_END'] else \
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
        ["trimmed/{sample}.R1.fastq.gz", "trimmed/{sample}.R2.fastq.gz"] if config['PAIR_END'] else \
            "trimmed/{sample}.fastq.gz"
    output:
        ["fastqc/details_trimmed/{sample}.R1_fastqc.html",
         "fastqc/details_trimmed/{sample}.R2_fastqc.html"] if config['PAIR_END'] else \
            "fastqc/details_trimmed/{sample}_fastqc.html",
        ["fastqc/details_trimmed/{sample}.R1_fastqc.zip",
         "fastqc/details_trimmed/{sample}.R2_fastqc.zip"] if config['PAIR_END'] else \
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


if config['ALIGNER'] == 'STAR':
    rule STAR_Index:
        input:
            fa=config['GENOME'],
            gtf=config["GTF"],
        output:
            config['INDEX'] + "/SAindex"
        params:
            index = config['INDEX']
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
            index = config['INDEX']
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
elif config['ALIGNER'] == 'HISAT2':
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
else:
    sys.exit("config['ALIGNER'] option not recognized")


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
        check="fastqc/CheckFile/CheckFile.{sample}.txt"
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

rule zip_bam_qc:
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
        check="fastqc/CheckFile/CheckFile.{sample}.txt"
    output:
        "bigWig/{sample}.{mode}.cpm.bw"
    params:
        "--minMappingQuality 20" if MODE == "strict" else " "
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
            if MODE == 'strict' else \
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
        "script/strandness_detection.py"


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
        "script/strandness_detection.ge.py"


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

if config['DESEQ2_ANALYSIS']:
    rule DESeq2:
        input:
            cnt=DESeq2_input(config),
            meta=config['META'],
            contrast=config["CONTRAST_DE"],
        output:
            "DESeq2/DESeq2.html",
            expand("DESeq2/rnk/{contrast}.rnk",contrast=CONTRASTS_DE)
        conda:
            "../envs/deseq2.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4000,
        params:
            rmd="'./DESeq2/DESeq2.Rmd'",
            fdr=config['MAX_FDR'],
            lfc=config['MIN_LFC'],
            independentFilter=config["independentFilter"],
            cooksCutoff=config["cooksCutoff"],
            blackSamples=config['blackSamples'] if 'blackSamples' in config else "",
            anno_tab=config['ANNO_TAB'],
            o="'DESeq2.html'"
        priority:
            100
        threads:
            1
        log:
            "DESeq2/DESeq2.log"
        benchmark:
            "DESeq2/DESeq2.benchmark"
        shell:
            #Rscript -e rmarkdown::render({params.rmd})
            'cp workflow/script/DESeq2.Rmd DESeq2; '
            'Rscript -e "rmarkdown::render( \
            {params.rmd}, \
            params=list( \
                max_fdr={params.fdr}, \
                min_lfc={params.lfc}, \
                cookscutoff={params.cooksCutoff}, \
                indfilter={params.independentFilter}, \
                countFile=\'../{input.cnt}\', \
                annoFile=\'{params.anno_tab}\', \
                metaFile=\'../{input.meta}\', \
                contrastFile=\'../{input.contrast}\', \
                blackSamples=\'{params.blackSamples}\' \
                ), \
            output_file={params.o})" > {log} 2>&1 ;'
            'D=DESeq2; rm -f $D/$D.zip && [ -d $D ] && zip -rq  $D/$D.zip $D/ >> {log} 2>&1;'


rule GSEA:
    """
    config[START] == RNK, then find config[RNKS], else, find DESeq2/rnk or CleanUpRNAseqDE/rnk
    config[RNKS]  must be put in ./meta/ folder, and has suffix of txt or xlsx
    """
    input:
        rnk=lambda wildcards: input_rnk_fname1(wildcards,config),
        db=lambda wildcards: os.path.join(config['GSEA_DB_PATH'],wildcards["db"])
    output:
        "gsea/{fname}/{db}.GseaPreranked/index.html"
    conda:
        "../envs/java11.yaml"  # test
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8100
    params:
        svg=config["GSEA_PLOT_SVG"],
        nplot=config["GSEA_NPLOTS"],
        label_db=lambda wildcards: wildcards["db"][:-12],
        rnk_flat_file=lambda wildcards: input_rnk_fname2(wildcards,config),
        gmt_fmted=lambda wildcards: "gsea/gene_sets/" + wildcards['fname'] + "." + wildcards['db']
    threads:
        1
    log:
        "gsea/log/{fname}/{db}.log"
    benchmark:
        "gsea/log/{fname}/{db}.benchmark"
    shell:
        """
        echo rnkfile: {input.rnk} > {log}
        echo rnkfile.cleaned: {params.rnk_flat_file} >> {log}
        echo label_db: {params.label_db} >> {log}
        echo fname: {wildcards.fname} >> {log}
        echo dbname: {wildcards.db} >> {log}

        mkdir -p gsea gsea/gene_sets gsea/{wildcards.fname} >> {log} 2>&1
        rm -rf gsea/{wildcards.fname}/{wildcards.db}.GseaPreranked.*/  # avoid confusion
        rm -rf gsea/{wildcards.fname}/{wildcards.db}.GseaPreranked/  # avoid dir structure mistake
        rm -rf gsea/{wildcards.fname}/error_{wildcards.db}*GseaPreranked*/  # avoid confusion of temp files in multiple run attempts

        python workflow/script/rnk_to_upper.py {input.rnk} >> {log} 2>&1;  # fix gene symbol error; standardize file format (xlsx, rnk, etc)
        python workflow/script/gmt_to_upper.py -f {input.db} 1> {params.gmt_fmted} 2> {log}  # gene symbols to upper
        
        if workflow/envs/GSEA_4.3.2/gsea-cli.sh GSEAPreranked \
        -gmx {params.gmt_fmted} -rnk {params.rnk_flat_file} -rpt_label {wildcards.db} \
        -norm meandiv -nperm 1000  -scoring_scheme classic \
        -create_svgs {params.svg} -make_sets true  -rnd_seed timestamp -zip_report false \
        -set_max 15000 -set_min 15 \
        -plot_top_x {params.nplot} -out ./gsea/{wildcards.fname} >> {log} 2>&1;
        then
            mv gsea/{wildcards.fname}/{wildcards.db}.GseaPreranked.*/ \
            gsea/{wildcards.fname}/{wildcards.db}.GseaPreranked/  >> {log} 2>&1;
        else
            mv gsea/{wildcards.fname}/*{wildcards.db}.GseaPreranked*/ \
            gsea/{wildcards.fname}/{wildcards.db}.GseaPreranked/  >> {log} 2>&1;
        fi
    
        cp workflow/envs/GSEA_ReadMe.html gsea/ >> {log} 2>&1;

        # # compress
        # cd gsea/{wildcards.fname}/
        # rm -f {wildcards.db}.GseaPreranked.zip
        # zip -rq {wildcards.db}.GseaPreranked.zip {wildcards.db}.GseaPreranked/ >> {log} 2>&1
        """

rule GSEA_compression:
    input:
        GSEA_OUTPUT(config)
    output:
        "gsea/{contrast}.tar.gz"
    resources:
        mem_mb=1000
    threads:
        4
    shell:
        "tar cf - gsea/{wildcards.contrast} | pigz -p {threads} > {output} "

rule GSEA_SingleBubblePlot:
    input:
        "gsea/{contrast}.tar.gz"
    output:
        touch('gsea_bubble/log/{contrast}.SingleBubblePlot.done')
    conda:
        "../envs/deseq2.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        'gsea_bubble/log/{contrast}.SingleBubblePlot.log'
    threads:
        1
    benchmark:
        'gsea_bubble/log/{contrast}.SingleBubblePlot.benchmark'
    shell:
        "Rscript workflow/script/gsea_bubble.R {input} {wildcards.contrast} &> {log}"

if config["GSEA_ANALYSIS"]:
    if config["START"] in ["FASTQ", "BAM", "COUNT"]:
        GSEA_compression_OUTPUT = expand("gsea/{contrast}.tar.gz",contrast=CONTRASTS_DE)
    else:
        GSEA_compression_OUTPUT = expand("gsea/{contrast}.tar.gz",contrast=config["RNKS"])

    rule GSEA_MultiBubblePlot:
        input:
            GSEA_compression_OUTPUT
        output:
            touch('gsea_bubble/log/MultiBubblePlot.done')
        conda:
            "../envs/deseq2.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4000
        priority: 100
        log:
            'gsea_bubble/log/MultiBubblePlot.log'
        threads:
            1
        benchmark:
            'gsea_bubble/log/MultiBubblePlot.benchmark'
        shell:
            "Rscript workflow/script/gsea_bubble.R {input} MultiBubblePlot &> {log}"

rule Read_Length_Detection:
    input:
        expand("bam_qc/stats/{sample}.stats.txt",sample=SAMPLES)
        if config['START'] in ["FASTQ", "BAM"]
        else "BAM file not provided. AND can't be generated from FASTQ, because FASTQ not provided"
    output:
        "meta/read_length.median.txt",
        "meta/read_length.max.txt"
    resources:
        mem_mb=1000,
    threads:
        1
    log:
        "meta/log/read_length.log"
    run:
        import re
        import statistics

        median_lens = []
        max_lens = []
        for f in input:
            p1 = re.compile("average length:\s*(\d*)")
            p2 = re.compile("maximum length:\s*(\d*)")
            for i, line in enumerate(open(f,"r")):
                for match1 in re.finditer(p1,line):
                    median_lens.append(match1.group(1))
                for match2 in re.finditer(p2,line):
                    max_lens.append(match2.group(1))
        median_lens = list(map(int,median_lens))
        max_lens = list(map(int,max_lens))
        print("average lengths detected from BAM:",median_lens,"\n")
        print("max lengths detected from BAM:",max_lens,"\n")
        if (len(set(max_lens)) < 1):
            sys.exit("read lengths detection from BAM failed")
        with open(output[0],"w") as out:
            out.write(str(int(statistics.median(median_lens))))
        with open(output[1],"w") as out:
            out.write(str(max(max_lens)))

if config['RMATS_ANALYSIS']:
    rule rMATS:
        #todo: update to docker: simple to install, faster, recommended on website
        input:
            b1=lambda wildcards: B1S[int(wildcards['ascn']) - 1],
            b2=lambda wildcards: B2S[int(wildcards['ascn']) - 1],
            gtf=config["GTF"],
            length_file="meta/read_length.median.txt",
            strand_file="meta/strandness.detected.txt",
        output:
            "rMATS.{ascn}/output/Results_JunctionCountsBased/SE.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/MXE.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A3SS.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A5SS.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/RI.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/SE.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/MXE.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A3SS.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A5SS.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/RI.MATS.JCEC.txt",
        conda:
            "../envs/rmats.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 10000,
        threads:
            4
        params:
            b1=lambda wildcards: B1[int(wildcards['ascn']) - 1],
            b2=lambda wildcards: B2[int(wildcards['ascn']) - 1],
            type="paired" if config['PAIR_END'] else "single",
            analysis="P" if config['PAIR_END'] else "U",
            length=get_read_length("meta/read_length.median.txt"),
            strandness=get_strandness("meta/strandness.detected.txt",config),
            MAX_FDR=config['MAX_FDR']
        log:
            "rMATS.{ascn}/rMATS.log"
        benchmark:
            "rMATS.{ascn}/rMATS.benchmark"
        shell:
            """
            rm -rf rMATS.{wildcards.ascn}
            mkdir -p rMATS.{wildcards.ascn}/
            mkdir -p rMATS.{wildcards.ascn}/params/
            echo {params.b1} > rMATS.{wildcards.ascn}/params/b1.txt
            echo {params.b2} > rMATS.{wildcards.ascn}/params/b2.txt

            rmats.py  \
            --b1 rMATS.{wildcards.ascn}/params/b1.txt \
            --b2  rMATS.{wildcards.ascn}/params/b2.txt \
            --gtf {input.gtf} \
            -t {params.type} \
            --readLength {params.length} \
            --variable-read-length \
            --libType {params.strandness} \
            --nthread {threads} \
            --tstat {threads} \
            --cstat 0.2 \
            --allow-clipping \
            --novelSS --mil 20 --mel 2000 \
            --od rMATS.{wildcards.ascn}/output/ \
            --tmp rMATS.{wildcards.ascn}/tmp/  > {log} 2>&1

            cp workflow/envs/rMATS_ReadMe.html rMATS.{wildcards.ascn}/output/  >>{log} 2>&1;
            rm -rf rMATS.{wildcards.ascn}/tmp/   >> {log} 2>&1;
            mkdir -p rMATS.{wildcards.ascn}/output/Results_JunctionCountsBased/ >> {log} 2>&1;
            mkdir -p rMATS.{wildcards.ascn}/output/Results_JunctionCountsAndExonCountsBased/  >> {log} 2>&1;
            mkdir -p rMATS.{wildcards.ascn}/output/fromGTF/ >> {log} 2>&1;
            mv rMATS.{wildcards.ascn}/output/*JCEC.txt rMATS.{wildcards.ascn}/output/Results_JunctionCountsAndExonCountsBased/ &>> {log}
            mv rMATS.{wildcards.ascn}/output/*JC.txt rMATS.{wildcards.ascn}/output/Results_JunctionCountsBased/ &>> {log}
            mv rMATS.{wildcards.ascn}/output/fromGTF*txt  rMATS.{wildcards.ascn}/output/fromGTF/ &>> {log}
            rm -f rMATS.{wildcards.ascn}/output/*raw.input* &>> {log}
            """

    rule rMATS_FLT:
        input:
            "rMATS.{ascn}/output/Results_JunctionCountsBased/SE.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/MXE.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A3SS.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A5SS.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/RI.MATS.JC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/SE.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/MXE.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A3SS.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A5SS.MATS.JCEC.txt",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/RI.MATS.JCEC.txt",
        output:
            "rMATS.{ascn}/output/Results_JunctionCountsBased/SE.MATS.JC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/MXE.MATS.JC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A3SS.MATS.JC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/A5SS.MATS.JC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsBased/RI.MATS.JC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/SE.MATS.JCEC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/MXE.MATS.JCEC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A3SS.MATS.JCEC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/A5SS.MATS.JCEC.sig.tsv",
            "rMATS.{ascn}/output/Results_JunctionCountsAndExonCountsBased/RI.MATS.JCEC.sig.tsv",
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4000
        threads:
            1
        run:
            for f in input:
                df = pd.read_table(f)
                out = df[df['PValue'] < 1]
                outname = re.sub(".txt$",".sig.tsv",f)
                out.to_csv(outname,sep="\t")
            import shutil

            shutil.make_archive('rMATS.' + str(wildcards.ascn),'zip', \
                root_dir="./",base_dir='rMATS.' + str(wildcards.ascn))
            shutil.move('rMATS.' + str(wildcards.ascn) + '.zip','rMATS.' + str(wildcards.ascn))

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
    input:
        bam="sorted_reads/{sample}.bam",
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
        "python workflow/script/dexseq_count.py -f bam -a {params.MAPQ} {params.readType} {params.strand} {DEXSeq_GFF} {input.bam} {output} &> {log}"
    # todo: maybe -r name (sort by name) is faster

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
            max_fdr = config['MAX_FDR'],
            min_lfc = config['MIN_LFC']
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

rule Genome_Faidx:
    input:
        config['GENOME']
    output:
        config['GENOME'] + '.fai'
    conda:
        "../envs/samtools.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    threads:
        1
    log:
        config['GENOME'] + '.fai' + '.log'
    shell:
        "samtools faidx {input}  &> {log}"


rule GATK_CreateSequenceDictionary:
    input:
        config['GENOME'],
    output:
        re.sub("fa$","dict",config['GENOME'])
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    conda:
        "../envs/gatk.yaml"  # test
    threads:
        1
    log:
        "log/GATK_CreateSequenceDictionary/GATK_CreateSequenceDictionary.log"
    benchmark:
        "log/GATK_CreateSequenceDictionary/GATK_CreateSequenceDictionary.benchmark"
    shell:
        "gatk CreateSequenceDictionary -R {input} &> {log}"


rule GATK_ASEReadCounter:
    input:
        genome=config['GENOME'],
        genome_idx=config['GENOME'] + '.fai',
        genome_ref=re.sub("fa$","dict",config['GENOME']),
        bam="sorted_reads/{sample}.bam",
        vcf=config['VCF']
    output:
        table="GATK_ASEReadCounter/{sample}.table"
    conda:
        "../envs/gatk.yaml"  # test
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    threads:
        1
    log:
        "log/GATK_ASEReadCounter/{sample}.log"
    benchmark:
        "log/GATK_ASEReadCounter/{sample}.benchmark"
    shell:
        "gatk ASEReadCounter -R {input.genome} -I {input.bam} -V {input.vcf} -O {output} \
        --min-depth-of-non-filtered-base 10 --min-mapping-quality 15 --min-base-quality 20 \
        &> {log};"  # todo: more thoughts on detailed parameters, e.g.    -drf DuplicateRead -U ALLOW_N_CIGAR_READS
        "D=GATK_ASEReadCounter; "
        "rm -f $D/$D.zip && [ -d $D ] && zip -rq  $D/$D.zip $D/ &>> {log};"


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


rule Create_DAG:
    input:
        "meta/configCheck.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    threads:
        1
    output:
        "Workflow_DAG.all.pdf",
        "rulegraph.pdf"
    log:
        "Workflow_DAG.all.pdf.log"
    shell:
        "snakemake --dag targets > dag 2> {log};"
        "cat dag|dot -Tpdf > Workflow_DAG.all.pdf 2>> {log};"
        "snakemake  --rulegraph targets > rulegraph; cat rulegraph| dot -Tpdf > rulegraph.pdf 2>> {log}"


rule reset:
    shell:
        """
        echo 'deleting files..'
        rm -rf lsf.log  meta/log/ gsea_bubble/ log/ workflow.log fastqc/ bam_qc/ trimmed/ mapped_reads/ sorted_reads/ bam_qc/ bigWig/ \
        feature_count/ fastq_salmon SalmonTE_output/ DESeq2/ salmon/ gsea/ gsea_compressed/ \
        GATK_ASEReadCounter/ DEXSeq_count/  DEXSeq/ rMATS.*/ CleanUpRNAseqQC/ CleanUpRNAseqDE/ \
        _STARgenome _STARtmp \
        feature_count_gene_level hisat2 stringtie \
        lsf.log Log.out nohup.out report.log report.html  dag Workflow_DAG.all.pdf \
        rulegraph  rulegraph.pdf report.log workflow.log log/ \
        meta/strandness.detected.txt  meta/decoder.txt meta/read_length.median.txt \
        meta/read_length.max.txt Workflow_DAG.all.pdf.log
        echo 'unlocking dir..'
        snakemake -j 1 --unlock
        """
