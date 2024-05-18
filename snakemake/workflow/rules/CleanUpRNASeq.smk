# CONFIG
GTF = config['GTF']
GENOME = config['GENOME']
# SALMON
TRANSCRIPTS = GTF + '.fa'
GENTROME = GTF + '.gentrome.fa'
DECOYS = GENTROME + '.decoys.txt'
# CleanUpRNAseq
META = config['META']
ENSDB = GTF + '.ensdb.sqlite'


rule gff_read:
    input:
        fasta=GENOME,
        annotation=GTF
    output:
        records=TRANSCRIPTS
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        GTF + '.fa.log'
    params:
        extra=""
    wrapper:
        "v3.10.2/bio/gffread"


rule salmon_decoy:
    input:
        transcriptome=TRANSCRIPTS,
        genome=GENOME
    output:
        gentrome=GENTROME,
        decoys=DECOYS
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        DECOYS + '.log'
    benchmark:
        DECOYS + '.benchmark'
    wrapper:
        "v3.10.2/bio/salmon/decoys"


rule salmon_index:
    input:
        sequences=GENTROME,
        decoys=DECOYS
    output:
        multiext(
            GENTROME + ".salmon_idx/",
            "complete_ref_lens.bin",
            "ctable.bin",
            "ctg_offsets.bin",
            "duplicate_clusters.tsv",
            "info.json",
            "mphf.bin",
            "pos.bin",
            "pre_indexing.log",
            "rank.bin",
            "refAccumLengths.bin",
            "ref_indexing.log",
            "reflengths.bin",
            "refseq.bin",
            "seq.bin",
            "versionInfo.json"
        ),
    log:
        GENTROME + ".salmon_idx/log",
    benchmark:
        GENTROME + ".salmon_idx/benchmark",
    # hg38-ensembl s       h:m:s   max_rss max_vms max_uss max_pss io_in   io_out  mean_load       cpu_time
    # 1378.6924       0:22:58 19463.66        42060.36        19451.28        19452.22        9.68    18896.98        446.52  6156.25
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    params:
        extra=""
    wrapper:
        "v3.10.2/bio/salmon/index"


if config["PAIR_END"]:
    rule salmon_quant_pe:
        input:
            r1="trimmed/{sample}.R1.fastq.gz",
            r2="trimmed/{sample}.R2.fastq.gz",
            index_flag=GENTROME + ".salmon_idx/complete_ref_lens.bin"
        output:
            quant="salmon/{sample}/quant.sf",
            lib="salmon/{sample}/lib_format_counts.json"
        log:
            "salmon/{sample}/log.txt",
        benchmark:
            "salmon/{sample}/benchmark.txt",
        params:
            # optional parameters
            index=GENTROME + ".salmon_idx/",
            libtype="A", #ISF
            extra="--seqBias --gcBias --posBias   --softclip  --softclipOverhangs",
            outdir="salmon/{sample}/"
        threads:
            16
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1500  # 18G for human
        conda:
            "../envs/salmon.yaml"  # docker: combinelab/salmon:latest
        shell:
            """
            salmon quant -i {params.index} -l {params.libtype} \
            -1 {input.r1} -2 {input.r2} -p {threads} {params.extra} \
            --validateMappings  -o {params.outdir} &> {log}
            """
else:
    rule salmon_quant_se:
        input:
            r="trimmed/{sample}.fastq.gz",
            index_flag=GENTROME + ".salmon_idx/complete_ref_lens.bin"
        output:
            quant="salmon/{sample}/quant.sf",
            lib="salmon/{sample}/lib_format_counts.json"
        log:
            "salmon/{sample}/log.txt",
        benchmark:
            "salmon/{sample}/benchmark.txt",
        params:
            # optional parameters
            index=GENTROME + ".salmon_idx/",
            libtype="U", #ISF
            extra="--seqBias --gcBias --posBias   --softclip  --softclipOverhangs",
            outdir="salmon/{sample}/"
        threads:
            16
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1500  # 18G for human
        conda:
            "../envs/salmon.yaml"  # docker: combinelab/salmon:latest
        shell:
            """
            salmon quant -i {params.index} -l {params.libtype} \
            -r {input.r} -p {threads} {params.extra} \
            --validateMappings  -o {params.outdir} &> {log}
            """


rule MakeEnsdb:
    input:
        gtf=GTF,
    output:
        ensdb=ENSDB  # ENSDB=GTF+'.ensdb.sqlite'
    log:
        ENSDB + ".log"
    conda:
        "../envs/cleanuprnaseq.yaml"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000
    script:
        "../script/make_ensdb.R"

rule MakeCleanUpMeta:
    input:
        META
    output:
        "CleanUpRNAseqQC/meta.cleanuprnaseq.csv"
    log:
        "CleanUpRNAseqQC/meta.cleanuprnaseq.csv.log"
    script:
        "../script/cleanupmakemeta.py"  # pandas

rule CleanUpRNAseqQC:
    input:
        meta="CleanUpRNAseqQC/meta.cleanuprnaseq.csv",
        bam=expand("mapped_reads/{sample}.bam",sample=SAMPLES),
        salmon=expand("salmon/{sample}/quant.sf",sample=SAMPLES),
        genome=GENOME,
        gtf=GTF,
        ensdb=ENSDB
    output:
        "CleanUpRNAseqQC/Fig7.PCA.showing.sample.variability.pdf",
        "CleanUpRNAseqQC/metadata.with.IR.rates.RDS"
    log:
        "CleanUpRNAseqQC/plots.log"
    threads:
        8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 6000
    conda:
        "../envs/cleanuprnaseq.yaml"
    script:
        "../script/cleanuprnaseq.R"


rule IR_DE:
    input:
        cnt=DESeq2_input(config),
        meta=config['META'],
        contrast=config["CONTRAST_DE"],
        ir_rate="CleanUpRNAseqQC/metadata.with.IR.rates.RDS"
    output:
        "CleanUpRNAseqDE/CleanUpRNAseqDE.html"
    conda:
        "../envs/deseq2.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
    params:
        rmd="'./CleanUpRNAseqDE/DESeq2.IR.Rmd'",
        fdr=MAX_FDR,
        lfc=MIN_LFC,
        independentFilter=config["independentFilter"],
        cooksCutoff=config["cooksCutoff"],
        blackSamples=config['blackSamples'] if 'blackSamples' in config else "",
        anno_tab=ANNO_TAB,
        o="'CleanUpRNAseqDE.html'"
    priority:
        100
    threads:
        1
    log:
        "CleanUpRNAseqDE/DESeq2.IR.log"
    benchmark:
        "CleanUpRNAseqDE/DESeq2.IR.benchmark"
    shell:
        'cp workflow/script/DESeq2.IR.Rmd CleanUpRNAseqDE; '
        
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
        
        'D=CleanUpRNAseqDE; rm -f $D/$D.zip && [ -d $D ] && zip -rq  $D/$D.zip $D/ >> {log} 2>&1;'


