# CONFIG
import pandas as pd
import sys

GTF = config['GTF']
GENOME = config['GENOME']
# SALMON
TRANSCRIPTS = GTF + '.fa'
GENTROME = GTF + '.gentrome.fa'
DECOYS = GENTROME + '.decoys.txt'
# CleanUpRNAseq
META = config['META']
ENSDB = GTF + '.ensdb.sqlite'
if config['PAIR_END']:
    LIBTYPES = ['ISF', 'ISR', 'IU']  # salmon strand libtype
else:
    LIBTYPES = ['SF', 'SR', 'U']  # salmon strand libtype


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
        extra="-g " + GENOME,
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
            index_flag=GENTROME + ".salmon_idx/complete_ref_lens.bin",
            check="meta/CheckFile/CheckFile.{sample}.txt"
        output:
            quant="salmon/{libtype}/{sample}/quant.sf",
            lib="salmon/{libtype}/{sample}/lib_format_counts.json"
        log:
            "salmon/{libtype}/{sample}/log.txt"
        benchmark:
            "salmon/{libtype}/{sample}/benchmark.txt"
        params:
            # optional parameters
            index=GENTROME + ".salmon_idx",
            extra="--seqBias --gcBias --posBias   --softclip  --softclipOverhangs",
            outdir="salmon/{libtype}/{sample}/"
        threads:
            16
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1500  # 18G for human
        conda:
            "../envs/salmon.yaml"  # docker: combinelab/salmon:latest
        shell:
            """
            salmon quant -i {params.index} -l {wildcards.libtype} \
            -1 {input.r1} -2 {input.r2} -p {threads} {params.extra} \
            --validateMappings  -o {params.outdir} &> {log}
            """
else:
    rule salmon_quant_se:
        input:
            r="trimmed/{sample}.fastq.gz",
            index_flag=GENTROME + ".salmon_idx/complete_ref_lens.bin",
            check="meta/CheckFile/CheckFile.{sample}.txt"
        output:
            quant="salmon/{libtype}/{sample}/quant.sf",
            lib="salmon/{libtype}/{sample}/lib_format_counts.json"
        log:
            "salmon/{libtype}/{sample}/log.txt"
        benchmark:
            "salmon/{libtype}/{sample}/benchmark.txt"
        params:
            # optional parameters
            index=GENTROME + ".salmon_idx",
            extra="--seqBias --gcBias --posBias   --softclip  --softclipOverhangs",
            outdir="salmon/{libtype}/{sample}/"
        threads:
            16
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1500  # 18G for human
        conda:
            "../envs/salmon.yaml"  # docker: combinelab/salmon:latest
        shell:
            """
            salmon quant -i {params.index} -l {wildcards.libtype} \
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
        META,
        "meta/strandness.detected.txt"
    output:
        "CleanUpRNAseqQC/meta.cleanuprnaseq.csv"
    log:
        "CleanUpRNAseqQC/meta.cleanuprnaseq.csv.log"
    script:
        "../script/cleanuprnaseq.makemeta.py"


def get_cleanuprnaseq_libtype(strandFile="meta/strandness.detected.txt"):
    if config["PAIR_END"]:
        book = {0: 'IU', 1: "ISF", 2: 'ISR'}  # PE
    else:
        book = {0: 'U', 1: "SF", 2: 'SR'}  # SE

    try:
        with open(strandFile,"r") as file:
            txt = file.readline()
        p = re.compile("counts.s(.)\.[liberal|strict]")
        res = p.search(txt)
        if res.group(1) is None:
            sys.exit("strandness detection wrong")
        libtype = book[int(res.group(1))]
        return (libtype)
    except FileNotFoundError:
        sys.stderr.write(strandFile + "will be found in real run, not in dry run\n")
        return ("placeholder")


def cleanuprnaseqqc_input(config):
    libtype = get_cleanuprnaseq_libtype()
    files = ["salmon/{}/".format(libtype) + sample + "/quant.sf" for sample in SAMPLES]
    return (files)


rule CleanUpRNAseqQC:
    input:
        meta="CleanUpRNAseqQC/meta.cleanuprnaseq.csv",
        bam=expand("mapped_reads/{sample}.bam",sample=SAMPLES),
        # salmon=cleanuprnaseqqc_input,  # PE/SE aware, strand aware
        salmon=expand("salmon/{libtype}/{sample}/quant.sf",libtype=LIBTYPES,sample=SAMPLES),
        genome=GENOME,
        gtf=GTF,
        ensdb=ENSDB,
        strand_detection="meta/strandness.detected.txt"
    output:
        "CleanUpRNAseqQC/Fig7.PCA.showing.sample.variability.pdf",
        "CleanUpRNAseqQC/metadata.with.IR.rates.RDS",
        "CleanUpRNAseqQC/Diagnostic.plots.objects.RDS"
    log:
        "CleanUpRNAseqQC/CleanUpRNAseqQC.log"
    benchmark:
        "CleanUpRNAseqQC/CleanUpRNAseqQC.benchmark"
    threads:
        16
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 500
    conda:
        "../envs/cleanuprnaseq.yaml"
    script:
        "../script/cleanuprnaseq.qc.R"


def get_salmon_files(config):  # meta.csv not possibly ready at workflow construction
    fname = "CleanUpRNAseqQC/meta.cleanuprnaseq.csv"
    try:
        df = pd.read_csv(fname)
    except FileNotFoundError:
        sys.stderr.write(fname + "will be found in real run, not in dry run\n")
        return (None)
    if 'salmon_quant_file' in df.columns:
        L = list(df['salmon_quant_file'])
    else:
        sys.exit(fname + " does not contain correct salmon_quant_file_columns\n")
    if 'salmon_quant_file_strand' in df.columns:
        L.extend(list(df['salmon_quant_file_strand']))
    if 'salmon_quant_file_reverse_strand' in df.columns:
        L.extend(list(df['salmon_quant_file_reverse_strand']))
    return (L)


rule CleanUpRNAseqCorrection:
    input:
        qc_rds="CleanUpRNAseqQC/Diagnostic.plots.objects.RDS",
        meta_rds="CleanUpRNAseqQC/metadata.with.IR.rates.RDS",
        meta_tab="CleanUpRNAseqQC/meta.cleanuprnaseq.csv",
        ensdb=ENSDB,
        salmon=expand("salmon/{libtype}/{sample}/quant.sf",sample=SAMPLES,libtype=LIBTYPES),# PE/SE aware, all LIBTYPES,
        # salmon = get_salmon_files,
        strandness="meta/strandness.detected.txt"
    output:
        "CleanUpRNAseqQC/global.corrected.count.csv"
    log:
        "CleanUpRNAseqQC/CleanUpRNAseqCorrection.log"
    benchmark:
        "CleanUpRNAseqQC/CleanUpRNAseqCorrection.benchmark"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    conda:
        "../envs/cleanuprnaseq.yaml"
    script:
        "../script/cleanuprnaseq.correction.R"


def deseq2_ir_df_is_enough(config):
    if config['START'] != 'FASTQ':
        return None
    output = True
    if config['META'].endswith('.csv'):
        meta = pd.read_csv(config['META'])
    elif config['META'].endswith('.xlsx'):
        meta = pd.read_excel(config['META'])
    elif config['META'].endswith('.txt'):
        meta = pd.read_table(config['META'])
    else:
        sys.exit(config['META'] + " not right suffix (format)")

    if config['CONTRAST_DE'].endswith('.csv'):
        contrast = pd.read_csv(config['CONTRAST_DE'])
    elif config['CONTRAST_DE'].endswith('.xlsx'):
        contrast = pd.read_excel(config['CONTRAST_DE'])
    elif config['CONTRAST_DE'].endswith('.txt'):
        contrast = pd.read_table(config['CONTRAST_DE'])
    else:
        sys.exit(config['CONTRAST_DE'] + " not right suffix (format)")

    samples = list(meta.iloc[:, 0])
    groups = list(meta.iloc[:, 1])
    batches = list(meta.iloc[:, 2])

    if len(set(samples)) != len(samples):
        sys.exit("In config['META'], duplicated samples exist")

    for j in range(contrast.shape[1]):
        contrast_groups_init = list(contrast.iloc[:, j])  # [N, P, N;P]
        contrast_groups = []
        for string in contrast_groups_init:
            if ";" in string:
                gs = string.split(";")
                contrast_groups.extend(gs)
            else:
                contrast_groups.append(string)
        relevant = [x in contrast_groups for x in groups]
        n_sample = sum(relevant)
        n_comparison = 1
        n_ir = 1
        batch_groups = [data for data, filter in zip(batches,relevant) if filter]
        n_batch = len(set(batch_groups)) - 1  # N_batch - 1 # 1 -> 0, 2 -> 1
        if n_sample - n_comparison - n_ir - n_batch <= 1:
            output = False
    return (output)


if config['START'] != 'RNK':
    if deseq2_ir_df_is_enough(config):
        rule DESeq2_IR:
            input:
                cnt=DESeq2_input(config),
                meta=config['META'],
                contrast=config["CONTRAST_DE"],
                ir_rate="CleanUpRNAseqQC/metadata.with.IR.rates.RDS"
            output:
                "CleanUpRNAseqDE/CleanUpRNAseqDE.html",
                expand("CleanUpRNAseqDE/rnk/{contrast}.rnk",contrast=CONTRASTS_DE)
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
    else:
        rule DESeq2_DNA_Corrected:
            input:
                cnt="CleanUpRNAseqQC/global.corrected.count.csv",
                meta=config['META'],
                contrast=config["CONTRAST_DE"],
            output:
                "CleanUpRNAseqDE/CleanUpRNAseqDE.html",
                expand("CleanUpRNAseqDE/rnk/{contrast}.rnk",contrast=CONTRASTS_DE)
            conda:
                "../envs/deseq2.yaml"
            resources:
                mem_mb=lambda wildcards, attempt: attempt * 4000,
            params:
                dir="./CleanUpRNAseqDE/",
                rmd="'./CleanUpRNAseqDE/DESeq2.Rmd'",
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
                "DESeq2/DESeq2.log"
            benchmark:
                "DESeq2/DESeq2.benchmark"
            shell:
                #Rscript -e rmarkdown::render({params.rmd})
                'cp workflow/script/DESeq2.Rmd {params.dir}; '
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
