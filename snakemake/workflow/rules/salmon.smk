GTF=config['GTF']
GENOME=config['GENOME']

TRANSCRIPTS=GTF+'.fa'
GENTROME=GTF+'.gentrome.fa'
DECOYS=GENTROME+'.decoys.txt'

rule gff_read:
    input:
        fasta=GENOME,
        annotation=GTF,
    output:
        records=TRANSCRIPTS,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        GTF+'.fa.log',
    params:
        extra="",
    wrapper:
        "v3.10.2/bio/gffread"

rule salmon_decoy:
    input:
        transcriptome=TRANSCRIPTS,
        genome=GENOME,
    output:
        gentrome=GENTROME,
        decoys=DECOYS,
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        DECOYS+'.log'
    benchmark:
        DECOYS+'.benchmark'
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
            "versionInfo.json",
        ),
    log:
        GENTROME + ".salmon_idx/log",
    benchmark:
        GENTROME + ".salmon_idx/benchmark",
        #hg38-ensembl s       h:m:s   max_rss max_vms max_uss max_pss io_in   io_out  mean_load       cpu_time
        #1378.6924       0:22:58 19463.66        42060.36        19451.28        19452.22        9.68    18896.98        446.52  6156.25
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    params:
        extra="",
    wrapper:
        "v3.10.2/bio/salmon/index"


if config["PAIR_END"]:
    rule salmon_quant_pe:
        input:
            r1="trimmed/{sample}.R1.fastq.gz",
            r2="trimmed/{sample}.R2.fastq.gz",
            index=GENTROME + ".salmon_idx/complete_ref_lens.bin",
        output:
            quant="salmon/{sample}/quant.sf",
            lib="salmon/{sample}/lib_format_counts.json",
        log:
            "salmon/{sample}/log.txt",
        benchmark:
            "salmon/{sample}/benchmark.txt",
        params:
            # optional parameters
            libtype="A",
            extra="",
        threads:
            4
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4000
        wrapper:
            "v3.10.2/bio/salmon/quant"
else:
    rule salmon_quant_se:
            input:
                r="trimmed/{sample}.fastq.gz",
                index=GENTROME + ".salmon_idx/complete_ref_lens.bin",
            output:
                quant="salmon/{sample}/quant.sf",
                lib="salmon/{sample}/lib_format_counts.json",
            log:
                "salmon/{sample}/log.txt",
            benchmark:
                "salmon/{sample}/benchmark",
            params:
                # optional parameters
                libtype="A",
                extra="",
            threads:
                4
            resources:
                mem_mb=lambda wildcards, attempt: attempt * 4000
            wrapper:
                "v3.10.2/bio/salmon/quant"