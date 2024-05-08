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
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    params:
        # optional parameters
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
