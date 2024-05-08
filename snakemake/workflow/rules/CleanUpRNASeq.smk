rule test_salmon_decoy:
    input:
        transcriptome="/home/rui.li-umw/genome/hg38_ensembl/Homo_sapiens.GRCh38.cdna.ncrna.fa.gz",
        genome=config['GENOME'],
    output:
        gentrome="salmon/decoy/gentrome.fasta.gz",
        decoys="salmon/decoy/decoys.txt",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        "salmon/log/decoys.log"
    benchmark:
        "salmon/log/decoys.benchmark"
    wrapper:
        "v3.10.2/bio/salmon/decoys"

rule salmon_index:
    input:
        sequences="salmon/decoy/gentrome.fasta.gz",
        decoys="salmon/decoy/decoys.txt"
    output:
        multiext(
            "salmon/transcriptome_index/",
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
        "log/salmon/transcriptome_index.log",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    params:
        # optional parameters
        extra="",
    wrapper:
        "v3.10.2/bio/salmon/index"

rule salmon_quant_reads_pe:
    input:
        r1="trimmed/{sample}.R1.fastq.gz"",
        r2="trimmed/{sample}.R2.fastq.gz",
        index="salmon/transcriptome_index",
    output:
        quant="salmon/{sample}/quant.sf",
        lib="salmon/{sample}/lib_format_counts.json",
    log:
        "log/salmon/{sample}.log",
    benchmark:
        "log/salmon/{sample}.benchmark",
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

#todo: SE

