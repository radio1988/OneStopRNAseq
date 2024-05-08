rule test_salmon_decoy:
    input:
        transcriptome="/home/rui.li-umw/genome/hg38_ensembl/Homo_sapiens.GRCh38.cdna.ncrna.fa",
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
        "salmon/transcriptome_index"
    log:
        "salmon/transcriptome_index.log",
    benchmark:
        "salmon/{sample}/benchmark.txt",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8000
    params:
        # optional parameters
        extra="",
    wrapper:
        "v3.10.2/bio/salmon/index"

rule salmon_quant_reads_pe:
    input:
        r1="trimmed/{sample}.R1.fastq.gz",
        r2="trimmed/{sample}.R2.fastq.gz",
        index="salmon/transcriptome_index",
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

#todo: SE

