GTF=config['GTF']
GENOME=config['GENOME']
META=config['META']
ENSDB=GTF+".ensdb.sqlite'

rule make_ensdb:
    input:
        gtf=GTF,
    output:
        ensdb=ENSDB
    script:
        "workflow/script/make_ensdb.R"

rule CleanUpMakeMeta:
    input:
        META
    output:
        "CleanUpRNASeqQC/meta.txt"
    script:
        "script/cleanupmakemeta.R"

rule CleanUpRNASeq:
    input:
        meta="CleanUpRNASeqQC/meta.txt", # todo
        bam=expand("mapped_reads/{sample}.bam", sample=SAMPLE),
        salmon=expand("salmon/{sample}/quant.sf", sample=SAMPLE),
        genome=GENOME,
        gtf=GTF,
        ensdb=ENSDB,
    output:
        "CleanUpRNASeqQC"
    script:
        "workflow/script/cleanuprnaseq.R"
