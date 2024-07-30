"""
ASEReadCounter
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
