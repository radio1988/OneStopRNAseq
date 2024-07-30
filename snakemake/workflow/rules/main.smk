def checkFileInput(wildcards):
    check="fastqc/CheckFile/CheckFile.{sample}.txt" if config['START'] == 'FASTQ' else 'Workflow_DAG.all.pdf'
    return check






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


rule Create_DAG:
    input:
        "config.yaml"
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
