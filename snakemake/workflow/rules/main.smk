def checkFileInput(wildcards):
    check = "fastqc/CheckFile/CheckFile.{sample}.txt" if config['START'] == 'FASTQ' else 'Workflow_DAG.all.pdf'
    return check


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
