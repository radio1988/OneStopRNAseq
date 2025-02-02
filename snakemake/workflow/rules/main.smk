localrules: Create_DAG, reset

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
        "workflow_full_DAG.pdf",
        "rulegraph.pdf"
    log:
        "workflow_full_DAG.pdf.log"
    shell:
        "snakemake --dag targets > dag 2> {log};"
        "cat dag|dot -Tpdf > workflow_full_DAG.pdf 2>> {log};"
        "snakemake  --rulegraph targets > rulegraph; cat rulegraph| dot -Tpdf > rulegraph.pdf 2>> {log}"


rule reset:
    run:
        shell(f"""
        echo 'deleting result and logs..'
        rm -rf lsf.log log/ meta/configCheck.log meta/configCheck.txt meta/log/ workflow.log 
        rm -rf dag rulegraph rulegraph.pdf workflow_full_DAG.pdf workflow_full_DAG.pdf.log
        rm -rf feature_count/ DESeq2/ gsea/ feature_count_gene_level/
        rm -rf trimmed/ fastqc/ mapped_reads/ sorted_reads/ bigWig/ bam_qc/
        rm -rf CleanUpRNAseqQC/ CleanUpRNAseqDE/
        rm -rf fastq_salmon/   SalmonTE_output/
        rm -rf salmon/ DEXSeq_count/  DEXSeq/ rMATS.*/
        rm -rf GATK_ASEReadCounter/   
        rm -rf _STARgenome _STARtmp 
        rm -rf hisat2 stringtie 
        rm -rf report.log report.html   workflow_full_DAG.pdf 
        rm -f meta/strandness.detected.txt  meta/decoder.txt meta/read_length.median.txt 
        rm -f meta/read_length.max.txt
                
        echo 'unlocking dir..'
        snakemake -j 1 --unlock
        """)

        if "MSHEET" in config and config["MSHEET"] and config["START"] == 'RNK':
            files = split_msheet_rnk_files(config)
            shell(f"rm -rf {' '.join(files)}")

