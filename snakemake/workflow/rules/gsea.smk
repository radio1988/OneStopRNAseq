"""
GSEA
GSEA_Bubble
"""
import pandas as pd
import os


rule GSEA:
    """
    config[START] == RNK, then find config[RNKS], else, find DESeq2/rnk or CleanUpRNAseqDE/rnk
    config[RNKS]  must be put in ./meta/ folder, and has suffix of txt or xlsx
    """
    input:
        rnk=lambda wildcards: input_rnk_fname1(wildcards,config),
        db=lambda wildcards: os.path.join(config['GSEA_DB_PATH'],wildcards["db"])
    output:
        html="gsea/{fname}/{db}.GseaPreranked/index.html",
        edb="gsea/{fname}/{db}.GseaPreranked/edb/results.edb"
    conda:
        "../envs/java11.yaml"  # test
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 8100
    params:
        svg=config["GSEA_PLOT_SVG"],
        nplot=config["GSEA_NPLOTS"],
        label_db=lambda wildcards: wildcards["db"][:-12],
        rnk_flat_file=lambda wildcards: input_rnk_fname2(wildcards,config),
        gmt_fmted=lambda wildcards: "gsea/gene_sets/" + wildcards['fname'] + "." + wildcards['db']
    threads:
        1
    log:
        "gsea/log/{fname}/{db}.log"
    benchmark:
        "gsea/log/{fname}/{db}.benchmark"
    shell:
        """
        echo rnkfile: {input.rnk} > {log}
        echo rnkfile.cleaned: {params.rnk_flat_file} >> {log}
        echo label_db: {params.label_db} >> {log}
        echo fname: {wildcards.fname} >> {log}
        echo dbname: {wildcards.db} >> {log}

        mkdir -p gsea gsea/gene_sets gsea/{wildcards.fname} >> {log} 2>&1
        rm -rf gsea/{wildcards.fname}/{wildcards.db}.GseaPreranked.*/  # avoid confusion
        rm -rf gsea/{wildcards.fname}/{wildcards.db}.GseaPreranked/  # avoid dir structure mistake
        rm -rf gsea/{wildcards.fname}/error_{wildcards.db}*GseaPreranked*/  # avoid confusion of temp files in multiple run attempts

        python workflow/script/rnk_to_upper.py {input.rnk} >> {log} 2>&1;  # fix gene symbol error; standardize file format (xlsx, rnk, etc)
        python workflow/script/gmt_to_upper.py -f {input.db} 1> {params.gmt_fmted} 2> {log}  # gene symbols to upper
        
        if workflow/envs/GSEA_4.3.2/gsea-cli.sh GSEAPreranked \
        -gmx {params.gmt_fmted} -rnk {params.rnk_flat_file} -rpt_label {wildcards.db} \
        -norm meandiv -nperm 1000  -scoring_scheme classic \
        -create_svgs {params.svg} -make_sets true  -rnd_seed timestamp -zip_report false \
        -set_max 15000 -set_min 15 \
        -plot_top_x {params.nplot} -out ./gsea/{wildcards.fname} >> {log} 2>&1;
        then
            mv gsea/{wildcards.fname}/{wildcards.db}.GseaPreranked.*/ \
            gsea/{wildcards.fname}/{wildcards.db}.GseaPreranked/  >> {log} 2>&1;
        else
            mv gsea/{wildcards.fname}/*{wildcards.db}.GseaPreranked*/ \
            gsea/{wildcards.fname}/{wildcards.db}.GseaPreranked/  >> {log} 2>&1;
        fi
    
        cp workflow/envs/GSEA_ReadMe.html gsea/ >> {log} 2>&1;

        # # compress
        # cd gsea/{wildcards.fname}/
        # rm -f {wildcards.db}.GseaPreranked.zip
        # zip -rq {wildcards.db}.GseaPreranked.zip {wildcards.db}.GseaPreranked/ >> {log} 2>&1
        """

rule GSEA_compression:
    input:
        GSEA_OUTPUT(config)
    output:
        "gsea/{contrast}.tar.gz"
    resources:
        mem_mb=1000
    threads:
        4
    benchmark:
        "gsea/log/{contrast}.tar.gz.benchmark"
    shell:
        "tar cf - -C gsea {wildcards.contrast} | pigz -p {threads} > {output} "

rule GSEA_SingleBubblePlot:
    input:
        "gsea/{contrast}.tar.gz"
    output:
        touch('gsea/gsea_bubble/log/{contrast}.SingleBubblePlot.done')
    conda:
        "../envs/deseq2.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000
    log:
        'gsea/gsea_bubble/log/{contrast}.SingleBubblePlot.log'
    threads:
        1
    benchmark:
        'gsea/gsea_bubble/log/{contrast}.SingleBubblePlot.benchmark'
    shell:
        "Rscript workflow/script/gsea_bubble.R {input} {wildcards.contrast} &> {log}"

if config["GSEA_ANALYSIS"]:
    if config["START"] in ["FASTQ", "BAM", "COUNT"]:
        GSEA_MultiBubblePlot_Input = lambda wildcards: expand(
            "gsea/{fname}/{db}.GseaPreranked/edb/results.edb",
            fname=DE_CONTRAST_NAMES,
            db=[wildcards.db]  # Use the single {db} wildcard
        )
    else:
        GSEA_MultiBubblePlot_Input = lambda wildcards: expand(
            "gsea/{fname}/{db}.GseaPreranked/edb/results.edb",
            fname=config["RNKS"],
            db=[wildcards.db]  # Use the single {db} wildcard
        )

    rule GSEA_MultiBubblePlot:
        input:
            GSEA_MultiBubblePlot_Input
        output:
            'gsea/gsea_bubble/{db}.pdf'
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4000
        priority: 100
        log:
            'gsea/gsea_bubble/log/MultiBubblePlot.{db}.log'
        threads:
            1
        benchmark:
            'gsea/gsea_bubble/log/MultiBubblePlot.{db}.benchmark'
        shell:
            "python workflow/script/gsea_bubble.py -edbs {input} -output {output} -alpha 0.05 -topn 1000 &> {log}"
