"""
GSEA
GSEA_Bubble
"""
import pandas as pd
import os

def split_msheet_rnk_file(config):
    if config['START'] == "RNK" and 'MSHEET' in config and config['MSHEET']:
        # check config
        if len(config['RNKS']) > 1:
            raise ValueError("If MSHEET is True, only one RNK file is allowed")

        msheet_fname = config['RNKS'][0]
        if not msheet_fname.endswith(".xlsx"):
            raise ValueError("If MSHEET is True, the RNK file must be xlsx")

        #split sheets
        os.makedirs("meta", exist_ok = True)  # Path.cwd is analysis root
        dfs = pd.read_excel(msheet_fname, sheet_name = None)
        rnk_file_names = []
        for sheet_name, sheet_df in dfs.items():
            comparison_name = sheet_df.columns[0]
            sheet_df.to_csv(f"meta/{comparison_name}.rnk.txt", sep = "\t", index = False)
            rnk_file_names.append(f"{comparison_name}.rnk.txt")

    print(rnk_file_names)
    return rnk_file_names


if config['START'] == "RNK" and 'MSHEET' in config and config['MSHEET']:
    rnk_file_names = split_msheet_rnk_file(config)
    config['RNKS'] = rnk_file_names  # only basename of rnk files


rule GSEA:
    """
    config[START] == RNK, then find config[RNKS], else, find DESeq2/rnk or CleanUpRNAseqDE/rnk
    config[RNKS]  must be put in ./meta/ folder, and has suffix of txt or xlsx
    """
    input:
        rnk=lambda wildcards: input_rnk_fname1(wildcards,config),
        db=lambda wildcards: os.path.join(config['GSEA_DB_PATH'],wildcards["db"])
    output:
        "gsea/{fname}/{db}.GseaPreranked/index.html"
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
        GSEA_compression_OUTPUT = expand("gsea/{contrast}.tar.gz",contrast=CONTRASTS_DE)
    else:
        GSEA_compression_OUTPUT = expand("gsea/{contrast}.tar.gz",contrast=config["RNKS"])

    rule GSEA_MultiBubblePlot:
        input:
            GSEA_compression_OUTPUT
        output:
            touch('gsea/gsea_bubble/log/MultiBubblePlot.done')
        conda:
            "../envs/deseq2.yaml"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 4000
        priority: 100
        log:
            'gsea/gsea_bubble/log/MultiBubblePlot.log'
        threads:
            1
        benchmark:
            'gsea/gsea_bubble/log/MultiBubblePlot.benchmark'
        shell:
            "Rscript workflow/script/gsea_bubble.R {input} MultiBubblePlot &> {log}"

    rule GSEA_Bubble_Compression:
        input:
            'gsea/gsea_bubble/log/MultiBubblePlot.done'
        output:
            'gsea/gsea_bubble.tar.gz'
        benchmark:
            'gsea/log/gsea_bubble.tar.gz.benchmark'
        resources:
            mem_mb=1000
        threads:
            4
        shell:
            "tar cf - -C gsea gsea_bubble | pigz -p {threads} > {output} "
