import pandas as pd
import sys
import re
import os
import math
import shutil


def check_config(config):
    """
    Check conflicts in config.yaml,
    call this in main snakefile
    sys.exit if not compatible
    """
    # INTRON and gDNA correction, only can choose one at most
    if config['START'] == 'FASTQ' and config["INTRON"] and config["CleanUpRNAseqCorrection"]:
        message = "INTRON mode and CleanUpRNAseqCorrection is incompatible.\n" + \
                  "gDNA correction is only possible for exon level rnaseq quantification"
        sys.exit(message)

    if config['ALIGNER'] != 'STAR' and config['ALIGNER'] != 'HISAT2':
        sys.exit("config['ALIGNER'] not STAR nor HISAT2")


# SPECIES and Analysis options
# Not now

def uncompress_gzip_genome_files(config):
    """
    # umcompress files (todo: make it a rule with checkpoints)
    """
    if config["GENOME"].endswith(".gz"):
        gunzip(config["GENOME"])
        config["GENOME"] = re.sub(".gz$","",config["GENOME"])

    if config["GTF"].endswith(".gz"):
        gunzip(config["GTF"])
        config["GTF"] = re.sub(".gz$","",config["GTF"])

    if config["VCF"].endswith(".gz"):
        gunzip(config["VCF"])
        config["VCF"] = re.sub(".gz$","",config["VCF"])

    return config


def check_fastq_size(config, SAMPLES):
    '''
    check size of all fastq files
    MAX_FILE_SIZE specified in config.yaml, unit GB
    if larger than MAX_FILE_SIZE, stop the workflow(mainly designed for website maintenance)
    '''
    if config['PAIR_END']:
        r1s = ["fastq/" + s + ".R1.fastq.gz" for s in SAMPLES]
        r2s = ["fastq/" + s + ".R2.fastq.gz" for s in SAMPLES]
        fastq_files = r1s + r2s
    else:
        fastq_files = ["fastq/" + s + ".fastq.gz" for s in SAMPLES]

    size = 0
    for f in fastq_files:
        size += os.path.getsize(f)  # will get true size even for softlinks

    if size / 1e9 > config['MAX_FASTQ_SIZE']:
        print('too large:',size / 1e9,'G detected')
        print('only',config['MAX_FASTQ_SIZE'],'G allowed')
        print(fastq_files)
        sys.exit('files too large')


def read_table(fname='meta/contrast.de.xlsx'):
    try:
        if fname.endswith(".txt"):
            df = pd.read_table(fname)
        elif fname.endswith(".csv"):
            df = pd.read_csv(fname)
        elif fname.endswith(".xlsx"):
            df = pd.read_excel(fname,engine='openpyxl')
        else:
            sys.exit("fname not xlsx nor txt")
    except:
        sys.exit(">>> Fatal Error: " + fname + " format error, it can't be read correctly." +
                 "\nPlease check if you saved txt file as xlsx or vice versa\n\n")
    return df


def gunzip(fname):
    import gzip
    outname = re.sub(".gz$","",fname)
    if outname == fname:
        sys.exit("config error: " + fname + " does not end with .gz")
    # todo: skip if uncompressed file exists (but how to check for previous uncompress integrity
    # todo: use bash gunzip to keep system timestamp of gz file
    with gzip.open(fname,'rb') as f_in:
        with open(outname,'wb') as f_out:
            shutil.copyfileobj(f_in,f_out)


### get CONTRASTS from meta/contrast.xlsx  ###
def get_contrast_fnames(fname):
    '''Get contrast name for DESeq2'''
    df = read_table(fname)
    CONTRASTS = []
    for j in range(df.shape[1]):
        c1 = df.iloc[0, j]
        c2 = df.iloc[1, j]
        c1 = c1.replace(" ","")
        c2 = c2.replace(" ","")
        c1 = re.sub(";$","",c1)  # remove extra ;
        c2 = re.sub(";$","",c2)
        c1 = re.sub(";",".",c1)  # remove extra ;
        c2 = re.sub(";",".",c2)
        name = c1 + "_vs_" + c2
        if len(name) > 100:
            CONTRASTS.append("contrast" + str(j + 1))
        else:
            CONTRASTS.append(c1 + "_vs_" + c2)
    # CONTRATS = ['KO_D8_vs_KO_D0', 'WT_D8.KO_D8_vs_WT_D0.KO_D0', 'KO_D0.KO_D2.KO_D8_vs_WT_D0.WT_D2.WT_D8']
    CONTRASTS = [x.replace("-",'.') for x in CONTRASTS]
    return CONTRASTS


def get_contrast_groups(fname):
    df2 = read_table(fname)
    C1S = [];
    C2S = []
    for j in range(df2.shape[1]):
        c1 = df2.iloc[0, j].strip()
        c2 = df2.iloc[1, j].strip()
        c1 = c1.replace(" ","")
        c2 = c2.replace(" ","")
        c1 = re.sub(";$","",c1)  # remove extra ;
        c2 = re.sub(";$","",c2)
        c1s = c1.split(";")
        c2s = c2.split(";")
        C1S.append(c1s)
        C2S.append(c2s)
    return [C1S, C2S]


def get_dict_from_meta(fname):
    df = read_table(fname)
    d = {}
    for i in range(df.shape[0]):
        sample = df.iloc[i, 0]
        group = df.iloc[i, 1]
        if group in d:
            d[group].append(sample)
        else:
            d[group] = []
            d[group].append(sample)
    return d


def swap_list_element(l, d):
    '''internal'''
    return [swap_list_element(x,d) if isinstance(x,list) else d.get(x,x) for x in l]


def flatten(l):
    '''internal'''
    flat_list = [item for sublist in l for item in sublist]
    return flat_list


def collapseL2(L3):
    '''internal'''
    L2 = []
    for l in L3:
        L2.append(flatten(l))
    return L2


def s2b(S):
    '''internal'''
    K = flatten(S)
    V = ["sorted_reads/" + x + ".bam" for x in K]
    return dict(zip(K,V))


def G2B_workflow(G, g2s):
    S = swap_list_element(G,g2s)
    S = collapseL2(S)
    s2b_dict = s2b(S)
    BS = swap_list_element(S,s2b_dict)
    B = [",".join(x) for x in BS]
    return BS, B


###  Input functions (functions to create input fnames for rules)  ###
def Merge_TE_and_Gene_input(config):
    if config['MODE'] == 'strict' and config['INTRON']:
        return "feature_count_gene_level/counts.strict.txt"
    if config['MODE'] == 'strict':
        return "feature_count/counts.strict.txt"
    if config['MODE'] == 'liberal' and config['INTRON']:
        return "feature_count_gene_level/counts.liberal.txt"
    if config['MODE'] == 'liberal':
        return "feature_count/counts.liberal.txt"
    raise Exception('config error')


def DESeq2_input(config):
    if config['START'] in ['COUNT']:
        return config['COUNT_FILE']

    if config['TE_ANALYSIS']:
        if config['INTRON']:
            return "feature_count_gene_level/TE_included.txt"
        else:
            return "feature_count/TE_included.txt"
    else:
        if config['INTRON']:
            folder = 'feature_count_gene_level'
        else:
            folder = 'feature_count'
        return folder + '/counts.' + config['MODE'] + '.txt'


def input_rnk_fname1(wildcards, config):
    if config['START'] == 'RNK':
        fname1 = 'meta/' + wildcards['fname']  # meta/test1.rnk.txt
    elif config['START'] == 'FASTQ' and config['CleanUpRNAseqCorrection']:
        fname1 = "CleanUpRNAseqDE/rnk/" + wildcards['fname'] + ".rnk"
    else:
        fname1 = "DESeq2/rnk/" + wildcards['fname'] + ".rnk"
    return fname1


def rnk_fname1_to_fname2(fname1):
    '''
    xxx.rnk -> xxx.rnk.txt
    xxx.rnk.xlsx -> xxx.rnk.txt
    xxx.rnk.txt -> xxx.rnk.txt
    internal
    '''
    import re
    import sys
    if fname1.endswith('.rnk.xlsx'):
        return re.sub('.rnk.xlsx$','.rnk.txt',fname1)
    elif fname1.endswith('.rnk'):
        return fname1 + ".txt"
    elif fname1.endswith('.rnk.txt'):
        return fname1
    else:
        sys.stderr.write(fname1)
        sys.exit('rnk file type not rnk, rnk.txt, or rnk.xlsx')
        return None


def input_rnk_fname2(wildcards, config):
    '''
    the corresponding flat rnk.txt file name (fname1) for corresponding fname1
    '''
    fname1 = input_rnk_fname1(wildcards,config)
    return rnk_fname1_to_fname2(fname1)


### Get parameters ###
RMATS_STRANDNESS = {0: 'fr-unstranded', 1: 'fr-firststrand', 2: 'fr-secondstrand'}


def get_strandness(strandFile="meta/strandness.detected.txt", config="config_dict"):
    try:
        with open(strandFile,"r") as file:
            txt = file.readline()
        p = re.compile("counts.s(.)\.[liberal|strict]")
        res = p.search(txt)
        if res.group(1) is None:
            sys.exit("strandness detection wrong")
        strand = (RMATS_STRANDNESS[int(res.group(1))])
        return strand
    except FileNotFoundError:
        sys.stderr.write(strandFile + "will be found in real run, not in dry run\n")
        return None


def get_strandness_for_dexseq(strandFile="meta/strandness.detected.txt"):
    book = {0: '-s no', 1: '-s yes', 2: '-s reverse'}
    try:
        with open(strandFile,"r") as file:
            txt = file.readline()
        p = re.compile("counts.s(.)\.[liberal|strict]")
        res = p.search(txt)
        if res.group(1) is None:
            sys.exit("strandness detection wrong")
        strand = book[int(res.group(1))]
        return strand
    except FileNotFoundError:
        sys.stderr.write(strandFile + "will be found in real run, not in dry run\n")
        return "Strand File not found"


def get_strandness_for_hisat2_PE(strandFile="meta/strandness.detected.txt"):
    book = {0: ' ', 1: '--rna-strandness FR', 2: '--rna-strandness RF'}
    try:
        with open(strandFile,"r") as file:
            txt = file.readline()
        p = re.compile("counts.s(.)\.[liberal|strict]")
        res = p.search(txt)
        if res.group(1) is None:
            sys.exit("strandness detection wrong")
        strand = book[int(res.group(1))]
        return strand
    except FileNotFoundError:
        sys.stderr.write(strandFile + "will be found in real run, not in dry run\n")
        return "Strand File not found"


def get_strandness_for_hisat2_SE(strandFile="meta/strandness.detected.txt"):
    book = {0: ' ', 1: '--rna-strandness F', 2: '--rna-strandness R'}
    try:
        with open(strandFile,"r") as file:
            txt = file.readline()
        p = re.compile("counts.s(.)\.[liberal|strict]")
        res = p.search(txt)
        if res.group(1) is None:
            sys.exit("strandness detection wrong")
        strand = book[int(res.group(1))]
        return strand
    except FileNotFoundError:
        sys.stderr.write(strandFile + "will be found in real run, not in dry run\n")
        return "Strand File not found"


def get_strandness_for_stringtie(strandFile="meta/strandness.detected.txt"):
    book = {0: ' ', 1: '--fr', 2: '--rf'}
    try:
        with open(strandFile,"r") as file:
            txt = file.readline()
        p = re.compile("counts.s(.)\.[liberal|strict]")
        res = p.search(txt)
        if res.group(1) is None:
            sys.exit("strandness detection wrong")
        strand = book[int(res.group(1))]
        return (strand)
    except FileNotFoundError:
        sys.stderr.write("meta/strandness.detected.txt will be found in real run, not in dry run\n")
        return "Strand File not found"


def get_read_length(lengthFile="meta/read_length.txt"):
    try:
        with open(lengthFile,'r') as file:
            txt = file.readline()
            length = int(float(txt))
        return (length)
    except FileNotFoundError:
        sys.stderr.write(lengthFile + "not found in dry run, will be found in real run\n")
        return None


def strand_detection_input(config):
    folder = 'feature_count_gene_level/' \
        if config['INTRON'] \
        else 'feature_count/'
    mode = config["MODE"]
    return ([folder + 'counts.s' + str(s) + '.' + config['MODE'] + '.txt.summary' \
             for s in config['STRAND']])


def mapped_bam_inputs(config, SAMPLES):
    if config['ALIGNER'] == 'STAR':
        return (["mapped_reads/{}.bam".format(s) for s in SAMPLES])
    elif config['ALIGNER'] == 'HISAT2':
        return (["hisat2/{}.bam".format(s) for s in SAMPLES])
    else:
        sys.exit("config['ALIGNER'] not recognized")


def mapped_bam_single_input(config):
    if config['ALIGNER'] == 'STAR':
        return "mapped_reads/{sample}.bam"
    elif config['ALIGNER'] == 'HISAT2':
        return "hisat2/{sample}.bam"
    else:
        sys.exit("config['ALIGNER'] not recognized")


def check_contrast_file(fname="meta/contrast.de.txt"):
    df = read_table(fname)
    if df.transpose().duplicated().any():
        raise ValueError(fname + " has duplicated contrasts, please fix and re-run")
    if df.shape[0] != 2:
        raise ValueError(
            fname +
            " does not have three rows:\n" +
            "1. contrast-name, e.g. KO_D8_vs_WT\n" +
            "2. A group names, e.g. KO_D8 \n" +
            "3. B group names, e.g. WT_D0; WT_D2; WT_D8\n" +
            "group names can be separated by `;` to combine GROUP_LABELS\n\n")


def check_meta_file(fname="meta/meta.txt"):
    df = read_table(fname)
    if df.duplicated().any():
        raise ValueError(fname + " has duplicated rows, please fix and re-run")
    if df.shape[1] != 3:
        raise ValueError(fname + " does not have three columns: SAMPLE_LABEL\tGROUP_LABEL\tBATCH\n\n")


def check_meta_data(config):
    if config['DESEQ2_ANALYSIS']:
        check_contrast_file(config['CONTRAST_DE'])
    if config['RMATS_ANALYSIS'] or config['DEXSEQ_ANALYSIS']:
        check_contrast_file(config['CONTRAST_AS'])
    if config['DESEQ2_ANALYSIS'] or config['RMATS_ANALYSIS'] or config['DEXSEQ_ANALYSIS']:
        check_meta_file(config['META'])


def GSEA_OUTPUT(config):
    L = []
    if config["GSEA_ANALYSIS"]:
        gsea_dbs = []
        for f in os.listdir(config['GSEA_DB_PATH']):
            if f.endswith('.gmt'):
                gsea_dbs.append(f)
        gsea_dbs = [os.path.basename(x) for x in gsea_dbs]
        if len(gsea_dbs) < 1:
            sys.exit("No gsea databases found in " + config['GSEA_DB_PATH'])
        if config["START"] in ["FASTQ", "BAM", "COUNT"]:
            L = expand("gsea/{contrast}/{db}.GseaPreranked/index.html",contrast=CONTRASTS_DE,db=gsea_dbs)
        else:
            L = expand("gsea/{contrast}/{db}.GseaPreranked/index.html",contrast=config["RNKS"],db=gsea_dbs)
    else:
        L = ["Workflow_DAG.all.pdf"]
    return L


def GSEACOMPRESS_OUTPUT(config):
    L = []
    if config["GSEA_ANALYSIS"]:
        if config["START"] in ["FASTQ", "BAM", "COUNT"]:
            L = expand("gsea/{contrast}.tar.gz",contrast=CONTRASTS_DE)
        else:
            L = expand("gsea/{contrast}.tar.gz",contrast=config["RNKS"])
    else:
        L = ["Workflow_DAG.all.pdf"]
    return L


def GSEA_SINGLEBUBBLE_OUTPUT(config):
    L = []
    if config["GSEA_ANALYSIS"]:
        if config["START"] in ["FASTQ", "BAM", "COUNT"]:
            L = expand("gsea/gsea_bubble/log/{contrast}.SingleBubblePlot.done",contrast=CONTRASTS_DE)
        else:
            L = expand("gsea/gsea_bubble/log/{contrast}.SingleBubblePlot.done",contrast=config["RNKS"])
    else:
        L = ["Workflow_DAG.all.pdf"]
    return L
