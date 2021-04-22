import pandas as pd
import sys
import re
import os
import math
import shutil


def gunzip(fname):
    import gzip
    outname = re.sub(".gz$", "", fname)
    if outname == fname:
        sys.exit("config error:", fname, "does not end with .gz")
    # todo: skip if uncompressed file exists (but how to check for previous uncompress integrity
    # todo: use bash gunzip to keep system timestamp of gz file
    with gzip.open(fname, 'rb') as f_in:
        with open(outname, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)


### get CONTRASTS from meta/contrast.xlsx  ###
def get_contrast_fnames (fname):
    '''Get contrast name for DESeq2'''
    df=pd.read_excel(fname, engine='openpyxl')
    CONTRASTS = []
    for j in range(df.shape[1]):
        c1 = df.iloc[0,j]
        c2 = df.iloc[1,j]
        c1 = c1.replace(" ", "")
        c2 = c2.replace(" ", "")
        c1 = re.sub(";$", "", c1)  # remove extra ;
        c2 = re.sub(";$", "", c2)
        c1 = re.sub(";", ".", c1)  # remove extra ;
        c2 = re.sub(";", ".", c2)
        name = c1 + "_vs_" + c2
        if len(name) > 200:
            CONTRASTS.append("contrast" + str(j+1))
        else:
            CONTRASTS.append(c1 + "_vs_" + c2)
    # CONTRATS = ['KO_D8_vs_KO_D0', 'WT_D8.KO_D8_vs_WT_D0.KO_D0', 'KO_D0.KO_D2.KO_D8_vs_WT_D0.WT_D2.WT_D8']
    return (CONTRASTS)

def get_contrast_groups (fname):
    df2=pd.read_excel(fname, engine='openpyxl')
    C1S = []; C2S = []
    for j in range(df2.shape[1]):
        c1 = df2.iloc[0, j].strip()
        c2 = df2.iloc[1, j].strip()
        c1 = c1.replace(" ", "")
        c2 = c2.replace(" ", "")
        c1 = re.sub(";$", "", c1)  # remove extra ;
        c2 = re.sub(";$", "", c2)
        c1s = c1.split(";")
        c2s = c2.split(";")
        C1S.append(c1s)
        C2S.append(c2s)
    return ([C1S, C2S])

def get_dict_from_meta (fname):
    df = pd.read_excel(fname, engine='openpyxl')
    d = {}
    for i in range(df.shape[0]):
        sample = df.iloc[i, 0]
        group = df.iloc[i, 1]
        if group in d:
            d[group].append(sample)
        else:
            d[group] = []
            d[group].append(sample)
    return (d)

def swap_list_element(l, d):
    '''internal'''
    return [swap_list_element(x, d) if isinstance(x, list) else d.get(x, x) for x in l]

def flatten(l):
    '''internal'''
    flat_list = [item for sublist in l for item in sublist]
    return (flat_list)

def collapseL2(L3):
    '''internal'''
    L2 = []
    for l in L3:
        L2.append(flatten(l))
    return (L2)

def s2b(S):
    '''internal'''
    K = flatten(S)
    V = ["sorted_reads/"+x+".bam" for x in K]
    return (dict(zip(K, V)))

def G2B_workflow(G, g2s):
    S = swap_list_element(G, g2s)
    S = collapseL2(S)
    s2b_dict = s2b(S)
    BS = swap_list_element(S, s2b_dict)
    B = [",".join(x) for x in BS]
    return (BS, B)

###  Input functions (functions to create input fnames for rules)  ###
def DESeq2_input(config):
    if config["START"] in ["FASTQ", "BAM"]:
        if config["MODE"] == "strict":
            if config["TE_ANALYSIS"] and config["START"] == "FASTQ":
                return("feature_count/TE_included.txt")
            else:
                return("feature_count/counts.strict.txt")
        elif config["MODE"] == "liberal":
            if config["TE_ANALYSIS"] and config["START"] == "FASTQ":
                return("feature_count/TE_included.txt")
            else:
                return("feature_count/counts.liberal.txt")
        else:
            raise Exception("START config not valid")
    elif config["START"] in ["COUNT"]:
            return(config["COUNT_FILE"])
    else:
        return("placeholder" )

def input_rnk_fname1(wildcards, config):
    '''
    fname1: the rnk file name as input for rule GSEA , 
    which is the output of upstream analysis
    '''
    if config['START'] == 'RNK':
        fname1 = 'meta/' + wildcards['fname']  
        # wildcards['fname']  = one of config['RNKS'], 
        # e.g. meta/wrong.rnk.xlsx, meta/name1.rnk.txt
    else:
        fname1 = "DESeq2/rnk/" + wildcards['fname'] + ".rnk"  
        # e.g. DESeq2/rnk/KO_D0.KO_D2.KO_D8_vs_WT_D0.WT_D2.WT_D8.rnk
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
        return re.sub('.rnk.xlsx$', '.rnk.txt', fname1)
    elif fname1.endswith('.rnk'):
        return fname1+".txt"
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
    fname1=input_rnk_fname1(wildcards, config)
    return rnk_fname1_to_fname2(fname1)


### Get parameters ###
RMATS_STRANDNESS= {0 : 'fr-unstranded', 1 : 'fr-firststrand',2 : 'fr-secondstrand'}
def get_strandness (strandFile="meta/strandness.detected.txt", config="config_dict"):
    try:
        with open(strandFile, "r") as file:
            txt = file.readline()
        p = re.compile("counts.s(.)\.[liberal|strict]")
        res = p.search(txt)
        if res.group(1) is None:
            sys.exit("strandness detection wrong")
        strand = (RMATS_STRANDNESS[int(res.group(1))])
        return (strand)
    except FileNotFoundError:
        sys.stderr.write(strandFile + "will be found in real run, not in dry run\n")
        return (None)

def get_strandness_for_dexseq (strandFile="meta/strandness.detected.txt"):
    book = {0 : '-s no', 1 : '-s yes', 2 : '-s reverse'}
    try:
        with open(strandFile, "r") as file:
            txt = file.readline()
        p = re.compile("counts.s(.)\.[liberal|strict]")
        res = p.search(txt)
        if res.group(1) is None:
            sys.exit("strandness detection wrong")
        strand = book[int(res.group(1))]
        return (strand)
    except FileNotFoundError:
        sys.stderr.write(strandFile + "will be found in real run, not in dry run\n")
        return (None)

def get_strandness_for_hisat2_PE (strandFile="meta/strandness.detected.txt"):
    book = {0 : ' ', 1 : '--rna-strandness FR', 2 : '--rna-strandness RF'}
    try:
        with open(strandFile, "r") as file:
            txt = file.readline()
        p = re.compile("counts.s(.)\.[liberal|strict]")
        res = p.search(txt)
        if res.group(1) is None:
            sys.exit("strandness detection wrong")
        strand = book[int(res.group(1))]
        return (strand)
    except FileNotFoundError:
        sys.stderr.write(strandFile + "will be found in real run, not in dry run\n")
        return (None)

def get_strandness_for_hisat2_SE (strandFile="meta/strandness.detected.txt"):
    book = {0 : ' ', 1 : '--rna-strandness F', 2 : '--rna-strandness R'}
    try:
        with open(strandFile, "r") as file:
            txt = file.readline()
        p = re.compile("counts.s(.)\.[liberal|strict]")
        res = p.search(txt)
        if res.group(1) is None:
            sys.exit("strandness detection wrong")
        strand = book[int(res.group(1))]
        return (strand)
    except FileNotFoundError:
        sys.stderr.write(strandFile + "will be found in real run, not in dry run\n")
        return (None)

def get_strandness_for_stringtie (strandFile="meta/strandness.detected.txt"):
    book = {0 : ' ', 1 : '--fr', 2 : '--rf'}
    try:
        with open(strandFile, "r") as file:
            txt = file.readline()
        p = re.compile("counts.s(.)\.[liberal|strict]")
        res = p.search(txt)
        if res.group(1) is None:
            sys.exit("strandness detection wrong")
        strand = book[int(res.group(1))]
        return (strand)
    except FileNotFoundError:
        sys.stderr.write("meta/strandness.detected.txt will be found in real run, not in dry run\n")
        return (None)

def read_length(lengthFile="meta/read_length.txt"):
    try:
        with open(lengthFile, 'r') as file:
            txt = file.readline()
        return (int(txt))
    except FileNotFoundError:
        sys.stderr.write(lengthFile + "not found in dry run, will be found in real run\n")
        return (None)
 
