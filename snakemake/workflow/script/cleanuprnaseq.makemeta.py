import sys
import pandas as pd
import re


def get_STRAND(strand_string):
    print(strand_string)
    p0 = re.compile("counts.s0\.[liberal|strict]")
    p1 = re.compile("counts.s1\.[liberal|strict]")
    p2 = re.compile("counts.s2\.[liberal|strict]")
    if re.search(p0, strand_string):
        STRAND = 'U'
    elif re.search(p1, strand_string):
        STRAND = 'F'
    elif re.search(p2, strand_string):
        STRAND = 'R'
    else:
        sys.exit("strand_string not recognized")
    return STRAND


log = open(snakemake.log[0], "w")
sys.stderr = log
sys.stdout = log

# get STRAND
with open(snakemake.input[1], "r") as file:  # meta/strandness.detected.txt
    strand_string = file.readline()
    strand_string = strand_string.rstrip("\n")
STRAND = get_STRAND(strand_string)

# output
fname = snakemake.input[0]  # config['META']
if fname.endswith('csv'):
    df = pd.read_csv(fname)
elif fname.endswith('xlsx'):
    df = pd.read_excel(fname)

df.columns = ['sample_name', 'group', 'batch']
df['BAM_file'] = 'mapped_reads/' + df['sample_name'] + '.bam'

if snakemake.config["PAIR_END"]:
    if STRAND == "U":
        LIBTYPE = "IU"
    elif STRAND == "F":
        LIBTYPE = "ISF"
        OPPO = "ISR"
    elif STRAND == "R":
        LIBTYPE = "ISR"
        OPPO = "ISF"
else:
    if STRAND == "U":
        LIBTYPE = "U"
    elif STRAND == "F":
        LIBTYPE = "SF"
        OPPO = "SR"
    elif STRAND == "R":
        LIBTYPE = "SR"
        OPPO = "SF"

df['salmon_quant_file'] = 'salmon/{}/'.format(LIBTYPE) + df['sample_name'] + '/quant.sf'

if STRAND == "F" or STRAND == "R":
    df['salmon_quant_file_strand'] = 'salmon/{}/'.format(LIBTYPE) + df['sample_name'] + '/quant.sf'
    df['salmon_quant_file_reverse_strand'] = 'salmon/{}/'.format(OPPO) + df['sample_name'] + '/quant.sf'
else:
    pass  # for U, no need for these two columns

df.to_csv(snakemake.output[0], index=False)

log.close()
