import pandas as pd 
import sys


# tab0 is exon_level, tab1 is gene_level
# usage: python get_max_of_featureCountsTable.py exon.featureCounts.txt gene_level.featureCounts.txt
# function: get max(exon-count, gene-count) for all counts in gene_level.featureCounts.txt (column 7:)
# then overwrite gene_level.featureCounts.txt

# file0="feature_count/counts.s0.strict.txt"
# file1="feature_count_gene_level/counts.s0.gene_level.strict.txt"

file0=sys.argv[0]
file1=sys.argv[1]

tab0 = pd.read_table(file0, comment='#')
tab1 = pd.read_table(file1, comment='#')

larger_exon_counts_idx = tab0.iloc[:, 7:].gt(tab1.iloc[:, 7:])

#tab1.iloc[:, 7:].where(~larger_exon_counts_idx, "bingo")
tab1.iloc[:, 7:] = tab1.iloc[:, 7:].where(~larger_exon_counts_idx, tab0.iloc[:, 7:])

print("replacing gene_count with max(gene_count, exon_count)")
tab1.to_csv(file1, sep="\t", index=False)
#larger_exon_counts_idx.to_csv("temp.csv")