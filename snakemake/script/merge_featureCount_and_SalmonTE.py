import pandas as pd
import sys
import os


print("usage: python merge_featureCount_and_SalmonTE.py gene_count.txt TE_EXPR.csv output.csv")
print("example: python script/merge_featureCount_and_SalmonTE.py \
feature_count/counts.strict.txt  SalmonTE_output/EXPR.csv  feature_count/TE_included.txt")
# input format rigit


def change_featureCount_colname(s):
    s = s.replace("mapped_reads/", "")
    s = s.replace(".bam", "")
    return (s)
    
def find_sf(path="../SalmonTE_output/"):
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(".sf"):
                 return(os.path.join(root, file))
                 
                 
# read featureCounts
genes = pd.read_csv(sys.argv[1], sep="\t", comment="#")
genes = genes.rename(columns=change_featureCount_colname)
genes_anno = genes.iloc[:, 0:6]
genes_num = genes.iloc[:,6:]
genes_num.sort_index(axis=1, inplace=True)


# read SalmonTE
te = pd.read_csv(sys.argv[2])
te_anno = te.iloc[:, 0:1]
te_num = te.iloc[:, 1:]
te_num.sort_index(axis=1, inplace=True)

empty = pd.DataFrame(index=range(0, te_num.shape[0]), columns=["Chr", "Start", "End", "Strand", "Length"])
sf_file = find_sf(os.path.dirname(sys.argv[2]))
sf = pd.read_csv(sf_file, sep="\t")
empty.loc[:, "Length"] = sf.loc[:, "Length"]
te = pd.concat([te_anno, empty, te_num], axis=1)
print("SalmonTE: ", te.columns)
print("featureCounts: ", genes.columns)
te.columns = genes.columns

# Merge and Output
merged = genes.append(te, ignore_index=True)
print("merged.shape", merged.shape)
# print(merged.head())
merged.to_csv(sys.argv[3], index=False, sep="\t", na_rep="NA")
