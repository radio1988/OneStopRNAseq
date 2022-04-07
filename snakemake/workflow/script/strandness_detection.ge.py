import pandas as pd 
import shutil

nums = []
names = []
for f in snakemake.input:
    print('reading:', f)
    df = pd.read_table(f, index_col=0)
    num =  sum(df.loc["Assigned", :])
    nums.append(num)
    names.append(f)

idx = nums.index(max(nums)) # assuming: if stranded, correctly stranded assignment > unstranded assignment, because of ambiguity issue
name = names[idx]
print("the correct strand is found in:", name)

# copy for DESeq2 input
shutil.copy(name.replace(".summary","") , snakemake.output[0])
print(snakemake.output[0], "created")

# write for later use 
with open("feature_count_gene_level/strandness.detected.txt","w") as output:
	output.write(name)
print("feature_count_gene_level/strandness.detected.txt created" )
