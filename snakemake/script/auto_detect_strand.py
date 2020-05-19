import pandas as pd 
import shutil
nums = []
names = []
for f in snakemake.input:
    print(f)
    df = pd.read_table(f, index_col=0)
    num =  sum(df.loc["Assigned", :])
    nums.append(num)
    names.append(f)

idx = nums.index(max(nums)) # assuming: if stranded, correctly stranded assignment > unstranded assignment, because of ambiguity issue
name = names[idx]
shutil.copy(name.replace(".summary","") , snakemake.output[0])