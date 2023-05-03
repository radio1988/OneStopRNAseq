import pandas as pd  
import sys
from osr import read_table
fname = sys.argv[1]
df=read_table(fname)
df.columns = ["sample.ID", "group.ID", "batch.ID"] 
df.to_csv("meta/decoder.txt", sep="\t", index=False)  
