import sys
import pandas as pd

fname = sys.argv[1]
max_pvalue = float(sys.argv[2])
outname = sys.argv[3]

df = pd.read_table(fname)
out = df[df['PValue'] < max_pvalue]
out.to_csv(outname, sep="\t") 
