import pandas as pd
import sys

print("usage: xlsx2txt.py name.txt")
print("output: name.txt.xlsx")

fname = sys.argv[1]
oname = fname + ".xlsx"
df = pd.read_table(fname)
print(df)
df.to_excel(oname, index=False)
