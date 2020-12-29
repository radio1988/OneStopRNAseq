import pandas as pd
import sys

print("usage: xlsx2txt.py name.xlsx")
print("output: name.xlsx.txt")

fname = sys.argv[1]
oname = fname + ".txt"
df = pd.read_excel(sys.argv[1])
print(df)
df.to_csv(oname, sep="\t", index=False, header=True)
