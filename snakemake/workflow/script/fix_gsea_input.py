import sys
import pandas as pd
import re

print("function: make all lower case letters in txt file upper case, remove all spaces/quotes in gene symbols, then save to txt file for GSEA")
print("usage: fix_gsea_input.py file.xlsx/file.txt/file.rnk")
print("path/file.txt or path/file.xlsx will be converted to path/file.rnk.txt")


fname = sys.argv[1]

over_write = False

if fname.endswith(".rnk.xlsx"):
    df = pd.read_excel(fname)
    outname = re.sub('.rnk.xlsx$', '.rnk.txt', fname)
    over_write = True
elif fname.endswith(".rnk.txt"):
    df = pd.read_table(fname, header=0)
    outname = fname
    over_write = False
elif fname.endswith(".rnk"):
    df = pd.read_table(fname, header=0)
    outname = fname+".txt"
    over_write = True
else:
    sys.exit("Error: File format not .rnk.xlsx, nor .rnk.txt, nor .rnk, Exit\n \
        Please note that '.rnk' suffix is necessary for xlsx and txt files, to indicate they are rnk files")

print("input:\n", df.head())

# lower gene name to upper
if any(df.iloc[:, 0].str.islower()):
    df = df.apply(lambda x: x.astype(str).str.upper())
    over_write = True

# remove space from gene symbols
if any(df.iloc[:, 0].str.find(" ") >= 0):
    df = df.apply(lambda x: x.astype(str).str.replace(" ", ""))
    over_write = True

if any(df.iloc[:, 0].str.find("\"") >= 0):
    df = df.apply(lambda x: x.astype(str).str.replace("\"", ""))
    over_write = True

if any(df.iloc[:, 0].str.find("\'") >= 0):
    df = df.apply(lambda x: x.astype(str).str.replace("\'", ""))
    over_write = True

# convert to upper 
if over_write:
    df.to_csv(outname, index=False, header=True, sep="\t")
    print("fixed spaces, lower cases, quotes and saved to ", outname)
    print("output:\n", df.head())
else:
    print("no lower case nor spaces found in file, not updating", fname)

