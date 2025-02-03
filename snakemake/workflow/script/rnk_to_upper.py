import sys
import pandas as pd
import re

print(
    "function: make all lower case letters in txt file upper case, remove all spaces/quotes in gene symbols, then save to txt file for GSEA")
print("usage: rnk_to_upper.py file.xlsx/file.txt/file.rnk")
print("path/file.txt or path/file.xlsx will be converted to path/file.rnk.txt")

fname = sys.argv[1]

over_write = False

if fname.endswith(".xlsx"):
    df = pd.read_excel(fname)
    outname = re.sub('.xlsx$', '.txt', fname)
    over_write = True
elif fname.endswith(".txt"):
    df = pd.read_table(fname, header=0)
    outname = fname
    over_write = False
elif fname.endswith(".rnk"):
    df = pd.read_table(fname, header=0)
    outname = fname + ".txt"
    over_write = True
else:
    sys.exit("Error: File format not .xlsx, nor .txt, nor .rnk, Exit\n")

print("input:\n", fname, "\n", df.head())


def containslower(string):
    if string:
        string = str(string)
        if any(char.islower() for char in string):
            return (True)
        else:
            return (False)
    else:
        return (False)


# header must contain #
if not df.columns[0].startswith("#"):
    df.rename(columns={df.columns[0]: "# " + df.columns[0].strip()}, inplace=True)
    over_write = True

# lower gene name to upper
if any(containslower(name) for name in df.iloc[:, 0]):
    df = df.apply(lambda x: x.astype(str).str.upper())
    over_write = True

# remove spaces, quotes from gene symbols
if any(df.iloc[:, 0].str.find(" ") >= 0):
    df = df.apply(lambda x: x.astype(str).str.replace(" ", ""))
    over_write = True

if any(df.iloc[:, 0].str.find("\"") >= 0):
    df = df.apply(lambda x: x.astype(str).str.replace("\"", ""))
    over_write = True

if any(df.iloc[:, 0].str.find("\'") >= 0):
    df = df.apply(lambda x: x.astype(str).str.replace("\'", ""))
    over_write = True

# output
if over_write:
    df.to_csv(outname, index=False, header=True, sep="\t")
    print("fixed spaces, lower cases, quotes and saved to ", outname)
    print("output:\n", outname, "\n", df.head())
else:
    print("no lower case nor spaces found in file, not updating", fname)
