"""
Keep only unique columns of a CSV file, the goal is to remove potential redundancies.

Usage: python uniq_csv.py input.csv output.csv
    
"""

import csv
import sys

filename = sys.argv[1]
outcsv = sys.argv[2]

col_dict = {}

with open(filename, 'r') as infile:
    reader = csv.reader(infile)
    header = next(reader)
    for head in header:
        col_dict[head] = []

with open(filename, 'r') as infile:
    reader = csv.reader(infile)
    header = next(reader)
    for row in reader:
        if not len(row) == 0:
            for i in range(len(col_dict.keys())):
                col_dict[list(col_dict.keys())[i]].append(row[i])
                
new_dict = {}
for key, value in col_dict.items():
    if value not in new_dict.values():
        new_dict[key] = value

with open(outcsv, "w") as f:
    header = ",".join(new_dict.keys()) + "\n"
    f.write(header)
    res_row_1 = ",".join([new_dict[key][0] for key in new_dict]) + "\n"
    res_row_2 = ",".join([new_dict[key][1] for key in new_dict]) + "\n"
    f.write(res_row_1)
    f.write(res_row_2)