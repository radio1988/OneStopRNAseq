# usage: 
# python fix_gtf_scaffold_names.py input.gtf name_mapping.txt
# name_mapping.txt: `new` `old`, in tsv format
# - output called modified.gtf will be created
import sys
import difflib
import json
print(sys.argv)

file_gtf = sys.argv[1]
file_mapping = sys.argv[2]

# read scaffold names
chr_names_gtf = set()
with open(file_gtf, 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue
        words = line.split("\t")
        chr_names_gtf.add(words[0])

# read mapping file
mapping = {}
with open(file_mapping, 'r') as file:
    for line in file:
        line = line.strip()
        if line.startswith('#'):
            continue
        words = line.split("\t")
        mapping[words[1]] = words[0] # format: new old

print("\n>>>Best Matches (filtered):\n")
for k,v in mapping.items():
    print(k,v)

with open('mapping.json', 'w') as file:
    json_str = json.dumps(mapping)
    file.write(json_str)

# fixing diff in gtf (leave genome untouched)
n = 0
with open(file_gtf, 'r') as file, open('modified.gtf', 'w') as outfile:
    for line in file:
        line0 = str(line)
        for a,b in mapping.items():
            line = line.replace(a,b)
        if line0 != line:
            n += 1
        outfile.write(line)

print(n, 'lines changed')
