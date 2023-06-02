# usage: 
# python fix_gtf_scaffold_names.py input.gtf input.genome.fasta
# - difference in scaffold names will be found
# - best match in genome based names will be found for each gtf unique names
# - output called modified.gtf will be created
import sys
import difflib
print(sys.argv)

# read scaffold names
file_gtf = sys.argv[1]
chr_names_gtf = set()
with open(file_gtf, 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue
        words = line.split("\t")
        chr_names_gtf.add(words[0])

file_genome = sys.argv[2]
chr_names_genome = set()
with open(file_genome, 'r') as file:
    for line in file:
        if line.startswith('>'):
            line = line.strip()
            word = line.replace('>','')
            chr_names_genome.add(word)

# find diff
to_fix =  chr_names_gtf - chr_names_genome
library = chr_names_genome - chr_names_gtf

best_matches = {}
best_match = ""
best_ratio = 0
for name1 in to_fix:
    for name2 in library:
        similarity_ratio = difflib.SequenceMatcher(None, name1, name2).ratio()   
        if similarity_ratio > best_ratio:
            best_ratio = similarity_ratio
            best_match = (name1, name2)
    best_matches[name1] = name2

print('Best Matches:', best_matches)

# fixing diff in gtf (leave genome untouched)
with open(file_gtf, 'r') as file, open('modified.gtf', 'w') as outfile:
    for line in file:
        if line.startswith('#'):
            outfile.write(line)
        line0 = str(line)
        for a,b in best_matches.items():
            line = line.replace(a,b)

        if line0 != line:
            print(line0, line)
        outfile.write(line)
