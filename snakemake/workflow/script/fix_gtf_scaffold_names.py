# usage: 
# python fix_gtf_scaffold_names.py input.gtf input.genome.fasta min_similarity
# - difference in scaffold names will be found
# - best match in genome based names will be found for each gtf unique names
# - filter if len(match)/len(gtf) > min_similarity
# - output called modified.gtf will be created
# e.g. python fix_gtf_scaffold_names.py gencode.v34.primary_assembly.annotation.gtf hg38.primary.fa 0.8 > fix_gtf_scaffold_names.py.log
import sys
import difflib
import json
print(sys.argv)

file_gtf = sys.argv[1]
file_genome = sys.argv[2]
min_similarity = float(sys.argv[3]) # e.g. 8 characters in GTF name matches in 10, 0.8

# read scaffold names
chr_names_gtf = set()
with open(file_gtf, 'r') as file:
    for line in file:
        if line.startswith('#'):
            continue
        words = line.split("\t")
        chr_names_gtf.add(words[0])

chr_names_genome = set()
with open(file_genome, 'r') as file:
    for line in file:
        if line.startswith('>'):
            line = line.strip()
            word = line.replace('>','')
            chr_names_genome.add(word)

# find diff
to_fix =  chr_names_gtf - chr_names_genome
with open('to_fix.txt', 'w') as output:
    for key in to_fix:
        output.write(key+"\n")

library = chr_names_genome - chr_names_gtf
with open('library.txt', 'w') as output:
    for key in library:
        output.write(key+"\n")

best_matches = {}
best_scores = []
for name1 in to_fix:
    best_match = ""
    best_ratio = 0
    for name2 in library:
        #name1_simp = name1.split('.')[0]
        match_size = difflib.SequenceMatcher(None, name1, name2).find_longest_match().size
        similarity_ratio = match_size / len(name1)
        if similarity_ratio > best_ratio:
            best_ratio = similarity_ratio
            best_match = name2
    if best_ratio >= min_similarity:
        best_matches[name1] = best_match
        library.remove(best_match)
    best_scores.append([name1, best_match, best_ratio])

print("\n>>>Best Matches (filtered):\n")
for k,v in best_matches.items():
    print(k,v)

print("\n>>>Best Similarities (not filtered):\n")
for s in best_scores:
    print(s)

with open('mapping.json', 'w') as file:
    json_str = json.dumps(best_matches)
    file.write(json_str)

# fixing diff in gtf (leave genome untouched)
with open(file_gtf, 'r') as file, open('modified.gtf', 'w') as outfile:
    for line in file:
        line0 = str(line)
        for a,b in best_matches.items():
            line = line.replace(a,b)
        if line0 != line:
            print(line0, line)
        outfile.write(line)
