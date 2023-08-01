import sys
import re

gene_ids = set()
pattern = r'gene_id\s+"([^"]+)";'

with open(sys.argv[1],'r') as gtf:
    for gtf_line in gtf:
        match = re.search(pattern, gtf_line)
        if match:
            gene_ids.add(match.group(1))
print(gene_ids)
print(len(gene_ids))


