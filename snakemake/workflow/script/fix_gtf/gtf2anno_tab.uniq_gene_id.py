import sys
import re

file = sys.argv[1]
outfile = file + '.txt'
print("Input: {}\nOutput: {}".format(file,outfile))


id2name = {}


with open(file) as f:
    while True:
        line = f.readline()
        if not line:
            break
        
        x = line.split("\t")

        if len(x) != 9:
            continue

        seqname,source,feature,start,end,score,strand,fname,attribute=x

#        if feature != "gene":
#            continue
        
        gene_name = re.search(r'gene_name "([^\"]+)"', attribute).groups()[0]
        gene_id = re.search(r'gene_id "([^\"]+)"', attribute).groups()[0]
        id2name[gene_id] = gene_name


o = open (outfile,'w')
for k,v in id2name.items():
    outline = "\t".join([k,v,'NA'])+"\n"
    o.write(outline)
o.close()

