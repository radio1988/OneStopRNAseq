import sys
import re

file = sys.argv[1]
outfile = file + '.txt'
print("Input: {}\nOutput: {}".format(file,outfile))


o = open (outfile,'w')
o.write("\t".join(['Gene','Name','Type','Chr','Strand','Start','End'])+"\n")

with open(file) as f:
    while True:
        line = f.readline()
        if not line:
            break
        
        x = line.split("\t")

        if len(x) != 9:
            continue

        seqname,source,feature,start,end,score,strand,fname,attribute=x

        if feature != "gene":
            continue
    
        gene_id = re.search(r'gene_id "([^\"]+)"', attribute).groups()[0]
        if re.search(r'gene_name "([^\"]+)"', attribute): 
            gene_name = re.search(r'gene_name "([^\"]+)"', attribute).groups()[0]
        else:
            gene_name = gene_id
        gene_type = re.search(r'gene_biotype "([^\"]+)"', attribute).groups()[0]

        outline = "\t".join([gene_id,gene_name,gene_type,seqname,strand,start,end])+"\n"
        o.write(outline)

o.close()

