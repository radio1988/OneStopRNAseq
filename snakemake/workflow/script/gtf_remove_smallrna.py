import sys
print('usage: xxx.py xxx.gtf smallRnaNames.txt')
print('output: nosmallrna.gtf')
fname = sys.argv[1]
dbname = sys.argv[2]
print('gtf file:', fname)
print('smallRNA file:', dbname)

with open(dbname, 'r') as db:
    lines = db.read()
    smallrnas = set(lines.splitlines())

output = open('nosmallrna.gtf', 'w') 
output2 = open('smallrna.gtf', 'w')
skipped, kept = 0, 0
with open(fname, 'r') as gtf:
    for l in gtf:
        if any(item in l for item in smallrnas):
            skipped += 1
            output2.write(l)
            continue
        kept += 1
        output.write(l)
        
print('skipped:', skipped, 'kept:', kept)            
