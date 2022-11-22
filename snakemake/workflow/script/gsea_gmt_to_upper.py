import argparse

parser = argparse.ArgumentParser(
    description = 'read gmt file, capitalize all gene names, does not touch gene set names, or source'
    )
parser.add_argument( '--fname', help = 'file name of the gmt file')
args = parser.parse_args()

F = open(args.fname, 'r')
Lines = F.readlines()
for line in Lines:
    line = line.rstrip()
    Eles = line.split("\t")
    Eles[2:] = [e.upper() for e in Eles[2:]]
    print("\t".join(Eles))

