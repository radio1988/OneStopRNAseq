"""
Read gmt file and a list of gene sets to keep,
and write a new gmt file with only the gene sets to keep.

e.g.
 python gmt_filter.py -gmt m5.go.bp.v2024.1.Mm.symbols.gmt \
 -gene_sets genesets.txt -out m5.go.bp.subset.gmt
"""


import sys
import os
import argparse
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-gmt", help="gmt file to filter")
    parser.add_argument("-gene_sets", help="file with gene sets to keep")
    parser.add_argument("-out", help="output gmt file")
    return parser.parse_args()

def read_gene_sets(gene_sets_file):
    with open(gene_sets_file) as f:
        gene_sets = set([line.strip() for line in f])
    return gene_sets

def filter_gmt(gmt_file, gene_sets, out_file):
    with open(gmt_file) as f:
        with open(out_file, "w") as out:
            for line in f:
                gene_set = line.split("\t")[0]
                if gene_set in gene_sets:
                    out.write(line)

def main():
    args = parse_args()
    gene_sets = read_gene_sets(args.gene_sets)
    filter_gmt(args.gmt, gene_sets, args.out)

if __name__ == "__main__":
    main()




