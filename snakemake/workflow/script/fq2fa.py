import argparse

argParser = argparse.ArgumentParser()
argParser.add_argument("-i", "--inputName", help="input fastq file name")
args = argParser.parse_args()

with open(args.inputName, 'r') as f:
    for i,l in enumerate(f.readlines()):
        l = l.rstrip()
        if i % 4 == 0:
            print(">" + l)
        if i % 4 == 1:
            print(l + "\n")
