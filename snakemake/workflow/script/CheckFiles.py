import os
import sys


def are_files_too_small(filenames):
    """
    if one file or more is too small, return False
    """
    Flag = True
    for filename in filenames:
        print(filename, os.path.getsize(filename), file=output)
        if os.path.getsize(filename) < 100:
            Flag = False
    return Flag


output = open(snakemake.output[0], "w")

file_list = snakemake.input
print(file_list, file=output)
if are_files_too_small(file_list):
    print("All trimmed.fastq.gz files are big enough", file=output)
else:
    print("At least one trimmed fastq.gz file is too small", file=output)
    sys.exit("At least one trimmed fastq.gz file is too small")

output.close()
