import os
import sys

output = open(snakemake.output[0], "w")


def count_of_small_files(filenames):
    """
    Input: list of filenames
    Return: count of small files less than 100 in size (int)
    """
    N = 0
    for filename in filenames:
        print(filename, os.path.getsize(filename), file=output)
        if os.path.getsize(filename) < 100:
            sys.stderr.write("file too small:" + os.path.getsize(filename) + "\n")
            N += 1
    return N


file_list = snakemake.input

N = count_of_small_files(file_list)

if N < 1:
    print("All files are big enough", file=output)
else:
    outstring = "{} files too small, workflow aborted!!!".format(N)
    sys.stderr.write(outstring + "\n\n")
    print(outstring, file=output)
    sys.exit(outstring)

output.close()
