import os
import sys

output = open(snakemake.output[0], "w")
log = open(snakemake.log[0], "w")


def count_of_small_files(filenames):
    """
    Input: list of filenames
    Return: count of small files less than 100 in size (int)
    """
    N = 0
    for filename in filenames:
        m = filename + str(os.path.getsize(filename))
        print(m, file=output)
        if os.path.getsize(filename) < 100:
            print("file too small: " + m, file=log)
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
log.close()
