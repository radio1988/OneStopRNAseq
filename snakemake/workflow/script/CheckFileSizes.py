import os, sys

output = open(snakemake.output[0], "w")
log = open(snakemake.log[0], "w")

filenames = snakemake.input

N = 0
for filename in filenames:
    message = filename + str(os.path.getsize(filename))
    print(message, file=output)
    if os.path.getsize(filename) < 100:
        print("file too small: " + message, file=log)
        N += 1

if N < 1:
    print("All files are big enough", file=output)
else:
    outstring = "{} files too small, workflow aborted!!!".format(N)
    print(outstring, file=output)
    print(outstring, file=log)
    sys.exit(outstring)

output.close()
log.close()
