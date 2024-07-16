import os, sys

output = open(snakemake.output[0], "w")
log = open(snakemake.log[0], "w")

filenames = snakemake.input

N = 0
for filename in filenames:
    message = filename + " " + str(os.path.getsize(filename)) + " bytes"
    # output
    print(message, file=output)
    # log
    if os.path.getsize(filename) < 100:
        print("file too small: " + message, file=log)
        N += 1
    else:
        print(message, file=log)

if N < 1:
    print("All files are big enough", file=output)
else:
    outstring = "{} files too small, aborted!!!".format(N)
    print(outstring, file=output)
    print(outstring, file=log)
    sys.exit(outstring)

output.close()
log.close()
