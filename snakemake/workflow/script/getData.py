"""
Given SRX, sampleLabel, and fastq folder path, perform srr download, md5 check, fastq-dump and renaming.
Usage:
python getData.py SRX_ID SampleLabel FastqFolder pair|single
"""

import sys
import time
import os
import subprocess


srx = sys.argv[1]
sampleName = sys.argv[2]
fastq_folder = sys.argv[3]
seqType      = sys.argv[4] # can be either "pair" or "single"

# Step1: retrieve SRX associating SRR:
jobCheckFile = fastq_folder + "/job_efecth_" + srx + ".txt"
efetchResFile = fastq_folder + "/res_efetch_" + srx + ".txt"
cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=500] -o " + jobCheckFile + " \"esearch -db sra -query " + srx + " | efetch -format runinfo | awk 'NR%2==0' | cut -d ',' -f 1 > " + efetchResFile + "\""
subprocess.call(cmd, shell=True)

efetchSuccess = 0
timer = 0
while (efetchSuccess == 0):
    jobCheckComplete = 0
    try:
        for line in open(jobCheckFile):
            if "Resource usage summary" in line:
                jobCheckComplete = 1
                break
    except:
        continue
    if jobCheckComplete:
        tooMany = 0
        for line in open(jobCheckFile):
            if "429 Too Many Requests" in line:
                tooMany = 1
                os.remove(jobCheckFile)
                subprocess.call(cmd, shell=True)
                break

        if tooMany:
            efetchSuccess = 0
        else:
            efetchSuccess = 1

        tem_1 = 0
        for line in open(efetchResFile):
            if len(line) > 5:
                tem_1 = len(line)

        if tem_1 < 5:
            try:
                os.remove(jobCheckFile)
            except:
                continue
            subprocess.call(cmd, shell=True)

        if tem_1 < 5:
            efetchSuccess = 0
        else:
            efetchSuccess = 1

    timer = timer + 10
    time.sleep(10)
    if timer > 10000:
        print("ERROR: " + srr + " efetch bsub failed!")
        exit()

srrList = [line.strip() for line in open(efetchResFile)]
# srrList = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True).decode("utf-8").split("\n")
srrList = [x for x in srrList if not x == ""]

errFile = fastq_folder + "/getDataPyInfo_" + srx + ".txt"

# testFile = fastq_folder + "/test_file.txt"
# # f = open(testFile, "a")
# # f.write("first" + srx + "\n")
# # f.close()

print("Downloading SRR records with prefetch for " + srx + " ...")
for srr in srrList:
    print("Processing " + srr + " ...")

    # f = open(testFile, "a")
    # f.write("now: " + srr + "\n")
    # f.close()


    jobCheckFile = fastq_folder + "/job_prefetch_" + srx + "_" + srr + ".txt"
    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=500] -o " + jobCheckFile + " \"prefetch -O " + fastq_folder + "/ " + srr + "\""
    # Can't exceeed -W 4:00 for short jobs.

    # f = open(testFile, "a")
    # f.write("before cmd" + "\n")
    # f.close()

    subprocess.call(cmd, shell=True)

    # f = open(testFile, "a")
    # f.write(cmd + "\n")
    # f.write("after cmd\n")
    # f.close()

    prefetchSuccess = 0
    timer = 0
    nthSubmit = 1
    while (prefetchSuccess == 0):
        jobCheckComplete = 0
        try:
            for line in open(jobCheckFile):
                if "Resource usage summary" in line:
                    jobCheckComplete = 1
                    break
        except:
            continue
        if jobCheckComplete:
            tooMany = 0
            for line in open(jobCheckFile):
                if "429 Too Many Requests" in line:
                    tooMany = 1
                    break
            if tooMany:
                f = open(errFile, "a")
                f.write("WARNING: too many requests for " + srx + " " + srr + ". Will re-bsub prefetch job, round: " + nthSubmit + "\n")
                f.close()
                nthSubmit = nthSubmit + 1
                os.remove(jobCheckFile)
                subprocess.call(cmd, shell=True)
            else:
                prefetchSuccess = 1
                f = open(errFile, "a")
                f.write("OKAY: prefetch job for " + srx + " " + srr + " is successful!\n")
                f.close()

        timer = timer + 30
        time.sleep(30)
        if timer > 10000:
            print("ERROR: " + srx + " " + srr + " prefetch bsub failed!")
            exit()

    # cmd = "cd " + fastq_folder + " && bsub -q short -n 1 -W 4:00 -R rusage[mem=1000] prefetch -O ./ " + srr
    fileName = fastq_folder + "/" + srr + ".sra"

    if not os.path.isfile(fileName):
        sleep    = 0
        download = 0
    else:
        print(srr + " downloaded!")
        download = 1

    # if not os.path.isfile(fileName):
    #     subprocess.call(cmd, shell=True)
    #     sleep    = 0
    #     download = 0
    # else:
    #     print(srr + " downloaded!")
    #     download = 1

    # check to see if downloading is complete or not:
    while (download == 0):
        print("  Sleeping " + str(sleep) + " till next check ...")
        if os.path.isfile(fileName):
            download = 1
            print(srr + " downloaded!")
        else:
            time.sleep(10)
            sleep = sleep + 10
        if sleep > 20000:
            print("ERROR: " + srr + " download failed!")
            f = open(errFile, "a")
            f.write("ERROR: " + srr + " download failed!\n")
            f.close()
            exit()

    # Check md5 with vdb-validate command:
    cmd = "cd " + fastq_folder + " && bsub -q short -n 1 -W 4:00 -R rusage[mem=1000] -o " + fastq_folder + "/" + srr + ".md5.log.txt" + " vdb-validate " + srr + ".sra"
    filename = fastq_folder + "/" + srr + ".md5.log.txt"
    # Check md5 check result:
    try:
        subprocess.call("rm " + filename, shell=True)
    except:
        continue
    subprocess.call(cmd, shell=True)

    md5_check = 0
    sleep     = 0
    while (md5_check == 0):
        print("Checking md5 integrity for " + srr + " ...")
        print(" Sleeping " + str(sleep) + " till next check ...")
        if os.path.isfile(filename):
            md5Complete = 0
            for line in open(filename):
                if ("The output (if any) is above this job summary." in line):
                    md5Complete = 1
            if md5Complete:
                md5_check = 1
                finished   = 0
                for line in open(filename):
                    if "is consistent" in line:
                        finished = 1
                if finished:
                    print("Done: md5check passed for " + srr)
                    break
                else:
                    print("ERROR: md5check failed for " + srr)
                    f = open(errFile, "a")
                    f.write("ERROR: md5check failed for " + srr + "!\n")
                    f.close()
                    exit()
        time.sleep(10)
        if sleep > 1000:
            print("ERROR: md5check not executed for " + srr + "!")
            f = open(errFile, "a")
            f.write("ERROR: md5check not executed for " + srr + "!\n")
            f.close()
            exit()
        sleep = sleep + 10

    # fastq-dump srr into fastq.gz files:
    # cmd = "cd " + fastq_folder + " && bsub -q short -n 1 -W 4:00 -R rusage[mem=20480] -o " + fastq_folder + "/" + srr + ".fastq_dump.log.txt" + " parallel-fastq-dump -s " + srr + ".sra -t 1 -O ./ --tmpdir ./ --split-files --gzip && rm " + srr + ".sra"
    cmd = "cd " + fastq_folder + " && bsub -q short -n 1 -W 4:00 -R rusage[mem=1000] -o " + fastq_folder + "/" + srr + ".fastq_dump.log.txt" + " parallel-fastq-dump -s " + srr + ".sra -t 1 -O ./ --tmpdir ./ --split-files --gzip"

    filename = fastq_folder + "/" + srr + ".fastq_dump.log.txt"
    try:
        subprocess.call("rm " + filename, shell=True)
    except:
        continue
    subprocess.call(cmd, shell=True)
    # print(cmd)

    dump_check = 0
    sleep     = 0
    while (dump_check == 0):
        print("Checking fastq_dump status for " + srr + " ...")
        print(" Sleeping " + str(sleep) + " till next check ...")
        if os.path.isfile(filename):
            fileComplete = 0
            for line in open(filename):
                if ("The output (if any) is above this job summary." in line):
                    fileComplete = 1
            if fileComplete:
                dump_check = 1
                finished   = 0
                for line in open(filename):
                    if "Successfully completed." in line:
                        finished = 1
                if finished:
                    print("Fastq.gz dumped for " + srr)
                    break
                else:
                    print("ERROR: fastq-dump failed for " + srr)
                    f = open(errFile, "a")
                    f.write("ERROR: fastq-dump failed for " + srr + "!\n")
                    f.close()
                    exit()
        time.sleep(30)
        if sleep > 10000:
            print("ERROR: fastq-dump not executed for " + srr + "!")
            f = open(errFile, "a")
            f.write("ERROR: fastq-dump not executed for " + srr + "!\n")
            f.close()
            exit()
        sleep = sleep + 30

# concatenate/rename fastq.gz files

## determine SE or PE:
# test = srrList[0]
# testR2 = fastq_folder + "/" + test + ".*2.fastq.gz"
# print(testR2)

# if PE:
# if os.path.isfile(testR2):
if (seqType == "pair"):
    print("Renaming PE fastq.gz files ...")
    temR1 = [item + "*1.fastq.gz" for item in srrList]
    tem_filesR1 = " ".join(temR1)
    cmdR1 = "cat " + tem_filesR1 + " > " + sampleName + ".R1.fastq.gz"

    temR2 = [item + "*2.fastq.gz" for item in srrList]
    tem_filesR2 = " ".join(temR2)
    cmdR2 = "cat " + tem_filesR2 + " > " + sampleName + ".R2.fastq.gz"

    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=20480] -o " + fastq_folder + "/" + srx + "_" + sampleName + ".R1.fastq.gz.log " + "\"" + "cd " + fastq_folder + " && " + cmdR1 + " && touch " + sampleName + ".R1.fastq.gz.done" + "\""
    subprocess.call(cmd, shell=True)
    print(cmd)

    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=20480] -o " + fastq_folder + "/" + srx + "_" + sampleName + ".R2.fastq.gz.log " + "\"" + "cd " + fastq_folder + " && " + cmdR2 + " && touch " + sampleName + ".R2.fastq.gz.done" + "\""
    subprocess.call(cmd, shell=True)
    print(cmd)
# if SE:
elif (seqType == "single"):
    print("Renaming SE fastq.gz files ...")
    temR1 = [item + "*1.fastq.gz" for item in srrList]
    tem_filesR1 = " ".join(temR1)
    cmdR1 = "cat " + tem_filesR1 + " > " + fastq_folder + "/" + sampleName + ".fastq.gz"

    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=20480] -o " + fastq_folder + "/" + srx + "_" + sampleName + ".fastq.gz.log " + "\"" + "cd " + fastq_folder + " && " + cmdR1 + " && touch " + sampleName + ".fastq.gz.done" + "\""

    print(cmd)
    subprocess.call(cmd, shell=True)

# clean intermediate data: check to see if the bsub cat finished or not first
sleep = 0
while (1):
    if seqType == "pair":
        if os.path.exists(fastq_folder + "/" + sampleName + ".R1.fastq.gz.done") and os.path.exists(fastq_folder + "/" + sampleName + ".R2.fastq.gz.done"):
            for srr in srrList:
                cmd = "cd " + fastq_folder + " && rm " + srr + ".sra && rm " + srr + "*.fastq.gz" # Note that if in bsub, the * must be \* like in temR1 =
                subprocess.call(cmd, shell=True)
            break
    elif seqType == "single":
        if os.path.exists(fastq_folder + "/" + sampleName + ".fastq.gz.done"):
            for srr in srrList:
                cmd = "cd " + fastq_folder + " && rm " + srr + ".sra && rm " + srr + "*.fastq.gz" # Note that if in bsub, the * must be \* like in temR1 =
                subprocess.call(cmd, shell=True)
            break
    sleep = sleep + 10
    time.sleep(10)
    if sleep > 1000:
        f = open(errFile, "a")
        f.write("WARNING: intemediate fastq.gz files not removed for " + srx + "!\n")
        f.close()
        exit()
        break

print("getData.py finished for " + srx + "!")
