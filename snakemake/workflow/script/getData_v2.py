"""
Given SRX, sampleLabel, and fastq folder path, perform srr download, md5 check, fastq-dump and renaming.
Usage:
python getData.py SRX_ID SampleLabel FastqFolder pair|single
"""

import sys
import time
import os
import subprocess
from datetime import datetime


srx = sys.argv[1]
sampleName = sys.argv[2]
fastq_folder = sys.argv[3]
seqType      = sys.argv[4] # can be either "pair" or "single"

# Step1: retrieve SRX associating SRR:
efetchLogFile = fastq_folder + "/efetch_log_" + srx + ".txt"
efetchResFile = fastq_folder + "/efetch_res_" + srx + ".txt"
# errorLogFile  = fastq_folder + "/getData_" + srx + ".txt"
logFile       = fastq_folder + "/getData_log_" + srx + ".txt"

cmd = "bsub -q short -n 1 -W 0:10 -R rusage[mem=500] -o " + efetchLogFile + " \"esearch -db sra -query " + srx + " | efetch -format runinfo | awk 'NR%2==0' | cut -d ',' -f 1 > " + efetchResFile + "\""
with open(logFile, "a") as f:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    f.write(current_time + " CMD efetch: ")
    f.write(cmd)
    f.write("\n")
subprocess.call(cmd, shell=True)

efetchDone = 0
timer = 0

while (efetchDone == 0):
    efetchLogComplete = 0
    efetchError = 0
    resubmit = 0

    if os.path.exists(efetchLogFile):
        for line in open(efetchLogFile):
            if "Resource usage summary" in line:
                efetchLogComplete = 1
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " efetch " + efetchLogFile + " complete.")
                    f.write("\n")
                # since we also check the res file, so the below check of the log file is redundant
                for line in open(efetchLogFile):
                    if "429 Too Many Requests" in line:
                        efetchError = 1
                        with open(logFile, "a") as f:
                            now = datetime.now()
                            current_time = now.strftime("%H:%M:%S")
                            f.write(current_time + " efetch " + efetchLogFile + " contains \"429 Too Many Requests\", need to resubmit job.")
                            resubmit = resubmit + 1
                            f.write("\n")
                        break

    if efetchLogComplete and not efetchError:
        if os.path.exists(efetchResFile):
            time.sleep(1) # make sure efetchResFile finishes saving to disk
            if not os.stat(efetchResFile).st_size: # if file empty:
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " efetch " + efetchResFile + " is done but is empty, need to resubmit job.")
                    f.write("\n")
                    resubmit = resubmit + 10
            else: # if not empty:
                efetchDone = 1
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " efetch " + efetchResFile + " is done and not empty. Move on.")
                    f.write("\n")

    if resubmit: # only resubmit once if job fails
        os.remove(efetchLogFile)
        os.remove(efetchResFile)
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " efetch resubmit code is: " + str(resubmit) + ", now resubmitting efetch job.")
            f.write("\n")
            subprocess.call(cmd, shell=True)

    timer = timer + 10
    time.sleep(10)

    if timer // 3000: # resubmit the job every 3000 seconds in case of the ssh jobs queue error.
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " efetch job reach 3000-second limit, resubmit job to hpc.")
            f.write("\n")
            subprocess.call(cmd, shell=True)

    if timer > 10000:
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " ERROR: efetch job time out.")
            f.write("\n")
            exit()

srrList = [line.strip() for line in open(efetchResFile)]
srrList = [x for x in srrList if not x == ""]

srrFile = " ".join(srrList)
with open(logFile, "a") as f:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    f.write(current_time + " " +  srx + " has the following srr files: " + srrFile + " moving to prefetch...")
    f.write("\n")

# print("Downloading SRR records with prefetch for " + srx + " ...")
for srr in srrList:
    with open(logFile, "a") as f:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        f.write(current_time + " prefetch for " + srx + ": " + srr + " ...")
        f.write("\n")

    prefetchLogFile = fastq_folder + "/prefetch_log_" + srx + "_" + srr + ".txt"
    prefetchSraFile = fastq_folder + "/" + srr + ".sra"
    # jobCheckFile = fastq_folder + "/job_prefetch_" + srx + "_" + srr + ".txt"
    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=500] -o " + prefetchLogFile + " \"prefetch -O " + fastq_folder + "/ " + srr + "\""

    with open(logFile, "a") as f:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        f.write(current_time + " CMD prefetch: ")
        f.write(cmd)
        f.write("\n")
    subprocess.call(cmd, shell=True)

    prefetchDone = 0
    timer = 0

    while (prefetchDone == 0):
        prefetchLogComplete = 0
        prefetchError = 0
        resubmit = 0

        if os.path.exists(prefetchLogFile):
            for line in open(prefetchLogFile):
                if "Resource usage summary" in line:
                    prefetchLogComplete = 1
                    with open(logFile, "a") as f:
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        f.write(current_time + " prefetch " + prefetchLogFile + " complete.")
                        f.write("\n")

        if prefetchLogComplete:
            time.sleep(10) # wait for 10 seconds for the .sra file to be ready.
            if not os.path.exists(prefetchSraFile):
                prefetchError = 1
                resubmit = resubmit + 1
            else:
                prefetchDone = 1
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " prefetch " + prefetchSraFile + " complete. Move on.")
                    f.write("\n")

            if resubmit: # only resubmit once if job fails
                os.remove(prefetchLogFile)
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " prefetch resubmit code is: " + str(resubmit) + ", now resubmitting prefetch job.")
                    f.write("\n")
                    subprocess.call(cmd_prefetch, shell=True)

            timer = timer + 10
            time.sleep(10)

            if timer // 9000 and not resubmit: # resubmit the job every 9000 seconds in case of the ssh jobs queue error.
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " prefetch job reach 9000-second limit, resubmit job to hpc.")
                    f.write("\n")
                    subprocess.call(cmd_prefetch, shell=True)

            if timer > 27000:
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " ERROR: prefetch job time out.")
                    f.write("\n")
                    exit()

    # Check md5 with vdb-validate command:
    cmd = "cd " + fastq_folder + " && bsub -q short -n 1 -W 0:30 -R rusage[mem=1000] -o " + fastq_folder + "/md5_" + srr + ".log.txt" + " vdb-validate " + srr + ".sra"
    md5LogFile = fastq_folder + "/md5_" + srr + ".log.txt"
    # Check md5 check result:
    if os.path.exists(md5LogFile):
        subprocess.call("rm " + md5LogFile, shell=True)

    with open(logFile, "a") as f:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        f.write(current_time + " md5check for " + srr + " ...")
        f.write("\n")
    subprocess.call(cmd, shell=True)

    md5Done = 0
    timer   = 0

    while (md5Done == 0):
        md5Complete = 0
        if os.path.exists(md5LogFile):
            for line in open(md5LogFile):
                if ("Resource usage summary" in line):
                    md5Complete = 1
                    break
        if md5Complete:
            status_code = 0
            for line in open(md5LogFile):
                if "is consistent" in line:
                    status_code = 1
                    break
            if status_code:
                md5Done = 1
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " md5check for " + srr + " passed! Move on.")
                    f.write("\n")
            else:
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " md5check for " + srr + " failed!")
                    f.write("\n")
                    exit()

        timer = timer + 10
        time.sleep(10)

        if timer // 1000: # resubmit the job every 1000 seconds in case of the ssh jobs queue error.
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " md5check job reach 1000-second limit, resubmit job to hpc.")
                f.write("\n")
                subprocess.call(cmd_prefetch, shell=True)

        if timer > 10000:
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " ERROR: md5check job time out.")
                f.write("\n")
                exit()

    # fastq-dump srr into fastq.gz files:
    cmd = "cd " + fastq_folder + " && bsub -q short -n 1 -W 4:00 -R rusage[mem=1000] -o " + fastq_folder + "/fastq_dump_" + srr + ".log.txt" + " parallel-fastq-dump -s " + srr + ".sra -t 1 -O ./ --tmpdir ./ --split-files --gzip"

    dumpLogFile = fastq_folder + "/fastq_dump_" + srr + ".log.txt"
    if os.path.exists(dumpLogFile):
        subprocess.call("rm " + dumpLogFile, shell=True)

    with open(logFile, "a") as f:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        f.write(current_time + " CMD: fastq_dump for " + srr + " ...")
        f.write(cmd)
        f.write("\n")
    subprocess.call(cmd, shell=True)

    dumpDone = 0
    timer    = 0
    while (dumpDone == 0):
        dumpLogComplete = 0
        dumpError       = 0
        resubmit        = 0
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " dumping, time passed: " + str(timer) + " ...")
            f.write("\n")

        if os.path.exists(dumpLogFile):
            for line in open(dumpLogFile):
                if ("Resource usage summary" in line):
                    dumpLogComplete = 1
            if dumpLogComplete:
                dumpError = 1
                for line in open(dumpLogFile):
                    if "Successfully completed." in line:
                        dumpDone = 1
                        dumpError = 0
                        break
                if dumpError:
                    with open(logFile, "a") as f:
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        f.write(current_time + " " + srr + " failed to dump. Need to resubmit job...")
                        f.write("\n")
                        resubmit = 1
                else:
                    with open(logFile, "a") as f:
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        f.write(current_time + " " + srr + " dumped. Move on.")
                        f.write("\n")
        if resubmit:
            os.remove(dumpLogFile)
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " fastq_dump resubmit code is: " + str(resubmit) + ", now resubmitting fastq_dump job.")
                f.write("\n")
                subprocess.call(cmd, shell=True)

        timer = timer + 30
        time.sleep(30)

        if (timer // 10000) and not resubmit: # resubmit the job every 10000 seconds in case of the ssh jobs queue error.
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " fastq_dump job reach 10000-second limit, resubmit job to hpc.")
                f.write("\n")
                subprocess.call(cmd, shell=True)

        if timer > 30000:
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " ERROR: fastq_dump job time out.")
                f.write("\n")
                exit()

# concatenate/rename fastq.gz files

## determine SE or PE:
# test = srrList[0]
# testR2 = fastq_folder + "/" + test + ".*2.fastq.gz"
# print(testR2)

# if PE:
# if os.path.isfile(testR2):
if (seqType == "pair"):
    with open(logFile, "a") as f:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        f.write(current_time + " Renaming PE fastq.gz files ...")
        f.write("\n")

    temR1 = [item + "*1.fastq.gz" for item in srrList]
    tem_filesR1 = " ".join(temR1)
    cmdR1 = "cat " + tem_filesR1 + " > " + sampleName + ".R1.fastq.gz"

    temR2 = [item + "*2.fastq.gz" for item in srrList]
    tem_filesR2 = " ".join(temR2)
    cmdR2 = "cat " + tem_filesR2 + " > " + sampleName + ".R2.fastq.gz"

    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=20480] -o " + fastq_folder + "/" + srx + "_" + sampleName + ".R1.fastq.gz.log " + "\"" + "cd " + fastq_folder + " && " + cmdR1 + " && touch " + sampleName + ".R1.fastq.gz.done" + "\""
    subprocess.call(cmd, shell=True)

    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=20480] -o " + fastq_folder + "/" + srx + "_" + sampleName + ".R2.fastq.gz.log " + "\"" + "cd " + fastq_folder + " && " + cmdR2 + " && touch " + sampleName + ".R2.fastq.gz.done" + "\""
    subprocess.call(cmd, shell=True)
# if SE:
elif (seqType == "single"):
    with open(logFile, "a") as f:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        f.write(current_time + " Renaming SE fastq.gz files ...")
        f.write("\n")

    temR1 = [item + "*1.fastq.gz" for item in srrList]
    tem_filesR1 = " ".join(temR1)
    cmdR1 = "cat " + tem_filesR1 + " > " + fastq_folder + "/" + sampleName + ".fastq.gz"

    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=20480] -o " + fastq_folder + "/" + srx + "_" + sampleName + ".fastq.gz.log " + "\"" + "cd " + fastq_folder + " && " + cmdR1 + " && touch " + sampleName + ".fastq.gz.done" + "\""
    subprocess.call(cmd, shell=True)

# clean intermediate data: check to see if the bsub cat finished or not first
timer = 0
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

    timer = timer + 10
    time.sleep(10)

    if timer // 3000: # resubmit the job every 3000 seconds in case of the ssh jobs queue error.
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " renaming job reach 3000-second limit, resubmit job to hpc.")
            f.write("\n")
            subprocess.call(cmd, shell=True)

    if timer > 10000:
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " ERROR: renaming job time out.")
            f.write("\n")
            exit()

with open(logFile, "a") as f:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    f.write(current_time + " getData job done for " + srx + "!")
    f.write("\n")
