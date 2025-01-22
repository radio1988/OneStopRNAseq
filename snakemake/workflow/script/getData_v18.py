"""
Given SRX, sampleLabel, and fastq folder path, perform srr download, md5 check, fastq-dump and renaming.
Usage:
python getData.py SRX_ID SampleLabel FastqFolder pair|single

v12 notes:
1. fix the .sra and .lite incompatibility and download to folder issue with "-o" flag when using prefetch.

V16 notes: based on v12
1. use long node for fastq dump
2. use 200GB for prefetch
3. sra subfolder issue already fixed before this version
"""

import sys
import time
import os
import subprocess
from datetime import datetime
from os.path import exists
import pathlib

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
cmd = "bsub -q short -n 1 -W 0:10 -R rusage[mem=500] -o " + efetchLogFile + " \"esearch -db sra -query " + srx + " | efetch -format runinfo | grep -v '^Run,' | awk 'NR==1' | cut -d ',' -f 1 | head -n 1 > " + efetchResFile + "\""
with open(logFile, "a") as f:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    f.write(current_time + " CMD efetch: \n")
    f.write(cmd + "\n")
subprocess.call(cmd, shell=True)

efetchDone = 0
timer = 0

# Examine for the efetch short job execution:
gap_efetch = 3000
gap_prefetch = 9000
gap_md5check = 1000
gap_fastq_dump = 50000
gap_cat_fastq = 1000
gap_rename = 3000

while (efetchDone == 0):
    efetchLogComplete = 0
    efetchError = 1
    resubmit = 0

    if os.path.exists(efetchLogFile):
        for line in open(efetchLogFile):
            if "The output (if any) is above this job summary." in line:
                efetchLogComplete = 1
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " efetch " + efetchLogFile + " complete.\n")
                # since we also check the res file, so the below check of the log file is redundant
                for line in open(efetchLogFile):
                    if "Successfully completed." in line:
                        efetchError = 0
                if efetchError == 1:
                    with open(logFile, "a") as f:
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        f.write(current_time + " efetch " + efetchLogFile + " didn't finish successfully, need to resubmit job.\n")
                        resubmit = resubmit + 1
                    break

    if efetchLogComplete and not efetchError:
        if os.path.exists(efetchResFile):
            time.sleep(1) # make sure efetchResFile finishes saving to disk
            efetchResFileContentList = [line for line in open(efetchResFile)] # keep a copy of efetchResFile so it won't be overwritten.
            # if not os.stat(efetchResFile).st_size: # if file empty:
            if len(efetchResFileContentList) < 1:
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " efetch " + efetchResFile + " is done but is empty, need to resubmit job.\n")
                    resubmit = resubmit + 1
            elif len(efetchResFileContentList[0]) < 3: # empty file could contain a single '\n'
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " efetch " + efetchResFile + " is done but is empty, need to resubmit job.\n")
                    resubmit = resubmit + 1
            else: # if not empty:
                efetchDone = 1
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " efetch " + efetchResFile + " is done and not empty. Move on.\n")

    if resubmit: # only resubmit once if job fails
        os.remove(efetchLogFile)
        os.remove(efetchResFile)
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " efetch resubmit code is: " + str(resubmit) + ", now resubmitting efetch job.\n")
            subprocess.call(cmd, shell=True)

    if os.path.exists(efetchLogFile): # Start timer only when log file exists meaning that bsub job starts executing
        timer = timer + 10
    time.sleep(10)

    if timer // gap_efetch: # resubmit the job every 3000 seconds in case of the ssh jobs queue error.
        gap_efetch = gap_efetch + gap_efetch
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " efetch job reach 3000-second limit, resubmit job to hpc.\n")
            os.remove(efetchLogFile)
            os.remove(efetchResFile)
            subprocess.call(cmd, shell=True)

    if timer > 10000:
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " ERROR: efetch job time out.\n")
            exit()


srrList = [line.strip() for line in efetchResFileContentList]
srrList = [x for x in srrList if not x == ""]

srrFile = " ".join(srrList)
with open(logFile, "a") as f:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    f.write(current_time + " " +  srx + " has the following srr files: " + srrFile + " moving to prefetch...\n")

# print("Downloading SRR records with prefetch for " + srx + " ...")
for srr in srrList:
    with open(logFile, "a") as f:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        f.write(current_time + " prefetch for " + srx + ": " + srr + " ...\n")

    prefetchLogFile = fastq_folder + "/prefetch_log_" + srx + "_" + srr + ".txt"
    prefetchSraFile = fastq_folder + "/" + srr + ".sra"
    # jobCheckFile = fastq_folder + "/job_prefetch_" + srx + "_" + srr + ".txt"
    cmd_clearPrefetchLogFile  = "rm " + prefetchLogFile # Clear potential log file before submit new job
    cmd_clearPrefetchPrevious = "rm " + fastq_folder + "/" + srr + ".sra*" # get rid of srx.sra* files otherwise may prevent new prefetch from executing.
    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=500] -o " + prefetchLogFile + " \"prefetch -f yes --max-size 200GB -o " + fastq_folder + "/" + srr + ".sra " + srr + "\""
    
    with open(logFile, "a") as f:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        f.write(current_time + " CMD prefetch: \n")
        f.write(cmd + "\n")
    subprocess.call(cmd_clearPrefetchLogFile, shell=True)
    subprocess.call(cmd_clearPrefetchPrevious, shell=True)
    subprocess.call(cmd, shell=True)


    prefetchDone = 0
    timer = 0

    while (prefetchDone == 0):
        prefetchLogComplete = 0
        prefetchError = 1
        resubmit = 0

        if os.path.exists(prefetchLogFile):
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " prefetch " + prefetchLogFile + " exists.\n")
            for line in open(prefetchLogFile):
                if "The output (if any) is above this job summary." in line:
                    prefetchLogComplete = 1
                    with open(logFile, "a") as f:
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        f.write(current_time + " prefetch " + prefetchLogFile + " complete.\n")
                    # since we also check the res file, so the below check of the log file is redundant
                    for line in open(prefetchLogFile):
                        if "Successfully completed." in line:
                            with open(logFile, "a") as f:
                                now = datetime.now()
                                current_time = now.strftime("%H:%M:%S")
                                f.write(current_time + " prefetch " + prefetchLogFile + " log file successfully completed.\n")
                                # move sra out of its containing folder if so:
                                srafiles = list(pathlib.Path(fastq_folder).glob("**/" + srr + ".sra"))
                                if len(srafiles) > 0:
                                    pathlib.Path(srafiles[0]).rename(fastq_folder + "/" + srr + ".sra")
                            prefetchError = 0
                            break
                    if prefetchError == 1:
                        with open(logFile, "a") as f:
                            now = datetime.now()
                            current_time = now.strftime("%H:%M:%S")
                            f.write(current_time + " prefetch " + prefetchLogFile + " didn't finish successfully, need to resubmit job.\n")
                            resubmit = resubmit + 1
                        break

        if prefetchLogComplete and not prefetchError:
            time.sleep(10) # wait for 10 seconds for the .sra file to be ready.
            if not os.path.exists(prefetchSraFile):
                prefetchError = 1
                resubmit = resubmit + 1
            else:
                prefetchDone = 1
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " prefetch " + prefetchSraFile + " complete. Move on.\n")

        if resubmit: # only resubmit once if job fails
            #os.remove(prefetchLogFile)
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " prefetch resubmit code is: " + str(resubmit) + ", now resubmitting prefetch job.\n")
                cmd_clearPrefetchLogFile  = "rm " + prefetchLogFile # Clear potential log file before submit new job
                cmd_clearPrefetchPrevious = "rm " + fastq_folder + "/" + srr + ".sra*" # get rid of srx.sra* files otherwise prefetch complians
                cmd_clearPrefetchFolder = "rm -r " + fastq_folder + "/" + srr # Occasionally, there will be a lingering SRR folder that prevent future prefetch downloading.

                subprocess.call(cmd_clearPrefetchLogFile, shell=True)
                subprocess.call(cmd_clearPrefetchPrevious, shell=True)
                subprocess.call(cmd_clearPrefetchFolder, shell=True)
                subprocess.call(cmd, shell=True)
        if os.path.exists(prefetchLogFile):
            timer = timer + 300
        time.sleep(60)

        if timer // gap_prefetch and not resubmit: # resubmit the job every 9000 seconds in case of the ssh jobs queue error.
            gap_prefetch = gap_prefetch + gap_prefetch
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " prefetch job reach 9000-second limit, resubmit job to hpc.\n")
                # os.remove(prefetchLogFile)
                os.rename(prefetchLogFile, prefetchLogFile + "_failed")
                subprocess.call(cmd, shell=True)

        if timer > 27000:
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " ERROR: prefetch job time out.")
                f.write("\n")
                exit()

    # Check md5 with vdb-validate command:
    cmd = "cd " + fastq_folder + " && bsub -q short -n 1 -W 4:00 -R rusage[mem=10000] -o " + fastq_folder + "/md5_" + srr + ".log.txt" + " vdb-validate " + srr + ".sra"
    md5LogFile = fastq_folder + "/md5_" + srr + ".log.txt"
    # Check md5 check result:
    if os.path.exists(md5LogFile):
        subprocess.call("rm " + md5LogFile, shell=True)

    with open(logFile, "a") as f:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        f.write(current_time + " md5check for " + srr + " ...\n")
    subprocess.call(cmd, shell=True)

    md5Done = 0
    timer   = 0

    while (md5Done == 0):
        md5LogComplete = 0
        md5Error = 1
        resubmit = 0

        if os.path.exists(md5LogFile):
            for line in open(md5LogFile):
                if "The output (if any) is above this job summary." in line:
                    md5LogComplete = 1
                    with open(logFile, "a") as f:
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        f.write(current_time + " md5check " + md5LogFile + " complete.\n")
                    for line in open(md5LogFile):
                        if "Successfully completed." in line:
                            md5Error = 0
                            break
                    if md5Error == 1:
                        with open(logFile, "a") as f:
                            now = datetime.now()
                            current_time = now.strftime("%H:%M:%S")
                            f.write(current_time + " md5check " + md5LogFile + " didn't finish successfully, need to resubmit job.\n")
                            resubmit = resubmit + 1
                        break

        if md5LogComplete and not md5Error:
            md5Done = 1
            status_code = 0
            for line in open(md5LogFile):
                if "is consistent" in line:
                    status_code = 1
                    break

            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                if status_code:
                    f.write(current_time + " md5check complete, is consistent, move on.\n")
                else:
                    f.write(current_time + " md5check complete, not consistent, quit.\n")

            if resubmit: # only resubmit once if job fails
                os.remove(md5LogFile)
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " md5check resubmit code is: " + str(resubmit) + ", now resubmitting md5check job.\n")
                    subprocess.call(cmd, shell=True)

            if timer // gap_md5check and not resubmit: # resubmit the job every 9000 seconds in case of the ssh jobs queue error.
                gap_md5check = gap_md5check + gap_md5check
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " md5check job reach 1000-second limit, resubmit job to hpc.\n")
                    os.remove(md5LogFile)
                    subprocess.call(cmd, shell=True)

            if timer > 5000:
                with open(logFile, "a") as f:
                    now = datetime.now()
                    current_time = now.strftime("%H:%M:%S")
                    f.write(current_time + " ERROR: md5check job time out.\n")
                    exit()

        if os.path.exists(md5LogFile):
            timer = timer + 10
        time.sleep(10)

        if timer // gap_md5check and not resubmit: # resubmit the job every 1000 seconds in case of the ssh jobs queue error.
            gap_md5check = gap_md5check + gap_md5check
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " md5check job reach 1000-second limit, resubmit job to hpc.\n")
                os.remove(md5LogFile)
                subprocess.call(cmd, shell=True)

        if timer > 10000:
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " ERROR: md5check job time out.\n")
                exit()

    # fastq-dump srr into fastq.gz files:
    # In case sra files saved under SRRXXX folder, move them out first:
    sra_tem = fastq_folder + "/" + srr + "/" + srr + ".sra"
    if exists(sra_tem):
        cmd_mv = "mv " + sra_tem + " " + fastq_folder
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " CMD: mv sra")
            f.write(cmd_mv + "\n")
        subprocess.call(cmd_mv, shell = True)

    cmd = "cd " + fastq_folder + " && bsub -q long -n 1 -W 40:00 -R rusage[mem=10000] -o " + fastq_folder + "/fastq_dump_" + srr + ".log.txt" + " parallel-fastq-dump -s " + srr + ".sra -t 1 -O ./ --tmpdir ./ --split-files --gzip"

    dumpLogFile = fastq_folder + "/fastq_dump_" + srr + ".log.txt"
    if os.path.exists(dumpLogFile):
        subprocess.call("rm " + dumpLogFile, shell=True)

    with open(logFile, "a") as f:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        f.write(current_time + " CMD: fastq_dump for " + srr + " ...")
        f.write(cmd + "\n")
    subprocess.call(cmd, shell=True)

    dumpDone = 0
    timer    = 0

    while (dumpDone == 0):
        dumpLogComplete = 0
        dumpError       = 1
        resubmit        = 0

        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " dumping, time passed: " + str(timer) + " ...\n")

        if os.path.exists(dumpLogFile):
            for line in open(dumpLogFile):
                if ("The output (if any) is above this job summary." in line):
                    dumpLogComplete = 1
                    with open(logFile, "a") as f:
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        f.write(current_time + " dump " + dumpLogFile + " complete.\n")
                    for line in open(dumpLogFile):
                        if "Successfully completed." in line:
                            dumpDone = 1
                            dumpError = 0
                            break
                    if dumpError == 1:
                        with open(logFile, "a") as f:
                            now = datetime.now()
                            current_time = now.strftime("%H:%M:%S")
                            f.write(current_time + " " + srr + " failed to dump. Need to resubmit job...\n")
                            resubmit = 1
                    else:
                        with open(logFile, "a") as f:
                            now = datetime.now()
                            current_time = now.strftime("%H:%M:%S")
                            f.write(current_time + " " + srr + " dumped. Move on.\n")
        if resubmit:
            os.remove(dumpLogFile)
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " fastq_dump resubmit code is: " + str(resubmit) + ", now resubmitting fastq_dump job.\n")
                subprocess.call(cmd, shell=True)

        if os.path.exists(dumpLogFile):
            timer = timer + 30
        time.sleep(30)

        if (timer // gap_fastq_dump) and not resubmit: # resubmit the job every 50000 seconds in case of the ssh jobs queue error.
            gap_fastq_dump = gap_fastq_dump + gap_fastq_dump
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " fastq_dump job reach 10000-second limit, resubmit job to hpc.\n")
                os.remove(dumpLogFile)
                subprocess.call(cmd, shell=True)

        if timer > gap_fastq_dump * 2:
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
        f.write(current_time + " Renaming PE fastq.gz files ...\n")

    temR1 = [item + "*1.fastq.gz" for item in srrList]
    tem_filesR1 = " ".join(temR1)
    cmdR1 = "cat " + tem_filesR1 + " > " + sampleName + ".R1.fastq.gz"

    temR2 = [item + "*2.fastq.gz" for item in srrList]
    tem_filesR2 = " ".join(temR2)
    cmdR2 = "cat " + tem_filesR2 + " > " + sampleName + ".R2.fastq.gz"

    # Note that don't create data_geo_ for R1, only do it after R2 finishes, otherwise data.txt will be created prematurely...
    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=20480] -o " + fastq_folder + "/" + srx + "_" + sampleName + ".R1.fastq.gz.log " + "\"" + "cd " + fastq_folder + " && " + cmdR1 + " && touch " + sampleName + ".R1.fastq.gz.done" + "\""
    subprocess.call(cmd, shell=True)

    catLogFile = fastq_folder + "/" + srx + "_" + sampleName + ".R1.fastq.gz.log"
    catDone = 0
    timer = 0
    while (catDone == 0):
        catLogComplete = 0
        catError = 1
        resubmit = 0

        if os.path.exists(catLogFile):
            for line in open(catLogFile):
                if "The output (if any) is above this job summary." in line:
                    catLogComplete = 1
                    with open(logFile, "a") as f:
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        f.write(current_time + " catFastq " + catLogFile + " complete.\n")

                    for line in open(catLogFile):
                        if "Successfully completed." in line:
                            catError = 0
                            break
                    if catError == 1:
                        with open(logFile, "a") as f:
                            now = datetime.now()
                            current_time = now.strftime("%H:%M:%S")
                            f.write(current_time + " catFastq " + catLogFile + " didn't finish successfully, need to resubmit job.")
                            resubmit = resubmit + 1
                            f.write("\n")
                        break

        if catLogComplete and not catError:
            catDone = 1
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " catFastq complete, move on.\n")
        if resubmit:
            os.remove(catLogFile)
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " catFastq resubmit code is: " + str(resubmit) + ", now resubmitting catFastq job.\n")
                subprocess.call(cmd, shell=True)

        if os.path.exists(catLogFile):
            timer = timer + 30
        time.sleep(30)

        if timer // gap_cat_fastq:
            gap_cat_fastq = gap_cat_fastq + gap_cat_fastq
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " catFastq job reach 1000-second limit, resubmit job to hpc.\n")
                os.remove(catLogFile)
                subprocess.call(cmd, shell=True)
        if timer > 10000:
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " ERROR: catFastq job time out.\n")
                exit()

    # print("R1 done ...")
    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=20480] -o " + fastq_folder + "/" + srx + "_" + sampleName + ".R2.fastq.gz.log " + "\"" + "cd " + fastq_folder + " && " + cmdR2 + " && touch " + sampleName + ".R2.fastq.gz.done" + " && touch " + fastq_folder + "/../log_ssh2/data_geo_" + srx + ".txt" + "\""
    # print("R2 calling ...")
    subprocess.call(cmd, shell=True)

    catLogFile = fastq_folder + "/" + srx + "_" + sampleName + ".R2.fastq.gz.log"
    catDone = 0
    timer = 0
    while (catDone == 0):
        catLogComplete = 0
        catError = 1
        resubmit = 0

        if os.path.exists(catLogFile):
            for line in open(catLogFile):
                if "The output (if any) is above this job summary." in line:
                    catLogComplete = 1
                    with open(logFile, "a") as f:
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        f.write(current_time + " catFastq " + catLogFile + " complete.\n")
                    for line in open(catLogFile):
                        if "Successfully completed." in line:
                            catError = 0
                            break
                    if catError == 1:
                        with open(logFile, "a") as f:
                            now = datetime.now()
                            current_time = now.strftime("%H:%M:%S")
                            f.write(current_time + " catFastq " + catLogFile + " didn't finish successfully, need to resubmit job.")
                            resubmit = resubmit + 1
                            f.write("\n")
                        break

        if catLogComplete and not catError:
            catDone = 1
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " catFastq complete, move on.\n")
        if resubmit:
            os.remove(catLogFile)
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " catFastq resubmit code is: " + str(resubmit) + ", now resubmitting catFastq job.\n")
                subprocess.call(cmd, shell=True)

        if os.path.exists(catLogFile):
            timer = timer + 30
        time.sleep(30)

        if timer // gap_cat_fastq:
            gap_cat_fastq = gap_cat_fastq + gap_cat_fastq
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " catFastq job reach 1000-second limit, resubmit job to hpc.\n")
                os.remove(catLogFile)
                subprocess.call(cmd, shell=True)
        if timer > 10000:
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " ERROR: catFastq job time out.\n")
                exit()

# if SE:
elif (seqType == "single"):
    with open(logFile, "a") as f:
        now = datetime.now()
        current_time = now.strftime("%H:%M:%S")
        f.write(current_time + " Renaming SE fastq.gz files ...\n")

    temR1 = [item + "*.fastq.gz" for item in srrList]
    tem_filesR1 = " ".join(temR1)
    cmdR1 = "cat " + tem_filesR1 + " > " + fastq_folder + "/" + sampleName + ".fastq.gz"

    cmd = "bsub -q short -n 1 -W 4:00 -R rusage[mem=20480] -o " + fastq_folder + "/" + srx + "_" + sampleName + ".fastq.gz.log " + "\"" + "cd " + fastq_folder + " && " + cmdR1 + " && touch " + sampleName + ".fastq.gz.done" + " && touch " + fastq_folder + "/../log_ssh2/data_geo_" + srx + ".txt" + "\""
    subprocess.call(cmd, shell=True)

    catLogFile = fastq_folder + "/" + srx + "_" + sampleName + ".fastq.gz.log"
    catDone = 0
    timer = 0
    while (catDone == 0):
        catLogComplete = 0
        catError = 1
        resubmit = 0

        if os.path.exists(catLogFile):
            for line in open(catLogFile):
                if "The output (if any) is above this job summary." in line:
                    catLogComplete = 1
                    with open(logFile, "a") as f:
                        now = datetime.now()
                        current_time = now.strftime("%H:%M:%S")
                        f.write(current_time + " catFastq " + catLogFile + " complete.\n")
                    for line in open(catLogFile):
                        if "Successfully completed." in line:
                            catError = 0
                            break
                    if catError == 1:
                        with open(logFile, "a") as f:
                            now = datetime.now()
                            current_time = now.strftime("%H:%M:%S")
                            f.write(current_time + " catFastq " + catLogFile + " didn't finish successfully, need to resubmit job.")
                            resubmit = resubmit + 1
                            f.write("\n")
                        break

        if catLogComplete and not catError:
            catDone = 1
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " catFastq complete, move on.\n")
        if resubmit:
            os.remove(catLogFile)
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " catFastq resubmit code is: " + str(resubmit) + ", now resubmitting catFastq job.\n")
                subprocess.call(cmd, shell=True)

        if os.path.exists(catLogFile):
            timer = timer + 30
        time.sleep(30)

        if timer // gap_cat_fastq:
            gap_cat_fastq = gap_cat_fastq + gap_cat_fastq
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " catFastq job reach 1000-second limit, resubmit job to hpc.\n")
                os.remove(catLogFile)
                subprocess.call(cmd, shell=True)
        if timer > 10000:
            with open(logFile, "a") as f:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                f.write(current_time + " ERROR: catFastq job time out.\n")
                exit()

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

    if timer // gap_rename: # resubmit the job every 3000 seconds in case of the ssh jobs queue error.
        gap_rename = gap_rename + gap_rename
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " renaming job reach 3000-second limit, resubmit job to hpc.\n")
            subprocess.call(cmd, shell=True)

    if timer > 10000:
        with open(logFile, "a") as f:
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            f.write(current_time + " ERROR: renaming job time out.\n")
            exit()

with open(logFile, "a") as f:
    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    f.write(current_time + " getData job done for " + srx + "!\n")
