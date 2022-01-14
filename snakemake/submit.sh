#!/bin/bash
# bsub -q long -W 144:00 -R rusage[mem=4000] 'bash submit.sh'

module purge
rm -f lsf.log
source activate osr-base > workflow.log  2>&1 

# first pass with short queue
snakemake -p -k --jobs 99 \
--use-conda --conda-prefix ~/anaconda3/envs \
--latency-wait 120 --ri --restart-times 0 \
--cluster 'bsub -q short -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -R select[rh=8] -W 4:00' >> workflow.log  2>&1 

snakemake -j 1 --report report.html > report.log  2>&1

# slower jobs with long queue
# low RAM jobs ran 3 times with triple RAM
snakemake -j 1 --unlock

snakemake -p -k --jobs 99 \
--use-conda --conda-prefix ~/anaconda3/envs \
--latency-wait 120 --ri --restart-times 2 \
--cluster 'bsub -q long -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -R select[rh=8] -W 140:00' >> workflow.log  2>&1 

snakemake -j 1 --report report.html > report.log  2>&1
