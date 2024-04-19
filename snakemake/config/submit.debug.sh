#!/bin/bash
# bsub -q long -W 144:00 -R rusage[mem=4000]  'bash submit.debug.sh'

rm -f lsf.log
source ~/anaconda3/etc/profile.d/conda.sh
conda activate osr-base > workflow.log  2>&1

snakemake -pk --jobs 99 --rerun-triggers mtime \
--keep-incomplete   \
--use-conda --conda-prefix ~/anaconda3/envs/osr_envs \
--nt --latency-wait 10 --ri --restart-times 0 \
--cluster 'bsub -q long -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 140:00' \
--cluster-cancel bkill \
>> workflow.log  2>&1

snakemake -j 1 --report report.html > report.log  2>&1
