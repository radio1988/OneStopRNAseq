#!/bin/bash
# bsub -q long -W 144:00 -R rusage[mem=4000]  'bash submit.sh'

source /home/rui.li-umw/anaconda3/etc/profile.d/conda.sh
conda activate snakemake > workflow.log  2>&1
mkdir -p log/lsf/

snakemake -p -k --jobs 99 \
--use-conda --conda-prefix ~/anaconda3/envs/osr_envs \
--ri --restart-times 0 --rerun-triggers mtime \
--cluster 'bsub -q long -N -o log/lsf/%J.lsf.txt -e log/lsf/%J.snakemake.txt -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 140:00' \
>> workflow.log  2>&1 


snakemake -p -k --jobs 99 \
--use-conda --conda-prefix ~/anaconda3/envs/osr_envs \
--ri --restart-times 2 --rerun-triggers mtime \
--cluster 'bsub -q long -N -o log/lsf/%J.lsf.txt -e log/lsf/%J.snakemake.txt -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 140:00' \
--cluster-cancel bkill \
>> workflow.log  2>&1 


snakemake -j 1 --report report.html > report.log  2>&1
