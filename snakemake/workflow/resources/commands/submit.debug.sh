#!/bin/bash
# bsub -q long -W 144:00 -R span[hosts=1] -R rusage[mem=2000] -n 2  'bash submit.debug.sh'

source ~/anaconda3/etc/profile.d/conda.sh
conda activate snakemake > workflow.log  2>&1

mkdir -p log
mkdir -p log/lsf/

snakemake -pk --jobs 999 --rerun-triggers mtime \
--keep-incomplete   --nt --ri \
--use-conda --conda-prefix ~/anaconda3/envs/osr_envs \
--restart-times 0 \
--cluster 'bsub -q short -N -o log/lsf/%J.lsf.txt -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 4:00' \
>> workflow.log  2>&1

snakemake -j 1 --report report.html > report.log  2>&1

tar cf - -C gsea gsea_bubble| pigz -p 2 > gsea/gsea_bubble.tar.gz
