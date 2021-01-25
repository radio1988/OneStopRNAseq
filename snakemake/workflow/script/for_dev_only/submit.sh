#!/bin/bash
# run with:
# nohup bash submit.sh &
# OR: bsub -q long -W 144:00 -R rusage[mem=4000] 'bash submit.sh'

module purge
source activate osr >> nohup.out  2>&1 

# fast jobs with short queue
snakemake -p -k --jobs 99 \
--latency-wait 120 \
--ri --restart-times 0 \
--cluster 'bsub -q short -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 4:00' >> nohup.out  2>&1 

# slower jobs with long queue
# low RAM jobs ran 3 times with triple RAM
snakemake -j 1 --unlock
snakemake -p -k --jobs 99 \
--latency-wait 120 \
--ri --restart-times 2 \
--cluster 'bsub -q long -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 48:00' >> nohup.out  2>&1 

snakemake -j 1 --report report.html > report.log  2>&1


# gsea compression (should be skipped)
# bash script/zip.sh
# [ -d 'gsea/' ] && tar cf - gsea/  | pigz -p 2 -f > gsea.tar.gz && mv gsea.tar.gz gsea

## Handy commands for development
# note: envmodules > singularity > conda > osr
# bkill -r rl44w -q short 0
# snakemake --unlock -j 1 
# snakemake -np -j 1 
# snakemake -nq -j 1 
# snakemake --export-cwl workflow.cwl

