#!/bin/bash
# run with:
# nohup bash submit.sh &
# OR: bsub -q long -W 144:00 -R rusage[mem=4000] 'bash submit.sh'

module purge
module load singularity/singularity-current > nohup.out   2>&1 
source activate osr >> nohup.out  2>&1 


# fast jobs with short queue
snakemake -p -k --jobs 16 \
--use-envmodules \
--use-singularity \
--use-conda  --conda-prefix "/project/umw_mccb/OneStopRNAseq/conda/" \
--latency-wait 120 \
--ri --restart-times 1 \
--cluster 'bsub -q short -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 4:00' >> nohup.out  2>&1 

# slower jobs with long queue
# low RAM jobs ran 3 times with triple RAM
snakemake -j 1 --unlock
snakemake -p -k --jobs 16 \
--use-envmodules \
--use-singularity \
--use-conda  --conda-prefix "/project/umw_mccb/OneStopRNAseq/conda/" \
--latency-wait 300 \
--ri --restart-times 2 \
--cluster 'bsub -q long -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 48:00' >> nohup.out  2>&1 


# note: envmodules > singularity > conda > osr

snakemake -j 1 --report report.html > report.log  2>&1

# bash script/zip.sh

# gsea compression (should be skipped)
# [ -d 'gsea/' ] && tar cf - gsea/  | pigz -p 2 -f > gsea.tar.gz && mv gsea.tar.gz gsea

## Handy commands for development
# bkill -r rl44w 0
# snakemake --unlock -j 1 
# snakemake -np -j 1 
# snakemake -nq -j 1 
# snakemake --export-cwl workflow.cwl
# mv DESeq2/ mapped_reads/ sorted_reads/ feature_count/ bigWig/ bam_qc/ fastqc/ trash

# From Kai
# cd let analysis_1
# source /home/rl44w/anaconda3/etc/profile.d/conda.sh
# conda activate osr
# snakemake --use-conda -k -p --jobs 999 --latency-wait 300 --cluster ’bsub -q short -o lsf.log -R ‘rusage[mem={params.mem_mb}]’ -n {threads} -R span[hosts=1] -W 4:00' > log.snakemake.txt 2>&1 && snakemake --report report.html
# touch /home/kh45w/project/umw_hira_goel/OneStopRNAseq/users/89/Example_study1_02.22.06-05.28.2020/log_ssh2/analysis1.txt

