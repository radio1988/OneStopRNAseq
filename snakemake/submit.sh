#!/bin/bash
# run with: nohup bash submit.sh &


# nohup snakemake --use-conda -k  --jobs 999 --latency-wait 604800 \
# --cluster 'bsub -q long -o lsf.log -R "rusage[mem={params.mem}]" -n {threads} -R span[hosts=1] -W 168:00' && snakemake --report report.html &


snakemake -p -k --jobs 999 --use-conda --latency-wait 300 \
--cluster 'bsub -q short -o lsf.log -R "rusage[mem={params.mem_mb}]" -n {threads} -R span[hosts=1] -W 4:00'

snakemake --report report.html 

# rm -rf lsf.log nohup.out meta/read_length.txt   meta/strandness.detected.txt log/ Workflow_DAG.all.svg     may*  DESeq2/ gsea/ report.html bigWig/ feature_count/ bam_qc/ rMATS*/

# bkill -r rl44w 0

# snakemake --unlock -j 1 

# snakemake -np -j 1 

# snakemake -nq -j 1 

# snakemake --export-cwl workflow.cwl

# mv DESeq2/ mapped_reads/ sorted_reads/ feature_count/ bigWig/ bam_qc/ fastqc/ trash


# Kai
# cd let analysis_1
# source /home/rl44w/anaconda3/etc/profile.d/conda.sh
# conda activate osr
# snakemake --use-conda -k -p --jobs 999 --latency-wait 300 --cluster ’bsub -q short -o lsf.log -R ‘rusage[mem={params.mem_mb}]’ -n {threads} -R span[hosts=1] -W 4:00' > log.snakemake.txt 2>&1 && snakemake --report report.html
# touch /home/kh45w/project/umw_hira_goel/OneStopRNAseq/users/89/Example_study1_02.22.06-05.28.2020/log_ssh2/analysis1.txt