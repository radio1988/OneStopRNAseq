nohup snakemake --use-conda -k  --jobs 999 --latency-wait 604800 \
--cluster 'bsub -q long -o lsf.log -R "rusage[mem={params.mem}]" -n {threads} -R span[hosts=1] -W 168:00' && snakemake --report report.html &

nohup snakemake -k --jobs 999 --latency-wait 1800 \
--cluster 'bsub -q short -o lsf.log -R "rusage[mem={params.mem_mb}]" -n {threads} -R span[hosts=1] -W 4:00' && snakemake --report report.html &

snakemake --export-cwl workflow.cwl

mv DESeq2/ mapped_reads/ sorted_reads/ feature_count/ bigWig/ bam_qc/ fastqc/ trash
