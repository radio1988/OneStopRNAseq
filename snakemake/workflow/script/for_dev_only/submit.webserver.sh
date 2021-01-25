cd $analysis_path

module purge
source /home/kh45w/src/anaconda3/etc/profile.d/conda.sh
conda activate /home/kh45w/umw_mccb/Kai/src/conda_env/osr_test >> nohup.out 2>&1

snakemake -j 1 --unlock
snakemake -p -k --jobs 16 --ri --restart-times 0 --latency-wait 120 --cluster 'bsub -q short -o lsf.log -R 'rusage[mem={resources.mem_mb}]' -n {threads} -R span[hosts=1] -W 4:00' > log.snakemake.txt 2>&1

snakemake -j 1 --unlock
snakemake -p -k --jobs 16 --ri --restart-times 0 --latency-wait 120 --cluster 'bsub -q short -o lsf.log -R 'rusage[mem={resources.mem_mb}]' -n {threads} -R span[hosts=1] -W 4:00' > log.snakemake.txt 2>&1

snakemake -j 1 --unlock
snakemake -p -k --jobs 16 --ri --restart-times 2 --latency-wait 120 --cluster 'bsub -q long -o lsf.log -R 'rusage[mem={resources.mem_mb}]' -n {threads} -R span[hosts=1] -W 48:00' > log.snakemake.txt 2>&1

snakemake -j 1 --report report.html > report.log 2>&1
touch /home/kh45w/project/umw_hira_goel/OneStopRNAseq/users/89/Example_study1_10.10.21-01.25.2021/log_ssh2/analysis1.txt
python /home/kh45w/umw_mccb/OneStopRNAseq/kai/script/logSnakemakeParser.py log.snakemake.txt > job_status.log.txt
