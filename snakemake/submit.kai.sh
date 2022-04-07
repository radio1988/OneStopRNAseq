cd /home/kh45w/project/umw_hira_goel/OneStopRNAseq/users/423/DIO3_07.22.24-03.19.2022/analysis_1

module purge
rm -f lsf.log
source /home/kh45w/src/anaconda3/etc/profile.d/conda.sh
conda activate /project/umw_mccb/OneStopRNAseq/conda/envs/osr-base > log.env.txt 2>&1

snakemake -p -k --jobs 16 --use-conda --conda-prefix /project/umw_mccb/OneStopRNAseq/conda/envs/conda_self_install --latency-wait 120 --ri --restart-times 2 --cluster 'bsub -q long -o lsf.log -R rusage[mem={resources.mem_mb}] -n {threads} -R span[hosts=1] -R select[rh=8] -W 140:00' >> log.snakemake.txt 2>&1
snakemake -j 1 --report report.html > report.log 2>&1
touch /home/kh45w/project/umw_hira_goel/OneStopRNAseq/users/423/DIO3_07.22.24-03.19.2022/log_ssh2/analysis1.txt
python /home/kh45w/project/umw_hira_goel/OneStopRNAseq/users/423/DIO3_07.22.24-03.19.2022/analysis_1/workflow/script/logSnakemakeParser.py log.snakemake.txt > job_status.log.txt


chmod -R 770 /home/kh45w/project/umw_hira_goel/OneStopRNAseq/users/423/DIO3_07.22.24-03.19.2022


