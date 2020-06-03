# Easy RNAseq Analysis with *oneStopRNAseq*


# Install

- install singularity
- install anaconda
- install osr env: `conda env create -n osr -f envs/env.yaml` 
- download and put hand_sandbox/ under envs/

# Essential parts
- Snakefile
- config.yaml
- envs/
- genome/
- meta/
- script/
- submit.sh
- README.md

# Running
- `mkdir analysis`
- ln -s ../snakemake/Snakefile 
- ln -s ../snakemake/config.yaml 
- ln -s ../snakemake/meta/
- ln -s ../snakemake/envs/
- ln -s ../snakemake/genome/
- ln -s ../snakemake/submit.sh 
- ln -s ../snakemake/fastq/
- cp -r ../snakemake/script/ .


- `source activate /home/rl44w/anaconda3/envs/osr`
- submit jobs 

	```
	kill previous running submission job
	snakemake --unlock -j 1 # if you have run this pipeline before and exitted or stopped
	nohup snakemake -k --jobs 999 --use-conda --latency-wait 300 \
	--cluster 'bsub -q short -o lsf.log -R "rusage[mem={params.mem_mb}]" -n {threads} -R span[hosts=1] -W 4:00' && snakemake --report report.html &
	```

## Start from FASTQ
- provide fastq files under fastq/
	- {sample}.R1.fastq.gz {sample}.R2.fastq.gz for PE reads
	- {sample}.fastq.gz for SE reads
- write sample names in config.yaml: 
	- {sample}
- write meta data
	- meta/meta.xlsx
	- meta/contrast.xlsx 
- set `START` flag in `config.yaml` as `FASTQ`

## Start from BAM

## Start from COUNT

## Start from 

## Caveats
- must copy script folder to let DESeq2 save results into the correct location


## todo:
- contrast.de.xlsx  contrast.as.xlsx
- use same contrast format for rMATS, at group level
- an easier way to link and run workflow: now script/ needs cp 
- handle gsea err

## Log:
- test success on Mac, Rui-HPCC-snakemake, Rui-HPCC-test_run
- test on u18 VM: 
	- test with osr (new)
		- star index output error 



