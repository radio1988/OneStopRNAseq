# Easy RNAseq Analysis with *oneStopRNAseq*


# Install

- install singularity
- install anaconda
- install osr env

# Essential parts
- Snakefile
- config.yaml
- envs/
- submit.sh
- README.md

# Running
- conda activate osr
- submit jobs

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