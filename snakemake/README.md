# Easy RNAseq Analysis with *oneStopRNAseq*


# Install (for any user on HPCC)

- load singularity: `module load singularity/singularity-current`
- install anaconda
- create a conda env called `osr`: `conda env create -n osr -f envs/env.yaml`
- download code and example data: `git clone git@github.com:radio1988/OneStopRNAseq.git`
- download hand_sandbox.simg and put softlink under `envs/`: `ln -s /home/rl44w/singularity/hand_sandbox.simg`


# Folder structure
- essential for execution of example workflow marked *


```
├── .gitignore
├── README.md
├── Snakefile *
├── config.yaml *
├── submit.sh *
├── LICENSE.md 
├── envs
│   ├── hand_sandbox.simg *
│   ├── env.yaml * 
│   ├── rmats.yaml *
│   ├── gsea_db *
│   ├── rMATS.3.2.5 *
│   ├── GSEA_4.0.3 *
├── scripts
│   ├── DESeq2.Rmd *
│   ├── strandness_detection.py *
├── meta
│   ├── contrast.as.xlsx *
│   ├── contrast.de.xlsx *
│   ├── meta.xlsx *
├── genome/mm10_chr19/
│   ├── mm10.chr19.fa *
│   ├── gencode.vM21.chr19.gtf *
```

# Running Example Dataset(on HPCC)
```
## Preps ##
mkdir analysis && cd analysis
snakemake=/project/umw_mccb/OneStopRNAseq/rui/test_run

ln -s $snakemake/Snakefile 
ln -s $snakemake/envs/
ln -s $snakemake/genome/
ln -s $snakemake/submit.sh 
cp $snakemake/config.yaml .
cp -r $snakemake/meta .
cp -r $snakemake/fastq .
cp -r $snakemake/script/ . # have to cp, softlink has problems (todo: fix)

## Running ##
source activate /home/rl44w/anaconda3/envs/osr
# If applicable, kill previous running submitted job, then 'snakemake --unlock -j 1'
nohup submit.sh &
```


## Start from FASTQ (For Kai's initial test)
- provide fastq files under fastq/
	- {sample}.R1.fastq.gz {sample}.R2.fastq.gz for PE reads
	- {sample}.fastq.gz for SE reads
- write sample names in config.yaml: 
	- {sample}
- write meta data
	- meta/meta.xlsx
	- meta/contrast.de.xlsx 
	- meta/contrast.as.xlsx 
- modify `config.yaml`
	- set `START: "FASTQ"`
	- set `STRAND: [0, 1, 2]` (for testing, June 03)
	- set `READ_LENGTH: 100` (for testing, June 03)
- Absolute Path
	- /project/umw_yong-xu_wang/singularity/hand_sandbox.simg

	
## Start from BAM

## Start from COUNT

## Start from RNK


# Development Note
## Workflow
- develop in HPCC: `/project/umw_mccb/OneStopRNAseq/rui/develop/`
	- mount folder to Mac for GUI
	- use u18 to update simg
	- git commit all function level changes
	- delete Mac version (Mac support for STAR is bad)
- deliver to Kai: `/project/umw_mccb/OneStopRNAseq/rui/delivery`


## Caveats and todo
- an easier way to link and run workflow: now script/ needs cp 
- handle gsea err
- Read Length detection from BAM 
- must copy script folder to let DESeq2 save results into the correct location
- must use --use-singularity and --use-conda; rmat conda env takes some time to init; 
- singularity 4GB to download
- Mac can't run STAR, evey tried singularity image build from u18
- √ osr and singularity(2.7.1) both have STAR, conda STAR(2.7.4) bad: now use conda

## Todo
- DESeq2 bug: no p-value histogram after the first one, just output as pdf?
- upgrade the latest GenCode for Human and Mouse, Add TE to anno file
- add custom genome and gtf
- rMATS in singularity to speed up overhead
- SalmonTE in github? in singularity? add installment code?

1. GSEA-DB upload
2. rMATS filter
3. sample names can't be numbers, must be strings





## Log:
- test from FQ, with `READ_LENGTH: 100`  `STRAND: [0, 1, 2]`
	- `--use-singularity --use-conda`
		- success on u18 VM: 
		- always fail on Mac: STAR zero mapping
		- fails on HPCC, unless removing gloabal singularity line in header
	- `--use-singularity` with gloabal singularity (less interference)
		- success on u18 VM: 
		- Mac: not working, singularity poor MAC support?
		- HPCC: works, except rMATS, not set up properlly 
		- 


