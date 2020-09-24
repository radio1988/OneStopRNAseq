# Easy RNAseq Analysis with *oneStopRNAseq*

- This is the backend of https://mccb.umassmed.edu/OneStopRNAseq/index.php.
- The installation guide will be modified for any Unix like workstation/serve in the future.

# Here is the step-by-step guide on how to install the backend pipeline on high performance computing cluster (HPCC). Please change the directory name accordingly.

- load singularity: `module load singularity/singularity-current`
- install anaconda
- create a conda env called `osr`: `conda env create -n osr -f envs/env.yaml`
- download code and example data: `git clone git@github.com:radio1988/OneStopRNAseq.git`
- download hand_sandbox.simg and put softlink under `envs/`: `ln -s /home/rl44w/singularity/hand_sandbox.simg`

# Here is how the folder structure should look like, which is essential for the successful execution of the example workflow.

```
├── README.md
├── Snakefile *
├── config.yaml *
├── submit.sh *
├── LICENSE.md 
├── envs *
├── scripts *
├── meta *
│   ├── contrast.as.xlsx *
│   ├── contrast.de.xlsx *
│   ├── meta.xlsx *
├── example_data/genome/mm10_chr19/
```

# Here are the commands on how to run the Example Dataset on HPCC.
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
cp -r $snakemake/script/ . # important to use cp instead of softlink


## submit jobs ##
nohup bash submit.sh &

# If necessary, please kill previously submitted jobs by typing 'snakemake --unlock -j 1'
```

## Here is how to modidfy the parameter setting in config.yaml with FASTQ input files.
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
	- set `STRAND: [0, 1, 2]`
	- set `READ_LENGTH: 100` 
	- etc.
- Set absolute Path
	- $path/hand_sandbox.simg

	
