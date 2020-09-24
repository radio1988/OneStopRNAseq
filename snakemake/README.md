# Easy RNAseq Analysis with *oneStopRNAseq*

- Currently, this is the backend of https://mccb.umassmed.edu/OneStopRNAseq/index.php and intended for internal use
- Installation will be simplified and make it suitable for download and usage in any Unix like workstation/serve in the future

# Install (for any user on HPCC)

- load singularity: `module load singularity/singularity-current`
- install anaconda
- create a conda env called `osr`: `conda env create -n osr -f envs/env.yaml`
- download code and example data: `git clone git@github.com:radio1988/OneStopRNAseq.git`
- download hand_sandbox.simg and put softlink under `envs/`: `ln -s /home/rl44w/singularity/hand_sandbox.simg`

# Folder structure
- essential for execution of example workflow marked *

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


## submit jobs ##
nohup bash submit.sh &

# If applicable, kill previous running submitted job, then 'snakemake --unlock -j 1'
```


## Start from FASTQ config.yaml settings
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
	- etc.
- Absolute Path
	- $path/hand_sandbox.simg

	
