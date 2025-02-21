# Easy RNAseq Analysis with *oneStopRNAseq*
- This is the backend of https://mccb.umassmed.edu/OneStopRNAseq/index.php
- Citation: Li R, Hu K, Liu H, Green MR, Zhu LJ. OneStopRNAseq: A Web Application for Comprehensive and Efficient Analyses of RNA-Seq Data. Genes (Basel). 2020 Oct 2;11(10):1165. doi: 10.3390/genes11101165. PMID: 33023248; PMCID: PMC7650687.
- Can be installed on Linux systems (tested on Ubuntu and CentOS, did not test on other flavors of Linux, not working on Mac yet)
- Would need knowledge in basic bash commands to use this workflow

# Installation
- install anaconda by following instructions on https://docs.anaconda.com/anaconda/install/  (installation tested in conda version conda 4.12.0)
- download OneStopRNASeq by `git clone https://github.com/radio1988/OneStopRNAseq.git`
- `cd OneStopRNAseq/snakemake`
- `conda env create -f workflow/envs/osr-base.yaml`  # create an conda env called 'osr-base'
- `conda activate osr-base`
- `Rscript -e 'install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz", repos=NULL, type="source");'` # install packages for QoRTs
- uncompress example dataset
    - `gunzip example_data/genome/mm10_chr19/mm10.chr19.fa.gz` 
    - `gunzip example_data/genome/mm10_chr19/gencode.vM25.primary_assembly.annotation.ch19.gtf.gz`
- the installation typically takes 15-30 mins

# Running OneStopRNASeq workflow on example datasets
## Start from FASTQ as input
### Preparations
```
mkdir fq_analysis && cd fq_analysis # create workdir
osr_path=$download_path/OneStopRNASeq/snakemake  # this is user specific, e.g. /home/user/git/OneStopRNAseq/snakemake
# put necessary files into workdir
ln -s $osr_path/meta
ln -s $osr_path/example_data
ln -s $osr_path/example_data/fastq_small/ fastq
cp $osr_path/config.fq_example.yaml config.yaml
mkdir -p workflow && cd workflow
ln -s $osr_path/workflow/envs
ln -s $osr_path/workflow/Snakefile
ln -s $osr_path/workflow/osr.py
rsync  -a $osr_path/workflow/script ./
```
### Running workflow
```
snakemake -j 1 -np   # quick test with a 'dry-run'
snakemake -j 2 -pk  --use-conda --conda-prefix ~/anaconda3/envs --ri # automatically download packages needed for the workflow, run the workflow on the example datase with two threads, takes around 30 min for the first run
```

# Quick start on snakemake advanced usage:
- If the workflow did not finish completely, try submitting the jobs again with `snakemake -j 2 -pk`
- If the workdir is locked in a second submission, please kill previously submitted jobs by typing 'snakemake --unlock -j 1', after you make sure the first submission has stopped running
- The workflow can be easily adapted to a LSF system by using `snakemake -pk --cluster 'bsub -q queue_name -o lsf.log -R "rusage[mem={resources.mem_mb}]" -n {threads} -R span[hosts=1] -W 24:00'` 

# Viewing results
- results are under subfolders in the workdir, e.g. DESeq2, gsea, rMATS.x, fastqc, bam_qc
- interpretation of results can be found here: https://mccb.umassmed.edu/OneStopRNAseq/documents/description_of_output_files.pdf 
- write up for method section in publication can be found here: https://mccb.umassmed.edu/OneStopRNAseq/documents/template_of_method_section.pdf 
