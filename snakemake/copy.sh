#!/bin/bash
rsync -aP /project/umw_mccb/OneStopRNAseq/rui/develop/OneStopRNAseq/snakemake/Snakefile .
rsync -aP  /project/umw_mccb/OneStopRNAseq/rui/develop/OneStopRNAseq/snakemake/script .
rsync -aP  /project/umw_mccb/OneStopRNAseq/rui/develop/OneStopRNAseq/snakemake/submit.sh .
rsync -aP  /project/umw_mccb/OneStopRNAseq/rui/develop/OneStopRNAseq/snakemake/meta .
ln -s /project/umw_mccb/OneStopRNAseq/rui/develop/OneStopRNAseq/snakemake/envs/
source activate osr
snakemake -j 1 reset
echo "you need to create config.yaml, fastq/"
