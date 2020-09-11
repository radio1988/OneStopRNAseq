#!/bin/bash
rsync -aP /project/umw_mccb/OneStopRNAseq/rui/develop/OneStopRNAseq/snakemake/Snakefile .
rsync -aP  /project/umw_mccb/OneStopRNAseq/rui/develop/OneStopRNAseq/snakemake/script .
rsync -aP  /project/umw_mccb/OneStopRNAseq/rui/develop/OneStopRNAseq/snakemake/submit.sh .
ln -s /project/umw_mccb/OneStopRNAseq/rui/develop/OneStopRNAseq/snakemake/envs/
source activate osr
echo "you need to create config.yaml, fastq/, meta/"
