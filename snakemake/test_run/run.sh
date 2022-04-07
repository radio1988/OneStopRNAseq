#!/user/bin/env bash

# first follow installation instructions on README.md
# modify ~/anaconda3/envs/osr_env2/ to the path you want to store conda envs
# then run this 

source activate osr-base

snakemake -p -k \
--use-conda --conda-prefix ~/anaconda3/envs/osr_env2/ \
--latency-wait 120 --ri --restart-times 0 \
-j 1

snakemake -j 1 --report report.html > report.log  2>&1
