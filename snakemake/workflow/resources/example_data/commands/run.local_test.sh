singularity exec osr-base_v0.sif snakemake -pk --ri --use-conda --cores 12 --resources mem_mb=36000 --nt  --rerun-trigger mtime
# --nt: no temp() used
# --rerun-trigger mtime: params,input,software-env,code not considered, only mtime considered
