# install anaconda
https://docs.anaconda.com/anaconda/install/

# create conda env called osr
conda clean --all
conda env create -n osr -f workflow/envs/osr.yaml

# activate osr
conda activate osr

# install R packages into R in osr
# for DESeq2 
Rscript -e "install.packages( c('BiocManager', 'ggrepel', 'gdata', 'plyr', 'pheatmap', 'PoiClaClu', 'gridExtra', 'rmarkdown', 'ashr'), repos='https://cloud.r-project.org')"
Rscript -e "BiocManager::install('scater')"
# for QoRTs
Rscript -e 'install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz", repos=NULL, type="source");'
pip install HTSeq

### test run ### 
snakemake -j 1 -pk



