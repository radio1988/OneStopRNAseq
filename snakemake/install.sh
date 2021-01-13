# currently, only Linux system tested
# if use assembly mode, need to have 160GB RAM, otherwise, 34GB is enough

# install anaconda
# Please follow instructions from https://docs.anaconda.com/anaconda/install/

echo '>>> cleanup old anaconda packages, if any'
conda clean --all

echo '>>> create conda env called osr'
conda env create -n osr -f workflow/envs/osr.yaml

echo '>>> activate osr'
conda activate osr

echo '>>> install R packages into R in osr'
echo '>> for DESeq2'
Rscript -e "install.packages( c('BiocManager', 'ggrepel', 'gdata', 'plyr', 'pheatmap', 'PoiClaClu', 'gridExtra', 'rmarkdown', 'ashr'), repos='https://cloud.r-project.org')"
Rscript -e "BiocManager::install('scater')"
echo '>> for QoRTs'
Rscript -e 'install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz", repos=NULL, type="source");'
echo '>> for DEXSeq'
pip install HTSeq

### test run ### 
snakemake -j 1 -pk



