# install anaconda
# https://docs.anaconda.com/anaconda/install/
# make sure conda version is 4.9.2

# create conda env called osr
conda clean --all
conda env create -n osr -f workflow/envs/env.yaml

# activate osr
conda activate osr

# install R packages into R in osr
# for DESeq2 
Rscript -e "install.packages( c('BiocManager', 'PoiClaClu', 'rmarkdown'), repos='https://cloud.r-project.org')"

# scater can be skipped if change R code, just normalize manually, calculateTPM, calculateFPKM
Rscript -e "BiocManager::install('scater')" # takes a long time


# for QoRTs
Rscript -e 'install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz", repos=NULL, type="source");'
#pip install HTSeq
# test3:  conda install -c bioconda htseq # seem to work

# For SalmonTE
#pip install docopt

# for rMATS
# conda install -c bioconda rmats=4.1.0

### test run ### 
# snakemake -j 1 -pkn
