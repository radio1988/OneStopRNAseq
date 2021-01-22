# install anaconda
# https://docs.anaconda.com/anaconda/install/
# tested under conda version is 4.9.2

# create conda env called osr
conda clean --all
conda env create -n osr -f workflow/envs/env.yaml

# activate osr
conda activate osr

# install R packages into R in osr
# for DESeq2
Rscript -e "install.packages( c('BiocManager', 'PoiClaClu', 'rmarkdown', 'gridExtra'), repos='https://cloud.r-project.org')"
# for QoRTs
Rscript -e 'install.packages("http://hartleys.github.io/QoRTs/QoRTs_STABLE.tar.gz", repos=NULL, type="source");'
