FROM ubuntu:latest
USER root

WORKDIR /app

RUN apt-get update && apt-get install -y curl
RUN apt-get -y install libcurl4-openssl-dev
RUN apt-get install -y wget 
RUN mkdir -p ~/miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
RUN chmod +x ~/miniconda3/miniconda.sh
RUN ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
ENV PATH=/root/miniconda3/bin:$PATH
RUN rm -rf ~/miniconda3/miniconda.sh
RUN ~/miniconda3/bin/conda init bash
RUN ~/miniconda3/bin/conda init zsh
RUN ~/miniconda3/bin/conda config --add channels conda-forge
RUN ~/miniconda3/bin/activate
RUN ~/miniconda3/bin/conda install -y r-base

# Installing R packages
RUN R -e "install.packages('BiocManager', repos = 'https://cloud.r-project.org')"
RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('ShortRead', force = TRUE)"
RUN R -e "BiocManager::install('dada2', force = TRUE)"
RUN R -e "BiocManager::install('ggplot2', force = TRUE)"
RUN R -e "BiocManager::install('phyloseq', force = TRUE)"
RUN R -e "BiocManager::install('Biostrings', force = TRUE)"
RUN R -e "install.packages('dplyr', repos = 'https://cloud.r-project.org')"
RUN R -e "BiocManager::install('easycsv', force = TRUE)"
RUN R -e "BiocManager::install('tidyverse', force = TRUE)"
#RUN apt-get update && apt-get install -y \
#    git \
#    libcurl4-openssl-dev \
#    libssl-dev \
#    libxml2-dev \
#    make
RUN ~/miniconda3/bin/conda install -y conda-forge::r-devtools
RUN R -e "devtools::install_github('mhahsler/rBLAST', force = TRUE)"