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

RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('DESeq2')"