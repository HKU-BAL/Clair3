FROM ubuntu:16.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/clair/bin:/opt/conda/bin:$PATH

# update ubuntu packages
RUN apt-get update --fix-missing && \
    yes|apt-get upgrade && \
    apt-get install -y \
        wget \
        bzip2 \
        make \
        gcc \
        vcftools && \
    rm -rf /bar/lib/apt/lists/*

WORKDIR /opt/clair
COPY . .

# install anaconda
RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-x86_64.sh && \
    bash Anaconda3-2019.10-Linux-x86_64.sh -b -p /opt/conda && \
    rm Anaconda3-2019.10-Linux-x86_64.sh

# create conda environment
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n clair-env -c bioconda -y clair
RUN echo "source activate clair-env" > ~/.bashrc
ENV PATH /opt/conda/envs/clair-env/bin:$PATH
RUN /bin/bash -c ". activate clair-env && \
    pypy3 -m ensurepip && \
    pypy3 -m pip install --no-cache-dir intervaltree blosc"
