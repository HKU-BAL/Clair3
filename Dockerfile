FROM ubuntu:16.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8 PATH=/opt/clair3/bin:/opt/conda/bin:$PATH

# update ubuntu packages
RUN apt-get update --fix-missing && \
    yes|apt-get upgrade && \
    apt-get install -y \
        wget \
        bzip2 \
        make \
        gcc \
        g++ && \
    rm -rf /bar/lib/apt/lists/*

WORKDIR /opt/clair3
COPY . .

# install anaconda
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda && \
    rm Miniconda3-latest-Linux-x86_64.sh

# create conda environment
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda create -n clair3 python=3.6 -y

ENV PATH /opt/conda/envs/clair3/bin:$PATH
ENV CONDA_DEFAULT_ENV clair3


RUN /bin/bash -c "source activate clair3" && \
    conda install -c conda-forge pypy3.6 -y && \
    pypy3 -m ensurepip && \
    pypy3 -m pip install intervaltree==3.0.2 && \
    pypy3 -m pip install python-Levenshtein==0.12.0 mpmath==1.2.1 && \
    pip install tensorflow==2.2.0 && \
    pip install intervaltree==3.0.2  tensorflow-addons==0.11.2 tables==3.6.1 python-Levenshtein==0.12.0 && \
    conda install -c anaconda pigz==2.4 -y && \
    conda install -c conda-forge parallel=20191122 zstd=1.4.4 -y && \
    conda install -c conda-forge -c bioconda samtools=1.10 -y && \
    conda install -c conda-forge -c bioconda whatshap=1.0 -y && \
    conda install -c conda-forge boost=1.67.0 -y && \
    echo "source activate clair3" > ~/.bashrc

#cd Clair3/preprocess/realign
#g++ -std=c++14 -O2 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c  realigner.cpp
#g++ -std=c++11   -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp
#/bin/bash -c "source activate clair3" &&
# docker rmi -f clair3
#docker build -f /autofs/bal33/zxzheng/nas2/Clair_private.2to3/Clair3/Dockerfile -t clair3 .
# docker run -it  -v /autofs/bal33/zxzheng:/autofs/bal33/zxzheng  -v /mnt/bal36/zxzheng/testData:/mnt/bal36/zxzheng/testData clair3 /bin/bash
