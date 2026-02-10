# syntax=docker/dockerfile:1.4
FROM condaforge/miniforge3:latest

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
WORKDIR /opt/bin

# Step 1 in README: create a conda env with bioinformatics/build dependencies.
RUN mamba create -n clair3_v2 \
      -c conda-forge \
      -c bioconda \
      python=3.11 \
      samtools \
      whatshap \
      parallel \
      zstd \
      xz \
      zlib \
      bzip2 \
      automake \
      make \
      gcc \
      gxx \
      binutils \
      curl \
      pigz \
      boost-cpp \
      -y && \
    mamba clean --all -y

ENV CONDA_PREFIX=/opt/conda/envs/clair3_v2
ENV CONDA_DEFAULT_ENV=clair3_v2
ENV PATH=${CONDA_PREFIX}/bin:/opt/bin:/opt/conda/bin:${PATH}

# Step 2 & Step 4 in README: install uv, PyTorch (CPU), and Python dependencies.
RUN pip install --no-cache-dir uv && \
    uv pip install --python ${CONDA_PREFIX}/bin/python \
      torch torchvision \
      --index-url https://download.pytorch.org/whl/cpu && \
    uv pip install --python ${CONDA_PREFIX}/bin/python \
      numpy \
      h5py \
      hdf5plugin \
      numexpr \
      tqdm \
      cffi \
      torchmetrics

# Clair3 source tree
COPY . /opt/bin

# Build Clair3 native components.
RUN make PREFIX=${CONDA_PREFIX} PYTHON=${CONDA_PREFIX}/bin/python && \
    cd /opt/bin/preprocess/realign && \
    g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp && \
    g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp && \
    rm -rf /opt/bin/samtools-* /opt/bin/longphase-*

# Step 5 in README: install pypy3.11 and mpmath for pypy.
RUN wget -q https://downloads.python.org/pypy/pypy3.11-v7.3.20-linux64.tar.bz2 && \
    tar -xjf pypy3.11-v7.3.20-linux64.tar.bz2 && \
    rm pypy3.11-v7.3.20-linux64.tar.bz2 && \
    ln -sf /opt/bin/pypy3.11-v7.3.20-linux64/bin/pypy3 ${CONDA_PREFIX}/bin/pypy3 && \
    ln -sf /opt/bin/pypy3.11-v7.3.20-linux64/bin/pypy3 ${CONDA_PREFIX}/bin/pypy && \
    pypy3 -m ensurepip && \
    pypy3 -m pip install --no-cache-dir mpmath==1.2.1

# Step 6 in README: recursively download pre-trained PyTorch models from a directory URL.
ARG CLAIR3_MODELS_URL=https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/
RUN set -eux; \
    mkdir -p /opt/models /tmp/clair3-models; \
    base_url="${CLAIR3_MODELS_URL%/}"; \
    wget -r -np -nH --cut-dirs=2 -R "index.html*" -P /tmp/clair3-models "${base_url}/"; \
    if [ -d /tmp/clair3-models/clair3_models_pytorch ]; then \
        cp -a /tmp/clair3-models/clair3_models_pytorch/. /opt/models/; \
    else \
        cp -a /tmp/clair3-models/. /opt/models/; \
    fi; \
    rm -rf /tmp/clair3-models
