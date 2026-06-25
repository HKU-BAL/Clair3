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
# Pass --build-arg MODEL_CACHE_BUST=$(date +%s) to force a re-download without --no-cache.
ARG CLAIR3_MODELS_URL=https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/
# Recent ONT-provided (Rerio-converted) models bundled on top of the HKU baseline.
ARG CLAIR3_RERIO_MODELS_URL=https://www.bio8.cs.hku.hk/clair3/clair3_models_rerio_pytorch/
ARG CLAIR3_BUNDLED_RERIO_MODELS="r1041_e82_400bps_sup_v520 r1041_e82_400bps_hac_v520 r1041_e82_400bps_hac_v600"
ARG MODEL_CACHE_BUST=0
RUN set -eux; \
    echo "MODEL_CACHE_BUST=${MODEL_CACHE_BUST}"; \
    mkdir -p /opt/models /tmp/clair3-models; \
    base_url="${CLAIR3_MODELS_URL%/}"; \
    wget -r -np -nH --cut-dirs=2 -R "index.html*" -P /tmp/clair3-models "${base_url}/"; \
    if [ -d /tmp/clair3-models/clair3_models_pytorch ]; then \
        cp -a /tmp/clair3-models/clair3_models_pytorch/. /opt/models/; \
    else \
        cp -a /tmp/clair3-models/. /opt/models/; \
    fi; \
    rm -rf /tmp/clair3-models; \
    rerio_url="${CLAIR3_RERIO_MODELS_URL%/}"; \
    for m in ${CLAIR3_BUNDLED_RERIO_MODELS}; do \
        echo "Bundling ONT/Rerio model: ${m}"; \
        rm -rf "/opt/models/${m}"; mkdir -p "/opt/models/${m}"; \
        wget -r -np -nd -R "index.html*" -P "/opt/models/${m}" "${rerio_url}/${m}/"; \
        test -f "/opt/models/${m}/pileup.pt" -a -f "/opt/models/${m}/full_alignment.pt"; \
    done
