# syntax=docker/dockerfile:1.4
FROM rust:1.88.0-bookworm AS build

ARG TARGETARCH
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake clang libclang-dev pkg-config curl ca-certificates unzip \
    autoconf automake libtool \
    libssl-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    python3 python3-dev python3-pip samtools tabix \
    && rm -rf /var/lib/apt/lists/*
RUN python3 -m pip install --no-cache-dir --break-system-packages cmake==3.31.6 && \
    cmake --version
WORKDIR /opt
RUN if [ "${TARGETARCH}" = "arm64" ]; then \
        python3 -m pip install --no-cache-dir --break-system-packages \
            --index-url https://download.pytorch.org/whl/cpu \
            --extra-index-url https://pypi.org/simple torch==2.13.0 && \
        ln -s "$(python3 -c 'import torch; print(torch.__path__[0])')" /opt/libtorch; \
    else \
        curl -fLsS -o /tmp/libtorch.zip \
            https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-2.13.0%2Bcpu.zip && \
        unzip -q /tmp/libtorch.zip && rm /tmp/libtorch.zip; \
    fi && \
    mkdir -p /opt/clair3-libs && cp -aL /opt/libtorch/lib/. /opt/clair3-libs/
ENV LIBTORCH=/opt/libtorch
WORKDIR /opt/clair3-src
COPY Cargo.toml Cargo.lock build.rs ./
COPY src ./src
COPY csrc ./csrc
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/opt/clair3-src/target \
    cargo build --release --locked && \
    mkdir -p /opt/clair3-bin && cp target/release/clair3 /opt/clair3-bin/
RUN if [ "${TARGETARCH}" = "arm64" ]; then \
        curl -fLsS -o /tmp/longphase.tar.gz \
            https://github.com/twolinin/longphase/archive/refs/tags/v1.7.3.tar.gz && \
        tar -xzf /tmp/longphase.tar.gz -C /tmp && \
        cd /tmp/longphase-1.7.3 && autoreconf -i && ./configure && make -j2 && \
        cp longphase /opt/longphase; \
    else \
        curl -fLsS -o /tmp/longphase.tar.xz \
            https://github.com/twolinin/longphase/releases/download/v1.7.3/longphase_linux-x64.tar.xz && \
        tar -xJf /tmp/longphase.tar.xz -C /tmp && \
        cp /tmp/longphase_linux-x64 /opt/longphase; \
    fi

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl wget libgomp1 libstdc++6 samtools tabix zstd lz4 \
    build-essential python3-dev python3-pip \
    && rm -rf /var/lib/apt/lists/* \
    && python3 -m pip install --no-cache-dir --break-system-packages whatshap==2.8
COPY --from=build /opt/clair3-libs /opt/clair3/lib/clair3
COPY --from=build /opt/clair3-bin/clair3 /opt/clair3/bin/clair3
COPY --from=build /opt/longphase /opt/clair3/bin/longphase
COPY run_clair3.sh run_clair3.py clair3.py /opt/clair3/
COPY scripts/clair3.sh /opt/clair3/scripts/clair3.sh
RUN chmod +x /opt/clair3/bin/clair3 /opt/clair3/bin/longphase /opt/clair3/run_clair3.sh /opt/clair3/run_clair3.py /opt/clair3/clair3.py
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
ENV LD_LIBRARY_PATH=/opt/clair3/lib/clair3 \
    PATH=/opt/clair3/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
WORKDIR /data
ENTRYPOINT ["/opt/clair3/run_clair3.sh"]
