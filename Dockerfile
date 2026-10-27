# syntax=docker/dockerfile:1.4
FROM rust:1.88.0-bookworm AS build

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential cmake clang libclang-dev pkg-config curl ca-certificates unzip \
    libssl-dev zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /opt
RUN curl -fLsS -o /tmp/libtorch.zip \
    https://download.pytorch.org/libtorch/cpu/libtorch-shared-with-deps-2.13.0%2Bcpu.zip && \
    unzip -q /tmp/libtorch.zip && rm /tmp/libtorch.zip
ENV LIBTORCH=/opt/libtorch
WORKDIR /opt/clair3-src
COPY Cargo.toml Cargo.lock build.rs ./
COPY src ./src
COPY csrc ./csrc
COPY third_party ./third_party
RUN --mount=type=cache,target=/usr/local/cargo/registry \
    --mount=type=cache,target=/opt/clair3-src/target \
    cargo build --release --locked

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates curl wget libgomp1 libstdc++6 samtools tabix zstd lz4 \
    && rm -rf /var/lib/apt/lists/*
COPY --from=build /opt/libtorch/lib /opt/clair3/lib/clair3
COPY --from=build /opt/clair3-src/target/release/clair3 /opt/clair3/bin/clair3
COPY run_clair3.sh run_clair3.py clair3.py /opt/clair3/
COPY scripts/clair3.sh /opt/clair3/scripts/clair3.sh
RUN chmod +x /opt/clair3/bin/clair3 /opt/clair3/run_clair3.sh /opt/clair3/run_clair3.py /opt/clair3/clair3.py && \
    curl -fLsS -o /tmp/longphase.tar.xz https://github.com/twolinin/longphase/releases/download/v1.7.3/longphase_linux-x64.tar.xz && \
    tar -xJf /tmp/longphase.tar.xz -C /opt/clair3/bin && \
    mv /opt/clair3/bin/longphase_linux-x64 /opt/clair3/bin/longphase && \
    chmod +x /opt/clair3/bin/longphase && rm /tmp/longphase.tar.xz
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
