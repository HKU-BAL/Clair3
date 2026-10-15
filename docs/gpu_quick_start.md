# Clair3 GPU quick start guide

Starting from v1.2, Clair3 natively supports NVIDIA GPU acceleration. Using a single GPU, Clair3 can complete an ONT WGS 30x whole-genome variant calling in ~20 minutes on a Linux server with 32 CPU threads and an NVIDIA GeForce RTX 4090.

The quickest way to run Clair3 on GPU is the pre-built Docker / Singularity image `hkubal/clair3:gpu` (built on CUDA 12.1, bundled with all pre-trained models).

For AMD Instinct GPUs (ROCm on Linux), the quickest path is the prebuilt Clair3 from Bioconda paired with ROCm PyTorch — see [Option 4](#option-4-amd-instinct-rocm-prebuilt-clair3-from-bioconda) below; nothing is compiled. To build the current source tree instead, use `Dockerfile.amd` — see [Option 5](#option-5-amd-instinct-rocm-build-from-source).

## Requirements

- NVIDIA driver ≥ 530.30.02 on the host.
  - CUDA 12.1 is chosen for broad driver compatibility; newer drivers (including those for RTX 50-series / Blackwell) are backward compatible and work with this image.
- For Docker: [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) installed on the host.
- For Singularity: `--nv` support (no NVIDIA Container Toolkit required).

Verify GPU passthrough works:

```shell
docker run --rm --gpus all hkubal/clair3:gpu nvidia-smi
```

Expected output shows your GPU:

```
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 580.82.09              Driver Version: 580.82.09      CUDA Version: 12.1     |
|-----------------------------------------+------------------------+----------------------+
|   0  NVIDIA GeForce RTX 4090        On  |   00000000:CA:00.0 Off |                  Off |
+-----------------------------------------+------------------------+----------------------+
```

## Option 1. Docker (NVIDIA GPU)

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input  (absolute path)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
MODEL_NAME="[YOUR_MODEL_NAME]"         # e.g. r1041_e82_400bps_sup_v500

docker run -it --gpus all \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:gpu \
  /opt/bin/run_clair3.sh \
    --bam_fn=${INPUT_DIR}/input.bam \
    --ref_fn=${INPUT_DIR}/ref.fa \
    --threads=${THREADS} \
    --platform=ont \                       ## {ont,hifi,ilmn}
    --model_path=/opt/models/${MODEL_NAME} \
    --output=${OUTPUT_DIR} \
    --use_gpu                              ## enable GPU-accelerated calling
```

**Notes**

- Absolute paths are required for `INPUT_DIR` and `OUTPUT_DIR`.
- Select specific GPUs with `--gpus '"device=0,1"'` (Docker) and `--device=cuda:0,1` (Clair3). By default all visible GPUs are used.
- `python3 /opt/bin/run_clair3.py` can replace `/opt/bin/run_clair3.sh` in the command above.

## Option 2. Singularity (NVIDIA GPU)

```bash
singularity pull docker://hkubal/clair3:gpu

singularity exec --nv --cleanenv --env TMPDIR=/tmp \
  -B ${INPUT_DIR},${OUTPUT_DIR} \
  clair3_gpu.sif \
  /opt/bin/run_clair3.sh \
    --bam_fn=${INPUT_DIR}/input.bam \
    --ref_fn=${INPUT_DIR}/ref.fa \
    --threads=${THREADS} \
    --platform=ont \                       ## {ont,hifi,ilmn}
    --model_path=/opt/models/${MODEL_NAME} \
    --output=${OUTPUT_DIR} \
    --use_gpu
```

**Notes**

- `--nv` injects the host NVIDIA driver and libraries into the container (equivalent of Docker's `--gpus all`); no NVIDIA Container Toolkit required.
- `--cleanenv --env TMPDIR=/tmp` avoids `parallel` failing when the host `TMPDIR` points to a path not visible inside the container.

## Option 3. Step-by-step (Conda with CUDA PyTorch)

Use this option when the Docker/Singularity image does not fit your environment (older driver, custom CUDA runtime, HPC cluster without container toolkit, etc.).

```bash
# Step 1: create conda environment
mamba create -n clair3_v2 -c conda-forge -c bioconda -y \
  python=3.11 samtools whatshap parallel \
  zstd xz zlib bzip2 automake make gcc gxx curl pigz
mamba activate clair3_v2
pip install uv

# Step 2: install PyTorch with CUDA support
# Pick the right CUDA version for your driver from https://pytorch.org/get-started/locally/
uv pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121

# Step 3: build Clair3
git clone https://github.com/HKU-BAL/Clair3.git
cd Clair3
export CLAIR3_PATH=$(pwd)
uv pip install numpy h5py hdf5plugin numexpr tqdm cffi torchmetrics
make PREFIX=${CONDA_PREFIX}

# Step 4: install PyPy (see main README for the full PyPy step)

# Step 5: run Clair3 with GPU
python3 ${CLAIR3_PATH}/run_clair3.py \
  --bam_fn=input.bam \
  --ref_fn=ref.fa \
  --threads=${THREADS} \
  --platform=ont \
  --model_path=${CLAIR3_PATH}/models/${MODEL_NAME} \
  --use_gpu \
  --device=cuda:0 \
  --output=${OUTPUT_DIR}
```

See the [main README](../README.md#option-4-step-by-step-conda) for the complete step-by-step instructions (PyPy install, model download, etc.).

## Option 4. AMD Instinct (ROCm): prebuilt Clair3 from Bioconda

The quickest way onto AMD Instinct. Installs the prebuilt Clair3 from Bioconda — `libclair3.so`, the pre-trained models, `longphase` and `run_clair3.sh` — and pairs it with ROCm PyTorch. Nothing is compiled, so the image builds in ~2 minutes.

This works because Clair3's CFFI layer is plain C and links nothing from PyTorch, so the prebuilt package runs unchanged against a ROCm torch. All GPU work goes through the `torch.cuda.*` API, which ROCm maps onto HIP, so `--device=cuda:0` and `CUDA_VISIBLE_DEVICES` behave exactly as on NVIDIA.

**Requirements**

- The `amdgpu` kernel driver on the host, with the GPU passed through at `docker run`. No system ROCm is needed — the wheels bundle their own ROCm SDK and the container never uses `/opt/rocm`.
- A Clair3 release whose GPU scheduler recognises AMD devices. If `--use_gpu` reports that no GPU is available, your Bioconda build predates it — use [Option 5](#option-5-amd-instinct-rocm-build-from-source) instead.

```dockerfile
FROM python:3.11-slim

ARG GFX_ARCH=gfx950
ARG ROCM_VERSION=10.0.0
ARG ROCM_INDEX_URL=https://stable.repo.amd.com/rocm/whl-next/

# g++ is required: MIOpen JIT-compiles the pileup kernels at runtime.
# The remaining libraries are what the prebuilt libclair3.so links against.
RUN apt-get update && apt-get install -y --no-install-recommends \
      ca-certificates curl bzip2 g++ libgomp1 \
      libcurl4 libssl3 libdeflate0 liblzma5 libbz2-1.0 zlib1g && \
    rm -rf /var/lib/apt/lists/*

# Prebuilt Clair3 from Bioconda; drop the CPU PyTorch it pulls in.
RUN curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest \
      | tar -xj -C /usr/local bin/micromamba && \
    micromamba create -y -p /opt/clair3 -c conda-forge -c bioconda clair3 && \
    micromamba remove -y -p /opt/clair3 --force pytorch pytorch-cpu libtorch && \
    micromamba clean -a -y

# ROCm PyTorch into a venv -- deliberately NOT the conda prefix, see the notes below.
RUN python -m venv /opt/venv && \
    /opt/venv/bin/pip install --no-cache-dir --upgrade pip && \
    /opt/venv/bin/pip install --no-cache-dir --index-url ${ROCM_INDEX_URL} \
      "torch[device-${GFX_ARCH}]" "rocm==${ROCM_VERSION}" && \
    /opt/venv/bin/pip install --no-cache-dir numpy cffi

# venv Python first; Clair3's scripts, models and tools come from the Bioconda prefix.
ENV PATH=/opt/venv/bin:/opt/clair3/bin:${PATH}
ENV PYTHONPATH=/opt/clair3/bin
```

```bash
docker build -t clair3:rocm-prebuilt .
```

Run with ROCm device passthrough (not `--gpus all`):

```bash
docker run -it \
  --device=/dev/kfd --device=/dev/dri --group-add=video \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  clair3:rocm-prebuilt \
  run_clair3.sh \
    --bam_fn=${INPUT_DIR}/input.bam \
    --ref_fn=${INPUT_DIR}/ref.fa \
    --threads=${THREADS} \
    --platform=ont \                       ## {ont,hifi,ilmn}
    --model_path=/opt/clair3/bin/models/${MODEL_NAME} \
    --output=${OUTPUT_DIR} \
    --use_gpu
```

**Notes**

- **Install ROCm PyTorch into a venv, not into the conda prefix.** The same ROCm 10 wheels that work in a venv fail inside a conda prefix, with `device kernel image is invalid` (`hipErrorInvalidImage`) on every kernel launch. That is why PyTorch lives in `/opt/venv` above while Clair3 itself stays in the Bioconda prefix.
- Target another GPU with `--build-arg GFX_ARCH=gfx942` (MI300X); the default `gfx950` covers MI350X/MI355X. Pass a comma-separated list to bundle several architectures.
- Models ship under `/opt/clair3/bin/models`; list them with `ls /opt/clair3/bin/models`.
- Select specific GPUs with `--device=cuda:0,1` (Clair3) or by exporting `CUDA_VISIBLE_DEVICES` / `HIP_VISIBLE_DEVICES` before the run.
- `--cap-add=SYS_PTRACE --security-opt seccomp=unconfined --ipc=host` are **not** required for variant calling. Add them only for profiling or multi-GPU workloads that need host shared memory; they weaken container isolation.
- `miopenStatusUnknownError` during calling means `g++` is missing from the image.
- You get the latest Bioconda release; for anything newer, build from source with [Option 5](#option-5-amd-instinct-rocm-build-from-source).

## Option 5. AMD Instinct (ROCm): build from source

Use this when you need the current source tree rather than the latest Bioconda release.

`Dockerfile.amd` uses a plain Ubuntu base (Python 3.12+) and installs PyTorch from TheRock's multi-arch ROCm index (`https://stable.repo.amd.com/rocm/whl-next/`, ROCm 10.0). The build defaults to `gfx950` (MI350X/MI355X); target another GPU with `--build-arg GFX_ARCH=<arch>` (e.g. `gfx942` for MI300X). The base image — and thus the Python version — is also a build arg: `--build-arg BASE_IMAGE=ubuntu:26.04` targets Python 3.14 (default `ubuntu:24.04` = 3.12), subject to matching TheRock wheels. For ROCm 10.1, switch the index to `https://nightly.repo.amd.com/rocm/whl-next/` and add `--pre`. PyTorch is installed into a plain venv rather than a conda prefix, for the reason given under Option 4.

**Requirements**

- The `amdgpu` kernel driver on the host, with the GPU passed through at `docker run`. No system ROCm is needed — the wheels bundle their own ROCm SDK.

```bash
# default (gfx950)
docker build -f Dockerfile.amd -t clair3:rocm .

# or target another architecture, e.g. MI300X
docker build -f Dockerfile.amd --build-arg GFX_ARCH=gfx942 -t clair3:rocm .
```

Run with ROCm device passthrough (not `--gpus all`):

```bash
docker run -it \
  --device=/dev/kfd --device=/dev/dri --group-add=video \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  clair3:rocm \
  /opt/bin/run_clair3.sh \
    --bam_fn=${INPUT_DIR}/input.bam \
    --ref_fn=${INPUT_DIR}/ref.fa \
    --threads=${THREADS} \
    --platform=ont \                       ## {ont,hifi,ilmn}
    --model_path=/opt/models/${MODEL_NAME} \
    --output=${OUTPUT_DIR} \
    --use_gpu
```

**Notes**

- Select specific GPUs with `--device=cuda:0,1` (Clair3) or by exporting `CUDA_VISIBLE_DEVICES` / `HIP_VISIBLE_DEVICES` before the run.
- `--cap-add=SYS_PTRACE --security-opt seccomp=unconfined --ipc=host` are **not** required for variant calling; add them only for profiling or multi-GPU workloads that need host shared memory.
- For a host (no-Docker) install on AMD Instinct, install PyTorch into a **venv** (Python 3.11+, not conda) from TheRock's ROCm index, e.g. `python3 -m venv venv && . venv/bin/activate && pip install --index-url https://stable.repo.amd.com/rocm/whl-next/ "torch[device-gfx950]"` (swap `gfx950` for your GPU arch), then build Clair3 with `make PREFIX=$VIRTUAL_ENV PYTHON=$VIRTUAL_ENV/bin/python`. `g++` must be present for MIOpen's runtime kernel compilation.


<!--
## Installation on Apple Silicon

Install brew from https://brew.sh/ or using command below:

```shell
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Then install Clair3 envrionment step by step:

```shell
#step 1: install the envrionment using brew
brew install gnu-getopt bash llvm micromamba pypy3 samtools

#step 2: install PyTorch and other dependencies using mamba
mamba create -n clair3 python=3.11 autoconf automake zlib libdeflate cffi parallel -y
mamba activate clair3
python -m pip install torch==2.2.* torchvision==0.17.* torchaudio==2.2.*

#step 3: build the dependecies
git clone https://github.com/HKU-BAL/Clair3.git && cd Clair3
make PREFIX=${CONDA_PREFIX}

# Then run clair3 like this afterward
export PATH="/opt/homebrew/opt/gnu-getopt/bin:$PATH"
python3 run_clair3.py \
  --bam_fn=input.bam \                 ## change your bam file name here
  --ref_fn=ref.fa \                    ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="${CONDA_PREFIX}/bin/models/${MODEL_NAME}" \
  --use_longphase_for_intermediate_phasing \
  --output=${OUTPUT_DIR}               ## output path prefix
```
-->
