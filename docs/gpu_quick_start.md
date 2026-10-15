<!-- Modifications Copyright (C) 2026 Advanced Micro Devices, Inc. All rights reserved. -->
# Clair3 GPU quick start guide

Starting from v1.2, Clair3 natively supports NVIDIA GPU acceleration. Using a single GPU, Clair3 can complete an ONT WGS 30x whole-genome variant calling in ~20 minutes on a Linux server with 32 CPU threads and an NVIDIA GeForce RTX 4090.

The quickest way to run Clair3 on GPU is the pre-built Docker / Singularity image `hkubal/clair3:gpu` (built on CUDA 12.1, bundled with all pre-trained models).

For AMD Instinct GPUs (ROCm on Linux), build the ROCm image from `Dockerfile.amd` — see [Option 4](#option-4-amd-instinct-rocm-on-linux) below; the same pre-trained models apply.

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

## Option 4. AMD Instinct (ROCm on Linux)

Clair3's GPU path runs unchanged on AMD Instinct GPUs through ROCm PyTorch: `--device=cuda:0`, `torch.cuda.*`, and `CUDA_VISIBLE_DEVICES` all map onto HIP, so no Clair3 options change.

Build the ROCm GPU image from the repository. `Dockerfile.amd` uses a plain Ubuntu base (Python 3.12+) and installs PyTorch from TheRock's multi-arch ROCm index (`https://stable.repo.amd.com/rocm/whl-next/`, ROCm 10.0); the wheels bundle their own ROCm SDK, so no system ROCm is required. The build defaults to `gfx950` (MI350X/MI355X); target another GPU with `--build-arg GFX_ARCH=<arch>` (e.g. `gfx942` for MI300X). The base image — and thus the Python version — is also a build arg: `--build-arg BASE_IMAGE=ubuntu:26.04` targets Python 3.14 (default `ubuntu:24.04` = 3.12), subject to matching TheRock wheels. For ROCm 10.1, switch the index to `https://nightly.repo.amd.com/rocm/whl-next/` and add `--pre`. PyTorch is installed into a plain venv (not conda) so that MIOpen's pileup RNN kernel JIT-compiles correctly on gfx9xx.

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
 --cap-add=SYS_PTRACE --security-opt seccomp=unconfined --ipc=host \
 -v ${INPUT_DIR}:${INPUT_DIR} \
 -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
 clair3:rocm \
 /opt/bin/run_clair3.sh \
 --bam_fn=${INPUT_DIR}/input.bam \
 --ref_fn=${INPUT_DIR}/ref.fa \
 --threads=${THREADS} \
 --platform=ont \ ## {ont,hifi,ilmn}
 --model_path=/opt/models/${MODEL_NAME} \
 --output=${OUTPUT_DIR} \
 --use_gpu
```

**Notes**

- Select specific GPUs with `--device=cuda:0,1` (Clair3) or by exporting `CUDA_VISIBLE_DEVICES` / `HIP_VISIBLE_DEVICES` before the run.
- For a host (no-Docker) install on AMD Instinct, install PyTorch into a **venv** (Python 3.12+, not conda) from TheRock's ROCm index, e.g. `python3 -m venv venv && . venv/bin/activate && pip install --index-url https://stable.repo.amd.com/rocm/whl-next/ "torch[device-gfx950]"` (swap `gfx950` for your GPU arch), then build Clair3 with `make PREFIX=$VIRTUAL_ENV PYTHON=$VIRTUAL_ENV/bin/python`. On ROCm 10 a conda-prefixed install breaks MIOpen's pileup RNN-kernel JIT; a venv avoids it.

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
