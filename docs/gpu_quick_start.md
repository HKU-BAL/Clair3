# Clair3 GPU quick start guide

Clair3 supports CUDA for calling and training on Linux, and Metal Performance Shaders (MPS) for calling on Apple Silicon. The implementation is native Rust; Python and TensorFlow are not required.

## NVIDIA Docker

Install the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html), then build the image:

```bash
docker build -f Dockerfile.gpu -t clair3-gpu .
docker run --rm --gpus all clair3-gpu nvidia-smi
```

Run a variant-calling job with an absolute input and output path:

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"
THREADS=8
MODEL_NAME="[YOUR_MODEL_NAME]"

docker run -it --gpus all \
  -u $(id -u):$(id -g) \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  clair3-gpu \
  /opt/clair3/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \
  --ref_fn=${INPUT_DIR}/ref.fa \
  --threads=${THREADS} \
  --platform="ont" \
  --model_path="/opt/models/${MODEL_NAME}" \
  --use_longphase_for_intermediate_phasing \
  --use_gpu \
  --device="cuda:0" \
  --output=${OUTPUT_DIR}
```

The image contains the Rust executable, LibTorch, Gloo, and LongPhase. It does not install Python, PyPy, or TensorFlow.

## Build locally on Linux

Install Rust 1.88.0 or newer, CMake, a C++20-capable C/C++ compiler, and samtools. Download a CUDA 12.9 LibTorch distribution and set `LIBTORCH` to its root directory. Build the native executable with CUDA support:

```bash
git clone https://github.com/HKU-BAL/Clair3.git
cd Clair3
export LIBTORCH=/path/to/libtorch
cargo build --release --locked --features cuda

clair3 run \
  --bam_fn=input.bam \
  --ref_fn=ref.fa \
  --threads=8 \
  --platform=ont \
  --model_path=/path/to/model \
  --use_gpu \
  --device="cuda:0" \
  --output=output
```

`--device` accepts a CUDA device such as `cuda:0`. The default when `--use_gpu` is enabled is `cuda:0`.

## Apple Silicon

Install Rust, CMake, a C++20-capable C/C++ compiler, and samtools. The default native build downloads the CPU/MPS LibTorch distribution; no Python installation is required. The downloaded Apple Silicon wheel requires macOS 14 or newer. Intel Macs need a compatible local LibTorch through `LIBTORCH`, because official PyTorch 2.13 macOS wheels are ARM64-only:

```bash
git clone https://github.com/HKU-BAL/Clair3.git
cd Clair3
cargo build --release --locked

clair3 run \
  --bam_fn=input.bam \
  --ref_fn=ref.fa \
  --threads=8 \
  --platform=ont \
  --model_path=/path/to/model \
  --device="mps" \
  --output=output
```

LongPhase can be built with `make longphase` using an OpenMP-capable compiler. Apple Clang does not provide this by default; install GCC and select it with `LONGPHASE_CC` and `LONGPHASE_CXX` if needed. WhatsHap is optional for ONT and HiFi; LongPhase is used when WhatsHap is unavailable. Illumina phasing uses WhatsHap, or can be disabled with `--no_phasing_for_fa`.

## Distributed training

Native training uses C10d and does not require `torchrun` or Python. Start one process per rank with the standard distributed environment variables:

```bash
RANK=0 WORLD_SIZE=2 LOCAL_RANK=0 MASTER_ADDR=127.0.0.1 MASTER_PORT=29500 \
  clair3 Train --bin_fn /path/to/bins --ochk_prefix /path/to/output/pileup --pileup --use_gpu
RANK=1 WORLD_SIZE=2 LOCAL_RANK=1 MASTER_ADDR=127.0.0.1 MASTER_PORT=29500 \
  clair3 Train --bin_fn /path/to/bins --ochk_prefix /path/to/output/pileup --pileup --use_gpu
```

Use an appropriate GPU index and distributed backend for the target system. Gloo is vendored and built automatically by the default feature set.
