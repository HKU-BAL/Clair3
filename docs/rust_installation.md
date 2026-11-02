# Installing Clair3-Rust

This guide installs the `clair3-rust` branch with native Rust calling and unmodified upstream third-party tools. CPU is the default backend. MPS and CUDA are optional; see [GPU quick start](gpu_quick_start.md).

## 1. Requirements

- Rust 1.88 or newer, Cargo, CMake 3.26 or newer, a C++20 compiler, curl, and Make. The Docker build installs CMake 3.31.6.
- Python 3 and pip for WhatsHap phasing. Python 3.11 is used in the tested environment. Python development headers are needed to build the native read selector.
- samtools 1.15.1 or newer, bgzip, and tabix.
- Network access for the initial download of Rust crates, LibTorch, Gloo, WhatsHap, and optional models.
- A model directory containing both `pileup.pt` and `full_alignment.pt`. TensorFlow checkpoints require separate conversion.

Core calling and training do not require Python, PyPy, or TensorFlow. Python is used by WhatsHap and by the macOS PyTorch-wheel build path.

## 2. Create an environment

Mamba or Conda can provide the compiler and bioinformatics tools:

```bash
mamba create -n clair3_v2 \
  -c conda-forge \
  -c bioconda \
  python=3.11 \
  pip \
  samtools \
  longphase \
  zstd \
  xz \
  cmake \
  make \
  clang \
  curl \
  -y
mamba activate clair3_v2
```

Install Rust if it is not already available:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
rustup toolchain install 1.88.0
rustup default 1.88.0
```

## 3. Check out the branch

```bash
git clone --branch clair3-rust --single-branch https://github.com/HKU-BAL/Clair3.git
cd Clair3
export CLAIR3_PATH=$(pwd)
```

## 4. Install upstream Python tools

```bash
make python_tools PYTHON=python3
```

This installs unmodified `whatshap==2.8`. On macOS, it also installs the official `torch==2.13.0` wheel when no local `LIBTORCH` is selected.

The phasing optimizations are enabled only for the tested WhatsHap 2.8 API/ABI. Other versions use the stock phasing path. The upstream DP solver is not replaced or modified.

## 5. Build and install Clair3

Install into the active Conda environment:

```bash
make clair3 PREFIX="${CONDA_PREFIX}" PYTHON=python3
```

This installs:

- `${CONDA_PREFIX}/bin/clair3`
- LibTorch libraries under `${CONDA_PREFIX}/lib/clair3`

The installed binary has a relative runtime-library path, so it can be moved together with that `lib/clair3` directory.

The default `make` target also builds LongPhase and installs the Python tools:

```bash
make PREFIX="${CONDA_PREFIX}" PYTHON=python3
```

LongPhase is optional when WhatsHap is installed. On macOS, building it from source requires an OpenMP-capable compiler; Apple Clang alone is insufficient. Select a suitable GCC installation if necessary:

```bash
make longphase PREFIX="${CONDA_PREFIX}" \
  LONGPHASE_CC=gcc-15 LONGPHASE_CXX=g++-15
```

### Use an existing LibTorch installation

Set `LIBTORCH` to a directory containing `include` and `lib`:

```bash
make clair3 PREFIX="${CONDA_PREFIX}" \
  LIBTORCH=/path/to/libtorch PYTHON=python3
```

If the distribution uses the old C++ ABI, also set `LIBTORCH_CXX11_ABI=0`. Official PyTorch 2.13 macOS wheels are ARM64-only, so Intel Macs need a compatible local LibTorch installation.

For a direct Cargo build with an official Python PyTorch wheel:

```bash
LIBTORCH_USE_PYTORCH=1 PYTHON=python3 cargo build --release --locked
```

On Linux, the default Cargo build can download LibTorch automatically. CPU and CUDA use separate Docker configurations; CUDA requires `--features cuda` for a local build.

## 6. Download models

Download all pre-trained PyTorch models into the installation prefix:

```bash
make models PREFIX="${CONDA_PREFIX}"
```

Alternatively, download a specific model from the [Clair3 model directory](https://www.bio8.cs.hku.hk/clair3/clair3_models_pytorch/). A usable model directory contains both:

```text
pileup.pt
full_alignment.pt
```

For example, a HiFi Revio model is typically located at `${CONDA_PREFIX}/bin/models/hifi_revio` after `make models`.

## 7. Verify the installation

```bash
clair3 --version
whatshap --version
samtools --version
```

Expected major versions in the tested installation are Clair3 2.0, WhatsHap 2.8, and samtools 1.15.1 or newer.

Continue with the [Clair3-Rust quick test](rust_quick_test.md).

## Linux build verification

Build and test the CPU image with a local Docker daemon:

```bash
docker build -f Dockerfile -t clair3-rust:local .
docker run --rm --entrypoint /opt/clair3/bin/clair3 clair3-rust:local --version
```

The image supports Linux x86_64 and ARM64. On ARM64, it uses the official PyTorch wheel and builds LongPhase from the unmodified upstream source; x86_64 retains the prebuilt LibTorch and LongPhase downloads. Both architectures were validated locally, with x86_64 tested under `linux/amd64` emulation on an ARM64 host. On another architecture, select x86_64 explicitly with `docker build --platform linux/amd64 ...` and `docker run --platform linux/amd64 ...`. Use the fixture and phasing commands in the [quick-test guide](rust_quick_test.md) to verify model inference and native phasing.

## Docker

The images install upstream WhatsHap and the required Python development headers for the native selector. They do not copy or patch vendored dependency sources. The tags below refer to locally built images; published Docker Hub images are not asserted to contain the `clair3-rust` branch.

Build a CPU image:

```bash
docker build -f Dockerfile -t clair3-rust:local .
```

Build a CUDA image on a host with the NVIDIA Container Toolkit:

```bash
docker build -f Dockerfile.gpu -t clair3-rust:cuda .
docker run --rm --gpus all --entrypoint nvidia-smi clair3-rust:cuda
```

Use `/opt/clair3/run_clair3.sh` as the container entry point. Pre-trained models are installed under `/opt/models`. Use absolute host paths when mounting input and output directories.

## Dependency policy

- WhatsHap comes from its official Python package and is used without source changes.
- The native selector is first-party Clair3 code. It builds against unmodified WhatsHap 2.8 headers downloaded from PyPI with SHA-256 `073ee7946d563e0125dbce58d9fc86a7a3fae377f27859dacb398309374b2cbe`. Linux builds include old- and new-C++-ABI variants, and the loader selects one using symbols exported by the installed core. Variant position, allele, and quality fields are deserialized in an explicit sequence for consistent behavior across compilers.
- `torch-sys` is the official crates.io package recorded in `Cargo.lock`.
- Gloo is downloaded from upstream commit `bcd1672ee07538123ea8f4fac76832efc58fb8ef`.
- LibTorch and LongPhase retain their upstream download paths.
- The CUDA image installs unmodified official NVIDIA cuSPARSE-Lt 0.8.1, NCCL 2.29.7, NVSHMEM 3.4.5, and CUPTI packages required by LibTorch 2.13, plus an unmodified upstream samtools 1.24 build. Custom CUDA library layouts can be supplied through `CLAIR3_CUDA_LIBRARIES` as a colon-separated list during the build.

No downloaded third-party source is patched by the build.

## Troubleshooting

- **Ubuntu samtools is too old:** Ubuntu 22.04 packages samtools 1.13, below the required 1.15.1. Use the Conda package or build a newer unmodified upstream release, such as samtools 1.24, with ncurses and libdeflate support. The local Docker image uses Debian's samtools 1.16.1 package.
- **Missing Python development headers:** install the package appropriate for the environment, such as `python3-dev` on Debian/Ubuntu, then rebuild. Without them, native selector acceleration is unavailable and the pipeline falls back to the installed implementation.
- **LibTorch not found:** select a complete distribution through `LIBTORCH`, or use the official Python wheel with `LIBTORCH_USE_PYTORCH=1`.
- **Runtime library errors:** keep `bin/clair3` and `lib/clair3` together, or set `LD_LIBRARY_PATH` on Linux or `DYLD_LIBRARY_PATH` on macOS to the installed `lib/clair3` directory.
- **WhatsHap is unavailable:** ONT and HiFi can use LongPhase. Illumina requires WhatsHap for phasing unless `--no_phasing_for_fa` is used.
- **Compiler errors while building LongPhase:** use an OpenMP-capable compiler, or use an existing LongPhase installation instead.
