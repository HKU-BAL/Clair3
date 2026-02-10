

# Clair3 GPU quick start guide

Starting from v1.2, Clair3 natively supports using GPU with CUDA or Apple Silicon to speed up calling. In Linux with CUDA, Clair3 automatically uses as many GPUs that are exposed to it. For **Apple Silicons** (we tested in M1/M2/M3), Clair3 uses PyTorch with Metal (MPS) for GPU acceleration. Using a single GPU or Apple Silicon, Clair3 can complete a ONT WGS 30x whole-genome variant calling in ~20 minutes, using either a Linux server with 32 CPU threads and an NVIDIA GeForce RTX 4090 GPU, or a Mac Studio with a 32-core M3 Ultra.

## Installation on Linux 

### Option 1.  NVIDIA Docker pre-built image

To use Clair3 in GPU-accelerated mode within Docker, install the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) to enable GPU support in Docker containers. Then, test if Docker can access your GPU by running:

```shell
docker run --rm --gpus all hkubal/clair3-gpu:latest nvidia-smi
```

Expected output should display your GPU information similar to:

```shell
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 555.42.02              Driver Version: 555.42.02      CUDA Version: 12.5     |
|-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA GeForce RTX 4090        On  |   00000000:CA:00.0 Off |                  Off |
|  0%   32C    P8             27W /  450W |       4MiB /  24564MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+
```

With a correct environment, you can run Clair3 GPU using a single command.

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
MODEL_NAME="[YOUR_MODEL_NAME]"         # e.g. r1041_e82_400bps_sup_v500

#--gpus is 1-index in docker
docker run -it --gpus all \
  -u $(id -u):$(id -g) \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3-gpu:latest \
  python3 /opt/bin/run_clair3.py \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/${MODEL_NAME}" \
  --use_longphase_for_intermediate_phasing \
  --use_gpu \                          ## use gpu for variant calling
  --device='cuda:0,1' \                ## select the GPU device (0-index) for calling, default: use all visual GPUs
  --output=${OUTPUT_DIR}               ## absolute output path prefix 
```

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR`. 

### Option 2.  Bioconda

```bash
# make sure channels are added in conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment named "clair3"
# install clair3 and cudatoolkit cudnn using one command in the machine with GPU 
conda create -n clair3 -c bioconda clair3 cudatoolkit cudnn python=3.9.0 -y
conda activate clair3

# run clair3 like this afterward
python3 run_clair3.py \
  --bam_fn=input.bam \                 ## change your bam file name here
  --ref_fn=ref.fa \                    ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="${CONDA_PREFIX}/bin/models/${MODEL_NAME}" \ 
  --use_gpu \                          ## use gpu for variant calling
  --use_longphase_for_intermediate_phasing \
  --device='cuda:0,1' \                ## select the GPU device for calling, default: use all Nvidia GPUs
  --output=${OUTPUT_DIR}               ## output path prefix 
```

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
