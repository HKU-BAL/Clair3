## Illumina Variant Calling Quick Demo
Here is a quick demo for the Illumina NGS variant calling using GIAB HG003 chromosome 20 data. We provide docker pre-built image (**with root privileges**) and  anaconda virtual environment (**without root privileges**) quick demo.

```bash
Platform:    Illumina
Sample:      GIAB HG003
Coverage:    ~39x
Reference:   GRCh38
Aligner:     BWA-MEM ALT-Aware
Region:      chr20:100000-300000
Instruments: NovaSeq
```

### Download data

```bash
# Parameters
PLATFORM='ilmn'
INPUT_DIR="${HOME}/clair3_illumina_quickDemo"
OUTPUT_DIR="${INPUT_DIR}/output"

mkdir -p ${INPUT_DIR}
mkdir -p ${OUTPUT_DIR}

# Download quick demo data
# GRCh38 Reference
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/illumina/GRCh38_chr20.fa
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/illumina/GRCh38_chr20.fa.fai
# BAM chr20:100000-300000
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/illumina/HG003_chr20_demo.bam
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/illumina/HG003_chr20_demo.bam.bai
# GIAB Truth VCF and BED
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/illumina/HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/illumina/HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz.tbi
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/illumina/HG003_GRCh38_chr20_v4.2.1_benchmark_noinconsistent.bed

REF="GRCh38_chr20.fa"
BAM="HG003_chr20_demo.bam"
BASELINE_VCF_FILE_PATH="HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz"
BASELINE_BED_FILE_PATH="HG003_GRCh38_chr20_v4.2.1_benchmark_noinconsistent.bed"
OUTPUT_VCF_FILE_PATH="merge_output.vcf.gz"

CONTIGS="chr20"
START_POS=100000
END_POS=300000
echo -e "${CONTIGS}\t${START_POS}\t${END_POS}" > ${INPUT_DIR}/quick_demo.bed
```

### Option 1. Docker pre-built image

```bash
THREADS=4
cd ${OUTPUT_DIR}

cd ${OUTPUT_DIR}
# Run Clair3 using one command
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/${BAM} \
  --ref_fn=${INPUT_DIR}/${REF} \
  --threads=${THREADS} \
  --platform=${PLATFORM} \
  --model_path="/opt/models/${PLATFORM}" \
  --output=${OUTPUT_DIR} \
  --bed_fn=${INPUT_DIR}/quick_demo.bed
```

**Run hap.py for benchmarking (optional)**

```bash
docker run \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
    ${INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
    ${OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
    -f "${INPUT_DIR}/${BASELINE_BED_FILE_PATH}" \
    -r "${INPUT_DIR}/${REF}" \
    -o "${OUTPUT_DIR}/happy" \
    -l ${CONTIGS}:${START_POS}-${END_POS} \
    --engine=vcfeval \
    --threads="${THREADS}" \
    --pass-only
```

**Hap.py Expected output:**

|   Type    | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1-Score |
| :-------: | :------: | :------: | :------: | :-----------: | :--------------: | :-------------: |
| **INDEL** |    59    |    0     |    0     |     1.000     |      1.000       |      1.000      |
|  **SNP**  |   386    |    16    |    1     |     0.960     |      0.997       |      0.978      |

Run all commands above:

```bash
cd ${HOME}
wget "http://www.bio8.cs.hku.hk/clair3/demo/clair3_ilmn_quick_demo.sh"
chmod +x clair3_ilmn_quick_demo.sh
./clair3_ilmn_quick_demo.sh
```

Check the results using `less ${HOME}/clair3_illumina_quickDemo/output/merge_output.vcf.gz`

### Option 2. Anaconda virtual environment

**You can run the quick demo below using Google Colab notebook:**

 [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/HKU-BAL/Clair3/blob/main/colab/clair3_illumina_quick_demo.ipynb)

##### Step 1. Install Clair3 and download pre-trained model, using [Installation - Option 3](https://github.com/HKU-BAL/Clair3#option-3-build-an-anaconda-virtual-environment)

##### Step 2. Install [Boost Graph Library](https://www.boost.org/doc/libs/1_65_1/libs/graph/doc/index.html) for Illumina realignment process

```bash
# Activate Clair3
conda activate clair3

# Install boost library
conda install -c conda-forge boost=1.67.0 -y
echo "Environment:" ${CONDA_PREFIX}
# Make sure in Clair3 directory
cd Clair3

cd preprocess/realign
g++ -std=c++14 -O1 -shared -fPIC -o realigner ssw_cpp.cpp ssw.c realigner.cpp
g++ -std=c++11 -shared -fPIC -o debruijn_graph -O3 debruijn_graph.cpp -I ${CONDA_PREFIX}/include -L ${CONDA_PREFIX}/lib

```

**Step 3. Run Clair3 without root privileges**

```bash
cd Clair3
./run_clair3.sh \
  --bam_fn=${INPUT_DIR}/${BAM} \
  --ref_fn=${INPUT_DIR}/${REF} \
  --threads=${THREADS} \
  --platform=${PLATFORM} \
  --model_path=`pwd`"/models/${PLATFORM}" \
  --output=${OUTPUT_DIR} \
  --bed_fn=${INPUT_DIR}/quick_demo.bed
```

**Run hap.py without root privileges for benchmarking (optional)**

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n happy-env -c bioconda hap.py -y
conda install -c bioconda rtg-tools -y
conda activate happy-env

# Benchmark using hap.py
hap.py \
    ${INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
    ${OUTPUT_DIR}/${OUTPUT_VCF_FILE_PATH} \
    -f "${INPUT_DIR}/${BASELINE_BED_FILE_PATH}" \
    -r "${INPUT_DIR}/${REF}" \
    -o "${OUTPUT_DIR}/happy" \
    -l ${CONTIGS}:${START_POS}-${END_POS} \
    --engine=vcfeval \
    --threads="${THREADS}" \
    --pass-only
```

**Hap.py Expected output:**

|   Type    | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1-Score |
| :-------: | :------: | :------: | :------: | :-----------: | :--------------: | :-------------: |
| **INDEL** |    59    |    0     |    0     |     1.000     |      1.000       |      1.000      |
|  **SNP**  |   386    |    16    |    1     |     0.960     |      0.997       |      0.978      |
