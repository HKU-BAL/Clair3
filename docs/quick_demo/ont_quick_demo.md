## ONT Variant Calling Quick Demo
Here is a quick demo for the Oxford Nanopore (ONT) variant calling using GIAB HG003 chromosome 20 data. We provide a Mamba/Conda virtual environment (**without root privileges**) and docker pre-built image (**with root privileges**) quick demo.

```bash
Platform:   ONT
Sample:     GIAB HG003
Coverage:   ~85x
Reference:  GRCh38_no_alt
Region:     chr20:100000-300000
Basecaller: Guppy 3.6.0
```

**Download data**

```bash
# Parameters
PLATFORM='ont'
INPUT_DIR="${HOME}/clair3_ont_quickDemo"
OUTPUT_DIR="${INPUT_DIR}/output"

mkdir -p ${INPUT_DIR}
mkdir -p ${OUTPUT_DIR}

# Download quick demo data
# GRCh38_no_alt reference
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/ont/GRCh38_no_alt_chr20.fa
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/ont/GRCh38_no_alt_chr20.fa.fai
# BAM chr20:100000-300000
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/ont/HG003_chr20_demo.bam
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/ont/HG003_chr20_demo.bam.bai
# GIAB Truth VCF and BED
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/ont/HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/ont/HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz.tbi
wget -P ${INPUT_DIR} http://www.bio8.cs.hku.hk/clair3/demo/quick_demo/ont/HG003_GRCh38_chr20_v4.2.1_benchmark_noinconsistent.bed

REF="GRCh38_no_alt_chr20.fa"
BAM="HG003_chr20_demo.bam"
BASELINE_VCF_FILE_PATH="HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz"
BASELINE_BED_FILE_PATH="HG003_GRCh38_chr20_v4.2.1_benchmark_noinconsistent.bed"
OUTPUT_VCF_FILE_PATH="merge_output.vcf.gz"

CONTIGS="chr20"
START_POS=100000
END_POS=300000
echo -e "${CONTIGS}\t${START_POS}\t${END_POS}" > ${INPUT_DIR}/quick_demo.bed
```

### Option 1. Mamba/Conda virtual environment

Install Clair3 following [Installation - Option 1](https://github.com/HKU-BAL/Clair3#option-1-build-an-environment-with-mambaconda).

**Run Clair3**

```bash
THREADS=4

cd Clair3
./run_clair3.sh \
  --bam_fn=${INPUT_DIR}/${BAM} \
  --ref_fn=${INPUT_DIR}/${REF} \
  --threads=${THREADS} \
  --platform=${PLATFORM} \
  --model_path=$(pwd)"/models/${PLATFORM}" \
  --output=${OUTPUT_DIR} \
  --bed_fn=${INPUT_DIR}/quick_demo.bed
```

**Note**: You can also use `python3 run_clair3.py` instead of `./run_clair3.sh` with the same arguments.

**Run hap.py for benchmarking (optional)**

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
| **INDEL** |    45    |    14    |    4     |     0.762     |      0.918       |      0.833      |
|  **SNP**  |   402    |    0     |    1     |     1.000     |      0.997       |      0.998      |

Check the results using `less ${HOME}/clair3_ont_quickDemo/output/merge_output.vcf.gz`

### Option 2. Docker pre-built image

```bash
THREADS=4

cd ${OUTPUT_DIR}
# Run Clair3 using one command
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:v2.0.0 \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/${BAM} \
  --ref_fn=${INPUT_DIR}/${REF} \
  --threads=${THREADS} \
  --platform=${PLATFORM} \
  --model_path="/opt/models/${PLATFORM}" \
  --output=${OUTPUT_DIR} \
  --bed_fn=${INPUT_DIR}/quick_demo.bed
```

**Note**: You can also use `python3 /opt/bin/run_clair3.py` instead of `/opt/bin/run_clair3.sh` in the Docker command above.

**Run hap.py for benchmarking (optional)**

```bash
# Run hap.py
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
| **INDEL** |    45    |    14    |    4     |     0.762     |      0.918       |      0.833      |
|  **SNP**  |   402    |    0     |    1     |     1.000     |      0.997       |      0.998      |
