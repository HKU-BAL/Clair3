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

Install Clair3 following the [Rust installation guide](../rust_installation.md).

**Run Clair3**

```bash
THREADS=4

cd Clair3
make models PREFIX="${CONDA_PREFIX}"
./run_clair3.sh \
  --bam_fn=${INPUT_DIR}/${BAM} \
  --ref_fn=${INPUT_DIR}/${REF} \
  --threads=${THREADS} \
  --platform=${PLATFORM} \
  --model_path="${CONDA_PREFIX}/bin/models/${PLATFORM}" \
  --output=${OUTPUT_DIR} \
  --bed_fn=${INPUT_DIR}/quick_demo.bed
```

**Note**: You can also use `./run_clair3.py` instead of `./run_clair3.sh` with the same arguments.

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

**Historical hap.py reference output (original release):**

These are historical reference metrics, not golden outputs for the Rust implementation; model versions and numerical backends can change calls.

|   Type    | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1-Score |
| :-------: | :------: | :------: | :------: | :-----------: | :--------------: | :-------------: |
| **INDEL** |    45    |    14    |    4     |     0.762     |      0.918       |      0.833      |
|  **SNP**  |   402    |    0     |    1     |     1.000     |      0.997       |      0.998      |

Check the results using `less ${HOME}/clair3_ont_quickDemo/output/merge_output.vcf.gz`

### Option 2. Local Docker image

Build the image from the `clair3-rust` repository root, then run the demo:

```bash
docker build -f Dockerfile -t clair3-rust:local .
THREADS=4

cd ${OUTPUT_DIR}
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  clair3-rust:local \
  --bam_fn=${INPUT_DIR}/${BAM} \
  --ref_fn=${INPUT_DIR}/${REF} \
  --threads=${THREADS} \
  --platform=${PLATFORM} \
  --model_path="/opt/models/${PLATFORM}" \
  --output=${OUTPUT_DIR} \
  --bed_fn=${INPUT_DIR}/quick_demo.bed
```

The image entry point invokes `run_clair3.sh` and supplies the `run` subcommand. To use the Python-named shell launcher instead, add `--entrypoint /opt/clair3/run_clair3.py` before the image name.

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

**Historical hap.py reference output (original release):**

These are the same historical reference metrics shown above, not golden outputs for the Rust implementation.

|   Type    | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1-Score |
| :-------: | :------: | :------: | :------: | :-----------: | :--------------: | :-------------: |
| **INDEL** |    45    |    14    |    4     |     0.762     |      0.918       |      0.833      |
|  **SNP**  |   402    |    0     |    1     |     1.000     |      0.997       |      0.998      |
