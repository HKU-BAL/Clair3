## Illumina Variant Calling Quick Demo
Here is a quick demo for the Illumina NGS variant calling using GIAB HG003 chromosome 20 data. We provide a Mamba/Conda virtual environment (**without root privileges**) and docker pre-built image (**with root privileges**) quick demo.

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

### Option 1. Mamba/Conda virtual environment

Install Clair3 following the [Rust installation guide](../rust_installation.md).

**Step 1. Native Illumina realignment**

Realignment is built into Clair3; no Boost installation or C++ compilation is required. These examples use `--no_phasing_for_fa` for a simpler self-contained run. Supply the installed upstream WhatsHap executable with `--whatshap` if you need phased full-alignment input.

**Step 2. Run Clair3**

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
  --bed_fn=${INPUT_DIR}/quick_demo.bed \
  --no_phasing_for_fa
```

**Note**: You can also use `bash run_clair3.py` instead of `./run_clair3.sh` with the same arguments.

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
| **INDEL** |    59    |    0     |    0     |     1.000     |      1.000       |      1.000      |
|  **SNP**  |   386    |    16    |    1     |     0.960     |      0.997       |      0.978      |

Check the results using `less ${HOME}/clair3_illumina_quickDemo/output/merge_output.vcf.gz`

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
  --bed_fn=${INPUT_DIR}/quick_demo.bed \
  --no_phasing_for_fa
```

The image entry point invokes `run_clair3.sh` and supplies the `run` subcommand. To use the Python-named shell launcher instead, add `--entrypoint /opt/clair3/run_clair3.py` before the image name.

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

**Historical hap.py reference output (original release):**

These are the same historical reference metrics shown above, not golden outputs for the Rust implementation.

|   Type    | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1-Score |
| :-------: | :------: | :------: | :------: | :-----------: | :--------------: | :-------------: |
| **INDEL** |    59    |    0     |    0     |     1.000     |      1.000       |      1.000      |
|  **SNP**  |   386    |    16    |    1     |     0.960     |      0.997       |      0.978      |
