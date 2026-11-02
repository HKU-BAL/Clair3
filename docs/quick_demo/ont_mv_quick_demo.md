## ONT Signal-Aware Variant Calling Quick Demo (Dwelling Time)
Here is a quick demo for the Oxford Nanopore (ONT) signal-aware variant calling using the dwelling time feature. This demo compares variant calling **with** and **without** dwelling time using GIAB HG003 chromosome 20 data. BAM files must contain `mv` (move table) tags for the dwelling time mode.

```bash
Platform:   ONT
Sample:     GIAB HG003
Coverage:   ~30x
Reference:  GRCh38_no_alt
Region:     chr20:100000-300000
Basecaller: Dorado (with --emit-moves)
Feature:    Dwelling time (signal-aware calling)
```

### Prerequisites

- Docker installed and running for the container commands
- Local Clair3 image: `clair3-rust:local`
- hap.py Docker image: `jmcdani20/hap.py:v0.3.12` for optional benchmarking
- BAM file must contain `mv` tags (produced by Dorado with `--emit-moves`)

Build the local image from the `clair3-rust` repository root:

```bash
docker build -f Dockerfile -t clair3-rust:local .
```

For a native Mamba/Conda run, follow the [Rust installation guide](../rust_installation.md) and install the models with `make models PREFIX="${CONDA_PREFIX}"`.

### Download data

```bash
# Parameters
PLATFORM='ont'
INPUT_DIR="${HOME}/clair3_ont_mv_quickDemo"
THREADS=4

mkdir -p ${INPUT_DIR}

# Download quick demo data
# GRCh38_no_alt reference
wget -P ${INPUT_DIR} https://www.bio8.cs.hku.hk/clair3/demo/ont_mv_quick_demo/GRCh38_no_alt_chr20.fa
wget -P ${INPUT_DIR} https://www.bio8.cs.hku.hk/clair3/demo/ont_mv_quick_demo/GRCh38_no_alt_chr20.fa.fai
# BAM chr20:100000-300000 (with mv tags)
wget -P ${INPUT_DIR} https://www.bio8.cs.hku.hk/clair3/demo/ont_mv_quick_demo/HG003_chr20_demo.bam
wget -P ${INPUT_DIR} https://www.bio8.cs.hku.hk/clair3/demo/ont_mv_quick_demo/HG003_chr20_demo.bam.bai
# GIAB Truth VCF and BED
wget -P ${INPUT_DIR} https://www.bio8.cs.hku.hk/clair3/demo/ont_mv_quick_demo/HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz
wget -P ${INPUT_DIR} https://www.bio8.cs.hku.hk/clair3/demo/ont_mv_quick_demo/HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz.tbi
wget -P ${INPUT_DIR} https://www.bio8.cs.hku.hk/clair3/demo/ont_mv_quick_demo/HG003_GRCh38_chr20_v4.2.1_benchmark_noinconsistent.bed
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

### Verify mv tags in your BAM

You can check whether your BAM file contains `mv` tags using samtools:

```bash
samtools view ${INPUT_DIR}/${BAM} | head -1 | tr '\t' '\n' | grep "^mv:"
```

If the output shows a line starting with `mv:B:c,`, your BAM file contains the required move table tags.

### Option 1. Run Clair3 with dwelling time

The key difference from standard ONT calling is the `--enable_dwell_time` flag and the dwelling-time-aware model `r1041_e82_400bps_hac_with_mv`.

**Docker:**

```bash
OUTPUT_DIR_DWELL="${INPUT_DIR}/output_dwell"
mkdir -p ${OUTPUT_DIR_DWELL}

docker run \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR_DWELL}:${OUTPUT_DIR_DWELL} \
  clair3-rust:local \
  --bam_fn=${INPUT_DIR}/${BAM} \
  --ref_fn=${INPUT_DIR}/${REF} \
  --threads=${THREADS} \
  --platform=${PLATFORM} \
  --model_path="/opt/models/r1041_e82_400bps_hac_with_mv" \
  --output=${OUTPUT_DIR_DWELL} \
  --bed_fn=${INPUT_DIR}/quick_demo.bed \
  --sample_name="HG003" \
  --enable_dwell_time
```

**Mamba/Conda:**

```bash
OUTPUT_DIR_DWELL="${INPUT_DIR}/output_dwell"
mkdir -p ${OUTPUT_DIR_DWELL}

cd Clair3
./run_clair3.sh \
  --bam_fn=${INPUT_DIR}/${BAM} \
  --ref_fn=${INPUT_DIR}/${REF} \
  --threads=${THREADS} \
  --platform=${PLATFORM} \
  --model_path="${CONDA_PREFIX}/bin/models/r1041_e82_400bps_hac_with_mv" \
  --output=${OUTPUT_DIR_DWELL} \
  --bed_fn=${INPUT_DIR}/quick_demo.bed \
  --sample_name="HG003" \
  --enable_dwell_time
```

**Note**: The local image entry point invokes `run_clair3.sh` and supplies the `run` subcommand. To use the Python-named shell launcher instead, add `--entrypoint /opt/clair3/run_clair3.py` before the image name. A native run can use `./run_clair3.py` instead of `./run_clair3.sh`. Neither shell launcher requires Python to execute.

**Note**: Dwell-time calling requires valid `mv` tags. Clair3 checks the first 50 alignments and reports an error if no valid move table is found. Re-run basecalling with Dorado `--emit-moves` when the tags are absent.

### Option 2. Run Clair3 without dwelling time

For comparison, run with the standard model `r1041_e82_400bps_hac_v500` without the `--enable_dwell_time` flag.

**Docker:**

```bash
OUTPUT_DIR_NO_DWELL="${INPUT_DIR}/output_no_dwell"
mkdir -p ${OUTPUT_DIR_NO_DWELL}

docker run \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR_NO_DWELL}:${OUTPUT_DIR_NO_DWELL} \
  clair3-rust:local \
  --bam_fn=${INPUT_DIR}/${BAM} \
  --ref_fn=${INPUT_DIR}/${REF} \
  --threads=${THREADS} \
  --platform=${PLATFORM} \
  --model_path="/opt/models/r1041_e82_400bps_hac_v500" \
  --output=${OUTPUT_DIR_NO_DWELL} \
  --bed_fn=${INPUT_DIR}/quick_demo.bed \
  --sample_name="HG003"
```

**Mamba/Conda:**

```bash
OUTPUT_DIR_NO_DWELL="${INPUT_DIR}/output_no_dwell"
mkdir -p ${OUTPUT_DIR_NO_DWELL}

cd Clair3
./run_clair3.sh \
  --bam_fn=${INPUT_DIR}/${BAM} \
  --ref_fn=${INPUT_DIR}/${REF} \
  --threads=${THREADS} \
  --platform=${PLATFORM} \
  --model_path="${CONDA_PREFIX}/bin/models/r1041_e82_400bps_hac_v500" \
  --output=${OUTPUT_DIR_NO_DWELL} \
  --bed_fn=${INPUT_DIR}/quick_demo.bed \
  --sample_name="HG003"
```

### Run hap.py for benchmarking

```bash
# Evaluate: with dwelling time
mkdir -p ${OUTPUT_DIR_DWELL}/happy
docker run \
  -v "${INPUT_DIR}":"${INPUT_DIR}" \
  -v "${OUTPUT_DIR_DWELL}":"${OUTPUT_DIR_DWELL}" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  ${INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
  ${OUTPUT_DIR_DWELL}/${OUTPUT_VCF_FILE_PATH} \
  -f "${INPUT_DIR}/${BASELINE_BED_FILE_PATH}" \
  -r "${INPUT_DIR}/${REF}" \
  -o "${OUTPUT_DIR_DWELL}/happy/happy" \
  -l ${CONTIGS}:${START_POS}-${END_POS} \
  --engine=vcfeval \
  --threads="${THREADS}" \
  --pass-only

# Evaluate: without dwelling time
mkdir -p ${OUTPUT_DIR_NO_DWELL}/happy
docker run \
  -v "${INPUT_DIR}":"${INPUT_DIR}" \
  -v "${OUTPUT_DIR_NO_DWELL}":"${OUTPUT_DIR_NO_DWELL}" \
  jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
  ${INPUT_DIR}/${BASELINE_VCF_FILE_PATH} \
  ${OUTPUT_DIR_NO_DWELL}/${OUTPUT_VCF_FILE_PATH} \
  -f "${INPUT_DIR}/${BASELINE_BED_FILE_PATH}" \
  -r "${INPUT_DIR}/${REF}" \
  -o "${OUTPUT_DIR_NO_DWELL}/happy/happy" \
  -l ${CONTIGS}:${START_POS}-${END_POS} \
  --engine=vcfeval \
  --threads="${THREADS}" \
  --pass-only
```

### Historical reference output

These tables record the original demonstration results. They are not golden outputs for the Rust implementation; model versions and numerical backends can change calls.

**With dwelling time** (model: `r1041_e82_400bps_hac_with_mv`, `--enable_dwell_time`):

|   Type    | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| :-------: | :------: | :------: | :------: | :-----------: | :--------------: | :-------------: |
| **INDEL** |    52    |    7     |    1     |    0.8814     |     0.9811       |     0.9286      |
|  **SNP**  |   402    |    0     |    1     |    1.0000     |     0.9975       |     0.9988      |

**Without dwelling time** (model: `r1041_e82_400bps_hac_v500`):

|   Type    | TRUTH.TP | TRUTH.FN | QUERY.FP | METRIC.Recall | METRIC.Precision | METRIC.F1_Score |
| :-------: | :------: | :------: | :------: | :-----------: | :--------------: | :-------------: |
| **INDEL** |    52    |    7     |    7     |    0.8814     |     0.8814       |     0.8814      |
|  **SNP**  |   402    |    0     |    1     |    1.0000     |     0.9975       |     0.9988      |

**Historical observation**: In this original comparison, enabling dwelling time improved INDEL precision from 0.88 to 0.98 and F1 score from 0.88 to 0.93, while SNP performance remained the same. False-positive INDEL calls decreased from 7 to 1.

For more details on the dwelling time feature, see [Dwelling Time Feature Documentation](../dwelling_time.md).
