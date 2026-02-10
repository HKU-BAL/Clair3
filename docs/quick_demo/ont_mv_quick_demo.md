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

- Docker installed and running
- Clair3 Docker image: `hkubal/clair3:v2.0.0`
- hap.py Docker image: `jmcdani20/hap.py:v0.3.12`
- BAM file must contain `mv` tags (produced by Dorado with `--emit-moves`)

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
  hkubal/clair3:v2.0.0 \
  /opt/bin/run_clair3.sh \
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
  --model_path=$(pwd)"/models/r1041_e82_400bps_hac_with_mv" \
  --output=${OUTPUT_DIR_DWELL} \
  --bed_fn=${INPUT_DIR}/quick_demo.bed \
  --sample_name="HG003" \
  --enable_dwell_time
```

**Note**: You can also use `python3 /opt/bin/run_clair3.py` (Docker) or `python3 run_clair3.py` (Conda) instead of `run_clair3.sh` with the same arguments.

**Note**: If the BAM file does not contain `mv` tags, Clair3 will still run but the dwell time channel will contain zero values, which may reduce accuracy. Ensure your BAM files were produced by Dorado with `--emit-moves`.

### Option 2. Run Clair3 without dwelling time

For comparison, run with the standard model `r1041_e82_400bps_hac_v500` without the `--enable_dwell_time` flag.

**Docker:**

```bash
OUTPUT_DIR_NO_DWELL="${INPUT_DIR}/output_no_dwell"
mkdir -p ${OUTPUT_DIR_NO_DWELL}

docker run \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR_NO_DWELL}:${OUTPUT_DIR_NO_DWELL} \
  hkubal/clair3:v2.0.0 \
  /opt/bin/run_clair3.sh \
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
  --model_path=$(pwd)"/models/r1041_e82_400bps_hac_v500" \
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

### Expected output

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

**Key observation**: With dwelling time enabled, INDEL precision improves significantly (0.98 vs 0.88) and F1 score improves (0.93 vs 0.88), while SNP performance remains the same. The dwelling time feature reduces false positive INDEL calls from 7 to 1.

For more details on the dwelling time feature, see [Dwelling Time Feature Documentation](../dwelling_time.md).
