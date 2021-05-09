# Clair3 Representation Unification

## Introduction

This document shows how Clair3 unifies representation between training materials and true variant set.

## Prerequisites

- Clair3 installed 
- GNU Parallel installed
- Pypy3 installed (although this submodule could run in Python 3, using pypy interpreter is recommended for better efficiency)

## Input/Output

**Input:**

- BAM: Indexed BAM  input
- REF: Indexed Reference input
- VCF: True variant VCF input
- BED: Confident bed input(optional)

**Ouput:**

- Unified VCF

## Contents

- [1. Setup variables](#1-setup-variables)
- [2.  Phase VCF file using WhatsHap](#2--phase-vcf-file-using-whatshap)
- [3.  Haplotag read alignment using WhatsHap](#3--haplotag-read-alignment-using-whatshap)
- [4.  Prepare true variant set and candidate input](#4--prepare-true-variant-set-and-candidate-input)
- [5.  Unify Representation for true variant set and candidate sites](#5--unify-representation-for-true-variant-set-and-candidate-sites)
- [6.  Merge and sort unified VCF output](#6--merge-and-sort-unified-vcf-output)
- [7.  Benchmark using unified VCF and true variant set (optional)](#7--benchmark-using-unified-vcf-and-true-variant-set-optional)

####  1. Setup variables

```bash
# Setup variables
CLAIR3="clair3.py"                               			   # clair3.py
PYPY="[PYPY_BIN_PATH]"                                         # e.g. pypy3
WHATSHAP="[WHATSHAP_BIN_PATH]"                                 # e.g. whatshap
PARALLEL="[PARALLEL_BIN_PATH]"                                 # e.g. parallel
TABIX="[TABIX_BIN_PATH]"                                       # e.g. tabix
SAMTOOLS="[SAMTOOLS_BIN_PATH]"                                 # e.g. samtools

# Input parameters
PLATFORM="[SEQUENCING_PLATFORM]"                               # e.g. {ont, hifi, ilmn}
VCF_FILE_PATH="[YOUR_VCF_FILE_PATH]"                           # e.g. hg003.vcf.gz
BAM_FILE_PATH="[YOUR_BAM_FILE_PATH]"                           # e.g. hg003.bam
REFERENCE_FILE_PATH="[YOUR_FASTA_FILE_PATH]"                   # e.g. hg003.fasta
BED_FILE_PATH="[YOUR_BED_FILE_PATH]"                           # e.g. hg003.bed
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER_PATH]"					       # e.g. output

# Chromosome prefix ("chr" if chromosome names have the "chr"-prefix)
CHR_PREFIX="chr"

# array of chromosomes (do not include "chr"-prefix)
CHR=(22)

# Number of threads to be used
THREADS=24

# The number of chucks to be divided into for parallel processing
chunk_num=15
CHUNK_LIST=`seq 1 ${chunk_num}`

# Minimum AF required for a candidate variant
MIN_AF=0.08

# Temporary working directory
CANDIDATE_DETAILS_PATH="${OUTPUT_DIR}/candidate_details"
SPLIT_BED_PATH="${OUTPUT_DIR}/split_beds"
VCF_OUTPUT_PATH="${OUTPUT_DIR}/vcf_output"
VAR_OUTPUT_PATH="${OUTPUT_DIR}/var"
PHASE_VCF_PATH="${OUTPUT_DIR}/phased_vcf"
PHASE_BAM_PATH="${OUTPUT_DIR}/phased_bam"

mkdir -p ${CANDIDATE_DETAILS_PATH}
mkdir -p ${SPLIT_BED_PATH}
mkdir -p ${VCF_OUTPUT_PATH}
mkdir -p ${VAR_OUTPUT_PATH}
mkdir -p ${PHASE_VCF_PATH}
mkdir -p ${PHASE_BAM_PATH}
```

#### 2.  Phase VCF file using WhatsHap

To apply representation unification,  using a phased read alignment is highly recommended in order to get more precious unified result.

```bash
cd ${OUTPUT_DIR}

# WhatsHap phasing vcf file if vcf file includes '|' in INFO tag
${WHATSHAP} unphase ${VCF_FILE_PATH} > ${OUTPUT_DIR}/INPUT.vcf.gz

# WhatsHap phase vcf file
${PARALLEL} --joblog ${PHASE_VCF_PATH}/phase.log -j${THREADS} \
"${WHATSHAP} phase \
    --output ${PHASE_VCF_PATH}/phased_{1}.vcf.gz \
    --reference ${REFERENCE_FILE_PATH} \
    --chromosome ${CHR_PREFIX}{1} \
    --ignore-read-groups \
    --distrust-genotypes \
    ${OUTPUT_DIR}/INPUT.vcf.gz \
    ${BAM_FILE_PATH}" ::: ${CHR[@]}
```

#### 3.  Haplotag read alignment using WhatsHap

```bash
# WhatsHap haplotags bam file
${PARALLEL} --joblog ${PHASE_BAM_PATH}/haplotag.log -j${THREADS} \
"${WHATSHAP} haplotag \
    --output ${PHASE_BAM_PATH}/{1}.bam \
    --reference ${REFERENCE_FILE_PATH} \
    --regions ${CHR_PREFIX}{1} \
    --ignore-read-groups \
    ${PHASE_VCF_PATH}/phased_{1}.vcf.gz \
    ${BAM_FILE_PATH}" ::: ${CHR[@]}

# Index the phased bam file using samtools
${PARALLEL} --joblog ${PHASE_BAM_PATH}/index.log -j ${THREADS} ${SAMTOOLS} index -@12 ${PHASE_BAM_PATH}/{1}.bam ::: ${CHR[@]}
```

#### 4.  Prepare true variant set and candidate input

```bash
# Split bed file regions according to the contig name and extend bed region
${PARALLEL} --joblog ${SPLIT_BED_PATH}/split_extend_bed.log -j${THREADS} \
"${PYPY} ${CLAIR3} SplitExtendBed \
    --bed_fn ${BED_FILE_PATH} \
    --output_fn ${SPLIT_BED_PATH}/{1} \
    --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]}

#Get true variant label information from VCF file
${PARALLEL} --joblog ${VAR_OUTPUT_PATH}/get_truth.log -j${THREADS} \
"${PYPY} ${CLAIR3} GetTruth \
    --vcf_fn ${PHASE_VCF_PATH}/phased_{1}.vcf.gz \
    --ctgName ${CHR_PREFIX}{1} \
    --var_fn ${VAR_OUTPUT_PATH}/var_{1}" ::: ${CHR[@]}

# Create candidate details for representation unification
${PARALLEL} --joblog ${CANDIDATE_DETAILS_PATH}/create_tensor.log -j${THREADS} \
"${PYPY} ${CLAIR3} CreateTensorFullAlignment \
    --bam_fn ${PHASE_BAM_PATH}/{1}.bam \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --indel_fn ${CANDIDATE_DETAILS_PATH}/{1}_{2} \
    --ctgName ${CHR_PREFIX}{1} \
    --samtools ${SAMTOOLS} \
    --min_af ${THRESHOLD} \
    --extend_bed ${SPLIT_BED_PATH}/{1} \
    --unify_repre_fn ${CANDIDATE_DETAILS_PATH}/{1}_{2} \
    --unify_repre \
    --phasing_info_in_bam \
    --bed_fn ${BED_FILE_PATH} \
    --chunk_id {2} \
    --chunk_num ${chunk_num}" ::: ${CHR[@]} ::: ${CHUNK_LIST[@]}
    
```

#### 5.  Unify Representation for true variant set and candidate sites

```bash
${PARALLEL} --joblog ${OUTPUT_DIR}/unify_repre.log -j${THREADS} \
"${PYPY} ${CLAIR3} UnifyRepresentation \
    --var_fn ${VAR_OUTPUT_PATH}/var_{1} \
    --candidate_details_fn_prefix ${CANDIDATE_DETAILS_PATH}/{1}_ \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --output_vcf_fn ${VCF_OUTPUT_PATH}/vcf_{1}_{2} \
    --chunk_id {2} \
    --chunk_num ${chunk_num} \
    --bed_fn ${BED_FILE_PATH} \
    --platform ${PLATFORM} \
    --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]} ::: ${CHUNK_LIST[@]} > ${OUTPUT_DIR}/RU.log
    
```

#### 6.  Merge and sort unified VCF output

```bash
cat ${FULL_ALIGNMENT_OUTPUT_PATH}${VCF_OUTPUT_PATH}/vcf_* | ${PYPY} ${CLAIR3} SortVcf --output_fn ${OUTPUT_DIR}/unified.vcf
bgzip -f ${OUTPUT_DIR}/unified.vcf
tabix -f -p vcf ${OUTPUT_DIR}/unified.vcf.gz

```

#### 7.  Benchmark using unified VCF and true variant set (optional)

```bash
# Install hap.py if not installed 
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -n happy-env -c bioconda hap.py
conda activate happy-env

# Benchmark using hap.py
hap.py \
${VCF_FILE_PATH} \
${OUTPUT_DIR}/unified.vcf.gz \
-o ${OUTPUT_DIR}/happy \
-r ${REFERENCE_FILE_PATH} \
-f ${BED_FILE_PATH} \
--threads ${THREADS} \
--engine=vcfeval \
-l "[YOUR_BENCHMARK_REGION]" # e.g. chr22

```
