# Train a model for Clair3 full-alignment calling

This document shows how to train or fine-tune a deep learning model for Clair3 full-alignment calling.  Compared with Clair3 [pileup](pileup_training.md) model. Full-alignment need much more time to train as the input as input size greatly increased. 

For full-alignment model, we also maintain multiple sample and multiple coverages training workflow. We divided all training materials into `SAMPLES`, `DEPTHS` and `CHR` tags, representing different training samples, different training coverages, and contig name, respectively.  All candidate variants are selected to create binary files and then were fed into training or fine-tune workflow.

## Prerequisites

- Clair3 installed
- GNU Parallel installed
- Sufficient hard-disk space
- Unified VCF file (recommended), check [here](representation_unification.md) for more details
- A powerful GPU

  - RTX Titan (tested)
  - GTX 2080 Ti (tested)

## Contents

* [I. Training data phasing and haplotaging](#i-training-data-phasing-and-haplotaging)
    - [1. Setup variables](#1-setup-variables)
    - [2.  Phase VCF file using WhatsHap](#2--phase-vcf-file-using-whatshap)
    - [3.  Haplotag read alignment using WhatsHap](#3--haplotag-read-alignment-using-whatshap)
    - [4.  Merge phased alignments](#4--merge-phased-alignments)
* [II. Build compressed binary files for full-alignment model training](#ii-build-compressed-binary-files-for-full-alignment-model-training)
    - [1. Setup variables](#1-setup-variables)
    - [2. Create temporary working folders for each submodule](#2-create-temporary-working-folders-for-each-submodule)
    - [3. Split and extend bed regions](#3-split-and-extend-bed-regions-using-the-splitextendbed-submodule)
    - [4. Create full-alignment tensor](#4-create-full-alignment-tensor-using-the-createtensorfullalignment-submodule)
    - [5. Get truth variants from unified VCF file](#5-get-truth-variants-from-unified-vcf-using-the-gettruth-submodule)
    - [6. Convert full-alignment tensor to compressed binary files](#6-convert-full-alignment-tensor-to-compressed-binary-file-using-the-tensor2bin-submodule)
* [III. Model training](#iii-model-training)
    - [1. full-alignment model training](#1-full-alignment-model-training)
    - [2. full-alignment model fine-tune using pre-trained model (optional)](#2-full-alignment-model-fine-tune-using-pre-trained-model-optional)

---

## I. Training data phasing and haplotaging

Full-alignment model integrates with phased alignment to achieve better performance, so phased alignment is highly recommended to use in training. 

**Hints**

> - If representation unification has applied, all phased alignment would be automatically generated in the `${OUTPUT_DIR}/phased_bam` folder, check [here](representation_unification.md#3--haplotag-read-alignment-using-whatshap) for more details.
> - WhatsHap `haplotag` submodule would occupy hard-disk usage of the same size of the input BAM.

#### 1. Setup variables

```bash
# Setup variables
CLAIR3="clair3.py"                               			   # clair3.py
WHATSHAP="[WHATSHAP_BIN_PATH]"                                 # e.g. whatshap
PARALLEL="[PARALLEL_BIN_PATH]"                                 # e.g. parallel
SAMTOOLS="[SAMTOOLS_BIN_PATH]"                                 # e.g. samtools

# Input parameters
PLATFORM="[SEQUENCING_PLATFORM]"                               # e.g. {ont, hifi, ilmn}
VCF_FILE_PATH="[YOUR_VCF_FILE_PATH]"                           # e.g. hg003.vcf.gz
BAM_FILE_PATH="[YOUR_BAM_FILE_PATH]"                           # e.g. hg003.bam
REFERENCE_FILE_PATH="[YOUR_FASTA_FILE_PATH]"                   # e.g. hg003.fasta
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER_PATH]"					       # e.g. output

# Temporary working directory
PHASE_VCF_PATH="${OUTPUT_DIR}/phased_vcf"
PHASE_BAM_PATH="${OUTPUT_DIR}/phased_bam"

mkdir -p ${PHASE_VCF_PATH}
mkdir -p ${PHASE_BAM_PATH}

# Name of training sample
SAMPLE="hg002"

# Chromosome prefix ("chr" if chromosome names have the "chr"-prefix)
CHR_PREFIX="chr"

# array of chromosomes (do not include "chr"-prefix)
CHR=(20)

# Number of threads to be used
THREADS=24

```

#### 2.  Phase VCF file using WhatsHap

```bash
cd ${OUTPUT_DIR}

# WhatsHap phasing vcf file if if vcf file includes '|' in INFO tag
${WHATSHAP} unphase ${VCF_FILE_PATH} > ${OUTPUT_DIR}/INPUT.vcf.gz

# WhatsHap phasing vcf file
${PARALLEL} --joblog ${PHASE_VCF_PATH}/phase.log -j${THREADS} \
"${WHATSHAP} phase \
    --output ${PHASE_VCF_PATH}/phased_{1}.vcf.gz \
    --reference ${REFERENCE_FILE_PATH} \
    --chromosome ${CHR_PREFIX}{1} \
    --ignore-read-groups \
    --distrust-genotypes \
    ${OUTPUT_DIR}/INPUT.vcf.gz \
    ${BAM_FILE_PATH}" ::: ${CHR[@]}
    
# Index the vcf file using tabix, which is neccesary for read haplotagging
${PARALLEL} -j ${THREADS} ${TABIX} -p vcf ${PHASE_VCF_PATH}/phased_{1}.vcf.gz ::: ${CHR[@]}

```

#### 3.  Haplotag read alignment using WhatsHap

```bash
## WhatsHap haplotags bam file
${PARALLEL} --joblog ${PHASE_BAM_PATH}/haplotag.log -j${THREADS} \
"${WHATSHAP} haplotag \
    --output ${PHASE_BAM_PATH}/{1}.bam \
    --reference ${REFERENCE_FILE_PATH} \
    --regions ${CHR_PREFIX}{1} \
    --ignore-read-groups \
    ${PHASE_VCF_PATH}/phased_{1}.vcf.gz \
    ${BAM_FILE_PATH}" ::: ${CHR[@]}

# Index the phased bam file using samtools
time ${PARALLEL} --joblog ${PHASE_BAM_PATH}/index.log -j ${THREADS} ${SAMTOOLS} index -@12 ${PHASE_BAM_PATH}/{1}.bam ::: ${CHR[@]}

```

#### 4.  Merge phased alignments

```bash
${SAMTOOLS} merge -@${THREADS} ${OUTPUT_DIR}/${SAMPLE}.haplotagged.bam ${PHASE_BAM_PATH}/*.bam

${SAMTOOLS} index -@${THREADS} ${OUTPUT_DIR}/${SAMPLE}.haplotagged.bam
```

## II. Build compressed binary files for full-alignment model training

**Hints**

> - The whole procedure was broken into blocks for better readability and error-tracing.
> - For each `parallel` command ran with the `--joblog` option, we can check the `Exitval` column from the job log output. If the column contains a non-zero value, it means error occurred; please try to rerun the block.
> - We suggest using absolute path everywhere.
> - We suggest using a unified VCF for true variant and non-variant labeling, check here for more details on how to generate a unified VCF for each training sample. 

This section shows how to build multiple compressed binary file for multiple samples with or without multiple coverages.

#### 1. Setup variables

```bash
# Setup variables
CLAIR3="clair3.py"                               			   # clair3.py
PYPY="[PYPY_BIN_PATH]"                                         # e.g. pypy3
PYTHON3="[PYTHON3_BIN_PATH]"                                   # e.g. python3
PARALLEL="[PARALLEL_BIN_PATH]"                                 # e.g. parallel
SAMTOOLS="[SAMTOOLS_BIN_PATH]"                                 # e.g. samtools

# Input parameters
PLATFORM="[SEQUENCING_PLATFORM]"                               # e.g. {ont, hifi, ilmn}
UNIFIED_VCF_FILE_PATH="[YOUR_VCF_FILE_PATH_ARRAY]"             # e.g. hg002.unified.vcf.gz
ALL_BAM_FILE_PATH="[YOUR_BAM_FILE_PATH_ARRAY]"                 # e.g. hg002.bam
DEPTHS="[YOUR_DEPTHS_OF_SAMPLES_ARRAY]"                        # e.g. 1000 (means no subsample)
ALL_REFERENCE_FILE_PATH="[YOUR_FASTA_FILE_PATH_ARRAY]"         # e.g. hg002.fasta
ALL_BED_FILE_PATH="[YOUR_BED_FILE_PATH_ARRAY]"                 # e.g. hg002.bed
ALL_SAMPLE="[YOUR_SAMPLE_NAME_ARRAY]"                          # e.g. hg002
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER_PATH]"					       # e.g. output_folder

# Each line represent one input bam with matched coverage with "DEPTH" array
ALL_PHASED_BAM_FILE_PATH=(
'hg002_1000.bam'
'hg002_800.bam'
'hg004_1000.bam'
)

# Each line represents subsample ratio to each sample, 1000 if no subsample applied
DEPTHS=(
1000
800
1000
)

# Each line represents one input sample name
ALL_SAMPLE=(
'hg002'
'hg002'
'hg004'
)

# Each line represents the reference file of each sample
ALL_REFERENCE_FILE_PATH=(
'GRch38.fa'
'GRch38.fa'
'GRch38.fa'
)

# Each line represents one BED region file for each sample
ALL_BED_FILE_PATH=(
'hg002.bed'
'hg002.bed'
'hg004.bed'
)

# Each line represents one unified VCF file for each sample
UNIFIED_VCF_FILE_PATH=(
'hg002_1000.unified.vcf.gz'
'hg002_800.unified.vcf.gz'
'hg004_1000.unified.vcf.gz'
)

# Chromosome prefix ("chr" if chromosome names have the "chr"-prefix)
CHR_PREFIX="chr"

# array of chromosomes (do not include "chr"-prefix) to training in all sample
CHR=(21 22)

# Number of threads to be used
THREADS=24

# Number of chucks to be divided into for parallel processing
chunk_num=15
CHUNK_LIST=`seq 1 ${chunk_num}`

# The number of chucks to be divided for bin file generation for parallel processing
bin_chunk_num=1
BIN_CHUNK_LIST=`seq 1 ${bin_chunk_num}`

# Minimum SNP and INDEL AF required for a candidate variant
MIN_AF=0.08

# Maximum non-variant ratio for full-alignment model training, for full-alignment model training, we use variant :non-variant = 1 : 1
MAXIMUM_NON_VARIANT_RATIO=1

```

#### 2. Create temporary working folders for each submodule

```bash
# Temporary working directory
DATASET_FOLDER_PATH="${OUTPUT_DIR}/build"
TENSOR_CANDIDATE_PATH="${DATASET_FOLDER_PATH}/tensor_can"
BINS_FOLDER_PATH="${DATASET_FOLDER_PATH}/bins"
CANDIDATE_DETAILS_PATH="${DATASET_FOLDER_PATH}/candidate_details"
SPLIT_BED_PATH="${DATASET_FOLDER_PATH}/split_beds"
VAR_OUTPUT_PATH="${DATASET_FOLDER_PATH}/var"

mkdir -p ${DATASET_FOLDER_PATH}
mkdir -p ${TENSOR_CANDIDATE_PATH}
mkdir -p ${BINS_FOLDER_PATH}
mkdir -p ${CANDIDATE_DETAILS_PATH}
mkdir -p ${SPLIT_BED_PATH}
mkdir -p ${VAR_OUTPUT_PATH}
```

#### 3. Split and extend bed regions using the `SplitExtendBed` submodule

```bash
cd ${OUTPUT_DIR}

# Split bed file regions according to the contig name and extend bed region
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/split_extend_bed.log -j${THREADS} \
"${PYPY} ${CLAIR3} SplitExtendBed \
    --bed_fn {4} \
    --output_fn ${SPLIT_BED_PATH}/{2}_{3}_{1} \
    --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_BED_FILE_PATH[@]}
```

#### 4. Create full-alignment tensor using the `CreateTensorFullAlignment` submodule

```bash

# Create full-alignment tensor for model training
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/create_tensor_full_alignment.log -j${THREADS} \
"${PYPY} ${CLAIR3} CreateTensorFullAlignment \
    --bam_fn {4} \
    --ref_fn {5} \
    --tensor_can_fn ${TENSOR_CANDIDATE_PATH}/tensor_can_{2}_{3}_{1}_{7} \
    --indel_fn ${CANDIDATE_DETAILS_PATH}/{2}_{3}_{1}_{7} \
    --ctgName ${CHR_PREFIX}{1} \
    --samtools ${SAMTOOLS} \
    --min_af ${MIN_AF} \
    --extend_bed ${SPLIT_BED_PATH}/{2}_{3}_{1} \
    --bed_fn {6} \
    --add_no_phasing_data_training \
    --phasing_info_in_bam \
    --chunk_id {7} \
    --platform ${PLATFORM} \
    --chunk_num ${chunk_num}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_PHASED_BAM_FILE_PATH[@]} :::+ ${ALL_REFERENCE_FILE_PATH[@]} :::+ ${ALL_BED_FILE_PATH[@]} ::: ${CHUNK_LIST[@]}

```

**Options**

 - `--zstd` : we recommended using [zstd](https://github.com/facebook/zstd) , an extremely fast and lossless compression tool to compress temporary tensor output, which provided much higher compression ratios than other compression tools.
 - `--phasing_info_in_bam` : we enable this option by default, which means we have phased a alignment using WhatsHap `haplotag`, phased alignment will have a `HP` tag in the alignment with haplotype information. 
 - `--add_no_phasing_data_training` :  full-alignment training uses phased alignment materials for training. To increase model robustness, we also add alignment without phasing information into training. For those candidates having phased information, we also generate another tensor disabling phasing channel for training.   

#### 5. Get truth variants from unified VCF using the `GetTruth` submodule

```bash
# Covert unified VCF file into simplified var file
${PARALLEL} --joblog ${VAR_OUTPUT_PATH}/get_truth.log -j${THREADS} \
"${PYPY} ${CLAIR3} GetTruth \
    --vcf_fn {4} \
    --ctgName ${CHR_PREFIX}{1} \
    --var_fn ${VAR_OUTPUT_PATH}/var_{2}_{3}_{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${UNIFIED_VCF_FILE_PATH[@]}

```

#### 6. Convert full-alignment tensor to compressed binary file using the `Tensor2Bin` submodule

```bash
# Convert full-alignment tensor to compressed bin
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/tensor2Bin.log -j${THREADS} \
"${PYTHON3} ${CLAIR3} Tensor2Bin \
    --tensor_fn ${TENSOR_CANDIDATE_PATH}/tensor_can_{2}_{3}_{1} \
    --var_fn ${VAR_OUTPUT_PATH}/var_{2}_{3}_{1} \
    --bin_fn ${BINS_FOLDER_PATH}/{2}_{3}_{1}_{4} \
    --chunk_id {4} \
    --chunk_num ${bin_chunk_num} \
    --allow_duplicate_chr_pos \
    --platform ${PLATFORM} \
    --shuffle \
    --maximum_non_variant_ratio ${MAXIMUM_NON_VARIANT_RATIO} \
    --candidate_details_fn_prefix ${CANDIDATE_DETAILS_PATH}/{2}_{3}_{1}_" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} ::: ${BIN_CHUNK_LIST[@]}

```

**Options**

 - `--allow_duplicate_chr_pos` : for multiple coverages training, this options are required to avoid replacing same variant sites from different coverages.
 - `--shuffle` :  as the input tensor are scanned in order of starting position, we shuffle the training data binary files with chunked iterator in advance to provide more variety. In the training process, we also apply index shuffling to reduce memory occupation.
 - `--maximum_non_variant_ratio` :  we set a maximum non-variant ratio (variant: non-variant = 1:1) for full-alignment model training, non-variants are randomly select from candidate set if exceeds the ratio,  otherwise, all non-variant will be selected for training. 

## III. Model training

#### 1. full-alignment model training 

```bash
# Full-alignment model training
MODEL_FOLDER_PATH="${OUTPUT_DIR}/train"
mkdir -p ${MODEL_FOLDER_PATH}

cd ${MODEL_FOLDER_PATH}
# Currently, we trained model in single GPU, we did not include multi-GPU and multi-worker workflow
export CUDA_VISIBLE_DEVICES="0"
${PYTHON3} ${CLAIR3} Train \
    --bin_fn ${BINS_FOLDER_PATH} \
    --ochk_prefix ${MODEL_FOLDER_PATH} \
    --add_indel_length True \
    --validation_dataset \
    --platform ${PLATFORM}
```

**Options**

 - `--pileup` : `pileup`  flag is the only flag to distinguish pileup model and full-alignment model.
 - `--add_indel_length` :  for full-alignment model, we currently enable two indel-length tasks.
 - `--validation_dataset`: we random select 10% from all candidate site as hold-out validation data, best-performance epoch from validation dataset are select for our default model.

#### 2. full-alignment model fine-tune using pre-trained model (optional)

```bash
# Full-alignment model fine-tune in new sample
MODEL_FOLDER_PATH="${OUTPUT_DIR}/train"
mkdir -p ${MODEL_FOLDER_PATH}

cd ${MODEL_FOLDER_PATH}

export CUDA_VISIBLE_DEVICES="0"
${PYTHON3} ${CLAIR3} Train \
    --bin_fn ${BINS_FOLDER_PATH} \
    --ochk_prefix ${MODEL_FOLDER_PATH} \
    --add_indel_length True \
    --validation_dataset \
    --platform ${PLATFORM} \
    --learning_rate 0.0001 \
    --chkpnt_fn "[YOUR_PRETRAINED_MODEL]"  ## use pre-trained full-alignment model weight here
```

We experimentally offer full-alignment model fine-tune using pre-trained Clair3 model, by using a smaller `learning_rate` and pre-trained checkpoint file `ochk_prefix`, We recommend to use a smaller learning rate 1e-4` to fine-tune pre-trained model.