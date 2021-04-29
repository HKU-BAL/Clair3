# Train a model for Clair3 pileup calling

This document shows how to train or fine-tune a deep learning model for Clair3 pileup calling.  We simplified the training workflow in [Clair](https://github.com/HKU-BAL/Clair/blob/master/docs/TRAIN.md), while still maintaining multiple sample and multiple coverages training workflow. We divided all training materials into `SAMPLES`, `DEPTHS` and `CHR` tags, which represent different training samples, different training coverages, and contig name, respectively.  All candidate variants are selected to create binary file and then were fed into training or fine-tune workflow.

## Prerequisites
- Clair3 installed
- GNU Parallel installed
- Sufficient hard-disk space
- Unified VCF file(recommended), check [here]() on how to generate unified VCF
- A powerful GPU
  
    - RTX Titan (tested)
    - GTX 2080 Ti (tested)
    - GTX 1080 Ti (tested)
    

## Contents

* [I. Training data subsamping](#i-training-data-subsamping-recommended)
* [II. Build compressed binary files](#ii-build-compressed-binary-files-for-pileup-model-training)
    - [1. Setup variables](#1-setup-variables)
    - [2. Create working directories ](#2-create-temporary-working-folders-for-each-submodule)
    - [3. Split and extend bed regions](#3-split-bed-regions-using-the-splitextendbed-submodule)
    - [4. Create pileup tensor](#4-create-pileup-tensor-using-the-createtensorpileup-submodule)
    - [5. Get truth variants from unified VCF file](#5-get-truth-variants-from-unified-vcf-using-the-gettruth-submodule)
    - [6. Convert pileup tensor to compressed binary file](#6-convert-pileup-tensor-to-compressed-binary-file-using-the-tensor2bin-submodule)
* [III. Model training](#iii-model-training)
    - [1. Pileup model training](#1-pileup-model-training)
    - [2. Pileup model fine-tune using pre-trained model](#2-pileup-model-fine-tune-using-pre-trained-model-optional)

---

## I. Training data subsamping (recommended)

To build a training dataset with multiple coverages, we need to create multiple subsampled BAM files from the original full-depth BAM file.

```bash
# Please make sure the provided bam file is sorted and samtools indexed
ALL_BAM_FILE_PATH=(
'hg002.bam'
'hg002.bam'
'hg002.bam'
)

# Each line represents one input sample name
ALL_SAMPLE=(
'hg002'
'hg002'
'hg002'
)

# Each line represents subsample ration to each sample
## FRAC values for 'samtools view -s INT.FRAC'
## please refer to samtools' documentation for further information
## here we set 90%, 60% and 30% of the full coverage as example
DEPTHS=(
900
600
300)

# Output folder to store all subsampled BAM
SUBSAMPLED_BAMS_FOLDER_PATH="[SUBSAMPLED_BAMS_FOLDER_PATH]"
mkdir -p ${SUBSAMPLED_BAMS_FOLDER_PATH}

# Other parameters
THREADS=12
PARALLEL='parallel'
SAMTOOLS='samtools'

# Subsample BAM to specific depth in parallel
${PARALLEL} -j${THREADS} "${SAMTOOLS} view -@12 -s {2}.{2} -b -o ${SUBSAMPLED_BAMS_FOLDER_PATH}/{2}_{1}.bam {3}" ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_BAM_FILE_PATH[@]}
${PARALLEL} -j${THREADS} "${SAMTOOLS} index ${SUBSAMPLED_BAMS_FOLDER_PATH}/{2}_{1}.bam" ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]}

# Add symbolic links for full coverage BAM and bai index
${PARALLEL} "ln -s {2} ${SUBSAMPLED_BAMS_FOLDER_PATH}/{1}_1000.bam" ::: ${ALL_SAMPLE[@]}  :::+ ${ALL_BAM_FILE_PATH[@]}
${PARALLEL} "ln -s {2}.bai ${SUBSAMPLED_BAMS_FOLDER_PATH}/{1}_1000.bam.bai" ::: ${ALL_SAMPLE[@]}  :::+ ${ALL_BAM_FILE_PATH[@]}

```

## II. Build compressed binary files for pileup model training

**Hints**

> - The whole procedure was break into blocks for better readability and error-tracing.
> - For each `parallel` command ran with the `--joblog` option, we can check the `Exitval` column from the job log output. If the column contains a non-zero value, it means error occurred, please try to rerun the block again.
> - We suggest to use absolute path everywhere.
> - We suggest to use a unified VCF for true variant and non-variant labeling, check here for more details on how to generate a unified VCF for each training sample. 

This section shows how to build multiple compressed binary file for multiple samples with or without multiple coverages.

#### 1. Setup variables
```bash
# Setup bin variables
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
## check the "Training data subsamping" section on how to apply BAM subsampling
ALL_BAM_FILE_PATH=(
'hg002_1000.bam'
'hg002_800.bam'
'hg004_1000.bam'
)

# Each line represents subsample ration to each sample, 1000 if no subsample applied
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
## pls note that we have excluded chr20 as hold-out set, so do not use chr20 for traning.
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
MIN_SNP_AF=0.08
MIN_INDEL_AF=0.15

# Maximum non-variant ratio for pileup model training, for pileup model training, we use variant:non-variant = 1:5
MAXIMUM_NON_VARIANT_RATIO=5

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
# Split bed file regions according to the contig name and extend bed region
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/split_extend_bed.log -j${THREADS} \
"${PYPY} ${CLAIR3} SplitExtendBed \
    --bed_fn {4} \
    --output_fn ${SPLIT_BED_PATH}/{2}_{3}_{1} \
    --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_BED_FILE_PATH[@]}
```

#### 4. Create pileup tensor using the `CreateTensorPileup` submodule

```bash
# Create pileup tensor for model training
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/create_tensor_pileup.log -j${THREADS} \
"${PYPY} ${CLAIR3} CreateTensorPileup \
    --bam_fn {4} \
    --ref_fn {5} \
    --tensor_can_fn ${TENSOR_CANDIDATE_PATH}/tensor_can_{2}_{3}_{1}_{7} \
    --indel_fn ${CANDIDATE_DETAILS_PATH}/{2}_{3}_{1}_{7} \
    --ctgName ${CHR_PREFIX}{1} \
    --samtools ${SAMTOOLS} \
    --snp_min_af ${MIN_SNP_AF} \
    --indel_min_af ${MIN_INDEL_AF} \
    --extend_bed ${SPLIT_BED_PATH}/{2}_{3}_{1} \
    --platform ${PLATFORM} \
    --bed_fn {6} \
    --chunk_id {7} \
    --chunk_num ${chunk_num}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_BAM_FILE_PATH[@]} :::+ ${ALL_REFERENCE_FILE_PATH[@]} :::+ ${ALL_BED_FILE_PATH[@]} ::: ${CHUNK_LIST[@]} |& tee  ${DATASET_FOLDER_PATH}/CTP.log

```

**Options**

 - `--zstd` : we recommended using [zstd](https://github.com/facebook/zstd) , an extremely fast and lossless compression tool to compress temporary tensor output, which provided much higher compression ratios compared with other compression tools.
 - `--max_depth` :  pileup input summarizes position-level read alignments, therefore depth information varies for various training materials. If need to re-train or fine-tune Clair3 pileup model, we recommend setting a maximum depth based on the maximum-coverage training materials accordingly.

#### 5. Get truth variants from unified VCF using the `GetTruth` submodule

```bash
${PARALLEL} --joblog ${VAR_OUTPUT_PATH}/get_truth.log -j${THREADS} \
"${PYPY} ${CLAIR3} GetTruth \
    --vcf_fn {4} \
    --ctgName ${CHR_PREFIX}{1} \
    --var_fn ${VAR_OUTPUT_PATH}/var_{2}_{3}_{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${UNIFIED_VCF_FILE_PATH[@]}

```

#### 6. Convert pileup tensor to compressed binary file using the `Tensor2Bin` submodule
```bash
# Convert pileup tensor to compressed bin
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/tensor2Bin.log -j${THREADS} \
"${PYTHON3} ${CLAIR3} Tensor2Bin \
    --tensor_fn ${TENSOR_CANDIDATE_PATH}/tensor_can_{2}_{3}_{1} \
    --var_fn ${VAR_OUTPUT_PATH}/var_{2}_{3}_{1} \
    --bin_fn ${BINS_FOLDER_PATH}/{2}_{3}_{1}_{4} \
    --chunk_id {4} \
    --chunk_num ${bin_chunk_num} \
    --pileup \
    --allow_duplicate_chr_pos \
    --platform ${PLATFORM} \
    --shuffle \
    --maximum_non_variant_ratio ${MAXIMUM_NON_VARIANT_RATIO} \
    --candidate_details_fn_prefix ${CANDIDATE_DETAILS_PATH}/{2}_{3}_{1}_" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} ::: ${BIN_CHUNK_LIST[@]}

```

**Options**

 - `--allow_duplicate_chr_pos` : for multiple coverages training, this options are required to avoid replace same variant sites from different coverages.
 - `--shuffle` :  as the input tensor are scanned in order of starting position, we shuffle the training data binary files with chunked iterator in advance to provide more variety. in the training process, we also apply index shuffling to reduce memory occupation.
 - `--maximum_non_variant_ratio` :  we set a maximum non-variant ratio (variant: non-variant = 1:5) for pileup model training, non-variants are randomly select from candidate set if exceeds the ratio,  otherwise, all non-variant will be selected for training. 

## III. Model training

#### 1. Pileup model training 

```bash
MODEL_FOLDER_PATH="${OUTPUT_DIR}/train"
mkdir -p ${MODEL_FOLDER_PATH}

cd ${MODEL_FOLDER_PATH}

# Currently, we trained model in single GPU, we did not include multi-GPU and multi-worker workflow here
export CUDA_VISIBLE_DEVICES="0"

${PYTHON3} ${CLAIR3} Train \
    --bin_fn ${BINS_FOLDER_PATH} \
    --ochk_prefix ${MODEL_FOLDER_PATH} \
    --pileup \
    --add_indel_length False \
    --validation_dataset \
    --platform ${PLATFORM}

```

**Options**

 - `--pileup` : `pileup`  flag is the only flag to distinguish pileup model and full-alignment model.
 - `--add_indel_length` :  for pileup model, we currently disabled two indel-length tasks.
 - `--validation_dataset`: we random select 10% from all candidate site as hold-out validation data, best-performance epoch from validation dataset are select for our default model.

#### 2. Pileup model fine-tune using pre-trained model (optional)

```bash
# Pileup model fine-tune in new sample

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
    --learning_rate 0.0005 \
    --chkpnt_fn "[YOUR_PRETRAINED_MODEL]"  ## use pre-trained full-alignment model weight here
```

We experimentally offer pileup model fine-tune using pre-trained Clair3 model, by using a smaller `learning_rate` and pre-trained checkpoint file `ochk_prefix`, We recommend to use a smaller learning rate `5e-4` to fine-tune pre-trained model.