# Train a model for Clair3 full-alignment calling (revision 1)

Compare to [revision 0](full_alignment_training.md), revision 1 runs Representation Unification only once on the full-depth of each sample. Revision 1 trains a multi-depth model faster, provides almost identical performance when depth is ≥ 30x. To train a model that maximizes the performance at depths ≤ 20x, please use [revision 0](full_alignment_training.md).

This document shows how to train and fine-tune a deep learning model for Clair3 full-alignment calling. For training a model for pileup calling, please check [here](pileup_training.md). Clair3 needs both a pileup model and a full-alignment model to work. Compared to Clair3's pileup model training, training a full-alignment model needs much longer time. The disk space requirement also increases significantly. The training materials are grouped according to sample, coverage, and chromosome. The groups are converted into tensor binaries. The binaries are much space-efficient and easier to process. As required, multiples tensor binaries can be used together for model training and fine-tuning. 

---

## Prerequisites

- Clair3 installed
- GNU Parallel installed
- Sufficient hard-disk space
- Truth VCF file after representation unification (check [here](https://github.com/HKU-BAL/Clair3/blob/main/docs/representation_unification.md) on how to generate unified VCF)
- A high-end GPU (have tested on RTX Titan, and RTX 2080Ti)

---

## Contents

* [I. Training data phasing and haplotaging](#i-training-data-phasing-and-haplotaging)
    - [1. Setup variables](#1-setup-variables)
    - [2.  Phase VCF file using WhatsHap](#2--phase-vcf-file-using-whatshap)
    - [3.  Haplotag read alignments using WhatsHap](#3--haplotag-read-alignments-using-whatshap)
* [II. Build compressed binary files for full-alignment model training](#ii-build-compressed-binary-files-for-full-alignment-model-training)
    - [1. Run Clair3 pileup model](#1-run-Clair3-pileup-model)
    - [2. Select low-quality pileup candidates](#2-Select-low-quality-pileup-candidates-using-the-SelectHetSnp-submodule)
    - [3. Split and extend bed regions](#3-split-and-extend-bed-regions-using-the-splitextendbed-submodule)
    - [4. Get truth variants from unified VCF file](#4-get-truth-variants-from-unified-vcf-using-the-gettruth-submodule)
    - [5. Create full-alignment tensor](#5-create-full-alignment-tensor-using-the-createtrainingtensor-submodule)
    - [6. Merge compressed binaries](#6-merge-compressed-binaries-using-the-mergebin-submodule)
* [III. Model training](#iii-model-training)
    - [1. full-alignment model training](#1-full-alignment-model-training)
    - [2. full-alignment model fine-tune using pre-trained model (optional)](#2-full-alignment-model-fine-tune-using-pre-trained-model-optional)

---

## I. Training data phasing and haplotaging

Full-alignment model utilizes phased alignment, phased alignments are required for training a full-alignment model. 

> - The whole procedure are breaking into blocks for better readability and error-tracing.
> - For each `parallel` command that run with the `--joblog` option, we can check the `Exitval` column from the job log output. If the column contains a non-zero value, it means error occurred; please rerun the failed block again.
> - We suggest using absolute path EVERYWHERE.
> - You can use a Truth VCF file without representation unification. You might want to do it only for testing because Clair3's performance would be significantly affected without representation unification.
> - If representation unification has applied, all phased alignment would be automatically generated in the `${OUTPUT_DIR}/phased_bam` folder, check [here](representation_unification.md#3--haplotag-read-alignment-using-whatshap) for more details.
> - WhatsHap `haplotag` submodule requires hard-disk space the same size as the input BAM.

#### 1. Setup variables

```bash
# Setup executable variables
CLAIR3_PATH="${CONDA_PREFIX}/bin"                    # clair3 installation path
CLAIR3="${CLAIR3_PATH}/clair3.py"                    # clair3.py
WHATSHAP="[WHATSHAP_BIN_PATH]"                       # e.g. whatshap
PARALLEL="[PARALLEL_BIN_PATH]"                       # e.g. parallel
SAMTOOLS="[SAMTOOLS_BIN_PATH]"                       # e.g. samtools
TABIX="[TABIX_BIN_PATH]"                             # e.g. tabix

# Input parameters
PLATFORM="[SEQUENCING_PLATFORM]"                     # e.g. {ont, hifi, ilmn}
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER_PATH]"               # e.g. output

ALL_UNPHASED_BAM_FILE_PATH=(
'hg002_1000.bam'
'hg002_800.bam'
'hg004_1000.bam'
)

# Each line represents a sample, a sample can be specified multiple times to allow downsampling
ALL_SAMPLE=(
'hg002'
'hg002'
'hg004'
)

# A downsampling numerator (1000 as denominator) for each sample in ALL_SAMPLE, 1000 means no downsampling, 800 means 80% (800/1000)
DEPTHS=(
1000
800
1000
)

# Reference genome file for each sample
ALL_REFERENCE_FILE_PATH=(
'GRCh38_no_alt.fasta'
'GRCh38_no_alt.fasta'
'GRCh38_no_alt.fasta'
)

# High-confident BED region file for each sample
ALL_BED_FILE_PATH=(
'HG002.bed'
'HG002.bed'
'HG004.bed'
)

# GIAB truth VCF files (without representation unification) for each sample
TRUTH_VCF_FILE_PATH=(
'HG002_GRCh38_v4.2.1.vcf.gz'
'HG002_GRCh38_v4.2.1.vcf.gz'
'HG004_GRCh38_v4.2.1.vcf.gz'
)

# Unified truth VCF file (with representation unification) for each sample
# For a same sample with multiple downsampling depths, only the unified truth VCF done at full depth is needed
UNIFIED_VCF_FILE_PATH=(
'hg002_1000.unified.vcf.gz'
'hg002_1000.unified.vcf.gz'
'hg004_1000.unified.vcf.gz'
)

# Chromosome prefix ("chr" if chromosome names have the "chr"-prefix)
CHR_PREFIX="chr"

# array of chromosomes (do not include "chr"-prefix) to training in all sample
CHR=(21 22)

# Number of threads to be used
THREADS=36
THREADS_LOW=$((${THREADS}*3/4))
if [[ ${THREADS_LOW} < 1 ]]; then THREADS_LOW=1; fi

# Number of chucks to be divided into for parallel processing
chunk_num=15
CHUNK_LIST=`seq 1 ${chunk_num}`

# Maximum non-variant ratio for full-alignment model training, for full-alignment model training, we use variant :non-variant = 1 : 1
MAXIMUM_NON_VARIANT_RATIO=1

# Temporary working directory
DATASET_FOLDER_PATH="${OUTPUT_DIR}/build"
TENSOR_CANDIDATE_PATH="${DATASET_FOLDER_PATH}/tensor_can"
BINS_FOLDER_PATH="${DATASET_FOLDER_PATH}/bins"
CANDIDATE_DETAILS_PATH="${DATASET_FOLDER_PATH}/candidate_details"
CANDIDATE_BED_PATH="${DATASET_FOLDER_PATH}/candidate_bed"
SPLIT_BED_PATH="${DATASET_FOLDER_PATH}/split_beds"
VAR_OUTPUT_PATH="${DATASET_FOLDER_PATH}/var"
PILEUP_OUTPUT_PATH="${OUTPUT_DIR}/pileup_output"
UNPHASED_TRUTH_VCF_PATH="${OUTPUT_DIR}/unphased_truth_vcf"
PHASE_VCF_PATH="${OUTPUT_DIR}/phased_vcf"
PHASE_BAM_PATH="${OUTPUT_DIR}/phased_bam"

mkdir -p ${DATASET_FOLDER_PATH}
mkdir -p ${TENSOR_CANDIDATE_PATH}
mkdir -p ${BINS_FOLDER_PATH}
mkdir -p ${CANDIDATE_DETAILS_PATH}
mkdir -p ${SPLIT_BED_PATH}
mkdir -p ${VAR_OUTPUT_PATH}
mkdir -p ${CANDIDATE_BED_PATH}
mkdir -p ${PILEUP_OUTPUT_PATH}
mkdir -p ${UNPHASED_TRUTH_VCF_PATH}
mkdir -p ${PHASE_VCF_PATH}
mkdir -p ${PHASE_BAM_PATH}
```

#### 2.  Phase VCF file using WhatsHap

```bash
cd ${OUTPUT_DIR}

# Remove the phasing information if the VCF input is already phased
${PARALLEL} -j${THREADS} "${WHATSHAP} unphase {3} > ${UNPHASED_TRUTH_VCF_PATH}/unphased_truth_{1}_{2}.vcf.gz" ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${TRUTH_VCF_FILE_PATH[@]}

# WhatsHap phasing
${PARALLEL} --joblog ${PHASE_VCF_PATH}/phase.log -j${THREADS} \
"${WHATSHAP} phase \
    --output ${PHASE_VCF_PATH}/phased_{2}_{3}_{1}.vcf.gz \
    --reference {5} \
    --chromosome ${CHR_PREFIX}{1} \
    --ignore-read-groups \
    --distrust-genotypes \
    ${UNPHASED_TRUTH_VCF_PATH}/unphased_truth_{2}_{3}.vcf.gz \
    {4}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_UNPHASED_BAM_FILE_PATH[@]} :::+ ${ALL_REFERENCE_FILE_PATH[@]} |& tee ${PHASE_VCF_PATH}/PHASE.log

# Index the phased VCF files using tabix, which is neccesary for read haplotagging
${PARALLEL} -j ${THREADS} ${TABIX} -p vcf ${PHASE_VCF_PATH}/phased_{2}_{3}_{1}.vcf.gz ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]}
```

#### 3.  Haplotag read alignments using WhatsHap

```bash
# WhatsHap haplotaging
${PARALLEL} --joblog ${PHASE_BAM_PATH}/haplotag.log -j${THREADS} \
"${WHATSHAP} haplotag \
    --output ${PHASE_BAM_PATH}/{2}_{3}_{1}.bam \
    --reference {5} \
    --regions ${CHR_PREFIX}{1} \
    --ignore-read-groups \
    ${PHASE_VCF_PATH}/phased_{2}_{3}_{1}.vcf.gz \
    {4}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_UNPHASED_BAM_FILE_PATH[@]} :::+ ${ALL_REFERENCE_FILE_PATH[@]} |& tee ${PHASE_VCF_PATH}/HAPLOTAG.log

# Index the phased bam files using samtools
${PARALLEL} --joblog ${PHASE_BAM_PATH}/index.log -j ${THREADS} ${SAMTOOLS} index -@12 ${PHASE_BAM_PATH}/{2}_{3}_{1}.bam ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]}

```



----

## II. Build compressed binary files for full-alignment model training

This section shows how to build multiple compressed tensor binary file for multiple samples with multiple coverages.

#### 1. Run Clair3 pileup model

```bash
# Call variants using Clair3‘s pileup model with the --pileup_only option
# Only select the candidates in the high-confident BED regions for model training (with --bed_fn)
${PARALLEL} -j1 ${CLAIR3_PATH}/run_clair3.sh \
  --bam_fn={3} \
  --ref_fn={4} \
  --threads=${THREADS} \
  --platform="ont" \
  --model_path="${CONDA_PREFIX}/bin/models/ont" \
  --output=${PILEUP_OUTPUT_PATH}/{1}_{2} \
  --bed_fn={5} \
  --pileup_only ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_UNPHASED_BAM_FILE_PATH[@]} :::+ ${ALL_REFERENCE_FILE_PATH[@]} :::+ ${ALL_BED_FILE_PATH[@]}
```

#### 2. Select low-quality pileup candidates using the `SelectHetSnp` submodule

```bash
# Select all pileup called variants (0/1, 1/1 and 1/2) and some pileup reference calls (0/0) for full-alignment model training
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/select_pileup_candidates.log -j${THREADS} \
"${PYPY} ${CLAIR3} SelectHetSnp \
--alt_fn ${PILEUP_OUTPUT_PATH}/{2}_{3}/pileup.vcf.gz \
--split_folder ${CANDIDATE_BED_PATH} \
--sampleName {2} \
--depth {3} \
--ref_pct_full 0.15 \
--var_pct_full 1.0 \
--chunk_num ${chunk_num} \
--phasing_info_in_bam \
--phase \
--ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]}
```

#### 3. Split and extend bed regions using the `SplitExtendBed` submodule

```bash
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/split_extend_bed.log -j${THREADS} \
"${PYPY} ${CLAIR3} SplitExtendBed \
    --bed_fn {4} \
    --output_fn ${SPLIT_BED_PATH}/{2}_{3}_{1} \
    --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_BED_FILE_PATH[@]}
```

#### 4. Get truth variants from unified VCF using the `GetTruth` submodule

```bash
# Convert an unified VCF file into a simplified var file
${PARALLEL} --joblog ${VAR_OUTPUT_PATH}/get_truth.log -j${THREADS} \
"${PYPY} ${CLAIR3} GetTruth \
    --vcf_fn {4} \
    --ctgName ${CHR_PREFIX}{1} \
    --var_fn ${VAR_OUTPUT_PATH}/var_{2}_{3}_{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${UNIFIED_VCF_FILE_PATH[@]}
```

#### 5. Create full-alignment tensor using the `CreateTrainingTensor` submodule

```bash
# Create full-alignment tensors for model training
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/create_tensor_full_alignment.log -j${THREADS_LOW} \
"${PYPY} ${CLAIR3} CreateTrainingTensor \
    --bam_fn ${PHASE_BAM_PATH}/{2}_{3}_{1}.bam \
    --ref_fn {5} \
    --var_fn ${VAR_OUTPUT_PATH}/var_{2}_{3}_{1} \
    --bin_fn ${TENSOR_CANDIDATE_PATH}/tensor_{2}_{3}_{1}_{7} \
    --ctgName ${CHR_PREFIX}{1} \
    --samtools ${SAMTOOLS} \
    --extend_bed ${SPLIT_BED_PATH}/{2}_{3}_{1} \
    --full_aln_regions ${CANDIDATE_BED_PATH}/{2}_{3}_{1}_{7} \
    --bed_fn {6} \
    --phasing_info_in_bam \
    --add_no_phasing_data_training \
    --allow_duplicate_chr_pos \
    --platform ${PLATFORM} \
    --shuffle \
    --maximum_non_variant_ratio ${MAXIMUM_NON_VARIANT_RATIO} \
    --chunk_id {7} \
    --chunk_num ${chunk_num}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_UNPHASED_BAM_FILE_PATH[@]} :::+ ${ALL_REFERENCE_FILE_PATH[@]} :::+ ${ALL_BED_FILE_PATH[@]} ::: ${CHUNK_LIST[@]}
```

**Options**

 - `--phasing_info_in_bam` : enabled by default, indicating the input BAM is phased using WhatsHap's `haplotag` module, and phased alignments are having a `HP` tag with haplotype details. 
 - `--allow_duplicate_chr_pos` : for multiple coverages training, this option is required to to allow different coverage training samples at the same variant site.
 - `--shuffle` :  as the input tensors are created in the order of starting positions, this option shuffles the training data in each chunk. During the training process, we also apply index reshuffling in each epoch.
 - `--maximum_non_variant_ratio` :  we set a maximum non-variant ratio (variant:non-variant = 1:1) for full-alignment model training, non-variants are randomly selected from the candidate set if the ratio is exceeded, or all non-variants will be used for training otherwise. 
 - `--add_no_phasing_data_training` : also include unphased alignments in additional to the phased alignments for training. We found including unphased alignments increased model robustness. 
 - `--full_aln_regions` : provide the pileup candidate regions to be included in full-alignment based calling. 

#### 6. Merge compressed binaries using the `MergeBin` submodule

```bash
# Merge compressed binaries
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/mergeBin.log -j${THREADS} \
"${PYTHON3} ${CLAIR3} MergeBin \
    ${TENSOR_CANDIDATE_PATH}/tensor_{2}_{3}_{1}_* \
    --platform ${PLATFORM} \
    --out_fn ${BINS_FOLDER_PATH}/bin_{2}_{3}_{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]}
```

----

## III. Model training

We provide two optional training mode:

​	**Option1**: Train  pileup model using new dataset, in this mode, we will use randomly initialized model weights and train the model until reaches max epochs(30) or converge.

​    **Option2**: Fine-tune pileup model using pre-trained parameters and choose a smaller learning rate for better converge in new dataset.

***We recommend using the fine-tune mode (option 2) for better robustness.***

#### 1. full-alignment model training 

```bash
# Full-alignment model training
MODEL_FOLDER_PATH="${OUTPUT_DIR}/train"
mkdir -p ${MODEL_FOLDER_PATH}

cd ${MODEL_FOLDER_PATH}

# A single GPU is used for model training
export CUDA_VISIBLE_DEVICES="0"
${PYTHON3} ${CLAIR3} Train \
    --bin_fn ${BINS_FOLDER_PATH} \
    --ochk_prefix ${MODEL_FOLDER_PATH}/full_alignment \
    --add_indel_length True \
    --random_validation \
    --platform ${PLATFORM}
```

**Options**

 - `--add_indel_length` :  enable or disable the two indel-length tasks. In the pre-trained models, the two tasks are enabled in full-alignment calling.
 - `--random_validation`: randomly holdout 10% from all candidate sites as validation data, the best-performing epoch on the validation data are selected as our pre-trained model.

#### 2. full-alignment model fine-tune using pre-trained model (optional)

```bash
# Full-alignment model fine-tuning using a new sample
MODEL_FOLDER_PATH="${OUTPUT_DIR}/train"
mkdir -p ${MODEL_FOLDER_PATH}

cd ${MODEL_FOLDER_PATH}

export CUDA_VISIBLE_DEVICES="0"
${PYTHON3} ${CLAIR3} Train \
    --bin_fn ${BINS_FOLDER_PATH} \
    --ochk_prefix ${MODEL_FOLDER_PATH}/full_alignment \
    --add_indel_length True \
    --random_validation \
    --platform ${PLATFORM} \
    --learning_rate 0.0001 \
    --chkpnt_fn "[YOUR_PRETRAINED_MODEL]"  ## use pre-trained full-alignment model here
```

We experimentally offer full-alignment model fine tuning using a pre-trained Clair3 full-alignment model, by using a smaller `learning_rate` and a pre-trained model `chkpnt_fn`. We recommend starting with a smaller learning rate such as `1e-4` to fine-tune a pre-trained full-alignment model.
