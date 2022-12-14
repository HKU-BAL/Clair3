# Train a model for Clair3 full-alignment calling (revision 0)

This document shows how to train and fine-tune a deep learning model for Clair3 full-alignment calling. For training a model for pileup calling, please check [here](pileup_training.md). Clair3 needs both a pileup model and a full-alignment model to work. Compared to Clair3's pileup model training, training a full-alignment model needs much longer time. The disk space requirement also increases significantly. The training materials are grouped according to sample, coverage, and chromosome. The groups are converted into tensor binaries. The binaries are much space-efficient and easier to process. As required, multiples tensor binaries can be used together for model training and fine-tuning. 

---

## Prerequisites

- Clair3 installed
- GNU Parallel installed
- Sufficient hard-disk space
- Truth VCF file after representation unification (check [here](https://github.com/HKU-BAL/Clair3/blob/main/docs/representation_unification.md) on how to generate unified VCF)
- A high-end GPU (have tested in RTX Titan, and RTX 2080Ti)

---

## Contents

* [I. Training data phasing and haplotaging](#i-training-data-phasing-and-haplotaging)
    - [1. Setup variables](#1-setup-variables)
    - [2.  Phase VCF file using WhatsHap](#2--phase-vcf-file-using-whatshap)
    - [3.  Haplotag read alignment susing WhatsHap](#3--haplotag-read-alignments-using-whatshap)
    - [4.  Merge phased alignments](#4--merge-phased-alignments)
* [II. Build compressed binary files for full-alignment model training](#ii-build-compressed-binary-files-for-full-alignment-model-training)
    - [1. Setup variables](#1-setup-variables)
    - [2. Create temporary working folders for each submodule](#2-create-temporary-working-folders-for-each-submodule)
    - [3. Split and extend bed regions](#3-split-and-extend-bed-regions-using-the-splitextendbed-submodule)
    - [4. Get truth variants from unified VCF file](#4-get-truth-variants-from-unified-vcf-using-the-gettruth-submodule)
    - [5. Create full-alignment tensor](#5-create-full-alignment-tensor-using-the-createtrainingtensor-submodule)
    - [6. Merge compressed binaries](#6-merge-compressed-binaries-using-the-mergebin-submodule)
* [III. Model training](#iii-model-training)
    - [1. full-alignment model training](#1-full-alignment-model-training)
    - [2. full-alignment model fine-tune using pre-trained model (optional)](#2-full-alignment-model-fine-tune-using-pre-trained-model-optional)

---

## I. Training data phasing and haplotaging

Full-alignment model utilizes phased alignment, phased alignments is required for training a full-alignment model. 

> - If representation unification has applied, all phased alignment would be automatically generated in the `${OUTPUT_DIR}/phased_bam` folder, check [here](representation_unification.md#3--haplotag-read-alignment-using-whatshap) for more details.
> - WhatsHap `haplotag` submodule requires hard-disk space the same size as the input BAM.

#### 1. Setup variables

```bash
# Setup executable variables
CLAIR3="clair3.py"                                   # clair3.py
WHATSHAP="[WHATSHAP_BIN_PATH]"                       # e.g. whatshap
PARALLEL="[PARALLEL_BIN_PATH]"                       # e.g. parallel
SAMTOOLS="[SAMTOOLS_BIN_PATH]"                       # e.g. samtools

# Input parameters
PLATFORM="[SEQUENCING_PLATFORM]"                     # e.g. {ont, hifi, ilmn}
VCF_FILE_PATH="[YOUR_VCF_FILE_PATH]"                 # e.g. hg003.vcf.gz
BAM_FILE_PATH="[YOUR_BAM_FILE_PATH]"                 # e.g. hg003.bam
REFERENCE_FILE_PATH="[YOUR_FASTA_FILE_PATH]"         # e.g. hg003.fasta
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER_PATH]"               # e.g. output

# Temporary working directories
PHASE_VCF_PATH="${OUTPUT_DIR}/phased_vcf"
PHASE_BAM_PATH="${OUTPUT_DIR}/phased_bam"

mkdir -p ${PHASE_VCF_PATH}
mkdir -p ${PHASE_BAM_PATH}

# Name of the training sample
SAMPLE="hg002"

# Chromosome prefix ("chr" if chromosome names have the "chr"-prefix)
CHR_PREFIX="chr"

# array of chromosomes (do not include "chr"-prefix)
CHR=(20)

# Number of threads to be used
THREADS=8
THREADS_LOW=$((${THREADS}*3/4))
if [[ ${THREADS_LOW} < 1 ]]; then THREADS_LOW=1; fi

```

#### 2.  Phase VCF file using WhatsHap

```bash
cd ${OUTPUT_DIR}

# Remove the phasing information if the VCF input is already phased
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

#### 3.  Haplotag read alignments using WhatsHap

```bash
## WhatsHap haplotags BAM file
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

----

## II. Build compressed binary files for full-alignment model training

> - The whole procedure are breaking into blocks for better readability and error-tracing.
> - For each `parallel` command that run with the `--joblog` option, we can check the `Exitval` column from the job log output. If the column contains a non-zero value, it means error occurred; please rerun the failed block again.
> - We suggest using absolute path EVERYWHERE.
> - You can use a Truth VCF file without representation unification. You might want to do it only for testing because Clair3's performance would be significantly affected without representation unification.

This section shows how to build multiple compressed tensor binary file for multiple samples either with or without multiple coverages.

#### 1. Setup variables

```bash
# Setup executable variables
CLAIR3="clair3.py"                                        # clair3.py
PYPY="[PYPY_BIN_PATH]"                                    # e.g. pypy3
PYTHON3="[PYTHON3_BIN_PATH]"                              # e.g. python3
PARALLEL="[PARALLEL_BIN_PATH]"                            # e.g. parallel
SAMTOOLS="[SAMTOOLS_BIN_PATH]"                            # e.g. samtools

# Input parameters
PLATFORM="[SEQUENCING_PLATFORM]"                          # e.g. {ont, hifi, ilmn}
UNIFIED_VCF_FILE_PATH="[YOUR_VCF_FILE_PATH_ARRAY]"        # e.g. hg002.unified.vcf.gz
ALL_BAM_FILE_PATH="[YOUR_BAM_FILE_PATH_ARRAY]"            # e.g. hg002.bam
DEPTHS="[YOUR_DEPTHS_OF_SAMPLES_ARRAY]"                   # e.g. 1000 (means no subsample)
ALL_REFERENCE_FILE_PATH="[YOUR_FASTA_FILE_PATH_ARRAY]"    # e.g. hg002.fasta
ALL_BED_FILE_PATH="[YOUR_BED_FILE_PATH_ARRAY]"            # e.g. hg002.bed
ALL_SAMPLE="[YOUR_SAMPLE_NAME_ARRAY]"                     # e.g. hg002
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER_PATH]"                    # e.g. output_folder

# Each line represent one input BAM with a matched coverage in the "DEPTH" array
## check the "Training data subsamping" section on how to apply BAM subsampling
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

# Each line represents one representation-unified VCF file for each sample
UNIFIED_VCF_FILE_PATH=(
'hg002_1000.unified.vcf.gz'
'hg002_800.unified.vcf.gz'
'hg004_1000.unified.vcf.gz'
)

# Chromosome prefix ("chr" if chromosome names have the "chr"-prefix)
CHR_PREFIX="chr"

# array of chromosomes (do not include the "chr"-prefix) to train in all sample
## pls note that in the pretrained Clair3 models, we have excluded chr20 as a hold-out set.
CHR=(21 22)

# Number of threads to be used
THREADS=8

# Number of chucks to be divided into for parallel processing
chunk_num=15
CHUNK_LIST=`seq 1 ${chunk_num}`

# The number of chucks to be divided for bin file generation for parallel processing
bin_chunk_num=1
BIN_CHUNK_LIST=`seq 1 ${bin_chunk_num}`

# Minimum SNP and INDEL AF required for a candidate variant
MIN_AF=0.08

# Maximum non-variant ratio for full-alignment model training, for full-alignment model training, we use variant:non-variant = 1:1
MAXIMUM_NON_VARIANT_RATIO=1

```

#### 2. Create temporary working folders for each submodule

```bash
# Temporary working directory
DATASET_FOLDER_PATH="${OUTPUT_DIR}/build"
TENSOR_CANDIDATE_PATH="${DATASET_FOLDER_PATH}/tensor_can"
BINS_FOLDER_PATH="${DATASET_FOLDER_PATH}/bins"
SPLIT_BED_PATH="${DATASET_FOLDER_PATH}/split_beds"
VAR_OUTPUT_PATH="${DATASET_FOLDER_PATH}/var"

mkdir -p ${DATASET_FOLDER_PATH}
mkdir -p ${TENSOR_CANDIDATE_PATH}
mkdir -p ${BINS_FOLDER_PATH}
mkdir -p ${SPLIT_BED_PATH}
mkdir -p ${VAR_OUTPUT_PATH}

```

#### 3. Split and extend bed regions using the `SplitExtendBed` submodule

```bash
cd ${OUTPUT_DIR}

# Split BED file regions according to the contig names and extend the bed regions
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/split_extend_bed.log -j${THREADS} \
"${PYPY} ${CLAIR3} SplitExtendBed \
    --bed_fn {4} \
    --output_fn ${SPLIT_BED_PATH}/{2}_{3}_{1} \
    --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_BED_FILE_PATH[@]}
    
```

#### 4. Get truth variants from unified VCF using the `GetTruth` submodule

```bash
# Covert unified VCF file into simplified var file
${PARALLEL} --joblog ${VAR_OUTPUT_PATH}/get_truth.log -j${THREADS} \
"${PYPY} ${CLAIR3} GetTruth \
    --vcf_fn {4} \
    --ctgName ${CHR_PREFIX}{1} \
    --var_fn ${VAR_OUTPUT_PATH}/var_{2}_{3}_{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${UNIFIED_VCF_FILE_PATH[@]}

```

#### 5. Create full-alignment tensor using the `CreateTrainingTensor` submodule

```bash
# Create full-alignment tensor for model training
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/create_tensor_full_alignment.log -j${THREADS_LOW} \
"${PYPY} ${CLAIR3} CreateTrainingTensor \
    --bam_fn {4} \
    --ref_fn {5} \
    --var_fn ${VAR_OUTPUT_PATH}/var_{2}_{3}_{1} \
    --bin_fn ${TENSOR_CANDIDATE_PATH}/tensor_{2}_{3}_{1}_{7} \
    --ctgName ${CHR_PREFIX}{1} \
    --samtools ${SAMTOOLS} \
    --min_af ${MIN_AF} \
    --extend_bed ${SPLIT_BED_PATH}/{2}_{3}_{1} \
    --bed_fn {6} \
    --phasing_info_in_bam \
    --add_no_phasing_data_training \
    --allow_duplicate_chr_pos \
    --platform ${PLATFORM} \
    --shuffle \
    --maximum_non_variant_ratio ${MAXIMUM_NON_VARIANT_RATIO} \
    --chunk_id {7} \
    --chunk_num ${chunk_num}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_PHASED_BAM_FILE_PATH[@]} :::+ ${ALL_REFERENCE_FILE_PATH[@]} :::+ ${ALL_BED_FILE_PATH[@]} ::: ${CHUNK_LIST[@]}

```

**Options**

 - `--phasing_info_in_bam` : enabled by default, indicating the input BAM is phased using WhatsHap's `haplotag` module, and phased alignments are having a `HP` tag with haplotype details. 
 - `--allow_duplicate_chr_pos` : for multiple coverages training, this option is required to to allow different coverage training samples at the same variant site.
 - `--shuffle` :  as the input tensors are created in the order of starting positions, this option shuffles the training data in each chunk. During the training process, we also apply index reshuffling in each epoch.
 - `--maximum_non_variant_ratio` :  we set a maximum non-variant ratio (variant:non-variant = 1:1) for full-alignment model training, non-variants are randomly selected from the candidate set if the ratio is exceeded, or all non-variants will be used for training otherwise. 
 - `--add_no_phasing_data_training` : also include unphased alignments in additional to the phased alignments for training. We found including unphased alignments increased model robustness.  

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
 - `--random_validation`: randomly holdout 10% from all candidate sites as validation data, the best-performing epoch on the validation data are selected as our pre-trained models.

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