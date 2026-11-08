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

# Number of chunks to be divided into for parallel processing
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
  --print_ref_calls \
  --model_path="${CONDA_PREFIX}/bin/models/ont" \
  --output=${PILEUP_OUTPUT_PATH}/{1}_{2} \
  --bed_fn={5} \
  --pileup_only ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_UNPHASED_BAM_FILE_PATH[@]} :::+ ${ALL_REFERENCE_FILE_PATH[@]} :::+ ${ALL_BED_FILE_PATH[@]}
```

#### 2. Select low-quality pileup candidates using the `SelectHetSnp` submodule

```bash
# Select all pileup called variants (0/1, 1/1 and 1/2) and some pileup reference calls (0/0) for full-alignment model training
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/select_pileup_candidates.log -j${THREADS} \
"${CLAIR3} SelectHetSnp \
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
"${CLAIR3} SplitExtendBed \
    --bed_fn {4} \
    --output_fn ${SPLIT_BED_PATH}/{2}_{3}_{1} \
    --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${ALL_BED_FILE_PATH[@]}
```

#### 4. Get truth variants from unified VCF using the `GetTruth` submodule

```bash
# Convert an unified VCF file into a simplified var file
${PARALLEL} --joblog ${VAR_OUTPUT_PATH}/get_truth.log -j${THREADS} \
"${CLAIR3} GetTruth \
    --vcf_fn {4} \
    --ctgName ${CHR_PREFIX}{1} \
    --var_fn ${VAR_OUTPUT_PATH}/var_{2}_{3}_{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]} :::+ ${UNIFIED_VCF_FILE_PATH[@]}
```

#### 5. Create full-alignment tensor using the `CreateTrainingTensor` submodule

```bash
# Create full-alignment tensors for model training
${PARALLEL} --joblog ${DATASET_FOLDER_PATH}/create_tensor_full_alignment.log -j${THREADS_LOW} \
"${CLAIR3} CreateTrainingTensor \
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
"${CLAIR3} MergeBin \
    ${TENSOR_CANDIDATE_PATH}/tensor_{2}_{3}_{1}_* \
    --platform ${PLATFORM} \
    --out_fn ${BINS_FOLDER_PATH}/bin_{2}_{3}_{1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} :::+ ${DEPTHS[@]}
```

----

## III. Model training

On `clair3-rust`, use the installed `clair3` executable or the `clair3.py` compatibility launcher. `--device=cpu` forces CPU training and `--device=cuda:0` explicitly selects the first CUDA device; otherwise, available CUDA is selected automatically. MPS training is not validated. Each epoch writes `<ochk_prefix>.NN.pt`, such as `full_alignment.01.pt` (or `fa.01.pt` in the quick script). `training.log` and `best_val_loss.pt` default to the current directory; `--log_fn` and `--best_checkpoint` can select other paths. Use either `--random_validation` or `--validation_fn` with an independent validation-bin directory, not both. Full-alignment training with indel-length tasks requires `--add_indel_length True`. For a complete native smoke workflow, see [Quick-testing Clair3-Rust](rust_quick_test.md#7-quick-model-training-workflow).

### Input correctness, independent validation, and restarts

Candidate decoding now advances the read query offset for every CIGAR insertion. Previously, a missing `Cigar::Ins` increment left subsequent SNP bases and inserted sequences at stale offsets, corrupting full-alignment candidate rows. A real indexed BAM regression with a `5M2I5M` read checks the insertion and following bases exactly. Representation-unification fixes for phased SNP/insertion/deletion matching and deterministic candidate ALT order are described in the [current validation report](rust_optimization.md#end-to-end-model-training-validation).

For `--validation_fn`, use genuinely held-out archives. Training-only `_20_` and `exclude_training_samples` filters no longer discard independent validation inputs, and an explicit directory with no complete validation batch is an error rather than a zero-loss success. In `scripts/run_quick_train_chr21.sh`, depth 222 is selected by `samtools view -s 222.222`; samtools `-U` writes the complementary, disjoint training QNAME set. The nominal depth-1000 training BAM therefore contains approximately 77.8% of the haplotagged reads, and additional training depths draw only from that complement. Bounded `BATCH_SIZE=1000` and `CHUNK_SIZE=200` defaults prevent empty batches in this smoke workflow and can be overridden for larger datasets.

`--chkpnt_fn` is a **weight-only restart**, not exact continuation: the trainer creates a fresh AdamW optimizer and epoch schedule, and does not restore optimizer, epoch, or RNG state. Controlled one-epoch tests with seed 29 exactly reproduced the fresh-run losses, validation losses, logs, and every checkpoint tensor, while a different initial checkpoint changed the result.

### Full-alignment correctness evidence and validation limitation

The deterministic original-versus-Rust fixture uses the same 40 labeled examples, initial checkpoint, seed 0, and one fixed 20-example batch. All four task losses, actual total loss, every updated checkpoint tensor, and reload outputs match with zero tolerance. Moving seeding after model construction/checkpoint loading preserved the original dropout stream; model-specific AdamW selection retains original full-alignment weight decay `1e-7` (and uses `1e-4` for pileup). Final full-alignment task losses were `[2.971264362335205, 1.2697056531906128, 3.38578462600708, 2.9342098236083984]`, with total-loss and maximum checkpoint differences exactly 0.

On the current complementary-read native data, three-epoch training/validation losses were `(9.195758, 0.650995)`, `(7.702460, 0.986795)`, and `(8.630635, 1.689801)`. CUDA likewise improved training loss overall while validation worsened: training `9.363139, 7.571423, 8.508289`; validation `0.504872, 0.900081, 1.585159`. **Full-alignment independent validation worsens during these bounded tests even though training improves overall; validation convergence or improvement is not claimed.** A zero-learning-rate probe also changed validation through training-mode state with unchanged weights, so constant validation was not a valid early-stopping substitute. The real-data test instead verifies the observed early stop at epoch 11 of 12 and the weight-only restart behavior.

Both native and original-Python checkpoints cross-load with finite `(2,90)` outputs on the tested platforms. Original CUDA DDP was not executed on the single-RTX-4090 host; the separate native two-process gradient-averaging test passed exactly. Platform timings, archive shapes, held-out chr20 F1, named baseline skips, and all residual differences are listed in the [current training validation report](rust_optimization.md#end-to-end-model-training-validation).

We provide two optional training mode:

​	**Option1**: Train  full-alignment model using new dataset, in this mode, we will use randomly initialized model weights and train the model until reaches max epochs(30) or converge.

​    **Option2**: Fine-tune full-alignment model using pre-trained parameters and choose a smaller learning rate for better converge in new dataset.

***We recommend using the fine-tune mode (option 2) for better robustness.***

#### 1. full-alignment model training 

```bash
# Full-alignment model training
MODEL_FOLDER_PATH="${OUTPUT_DIR}/train"
mkdir -p ${MODEL_FOLDER_PATH}

cd ${MODEL_FOLDER_PATH}

# A single GPU is used for model training
export CUDA_VISIBLE_DEVICES="0"
${CLAIR3} Train \
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
${CLAIR3} Train \
    --bin_fn ${BINS_FOLDER_PATH} \
    --ochk_prefix ${MODEL_FOLDER_PATH}/full_alignment \
    --add_indel_length True \
    --random_validation \
    --platform ${PLATFORM} \
    --learning_rate 0.0001 \
    --chkpnt_fn "[YOUR_PRETRAINED_MODEL]"  ## use pre-trained full-alignment model here
```

We experimentally offer full-alignment model fine tuning using a pre-trained Clair3 full-alignment model, by using a smaller `learning_rate` and a pre-trained model `chkpnt_fn`. We recommend starting with a smaller learning rate such as `1e-4` to fine-tune a pre-trained full-alignment model.
