#!/bin/bash

################################################################################
# Unified Training Pipeline for Clair3 (Multi-Sample Support)
# Pipeline: UnifyRepresentation -> Pileup Training -> Full-Alignment Training
# Supports constructing bins from multiple samples and training on the combined dataset.
################################################################################

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit if any command in a pipeline fails

################################################################################
# Control Variables - Set to 1 to run, 0 to skip
################################################################################
RUN_RU=1                    # Run UnifyRepresentation
RUN_DOWNSAMPLE=1            # Run downsampling (optional)
BUILD_PILEUP=1              # Build pileup dataset
TRAIN_PILEUP=1              # Train pileup model
BUILD_FULL_ALIGNMENT=0      # Build full-alignment dataset
TRAIN_FULL_ALIGNMENT=0      # Train full-alignment model
RUN_EVALUATION=0            # Run hap.py evaluation (per sample)

# Training configuration
TRAIN_NAME="mv_table"
USE_GPU=1
GPU_DEVICE="0"

################################################################################
# Sample Configuration
################################################################################
# Define your samples here. Arrays must be of equal length.
SAMPLES=(
    "hg001"
    "hg005"
    # "hg003" 
    # "hg004"
)

BAM_FILE_PATHS=(
    "/autofs/bal19/xyu/ont_open_data/HG001/hac_calls/HG001_hac_merged.bam"
    "/autofs/bal19/xyu/ont_open_data/HG005/hac_calls/HG005_hac_merged.bam"
    # "/path/to/hg003.bam"
    # "/path/to/hg004.bam"
)

VCF_FILE_PATHS=(
    "/autofs/bal33/zxzheng/data/vcf/hg38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    "/autofs/bal33/zxzheng/data/vcf/hg38/HG005_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    # "/path/to/hg003_benchmark.vcf.gz"
    # "/path/to/hg004_benchmark.vcf.gz"
)

BED_FILE_PATHS=(
    "/autofs/bal33/zxzheng/data/vcf/hg38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed"
    "/autofs/bal33/zxzheng/data/vcf/hg38/HG005_GRCh38_1_22_v4.2.1_benchmark.bed"
    # "/path/to/hg003_benchmark.bed"
    # "/path/to/hg004_benchmark.bed"
)

# Common Reference
REFERENCE_FILE_PATH="/autofs/bal19/zxzheng/testData/hifi/data/GRCh38_no_alt_analysis_set.fasta"

# Working directory (Base path)
INPUT_PATH="/autofs/nas5/xyu2/projects/Clair3_mv/hac_multisample"
mkdir -p ${INPUT_PATH}

################################################################################
# Environment Setup
################################################################################
echo "========================================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting Unified Training Pipeline (Multi-Sample)"
echo "========================================================================"

export JAVA_HOME=/autofs/bal33/zxzheng/software/java/jdk-17.0.2/
export PATH=$JAVA_HOME/bin:$PATH
export PERL5LIB=/autofs/bal33/zxzheng/nas2/software/vcftools_0.1.13/perl

################################################################################
# Tool Paths
################################################################################
CLAIR3_PATH="/autofs/nas5/xyu2/softwares/Clair3_no_change_v1.2"
CLAIR3="${CLAIR3_PATH}/clair3.py"
PYPY='/autofs/bal33/zxzheng/env/conda/envs/clair-somatic/bin/pypy3'
PYTHON3='/autofs/bal19/zxzheng/env/conda/envs/mamba/envs/somatic/bin/python3'
WHATSHAP='/autofs/bal33/zxzheng/env/conda/envs/clair-somatic/bin/whatshap'
PARALLEL='/autofs/bal33/zxzheng/env/conda/envs/clair-somatic/bin/parallel'
TABIX='/autofs/bal33/zxzheng/env/conda/envs/clair-somatic/bin/tabix'
SAMTOOLS='/autofs/bal33/zxzheng/env/conda/envs/clair-somatic/bin/samtools'
MOSDEPTH="/autofs/bal33/zxzheng/software/tool_bin/mosdepth"
VCF_SORT='/autofs/bal33/zxzheng/nas2/software/vcftools_0.1.13/bin/vcf-sort'
VCF_CONCAT='/autofs/bal33/zxzheng/nas2/software/vcftools_0.1.13/bin/vcf-concat'
RTG="/autofs/bal33/zxzheng/software/rtg-tools-3.9.1-6dde278/rtg"
HAP="/autofs/bal33/zxzheng/env/miniconda2/envs/singularity-env/bin/singularity exec -B /autofs/bal19/zxzheng,/autofs/bal33/zxzheng /autofs/bal19/zxzheng/singularity/hap.py_v0.3.12.sif /opt/hap.py/bin/hap.py"

################################################################################
# Processing Parameters
################################################################################
PLATFORM='ont'
CHR_PREFIX="chr"
CHR=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 21 22)
THREADS=80
THREADS_LOW=$((${THREADS}*3/4))
if [[ ${THREADS_LOW} < 1 ]]; then THREADS_LOW=1; fi

chunk_num=15
CHUNK_LIST=`seq 1 ${chunk_num}`

# AF thresholds
MIN_AF=0.08                     # For UnifyRepresentation
MIN_SNP_AF=0.07                 # For Pileup
MIN_INDEL_AF=0.08               # For Pileup
MIN_SNP_AF_FA=0.08              # For Full-alignment
MIN_INDEL_AF_FA=0.15            # For Full-alignment

# Training parameters
MAXIMUM_NON_VARIANT_RATIO_PILEUP=5
MAXIMUM_NON_VARIANT_RATIO_FA=3

################################################################################
# Helper Functions
################################################################################

check_file_exists() {
    if [ ! -f "$1" ]; then
        echo "[ERROR] Required file not found: $1"
        exit 1
    fi
}

check_dir_exists() {
    if [ ! -d "$1" ]; then
        echo "[ERROR] Required directory not found: $1"
        exit 1
    fi
}

log_step() {
    echo ""
    echo "========================================================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "========================================================================"
}

check_step_completion() {
    local marker_file="$1"
    local step_name="$2"
    
    if [ -f "${marker_file}" ]; then
        echo "[INFO] Step '${step_name}' already completed. Skipping..."
        return 0
    else
        return 1
    fi
}

mark_step_complete() {
    local marker_file="$1"
    local step_name="$2"
    
    touch "${marker_file}"
    echo "[INFO] Step '${step_name}' completed successfully."
}

################################################################################
# PHASE 1: Per-Sample Processing (RU, Downsample, Pileup Build)
################################################################################

for idx in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$idx]}"
    BAM_FILE_PATH="${BAM_FILE_PATHS[$idx]}"
    VCF_FILE_PATH="${VCF_FILE_PATHS[$idx]}"
    BED_FILE_PATH="${BED_FILE_PATHS[$idx]}"
    SAMPLE_UPPER=${SAMPLE^^}
    SAMPLE_NAME="${SAMPLE}"
    
    log_step "Processing Sample: ${SAMPLE} (${idx})"
    
    # Setup sample directory
    SAMPLE_DIR="${INPUT_PATH}/${SAMPLE}"
    mkdir -p "${SAMPLE_DIR}"
    
    # --- Depth Calculation ---
    echo "[INFO] Calculating depth information for ${SAMPLE}..."
    if [ ! -f ${SAMPLE_DIR}/${SAMPLE}.mosdepth.summary.txt ]; then
        echo "[INFO] Running mosdepth..."
        ${MOSDEPTH} -t 30 -n -x --quantize 0:15:150 ${SAMPLE_DIR}/${SAMPLE} ${BAM_FILE_PATH}
    fi
    
    tail -n 1 ${SAMPLE_DIR}/${SAMPLE}.mosdepth.summary.txt | cut -f 4 | xargs touch
    max_depth=$(tail -n 1 "${SAMPLE_DIR}/${SAMPLE}.mosdepth.summary.txt" | cut -f 4)
    echo "[INFO] RAW_DEPTH: ${max_depth}"
    
    # Calculate downsampling ratios for this sample
    EXPECT_DEPTHS=(60 50 40 30 20 10)
    DEPTHS=(1000)  # 1000 means no downsampling
    for depth in "${EXPECT_DEPTHS[@]}"; do
        pro=$(echo "scale=3; $depth / $max_depth * 1000" | bc)
        formatted_pro=$(printf "%03d" $(echo "$pro" | awk '{printf("%d\n", $1 + 0.5)}'))
        # echo "[INFO] Depth $depth -> Ratio $formatted_pro"
        DEPTHS+=("$formatted_pro")
    done
    echo "[INFO] SUBSAMPLE_RATIOS for ${SAMPLE}: ${DEPTHS[@]}"
    
    # --- STAGE 1: UnifyRepresentation (RU) ---
    if [ $RUN_RU -eq 1 ]; then
        RU_OUTPUT_DIR=${INPUT_PATH}/unifyrepre/${SAMPLE}
        PHASE_VCF_PATH="${RU_OUTPUT_DIR}/phased_vcf"
        PHASE_BAM_PATH="${RU_OUTPUT_DIR}/phased_bam"
        SPLIT_BED_PATH="${RU_OUTPUT_DIR}/split_beds"
        VAR_OUTPUT_PATH="${RU_OUTPUT_DIR}/var"
        VCF_OUTPUT_PATH="${RU_OUTPUT_DIR}/vcf_output"
        
        mkdir -p ${RU_OUTPUT_DIR} ${PHASE_VCF_PATH} ${PHASE_BAM_PATH} ${SPLIT_BED_PATH} ${VAR_OUTPUT_PATH} ${VCF_OUTPUT_PATH}
        cd ${RU_OUTPUT_DIR}
        
        # 1.1 Unphase
        if ! check_step_completion "${RU_OUTPUT_DIR}/.unphase_done" "Whatshap unphase"; then
            ${WHATSHAP} unphase ${VCF_FILE_PATH} > ${RU_OUTPUT_DIR}/INPUT.vcf.gz
            mark_step_complete "${RU_OUTPUT_DIR}/.unphase_done" "Whatshap unphase"
        fi
        
        # 1.2 Phase
        if ! check_step_completion "${RU_OUTPUT_DIR}/.phase_done" "Whatshap phase"; then
            ${PARALLEL} --joblog ${PHASE_VCF_PATH}/phase.log -j${THREADS} \
            "${WHATSHAP} phase --output ${PHASE_VCF_PATH}/phased_{1}.vcf.gz --reference ${REFERENCE_FILE_PATH} --chromosome ${CHR_PREFIX}{1} --ignore-read-groups --distrust-genotypes ${RU_OUTPUT_DIR}/INPUT.vcf.gz ${BAM_FILE_PATH}" ::: ${CHR[@]}
            mark_step_complete "${RU_OUTPUT_DIR}/.phase_done" "Whatshap phase"
        fi
        
        # 1.3 Index VCF
        if ! check_step_completion "${RU_OUTPUT_DIR}/.index_vcf_done" "Index VCF"; then
            ${PARALLEL} -j ${THREADS} ${TABIX} -p vcf ${PHASE_VCF_PATH}/phased_{1}.vcf.gz ::: ${CHR[@]}
            mark_step_complete "${RU_OUTPUT_DIR}/.index_vcf_done" "Index VCF"
        fi
        
        # 1.4 Haplotag
        if ! check_step_completion "${RU_OUTPUT_DIR}/.haplotag_done" "Whatshap haplotag"; then
            ${PARALLEL} --joblog ${PHASE_BAM_PATH}/haplotag.log -j${THREADS} \
            "${WHATSHAP} haplotag --output ${PHASE_BAM_PATH}/${SAMPLE_NAME}_{1}.bam --reference ${REFERENCE_FILE_PATH} --regions ${CHR_PREFIX}{1} --ignore-read-groups ${PHASE_VCF_PATH}/phased_{1}.vcf.gz ${BAM_FILE_PATH}" ::: ${CHR[@]}
            mark_step_complete "${RU_OUTPUT_DIR}/.haplotag_done" "Whatshap haplotag"
        fi
        
        # 1.5 Index BAM
        if ! check_step_completion "${RU_OUTPUT_DIR}/.index_bam_done" "Index BAM"; then
            time ${PARALLEL} --joblog ${PHASE_BAM_PATH}/index.log -j ${THREADS} ${SAMTOOLS} index -@12 ${PHASE_BAM_PATH}/${SAMPLE_NAME}_{1}.bam ::: ${CHR[@]}
            mark_step_complete "${RU_OUTPUT_DIR}/.index_bam_done" "Index BAM"
        fi
        
        # 1.6 Split extend bed
        if ! check_step_completion "${RU_OUTPUT_DIR}/.split_bed_done" "Split extend bed"; then
            ${PARALLEL} --joblog ${SPLIT_BED_PATH}/split_extend_bed.log -j${THREADS} \
            "${PYPY} ${CLAIR3} SplitExtendBed --bed_fn ${BED_FILE_PATH} --output_fn ${SPLIT_BED_PATH}/{1} --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]}
            mark_step_complete "${RU_OUTPUT_DIR}/.split_bed_done" "Split extend bed"
        fi
        
        # 1.7 Get truth
        if ! check_step_completion "${RU_OUTPUT_DIR}/.get_truth_done" "Get truth"; then
            ${PARALLEL} --joblog ${VAR_OUTPUT_PATH}/get_truth.log -j${THREADS} \
            "${PYPY} ${CLAIR3} GetTruth --vcf_fn ${VCF_FILE_PATH} --ctgName ${CHR_PREFIX}{1} --var_fn ${VAR_OUTPUT_PATH}/var_{1}" ::: ${CHR[@]}
            mark_step_complete "${RU_OUTPUT_DIR}/.get_truth_done" "Get truth"
        fi
        
        # 1.8 UnifyRepresentation
        if ! check_step_completion "${RU_OUTPUT_DIR}/.unify_repre_done" "UnifyRepresentation"; then
            time ${PARALLEL} --joblog ${RU_OUTPUT_DIR}/unify_repre.log -j15 \
            "${PYPY} ${CLAIR3} UnifyRepresentation --bam_fn ${PHASE_BAM_PATH}/${SAMPLE_NAME}_{1}.bam --var_fn ${VAR_OUTPUT_PATH}/var_{1} --ref_fn ${REFERENCE_FILE_PATH} --bed_fn ${BED_FILE_PATH} --extend_bed ${SPLIT_BED_PATH}/{1} --output_vcf_fn ${VCF_OUTPUT_PATH}/vcf_{1}_{2} --min_af ${MIN_AF} --chunk_id {2} --chunk_num ${chunk_num} --platform ${PLATFORM} --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]} ::: ${CHUNK_LIST[@]}
            mark_step_complete "${RU_OUTPUT_DIR}/.unify_repre_done" "UnifyRepresentation"
        fi
        
        # 1.9 Sort VCF
        if ! check_step_completion "${RU_OUTPUT_DIR}/.sort_vcf_done" "Sort VCF"; then
            cat ${VCF_OUTPUT_PATH}/vcf_* | ${PYPY} ${CLAIR3} SortVcf --output_fn ${RU_OUTPUT_DIR}/unified.vcf
            bgzip -f ${RU_OUTPUT_DIR}/unified.vcf
            tabix -f -p vcf ${RU_OUTPUT_DIR}/unified.vcf.gz
            mark_step_complete "${RU_OUTPUT_DIR}/.sort_vcf_done" "Sort VCF"
        fi
    else
        RU_OUTPUT_DIR=${INPUT_PATH}/unifyrepre/${SAMPLE}
        PHASE_BAM_PATH="${RU_OUTPUT_DIR}/phased_bam"
    fi

    # --- STAGE 1.5: Downsample ---
    if [ $RUN_DOWNSAMPLE -eq 1 ]; then
        if ! check_step_completion "${RU_OUTPUT_DIR}/.downsample_done" "Downsampling"; then
            # Links for 1000
            ${PARALLEL} ln -sf ${PHASE_BAM_PATH}/${SAMPLE_NAME}_{1}.bam ${PHASE_BAM_PATH}/${SAMPLE_NAME}_1000_{1}.bam ::: ${CHR[@]}
            ${PARALLEL} ln -sf ${PHASE_BAM_PATH}/${SAMPLE_NAME}_{1}.bam.bai ${PHASE_BAM_PATH}/${SAMPLE_NAME}_1000_{1}.bam.bai ::: ${CHR[@]}
            
            # Downsample
            ${PARALLEL} -j18 --joblog ${RU_OUTPUT_DIR}/parallel_downsample.log \
            "samtools view -@12 -s {2}.{2} -b -o ${PHASE_BAM_PATH}/${SAMPLE_NAME}_{2}_{1}.bam ${PHASE_BAM_PATH}/${SAMPLE_NAME}_1000_{1}.bam ${CHR_PREFIX}{1}" ::: ${CHR[@]} ::: ${DEPTHS[@]:1}
            
            # Index
            ${PARALLEL} -j24 "samtools index -@16 ${PHASE_BAM_PATH}/${SAMPLE_NAME}_{2}_{1}.bam" ::: ${CHR[@]} ::: ${DEPTHS[@]}
            
            mark_step_complete "${RU_OUTPUT_DIR}/.downsample_done" "Downsampling"
        fi
    fi

    # --- STAGE 2.1: Pileup Build ---
    if [ $BUILD_PILEUP -eq 1 ]; then
        PILEUP_OUTPUT_DIR="${INPUT_PATH}/pileup/${SAMPLE}"
        PILEUP_DATASET_PATH="${PILEUP_OUTPUT_DIR}/build"
        PILEUP_TENSOR_PATH="${PILEUP_DATASET_PATH}/tensor_can"
        PILEUP_BINS_PATH="${PILEUP_DATASET_PATH}/bins"
        PILEUP_SPLIT_BED_PATH="${PILEUP_DATASET_PATH}/split_beds"
        PILEUP_VAR_PATH="${PILEUP_DATASET_PATH}/var"
        
        mkdir -p ${PILEUP_DATASET_PATH}/tmp ${PILEUP_TENSOR_PATH} ${PILEUP_BINS_PATH} ${PILEUP_SPLIT_BED_PATH} ${PILEUP_VAR_PATH}
        
        UNIFIED_VCF_FILE_PATH="${RU_OUTPUT_DIR}/unified.vcf.gz"
        
        # Split Bed
        if ! check_step_completion "${PILEUP_DATASET_PATH}/.split_bed_done" "Pileup split bed"; then
            ${PARALLEL} --joblog ${PILEUP_DATASET_PATH}/parallel_1_split_extend_bed.log -j${THREADS} \
            "${PYPY} ${CLAIR3} SplitExtendBed --bed_fn ${BED_FILE_PATH} --output_fn ${PILEUP_SPLIT_BED_PATH}/${SAMPLE_NAME}_{1} --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]}
            mark_step_complete "${PILEUP_DATASET_PATH}/.split_bed_done" "Pileup split bed"
        fi
        
        # Get Truth
        if ! check_step_completion "${PILEUP_DATASET_PATH}/.get_truth_done" "Pileup get truth"; then
            ${PARALLEL} --joblog ${PILEUP_DATASET_PATH}/parallel_2_get_truth.log -j${THREADS} \
            "${PYPY} ${CLAIR3} GetTruth --vcf_fn ${UNIFIED_VCF_FILE_PATH} --ctgName ${CHR_PREFIX}{1} --truth_vcf_fn ${VCF_FILE_PATH} --var_fn ${PILEUP_VAR_PATH}/var_${SAMPLE_NAME}_1000_{1}" ::: ${CHR[@]}
            mark_step_complete "${PILEUP_DATASET_PATH}/.get_truth_done" "Pileup get truth"
        fi
        
        # Create Tensor
        if ! check_step_completion "${PILEUP_DATASET_PATH}/.create_tensor_done" "Pileup create tensor"; then
            time ${PARALLEL} --joblog ${PILEUP_DATASET_PATH}/parallel_3_create_tensor.log -j${THREADS_LOW} \
            "${PYPY} ${CLAIR3} CreateTrainingTensor --bam_fn ${PHASE_BAM_PATH}/${SAMPLE_NAME}_{2}_{1}.bam --ref_fn ${REFERENCE_FILE_PATH} --var_fn ${PILEUP_VAR_PATH}/var_${SAMPLE_NAME}_1000_{1} --bin_fn ${PILEUP_TENSOR_PATH}/tensor_${SAMPLE_NAME}_{2}_{1}_{3} --ctgName ${CHR_PREFIX}{1} --samtools ${SAMTOOLS} --extend_bed ${PILEUP_SPLIT_BED_PATH}/${SAMPLE_NAME}_{1} --bed_fn ${BED_FILE_PATH} --allow_duplicate_chr_pos --platform ${PLATFORM} --python ${PYTHON3} --shuffle --pileup --maximum_non_variant_ratio ${MAXIMUM_NON_VARIANT_RATIO_PILEUP} --chunk_id {3} --chunk_num ${chunk_num} --pypy ${PYPY}" ::: ${CHR[@]} ::: ${DEPTHS[@]} ::: ${CHUNK_LIST[@]}
            mark_step_complete "${PILEUP_DATASET_PATH}/.create_tensor_done" "Pileup create tensor"
        fi
        
        # Merge Bin
        if ! check_step_completion "${PILEUP_DATASET_PATH}/.merge_bin_done" "Pileup merge bin"; then
            time ${PARALLEL} --joblog ${PILEUP_DATASET_PATH}/parallel_4_mergeBin.log -j${THREADS} \
            "${PYTHON3} ${CLAIR3} MergeBin ${PILEUP_TENSOR_PATH}/tensor_${SAMPLE_NAME}_{1}_* --pileup --out_fn ${PILEUP_BINS_PATH}/bin_${SAMPLE_NAME}_{1}" ::: ${DEPTHS[@]}
            mark_step_complete "${PILEUP_DATASET_PATH}/.merge_bin_done" "Pileup merge bin"
        fi
    fi
done

################################################################################
# PHASE 2: Pileup Training (Combined)
################################################################################

if [ $BUILD_PILEUP -eq 1 ] || [ $TRAIN_PILEUP -eq 1 ]; then
    PILEUP_ALL_BINS_PATH="${INPUT_PATH}/pileup/all/build/bins"
    mkdir -p ${PILEUP_ALL_BINS_PATH}
    
    if [ $BUILD_PILEUP -eq 1 ]; then
        echo "[INFO] Linking pileup bins from all samples..."
        for idx in "${!SAMPLES[@]}"; do
            SAMPLE="${SAMPLES[$idx]}"
            SAMPLE_BINS_PATH="${INPUT_PATH}/pileup/${SAMPLE}/build/bins"
            if [ -d "${SAMPLE_BINS_PATH}" ]; then
                ln -sf ${SAMPLE_BINS_PATH}/* ${PILEUP_ALL_BINS_PATH}/
            fi
        done
    fi
    
    if [ $TRAIN_PILEUP -eq 1 ]; then
        log_step "STAGE 2.2: Training Pileup Model (All Samples)"
        PILEUP_MODEL_PATH="${INPUT_PATH}/pileup/all/train/${TRAIN_NAME}/${TRAIN_NAME}"
        mkdir -p ${PILEUP_MODEL_PATH}
        cd ${PILEUP_MODEL_PATH}
        
        if ! check_step_completion "${PILEUP_MODEL_PATH}/.train_done" "Pileup training"; then
            if [ $USE_GPU -eq 1 ]; then export CUDA_VISIBLE_DEVICES="${GPU_DEVICE}"; else export CUDA_VISIBLE_DEVICES=""; fi
            
            ${PYTHON3} ${CLAIR3} Train --bin_fn ${PILEUP_ALL_BINS_PATH} --ochk_prefix ${PILEUP_MODEL_PATH} --add_indel_length False --random_validation --pileup --platform ${PLATFORM} |& tee ${PILEUP_MODEL_PATH}/log.txt
            mark_step_complete "${PILEUP_MODEL_PATH}/.train_done" "Pileup training"
        fi
    fi
else
    PILEUP_MODEL_PATH="${INPUT_PATH}/pileup/all/train/${TRAIN_NAME}/${TRAIN_NAME}"
fi

################################################################################
# PHASE 3: Per-Sample Processing (FA Build)
################################################################################

for idx in "${!SAMPLES[@]}"; do
    SAMPLE="${SAMPLES[$idx]}"
    BAM_FILE_PATH="${BAM_FILE_PATHS[$idx]}"
    VCF_FILE_PATH="${VCF_FILE_PATHS[$idx]}"
    BED_FILE_PATH="${BED_FILE_PATHS[$idx]}"
    SAMPLE_UPPER=${SAMPLE^^}
    SAMPLE_NAME="${SAMPLE}"
    
    # Re-calculate depths
    SAMPLE_DIR="${INPUT_PATH}/${SAMPLE}"
    max_depth=$(tail -n 1 "${SAMPLE_DIR}/${SAMPLE}.mosdepth.summary.txt" | cut -f 4)
    EXPECT_DEPTHS=(60 50 40 30 20 10)
    DEPTHS=(1000)
    for depth in "${EXPECT_DEPTHS[@]}"; do
        pro=$(echo "scale=3; $depth / $max_depth * 1000" | bc)
        formatted_pro=$(printf "%03d" $(echo "$pro" | awk '{printf("%d\n", $1 + 0.5)}'))
        DEPTHS+=("$formatted_pro")
    done
    
    RU_OUTPUT_DIR=${INPUT_PATH}/unifyrepre/${SAMPLE}
    PHASE_BAM_PATH="${RU_OUTPUT_DIR}/phased_bam"
    UNIFIED_VCF_FILE_PATH="${RU_OUTPUT_DIR}/unified.vcf.gz"

    # --- STAGE 3.1: Full-Alignment Build ---
    if [ $BUILD_FULL_ALIGNMENT -eq 1 ]; then
        log_step "Building Full-Alignment Dataset for ${SAMPLE}"
        
        FA_OUTPUT_DIR="${INPUT_PATH}/full_alignment/${SAMPLE}"
        FA_DATASET_PATH="${FA_OUTPUT_DIR}/build"
        FA_TENSOR_PATH="${FA_DATASET_PATH}/tensor_can"
        FA_BINS_PATH="${FA_DATASET_PATH}/bins"
        FA_SPLIT_BED_PATH="${FA_DATASET_PATH}/split_beds"
        FA_VAR_PATH="${FA_DATASET_PATH}/var"
        FA_CANDIDATE_BED_PATH="${FA_DATASET_PATH}/candidate_bed"
        FA_PILEUP_OUTPUT_PATH="${FA_DATASET_PATH}/pileup_output"
        
        mkdir -p ${FA_DATASET_PATH}/tmp ${FA_TENSOR_PATH} ${FA_BINS_PATH} ${FA_SPLIT_BED_PATH} ${FA_VAR_PATH} ${FA_CANDIDATE_BED_PATH} ${FA_PILEUP_OUTPUT_PATH}
        
        PILEUP_MODEL_FOR_CALLING="${INPUT_PATH}/pileup/all/train/${TRAIN_NAME}/${TRAIN_NAME}/best_val_loss/variables/"
        if [ ! -d "${PILEUP_MODEL_FOR_CALLING}" ]; then
            echo "[ERROR] Pileup model not found at ${PILEUP_MODEL_FOR_CALLING}"
            exit 1
        fi
        
        # 3.1 Pileup calling
        if ! check_step_completion "${FA_DATASET_PATH}/.pileup_calling_done" "Pileup calling"; then
             ${PARALLEL} -j3 --joblog ${FA_PILEUP_OUTPUT_PATH}/pileup_calling.log \
            "docker run --rm --gpus all -u $(id -u):$(id -g) -v /autofs/nas5/xyu2/:/autofs/nas5/xyu2/ -v /autofs/bal19/zxzheng/:/autofs/bal19/zxzheng/ -v /autofs/bal33/zxzheng/:/autofs/bal33/zxzheng/ hkubal/clair3-gpu:latest /opt/bin/run_clair3.sh --bam_fn=${PHASE_BAM_PATH}/${SAMPLE_NAME}_{2}_{1}.bam --ref_fn=${REFERENCE_FILE_PATH} --threads=30 --platform='ont' --model_path=${PILEUP_MODEL_FOR_CALLING} --pileup_model_prefix=variables --fa_model_prefix=variables --output=${FA_PILEUP_OUTPUT_PATH}/${SAMPLE_NAME}_{2}_{1} --bed_fn=${BED_FILE_PATH} --snp_min_af=${MIN_SNP_AF_FA} --indel_min_af=${MIN_INDEL_AF_FA} --print_ref_calls --pileup_only" ::: ${CHR[@]} ::: ${DEPTHS[@]}
            mark_step_complete "${FA_DATASET_PATH}/.pileup_calling_done" "Pileup calling"
        fi
        
        # 3.2 Concat VCF
        if ! check_step_completion "${FA_DATASET_PATH}/.concat_vcf_done" "Concatenate VCF"; then
            ${PARALLEL} mkdir -p ${FA_PILEUP_OUTPUT_PATH}/${SAMPLE}_{1} ::: ${DEPTHS[@]}
            ${PARALLEL} -j10 --joblog ${FA_PILEUP_OUTPUT_PATH}/vcf_concat.log \
            "${VCF_CONCAT} ${FA_PILEUP_OUTPUT_PATH}/${SAMPLE_NAME}_{1}_*/pileup.vcf.gz | ${VCF_SORT} > ${FA_PILEUP_OUTPUT_PATH}/${SAMPLE}_{1}/pileup.vcf" ::: ${DEPTHS[@]}
            ${PARALLEL} -j10 bgzip -f ${FA_PILEUP_OUTPUT_PATH}/${SAMPLE}_{1}/pileup.vcf ::: ${DEPTHS[@]}
            ${PARALLEL} -j10 tabix -p vcf -f ${FA_PILEUP_OUTPUT_PATH}/${SAMPLE}_{1}/pileup.vcf.gz ::: ${DEPTHS[@]}
            mark_step_complete "${FA_DATASET_PATH}/.concat_vcf_done" "Concatenate VCF"
        fi
        
        # 3.3 Select candidates
        if ! check_step_completion "${FA_DATASET_PATH}/.select_candidates_done" "Select candidates"; then
            ${PARALLEL} --joblog ${FA_DATASET_PATH}/select_pileup_candidates.log -j${THREADS} \
            "${PYPY} ${CLAIR3} SelectHetSnp --alt_fn ${FA_PILEUP_OUTPUT_PATH}/${SAMPLE}_{2}/pileup.vcf.gz --split_folder ${FA_CANDIDATE_BED_PATH} --sampleName ${SAMPLE} --depth {2} --ref_pct_full 0.15 --var_pct_full 1.0 --chunk_num ${chunk_num} --phasing_info_in_bam --phase --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]} ::: ${DEPTHS[@]}
            mark_step_complete "${FA_DATASET_PATH}/.select_candidates_done" "Select candidates"
        fi
        
        # 3.4 Split extend bed
        if ! check_step_completion "${FA_DATASET_PATH}/.split_bed_done" "FA split bed"; then
            ${PARALLEL} --joblog ${FA_DATASET_PATH}/parallel_1_split_extend_bed.log -j${THREADS} \
            "${PYPY} ${CLAIR3} SplitExtendBed --bed_fn ${BED_FILE_PATH} --output_fn ${FA_SPLIT_BED_PATH}/${SAMPLE_NAME}_{1} --ctgName ${CHR_PREFIX}{1}" ::: ${CHR[@]}
            mark_step_complete "${FA_DATASET_PATH}/.split_bed_done" "FA split bed"
        fi
        
        # 3.5 Get truth
        if ! check_step_completion "${FA_DATASET_PATH}/.get_truth_done" "FA get truth"; then
            ${PARALLEL} --joblog ${FA_DATASET_PATH}/parallel_2_get_truth.log -j${THREADS} \
            "${PYPY} ${CLAIR3} GetTruth --vcf_fn ${UNIFIED_VCF_FILE_PATH} --ctgName ${CHR_PREFIX}{1} --truth_vcf_fn ${VCF_FILE_PATH} --var_fn ${FA_VAR_PATH}/var_${SAMPLE_NAME}_1000_{1}" ::: ${CHR[@]}
            mark_step_complete "${FA_DATASET_PATH}/.get_truth_done" "FA get truth"
        fi
        
        # 3.6 Create Tensor
        if ! check_step_completion "${FA_DATASET_PATH}/.create_tensor_done" "FA create tensor"; then
            time ${PARALLEL} --joblog ${FA_DATASET_PATH}/parallel_3_create_tensor.log --resume -j${THREADS_LOW} \
            "${PYPY} ${CLAIR3} CreateTrainingTensor --bam_fn ${PHASE_BAM_PATH}/${SAMPLE_NAME}_{2}_{1}.bam --ref_fn ${REFERENCE_FILE_PATH} --var_fn ${FA_VAR_PATH}/var_${SAMPLE_NAME}_1000_{1} --bin_fn ${FA_TENSOR_PATH}/tensor_${SAMPLE_NAME}_{2}_{1}_{3} --ctgName ${CHR_PREFIX}{1} --samtools ${SAMTOOLS} --extend_bed ${FA_SPLIT_BED_PATH}/${SAMPLE_NAME}_{1} --full_aln_regions ${FA_CANDIDATE_BED_PATH}/${SAMPLE}_{2}_{1}_{3} --bed_fn ${BED_FILE_PATH} --phasing_info_in_bam --add_no_phasing_data_training --allow_duplicate_chr_pos --platform ${PLATFORM} --python ${PYTHON3} --pypy ${PYPY} --shuffle --maximum_non_variant_ratio ${MAXIMUM_NON_VARIANT_RATIO_FA} --chunk_id {3} --chunk_num ${chunk_num} --enable_dwell_time" ::: ${CHR[@]} ::: ${DEPTHS[@]} ::: ${CHUNK_LIST[@]}
            mark_step_complete "${FA_DATASET_PATH}/.create_tensor_done" "FA create tensor"
        fi
        
        # 3.7 Merge Bin
        if ! check_step_completion "${FA_DATASET_PATH}/.merge_bin_done" "FA merge bin"; then
            ${PARALLEL} --joblog ${FA_DATASET_PATH}/parallel_4_mergeBin.log -j${THREADS} \
            "${PYTHON3} ${CLAIR3} MergeBin ${FA_TENSOR_PATH}/tensor_${SAMPLE_NAME}_{1}_* --out_fn ${FA_BINS_PATH}/bin_${SAMPLE_NAME}_{1} --enable_dwell_time" ::: ${DEPTHS[@]}
            mark_step_complete "${FA_DATASET_PATH}/.merge_bin_done" "FA merge bin"
        fi
    fi
done

################################################################################
# PHASE 4: Full-Alignment Training (Combined)
################################################################################

if [ $BUILD_FULL_ALIGNMENT -eq 1 ] || [ $TRAIN_FULL_ALIGNMENT -eq 1 ]; then
    FA_ALL_BINS_PATH="${INPUT_PATH}/full_alignment/all/build/bins"
    mkdir -p ${FA_ALL_BINS_PATH}
    
    if [ $BUILD_FULL_ALIGNMENT -eq 1 ]; then
        echo "[INFO] Linking full-alignment bins from all samples..."
        for idx in "${!SAMPLES[@]}"; do
            SAMPLE="${SAMPLES[$idx]}"
            SAMPLE_BINS_PATH="${INPUT_PATH}/full_alignment/${SAMPLE}/build/bins"
            if [ -d "${SAMPLE_BINS_PATH}" ]; then
                ln -sf ${SAMPLE_BINS_PATH}/* ${FA_ALL_BINS_PATH}/
            fi
        done
    fi
    
    if [ $TRAIN_FULL_ALIGNMENT -eq 1 ]; then
        log_step "STAGE 3.2: Training Full-Alignment Model (All Samples)"
        FA_MODEL_PATH="${INPUT_PATH}/full_alignment/all/train/${TRAIN_NAME}/${TRAIN_NAME}"
        mkdir -p ${FA_MODEL_PATH}
        cd ${FA_MODEL_PATH}
        
        if ! check_step_completion "${FA_MODEL_PATH}/.train_done" "FA training"; then
            if [ $USE_GPU -eq 1 ]; then export CUDA_VISIBLE_DEVICES="${GPU_DEVICE}"; else export CUDA_VISIBLE_DEVICES=""; fi
            
            ${PYTHON3} ${CLAIR3} Train --bin_fn ${FA_ALL_BINS_PATH} --ochk_prefix ${FA_MODEL_PATH} --add_indel_length True --random_validation --platform ${PLATFORM} --enable_dwell_time |& tee ${FA_MODEL_PATH}/log.txt
            mark_step_complete "${FA_MODEL_PATH}/.train_done" "FA training"
        fi
    fi
fi

echo "========================================================================"
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Multi-Sample Unified Training Pipeline Completed!"
echo "========================================================================"

