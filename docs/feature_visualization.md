# Clair3 Feature Visualization



## Introduction

This document shows how to visualize Clair3 input tensor and  feature maps. For input tensor, we encode the BAM and REF files to get full-alignment  input tensor. For feature maps visualization, we generated heatmaps of importance assigned by Clair3 to individual pixels in its input using [guided backpropagation](https://arxiv.org/abs/1311.2901), which is a technique that uses model gradients to assign importance scores to pixels. 

----

## Prerequisites

- Clair3 installed 
- Matplotlib installed

----

## Input/Output

**Input:**

- BAM: Indexed BAM input
- REF: Indexed Reference input
- CHR: The sequence name
- POS: a  genomic position to visualize 

**Ouput:**

- Output feature maps and tensor visualization

####  1. Setup variables

```bash
# Setup variables
PLATFORM='ont'
BAM_FILE_PATH="input.bam"
REF_FILE_PATH="input.fasta"
CHR="chr1"
POS=14160735
OUTPUT_PATH="vis_output"
SAMTOOLS='samtools'
PYPY="pypy3"
PYTHON="python3"
CLAIR3='Clair3/clair3.py'
```

####  2. Download pre-trained models

```bash
mkdir models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz 
tar -zxvf clair3_models.tar.gz -C ./models
CLIAR3_MODEL_PATH=`pwd`"/models/r941_prom_sup_g5014"
```

#### 3.  Generate input binary for visualization

```bash
mkdir -p ${OUTPUT_PATH}
touch ${OUTPUT_PATH}/var
START_POS=`expr ${POS} - 16`
END_POS=`expr ${POS} + 18`
echo -e "${CHR}\t${START_POS}\t${END_POS}" > ${OUTPUT_PATH}/input.bed

${PYPY} ${CLAIR3} CreateTrainingTensor \
    --bam_fn ${BAM_FILE_PATH} \
    --ref_fn ${REF_FILE_PATH} \
    --var_fn ${OUTPUT_PATH}/var \
    --bin_fn ${OUTPUT_PATH}/output_bin \
    --ctgName ${CHR} \
    --extend_bed ${OUTPUT_PATH}/input.bed \
    --bed_fn ${OUTPUT_PATH}/input.bed \
    --samtools ${SAMTOOLS} \
    --full_aln_regions ${OUTPUT_PATH}/input.bed \
    --phasing_info_in_bam \
    --allow_duplicate_chr_pos \
    --platform ${PLATFORM} \
    --shuffle \
    --chunk_id 0 \
    --chunk_num 1
```

#### 4.  Output visualization figures

```bash
${PYTHON} ${CLAIR3} VisualizeFeature \
	--platform ${PLATFORM} \
    --model_path ${CLIAR3_MODEL_PATH} \
    --input_bin ${OUTPUT_PATH}/output_bin \
    --output_dir ${OUTPUT_PATH}
```

Check the feature maps and tensor figure in `${OUTPUT_PATH}/guided_backpropagation_feature_maps.png` and `${OUTPUT_PATH}/tensor.png`.

Here is a [jupyter notebook](http://www.bio8.cs.hku.hk/clair3/scripts_answering_reviewers_comments/rex_feature_visualization/) visualization example for in four categories, SNP, Insertion, Deletion, and Non-variant site of ONT data.

