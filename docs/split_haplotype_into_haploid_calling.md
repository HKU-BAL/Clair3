# Run two rounds of Clair3 to first split reads into haplotypes and then do variant calling on single haplotypes

## Introduction

This document provides a step-by-step guide to first split reads into haplotypes according to phased variants and then do variant calling on single haplotypes.

## Prerequisites

- Clair3 installed 
- GNU Parallel installed
- Pypy3 installed

----

## Input/Output

**Input:**

- BAM: Indexed BAM input
- REF: Indexed Reference input

**Ouput:**

- Two output VCF files of from two haplotagged BAM files for downstream analysis

----

####  1. Run Clair3 to acquire the phased output VCF

```bash
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  python3 /opt/bin/run_clair3.py \
  --bam_fn=${INPUT_DIR}/input.bam \    
  --ref_fn=${INPUT_DIR}/ref.fa \      
  --threads=${THREADS} \               
  --platform="ont" \                   
  --model_path="/opt/models/${MODEL_NAME}" \
  --use_whatshap_for_final_output_phasing \ ## enable phased output
  --sample_name=${SAMPLE} \
  --output=${OUTPUT_DIR}
```

#### 2. Use WhatsHap to haplotag and split BAM into two haplotagged BAM files
```bash
cd ${OUTPUT_DIR}

# haplotag read alignment using WhatsHap
whatshap haplotag \
  --output-threads ${THREADS} \
  -o ${OUTPUT_DIR}/output_haplotagged.bam \
  --reference ${INPUT_DIR}/ref.fa \
  --ignore-read-groups \
  --skip-missing-contigs \
  --output-haplotag-list ${OUTPUT_DIR}/${SAMPLE}.split.tsv \
  ${OUTPUT_DIR}/${SAMPLE}.phased.vcf.gz ${INPUT_DIR}/input.bam

# Split haplotagged BAM file into two haplotagged BAM files
whatshap split \
  --output-h1 ${OUTPUT_DIR}/${SAMPLE}.h1.bam \
  --output-h2 ${OUTPUT_DIR}/${SAMPLE}.h2.bam \
  ${OUTPUT_DIR}/output_haplotagged.bam \
  ${OUTPUT_DIR}/${SAMPLE}.split.tsv
 
#Index two haplotagged BAM files
samtools index ${OUTPUT_DIR}/${SAMPLE}.h1.bam
samtools index ${OUTPUT_DIR}/${SAMPLE}.h2.bam
```

#### 3.  Run Clair3 on two haplotagged BAMs

```bash
# run Clair3 on haplotagged BAM file 1
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  python3 /opt/bin/run_clair3.py \
  --bam_fn=${OUTPUT_DIR}/${SAMPLE}.h1.bam \    
  --ref_fn=${INPUT_DIR}/ref.fa \      
  --threads=${THREADS} \               
  --platform="ont" \                   
  --model_path="/opt/models/${MODEL_NAME}" \
  --sample_name=${SAMPLE}_h1 \
  --output=${OUTPUT_DIR}_2
  
# run Clair3 on haplotagged BAM file 2
docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  python3 /opt/bin/run_clair3.py \
  --bam_fn=${OUTPUT_DIR}/${SAMPLE}.h1.bam \    
  --ref_fn=${INPUT_DIR}/ref.fa \      
  --threads=${THREADS} \               
  --platform="ont" \                   
  --model_path="/opt/models/${MODEL_NAME}" \
  --sample_name=${SAMPLE}_h2 \
  --output=${OUTPUT_DIR}_h2
```



