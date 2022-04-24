<div align="center">
  <a href="https://en.wiktionary.org/wiki/%E7%9C%BC" target="_blank">
    <img src="docs/images/clair3_logo.png" width = "110" height = "90" alt="Clair3">
  </a>
</div>

# Clair3 - Symphonizing pileup and full-alignment for high-performance long-read variant calling

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)  [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/clair3/README.html)

Contact: Ruibang Luo, Zhenxian Zheng  
Email: rbluo@cs.hku.hk, zxzheng@cs.hku.hk  

----

## Introduction

Clair3 is a germline small variant caller for long-reads. Clair3 makes the best of two major method categories: pileup calling handles most variant candidates with speed, and full-alignment tackles complicated candidates to maximize precision and recall. Clair3 runs fast and has superior performance, especially at lower coverage. Clair3 is simple and modular for easy deployment and integration.

Clair3 is the 3<sup>rd</sup> generation of [Clair](https://github.com/HKU-BAL/Clair) (the 2<sup>nd</sup>) and [Clairvoyante](https://github.com/aquaskyline/Clairvoyante) (the 1<sup>st</sup>).

A short preprint describing Clair3's algorithms and results is at [bioRxiv](https://www.biorxiv.org/content/10.1101/2021.12.29.474431v1).

----

## Contents

* [Introduction](#introduction)
* [Latest Updates](#latest-updates)
* [Pre-trained Models](#pre-trained-models)
  * [Guppy5,6 Model](docs/guppy5_20220113.md)
  * [R10.4 with the Kit 12 chemistry (Q20) Models](#ont-provided-models)
  * [Guppy3,4 Model](#pre-trained-models)
  * [Guppy2 Model](docs/guppy2.md)
* [What's New in Clair3](#whats-new-in-clair3)
* [Installation](#installation)
  + [Option 1. Docker pre-built image](#option-1--docker-pre-built-image)
  + [Option 2. Singularity](#option-2-singularity)
  + [Option 3. Bioconda](#option-3--bioconda)
  + [Option 4. Build an anaconda virtual environment](#option-4-build-an-anaconda-virtual-environment)
  + [Option 5. Docker Dockerfile](#option-5-docker-dockerfile)
* [Quick Demo](#quick-demo)
* [Usage](#usage)
* [Folder Structure and Submodule Descriptions](#folder-structure-and-submodule-descriptions)
* [Training Data](docs/training_data.md)
* [VCF/GVCF Output Formats](#vcfgvcf-output-formats)
* [Pileup Model Training](docs/pileup_training.md)
* [Full-Alignment Model Training](docs/full_alignment_training_r1.md)
* [Representation Unification](docs/representation_unification.md)
* [Visualization](docs)
  * [Model Input](docs/model_input_visualization.md)
  * [Representation Unification](docs/representation_unification_visualization.md)

----

## Latest Updates

*v0.1-r11 minor 2 (Apr 16)* : 1. fixed a bug in GVCF output that occasionally caused missing of non-variant positions at chunk boundaries. 2. fixed a bug in GVCF output that consumes too much memory for caching, now GVCF output mode takes amount of memory similar to VCF ([#88](https://github.com/HKU-BAL/Clair3/issues/88)).

*v0.1-r11 (Apr 4)* : 1. Variant calling ~2.5x faster than `v0.1-r10` tested with ONT Q20 data, with feature generation in both pileup and full-alignment now implemented in C (co-contributors @[cjw85](https://github.com/cjw85), @[ftostevin-ont](https://github.com/ftostevin-ont), @[EpiSlim](https://github.com/EpiSlim)). 2. Added the lightning-fast [longphase](https://github.com/twolinin/longphase) as an option for phasing. Enable using `longphase` with option `--longphase_for_phasing`. New option disabled by default to align with the default behavior of the previous versions, but we recommend enable when calling human variants with ≥20x long-reads). 3. Added `--min_coverage` and `--min_mq` options ([#83](https://github.com/HKU-BAL/Clair3/issues/83)). 4. Added `--min_contig_size` option to skip calling variants in short contigs when using genome assembly as input. 4. Reads haplotagging after phasing before full-alignment calling now integrated into full-alignment calling to avoid generating an intermediate BAM file. 5. Supported .`csi` BAM index for large references ([#90](https://github.com/HKU-BAL/Clair3/issues/90)). For more speedup details, please check [Notes on r11](docs/v0.1_r11_speedup.md).

*v0.1-r10 (Jan 13, 2022)* : 1. Added a new ONT Guppy5 model  (`r941_prom_sup_g5014`). Click [here](docs/guppy5_20220113.md) for some benchmarking results. This `sup` model is also applicable to reads called using the `hac` and `fast` mode. The old `r941_prom_sup_g506` model that was fine-tuned from the Guppy3,4 model is obsoleted. 2. Added `--var_pct_phasing` option to control the percentage of top ranked heterozygous pile-up variants used for WhatsHap phasing.

*v0.1-r9 (Dec 1)* : Added the `--enable_long_indel` option to output indel variant calls >50bp ([#64](https://github.com/HKU-BAL/Clair3/issues/64)), Click [here](https://github.com/HKU-BAL/Clair3/blob/main/docs/indel_gt50_performance.md) to see more benchmarking results.

*v0.1-r8 (Nov 11)* : 1. Added the `--enable_phasing` option that adds a step after Clair3 calling to output variants phased by WhatsHap ([#63](https://github.com/HKU-BAL/Clair3/issues/63)). 2. Fixed unexpected program termination on successful runs.

*v0.1-r7 (Oct 18)* : 1. Increased `var_pct_full` in ONT mode from 0.3 to 0.7. Indel F1-score increased ~0.2%, but took ~30 minutes longer to finish calling a ~50x ONT dataset. 2. Expand fall through to next most likely variant if network prediction has insufficient read coverage ([#53](https://github.com/HKU-BAL/Clair3/pull/53) commit 09a7d185, contributor @[ftostevin-ont](https://github.com/ftostevin-ont)), accuracy improved on complex Indels. 3. Streamized pileup and full-alignment training workflows. Reduce diskspace demand in model training ([#55](https://github.com/HKU-BAL/Clair3/pull/55) commit 09a7d185, contributor @[ftostevin-ont](https://github.com/ftostevin-ont)). 4.  Added `mini_epochs` option in Train.py, performance  slightly improved in training a model for ONT Q20 data using mini-epochs([#60](https://github.com/HKU-BAL/Clair3/pull/60), contributor @[ftostevin-ont](https://github.com/ftostevin-ont)). 5. Massively reduced disk space demand when outputting GVCF. Now compressing GVCF intermediate files with lz4, five times smaller with little speed penalty. 6. Added `--remove_intermediate_dir`to remove intermediate files as soon as no longer needed ([#48](https://github.com/HKU-BAL/Clair3/issues/48)). 7. Renamed ONT pre-trained models with [Medaka](https://github.com/nanoporetech/medaka/blob/master/medaka/options.py#L22)'s naming convention. 8. Fixed training data spilling over to validation data ([#57](https://github.com/HKU-BAL/Clair3/issues/57)).

*ONT-provided Models (Sep 23)*: ONT also provides Clair3 models for specific chemistries and basecallers through [Rerio](https://github.com/nanoporetech/rerio).

*v0.1-r6 (Sep 4)* : 1. Reduced memory footprint at the `SortVcf` stage([#45](https://github.com/HKU-BAL/Clair3/issues/45)). 2. Reduced `ulimit -n` (number of files simultaneously opened) requirement ([#45](https://github.com/HKU-BAL/Clair3/issues/45), [#47](https://github.com/HKU-BAL/Clair3/issues/47)). 3. Added Clair3-Illumina package in bioconda([#42](https://github.com/HKU-BAL/Clair3/issues/42)).

*v0.1-r5 (July 19)* : 1. Modified data generator in model training to avoid memory exhaustion and unexpected segmentation fault by Tensorflow (contributor @[ftostevin-ont](https://github.com/ftostevin-ont) ). 2. Simplified dockerfile workflow to reuse container caching (contributor @[amblina](https://github.com/amblina)). 3. Fixed ALT output for reference calls (contributor @[wdecoster](https://github.com/wdecoster)). 4. Fixed a bug in multi-allelic AF computation (AF of [ACGT]Del variants was wrong before r5). 5. Added AD tag to the GVCF output. 6. Added the `--call_snp_only` option to only call SNP only ([#40](https://github.com/HKU-BAL/Clair3/issues/40)). 7. Added pileup and full-alignment output validity check to avoid workflow crashing ([#32](https://github.com/HKU-BAL/Clair3/issues/32), [#38](https://github.com/HKU-BAL/Clair3/issues/38)).

*v0.1-r4 (June 28)* : 1. Install via [bioconda](https://github.com/HKU-BAL/Clair3#option-3--bioconda). 2. Added an ONT Guppy2 model to the images (`ont_guppy2`). Click [here](https://github.com/HKU-BAL/Clair3/blob/main/docs/guppy2.md) for more benchmarking results. **The results show you have to use the Guppy2 model for Guppy2 or earlier data**. 3. Added [google colab notebooks](https://github.com/HKU-BAL/Clair3/blob/main/colab) for quick demo. 4. Fixed a bug when there are too few variant candidates ([#28](https://github.com/HKU-BAL/Clair3/issues/28)).

*v0.1-r3 (June 9)* : 1. Added `ulimit -u` (max user processes) check (lowers the `THREADS` if the resource is insufficient) and automatic retries on failed jobs ([#20](https://github.com/HKU-BAL/Clair3/issues/20), [#23](https://github.com/HKU-BAL/Clair3/issues/23), [#24](https://github.com/HKU-BAL/Clair3/issues/24)). 2. Added an ONT Guppy5 model to the images (`ont_guppy5`). Click [here](docs/guppy5.md) for more benchmarks on the Guppy5 model and data.

*v0.1-r2 (May 23)* : 1. Fixed BED file out of range error ([#12](https://github.com/HKU-BAL/Clair3/issues/12)). 2. Added support for both `.bam.bai` and `.bai` BAM index filename ([#10](https://github.com/HKU-BAL/Clair3/issues/10)). 3. Added some boundary checks on inputs. 4. Added version checks on required packages and utilities. 5. Increased pipeline robusity.

*v0.1-r1 (May 18)* : 1. Support relative path in Conda, but Docker and Singularity still require absolute path ([#5](https://github.com/HKU-BAL/Clair3/issues/5)). 2. Fix `taskset` CPU-core visibility and provide a Singularity image ([#6](https://github.com/HKU-BAL/Clair3/issues/6)).

*v0.1 (May 17, 2021)*: Initial release.

---

## Pre-trained Models

### HKU-provided Models

Download models from [here](http://www.bio8.cs.hku.hk/clair3/clair3_models/) or click on the links below.

In a docker installation, models are in `/opt/models/`. In a bioconda installation, models are in `{CONDA_PREFIX}/bin/models/`.

|           Model name           |  Platform   |                       Training samples                       | Included in the bioconda package | Included in the docker image |   Date   |  Basecaller  | File                                |                             Link                             |
| :----------------------------: | :---------: | :----------------------------------------------------------: | -------------------------------- | :--------------------------: | :------: | :----------: | ----------------------------------- | :----------------------------------------------------------: |
|      r941_prom_sup_g5014       |     ONT     |                    HG002,4,5 (Guppy5_sup)                    | Yes                              |             Yes              | 20220112 |  Guppy5 sup  | r941_prom_sup_g5014.tar.gz          | [Download](http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_sup_g5014.tar.gz) |
|    r941_prom_hac_g360+g422     |     ONT     |                         HG001,2,4,5                          | Yes                              |             Yes              | 20210517 | Guppy3,4 hac | r941_prom_hac_g360+g422.tar.gz      | [Download](http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_hac_g360+g422.tar.gz) |
|  r941_prom_hac_g360+g422_1235  |     ONT     |                         HG001,2,3,5                          |                                  |                              | 20210517 | Guppy3,4 hac | r941_prom_hac_g360+g422_1235.tar.gz | [Download](http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_hac_g360+g422_1235.tar.gz) |
|       r941_prom_hac_g238       |     ONT     |                         HG001,2,3,4                          |                                  |             Yes              | 20210627 |    Guppy2    | r941_prom_hac_g238.tar.gz           | [Download](http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_hac_g238.tar.gz) |
| ~~r941_prom_sup_g506~~ |     ONT     | Base model: HG001,2,4,5 (Guppy3,4) <br>Fine-tuning data: HG002 (Guppy5_sup) |                                  |                              | 20210609 |  Guppy5 sup  | r941_prom_sup_g506.tar.gz           | [Download](http://www.bio8.cs.hku.hk/clair3/clair3_models/r941_prom_sup_g506.tar.gz) |
|              hifi              | PacBio HiFi |                         HG001,2,4,5                          | Yes                              |             Yes              | 20210517 |      NA      | hifi.tar.gz                         | [Download](http://www.bio8.cs.hku.hk/clair3/clair3_models/hifi.tar.gz) |
|              ilmn              |  Illumina   |                         HG001,2,4,5                          | Yes                              |             Yes              | 20210517 |      NA      | ilmn.tar.gz                         | [Download](http://www.bio8.cs.hku.hk/clair3/clair3_models/ilmn.tar.gz) |

### ONT-provided Models

ONT provides models for some latest or specific chemistries and basecallers through [Rerio](https://github.com/nanoporetech/rerio). These models are tested and supported by the ONT developers. Avaiable model in Rerio including:

|       Config       | Chemistry  | Guppy basecaller |
| :----------------: | :-------:  | :--------------: |
| r104_e81_sup_g5015 | R10.4 E8.1 | v5.0.15 SUP      |
| r104_e81_hac_g5015 | R10.4 E8.1 | v5.0.15 HAC      |

----

## What's New in Clair3

* **New Architecture.** Clair3 integrates both pileup  (summarized alignment statistics) model and full-alignment model for variant calling. While a pileup model determines the result of a majority of variant candidates, candidates with uncertain results are further processed with a more computational-intensive haplotype-resolved full-alignment model.  
* **Improved Performance.** Using HG003 85-fold coverage ONT data from PrecisionFDA for benchmarking, Clair3 achieved 99.69% SNP F1-score and 80.58% Indel F1-score. Compare to Clair, Clair3 reduced SNP errors by **~78%**,  and Indel errors by **~48%**.  
* **High Efficiency.** Using 36 CPU cores,
  * Clair3 takes ~8 hours to process 50-fold WGS ONT data (~4x faster than PEPPER (r0.4) and ~14x faster than Medaka (v1.3.2)). Memory consumption of Clair3 is capped at 1 GB per CPU thread,  which is roughly five times lower than Clair. 
  * Clair3 takes ~2 hours to process 35-fold WGS PacBio HiFi data (13x faster than DeepVariant (v1.1.0)).
* **Using data from newer basecallers.**  Clair3 models were trained using data from Guppy version 3.6.0 and 4.2.2, please check [Training Data](docs/training_data.md) for details and links.  
* **GVCF Support.**  Clair3 can output GVCF using the ```--gvcf``` option, enabling downstream joint-sample genotyping and cohort merging. 

----

## Quick Demo

*   Oxford Nanopore (ONT) data, see [ONT Quick Demo](docs/quick_demo/ont_quick_demo.md).
*   PacBio HiFi data, see [PaBio HiFi Quick Demo](docs/quick_demo/pacbio_hifi_quick_demo.md).
*   Illumina NGS data, see [Illumina Quick Demo](docs/quick_demo/illumina_quick_demo.md).

**Run Clair3 ONT quick demo**: 

- **(Option 1) using Google Colab notebook:**

   [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/HKU-BAL/Clair3/blob/main/colab/clair3_ont_quick_demo.ipynb)

- **(Option 2) using pre-built docker image:**

```bash
cd ${HOME}
wget "http://www.bio8.cs.hku.hk/clair3/demo/clair3_ont_quick_demo.sh"
chmod +x clair3_ont_quick_demo.sh
./clair3_ont_quick_demo.sh
```

Check the results using `less ${HOME}/clair3_ont_quickDemo/output/merge_output.vcf.gz`

----

## Installation

### Option 1.  Docker pre-built image

A pre-built docker image is available [here](https://hub.docker.com/r/hkubal/clair3). With it you can run Clair3 using a single command.

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR`. 

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
MODEL_NAME="[YOUR_MODEL_NAME]"         # e.g. r941_prom_hac_g360+g422

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/${MODEL_NAME}" \
  --output=${OUTPUT_DIR}               ## absolute output path prefix 
```

Check [Usage](#Usage) for more options.

### Option 2. Singularity

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR`. 

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
MODEL_NAME="[YOUR_MODEL_NAME]"         # e.g. r941_prom_hac_g360+g422

conda config --add channels defaults
conda create -n singularity-env -c conda-forge singularity -y
conda activate singularity-env

# singularity pull docker pre-built image
singularity pull docker://hkubal/clair3:latest

# run clair3 like this afterward
singularity exec clair3_latest.sif \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/${MODEL_NAME}" \
  --output=${OUTPUT_DIR}               ## absolute output path prefix
```

### Option 3.  Bioconda

*For using Clair3 with Illumina data, install [clair3-illumina](https://anaconda.org/bioconda/clair3-illumina) package in bioconda channel instead.*

```bash
# make sure channels are added in conda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# create conda environment named "clair3"
# replace clair3 by clair3-illumina for using illumina data
conda create -n clair3 -c bioconda clair3 python=3.6.10 -y
conda activate clair3

# run clair3 like this afterward
MODEL_NAME="[YOUR_MODEL_NAME]"         # e.g. r941_prom_hac_g360+g422

run_clair3.sh \
  --bam_fn=input.bam \                 ## change your bam file name here
  --ref_fn=ref.fa \                    ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="${CONDA_PREFIX}/bin/models/${MODEL_NAME}" \ 
  --output=${OUTPUT_DIR}               ## output path prefix 
```

Check [Usage](#Usage) for more options. [Pre-trained models](#pre-trained-models) are already included in the bioconda package.

### Option 4. Build an anaconda virtual environment

**Anaconda install**:

Please install anaconda using the official [guide](https://docs.anaconda.com/anaconda/install) or using the commands below:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./Miniconda3-latest-Linux-x86_64.sh 
./Miniconda3-latest-Linux-x86_64.sh
```

**Install Clair3 using anaconda step by step:**

*For using Clair3 on Illumina data, additional installation steps after the following steps are mandatory. Please follow this [guide](https://github.com/HKU-BAL/Clair3/blob/main/docs/quick_demo/illumina_quick_demo.md#step-2-install-boost-graph-library-for-illumina-realignment-process) for the additional steps.*

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. ./input
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. ./output
THREADS="[MAXIMUM_THREADS]"            # e.g. 8

# create and activate an environment named clair3
conda create -n clair3 python=3.6.10 -y
source activate clair3

# install pypy and packages in the environemnt
conda install -c conda-forge pypy3.6 -y
pypy3 -m ensurepip
pypy3 -m pip install mpmath==1.2.1

# install python packages in environment
pip3 install tensorflow==2.2.0
pip3 install tensorflow-addons==0.11.2 tables==3.6.1
conda install -c anaconda pigz==2.4 cffi==1.14.4 -y
conda install -c conda-forge parallel=20191122 zstd=1.4.4 -y
conda install -c conda-forge -c bioconda samtools=1.10 -y
conda install -c conda-forge -c bioconda whatshap=1.0 -y
conda install -c conda-forge xz zlib bzip2 automake curl -y
    
# clone Clair3
git clone https://github.com/HKU-BAL/Clair3.git
cd Clair3

# compile samtools, longphase and cffi library for c implement
# after building, longphase binary is in `Clair3` folder
source activate clair3 && make PREFIX=${CONDA_PREFIX}

# download pre-trained models
mkdir models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz 
tar -zxvf clair3_models.tar.gz -C ./models

# run clair3
MODEL_NAME="[YOUR_MODEL_NAME]"         # e.g. r941_prom_hac_g360+g422
./run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path=`pwd`"/models/${MODEL_NAME}" \
  --output=${OUTPUT_DIR}               ## output path prefix
```



### Option 5. Docker Dockerfile

This is the same as option 1 except that you are building a docker image yourself. Please refer to option 1 for usage. 

```bash
# clone Clair3
git clone https://github.com/hku-bal/Clair3.git
cd Clair3

# build a docker image named hkubal/clair3:latest
# might require docker authentication to build docker image 
docker build -f ./Dockerfile -t hkubal/clair3:latest .

# run clair3 docker image like option 1
docker run -it hkubal/clair3:latest /opt/bin/run_clair3.sh --help
```

----

## Usage

### General Usage

**Caution**:  Use `=value` for optional parameters, e.g. `--bed_fn=fn.bed` instead of `--bed_fn fn.bed`.

```bash
./run_clair3.sh \
  --bam_fn=${BAM} \
  --ref_fn=${REF} \
  --threads=${THREADS} \  		     
  --platform="ont" \               ## options: {ont,hifi,ilmn}
  --model_path=${MODEL_PREFIX} \   ## absolute model path prefix
  --output=${OUTPUT_DIR}           ## absolute output path prefix
## pileup output file: ${OUTPUT_DIR}/pileup.vcf.gz
## full-alignment output file: ${OUTPUT_DIR}/full_alignment.vcf.gz
## Clair3 final output file: ${OUTPUT_DIR}/merge_output.vcf.gz
```

### Options

**Required parameters:**

```bash
  -b, --bam_fn=FILE             BAM file input. The input file must be samtools indexed.
  -f, --ref_fn=FILE             FASTA reference file input. The input file must be samtools indexed.
  -m, --model_path=STR          The folder path containing a Clair3 model (requiring six files in the folder, including pileup.data-00000-of-00002, pileup.data-00001-of-00002 pileup.index, full_alignment.data-00000-of-00002, full_alignment.data-00001-of-00002  and full_alignment.index).
  -t, --threads=INT             Max threads to be used. The full genome will be divided into small chunks for parallel processing. Each chunk will use 4 threads. The chunks being processed simultaneously is ceil($threads/4)*3. 3 is the overloading factor.
  -p, --platform=STR            Select the sequencing platform of the input. Possible options: {ont,hifi,ilmn}.
  -o, --output=PATH             VCF/GVCF output directory.
```

**Other parameters:**

 **Caution**:  Use `=value` for optional parameters, e.g., `--bed_fn=fn.bed` instead of `--bed_fn fn.bed`

```bash
      --bed_fn=FILE             Call variants only in the provided bed regions.
      --vcf_fn=FILE             Candidate sites VCF file input, variants will only be called at the sites in the VCF file if provided.
      --ctg_name=STR            The name of the sequence to be processed.
      --sample_name=STR         Define the sample name to be shown in the VCF file.
      --qual=INT                If set, variants with >$qual will be marked PASS, or LowQual otherwise.
      --samtools=STR            Path of samtools, samtools version >= 1.10 is required.
      --python=STR              Path of python, python3 >= 3.6 is required.
      --pypy=STR                Path of pypy3, pypy3 >= 3.6 is required.
      --parallel=STR            Path of parallel, parallel >= 20191122 is required.
      --whatshap=STR            Path of whatshap, whatshap >= 1.0 is required.
      --longphase=STR           Path of longphase, longphase >= 1.0 is required.
      --chunk_size=INT          The size of each chuck for parallel processing, default: 5Mbp.
      --pileup_only             Use the pileup model only when calling, default: disable.
      --print_ref_calls         Show reference calls (0/0) in vcf file, default: disable.
      --include_all_ctgs        Call variants on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}, default: disable.
      --gvcf                    Enable GVCF output, default: disable.
      --enable_phasing          Output phased variants using whatshap, default: disable.
      --longphase_for_phasing   Use longphase for phasing, default: enable.
      --disable_c_impl          Disable C implement with cffi for pileup and full-alignment create tensor, default: enable.
      --remove_intermediate_dir Remove intermediate directory, including intermediate phased BAM, pileup and full-alignment results. default: disable.
      --snp_min_af=FLOAT        Minimum SNP AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.08,hifi:0.08,ilmn:0.08.
      --indel_min_af=FLOAT      Minimum INDEL AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.15,hifi:0.08,ilmn:0.08.
      --var_pct_full=FLOAT      EXPERIMENTAL: Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for full-alignment mode calling, default: 0.3.
      --ref_pct_full=FLOAT      EXPERIMENTAL: Specify an expected percentage of low quality 0/0 variants called in the pileup mode for full-alignment mode calling, default: 0.3 for ilmn and hifi, 0.1 for ont.
      --var_pct_phasing=FLOAT   EXPERIMENTAL: Specify an expected percentage of high quality 0/1 variants used in WhatsHap phasing, default: 0.8 for ont guppy5 and 0.7 for other platforms.
      --pileup_model_prefix=STR EXPERIMENTAL: Model prefix in pileup calling, including $prefix.data-00000-of-00002, $prefix.data-00001-of-00002 $prefix.index. default: pileup.
      --fa_model_prefix=STR     EXPERIMENTAL: Model prefix in full-alignment calling, including $prefix.data-00000-of-00002, $prefix.data-00001-of-00002 $prefix.index, default: full_alignment.
      --min_mq=INT              EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered, default: 5.
      --min_coverage=INT        EXPERIMENTAL: Minimum coverage required to call a variant, default: 2.
      --min_contig_size=INT     EXPERIMENTAL: If set, contigs with contig size<$min_contig_size are filtered, default: 0.
      --fast_mode               EXPERIMENTAL: Skip variant candidates with AF <= 0.15, default: disable.
      --haploid_precise         EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant, default: disable.
      --haploid_sensitive       EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant, default: disable.
      --no_phasing_for_fa       EXPERIMENTAL: Call variants without whatshap phasing in full alignment calling, default: disable.
      --call_snp_only           EXPERIMENTAL: Call candidates pass SNP minimum AF only, ignore Indel candidates, default: disable.
      --enable_long_indel       EXPERIMENTAL: Call long Indel variants(>50 bp), default: disable.
```

#### Call variants in a chromosome

```bash
CONTIGS_LIST="[YOUR_CONTIGS_LIST]"     # e.g "chr21" or "chr21,chr22"
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input  (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
MODEL_NAME="[YOUR_MODEL_NAME]"         # e.g. r941_prom_hac_g360+g422

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/${MODEL_NAME}" \
  --output=${OUTPUT_DIR} \             ## absolute output path prefix
  --ctg_name=${CONTIGS_LIST}
```

#### Call variants at known variant sites

```bash
KNOWN_VARIANTS_VCF="[YOUR_VCF_PATH]"   # e.g. /home/user1/known_variants.vcf.gz (absolute path needed)
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
MODEL_NAME="[YOUR_MODEL_NAME]"         # e.g. r941_prom_hac_g360+g422

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/${MODEL_NAME}" \
  --output=${OUTPUT_DIR} \             ## absolute output path prefix
  --vcf_fn=${KNOWN_VARIANTS_VCF}
```

#### Call variants at specific sites or bed regions

We highly recommended using BED file to define the regions of interest like:

```shell
# define 0-based "ctg start end" if at specific sites
CONTIGS="[YOUR_CONTIGS_NAME]"          # e.g. chr22
START_POS="[YOUR_START_POS]"           # e.g. 0
END_POS="[YOUR_END_POS]"               # e.g 10000
echo -e "${CONTIGS}\t${START_POS}\t${END_POS}" > /home/user1/tmp.bed ## change directory accordingly
```

Then run Clair3 like this:

```bash
BED_FILE_PATH="[YOUR_BED_FILE]"        # e.g. /home/user1/tmp.bed (absolute path needed)
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
MODEL_NAME="[YOUR_MODEL_NAME]"         # e.g. r941_prom_hac_g360+g422

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/${MODEL_NAME}" \
  --output=${OUTPUT_DIR} \             ## absolute output path prefix
  --bed_fn=${BED_FILE_PATH}
```

#### Call variants in non-diploid organisms (Haploid calling)

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
MODEL_NAME="[YOUR_MODEL_NAME]"         # e.g. r941_prom_hac_g360+g422

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference file name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/${MODEL_NAME}" \
  --output=${OUTPUT_DIR} \
  --no_phasing_for_fa \                ## disable phasing for full-alignment
  --include_all_ctgs \                 ## call variants on all contigs in the reference fasta
  --haploid_precise                    ## optional(enable --haploid_precise or --haploid_sensitive) for haploid calling
```

----

## Folder Structure and Submodule Descriptions

Submodules in __`clair3/`__ are for variant calling and model training. Submodules in __`preprocess`__ are for data preparation.

*For all the submodules listed below, you can use `-h` or `--help` for available options.*

`clair3/` | Note: submodules under this folder are pypy incompatible, please run using python
---: | ---
`CallVariants` | Call variants using a trained model and tensors of candidate variants.
`CallVarBam` | Call variants using a trained model and a BAM file.
`Train` | Training a model using the `RectifiedAdam` optimizer. We also use the `Lookahead` optimizer to adjust the `RectifiedAdam` parameters dynamically. The initial learning rate is `1e-3` with `0.1` learning rate warm-up. Input a binary containing tensors created by `Tensor2Bin`. 



`preprocess/` | Note: submodules under this folder is Pypy compatible unless specified.
---: | ---
`CheckEnvs`| Check the environment and  validity of the input variables, preprocess the BED input if necessary, `--chunk_size` sets the chuck size to be processed per parallel job. 
`CreateTensorPileup`| Generate variant candidate tensors in pileup format for training or calling. 
`CreateTensorFullAlignment`| Generate variant candidate tensors in phased full-alignment format for training or calling. 
`GetTruth`| Extract the variants from a truth VCF. Input: VCF; Reference FASTA if the VCF contains asterisks in ALT field.
`MergeVcf` | Merge pileup and full-alignment VCF/GVCF.
`RealignReads` | Reads local realignment for Illumina platform.
`SelectCandidates`| Select pileup candidates for full-alignment calling.
`SelectHetSnp` | Select heterozygous SNP candidates for whatshap phasing.
`SelectQual` | Select a quality cutoff using the pileup calling results. Variants below the cutoff are included in phasing and full-alignment calling. 
`SortVcf` | Sort VCF file. 
`SplitExtendBed` | Split BED file regions according to the contig names and extend bed region by 33bp by default for variant calling. 
`UnifyRepresentation` | Representation unification between candidate sites and true variants. 
`MergeBin` | Combine tensor binaries into a single file. 
`CreateTrainingTensor` | Create tensor binaries for pileup or full-alignment training. 
`Tensor2Bin` | Combine the variant and non-variant tensors and convert them to a binary, using `blosc:lz4hc` meta-compressor, the overall training memory is 10~15G (pypy incompatible). 

----

## Training Data

Clair3 trained both its pileup and full-alignment models using four GIAB samples (HG001, HG002, HG004 and HG005), excluded HG003. On ONT, we also trained a model using HG001, 2, 3, and 5, excluded HG004. All models were trained with chr20 excluded (including only chr1-19, 21, 22). 

|  Platform   |   Reference   |      Aligner      | Training samples |
| :---------: | :-----------: | :---------------: | :--------------: |
|     ONT     | GRCh38_no_alt |     minimap2      | HG001,2,(3\|4),5 |
| PacBio HiFi | GRCh38_no_alt |       pbmm2       |   HG001,2,4,5    |
|  Illumina   |    GRCh38     | BWA-MEM/NovoAlign |   HG001,2,4,5    |

Please find more details about the training data and links at [Training Data](docs/training_data.md).

----

## VCF/GVCF Output Formats

Clair3 supports both VCF and GVCF output formats. Clair3 uses VCF version 4.2 specifications. Specifically, Clair3 adds a `P` INFO tag to the results called using a pileup model, and a `F` INFO tag to the results called using a full-alignment model.

Clair3 outputs a GATK-compatible GVCF format that passes GATK's `ValidateVariants` module. Different from DeepVariant that uses `<*>` to represent any possible alternative allele, Clair3 uses `<NON_REF>`, the same as GATK.
