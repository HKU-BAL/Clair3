<div align="center">
  <a href="https://en.wiktionary.org/wiki/%E7%9C%BC" target="_blank">
    <img src="docs/images/clair3_logo.png" width = "110" height = "90" alt="Clair3">
  </a>
</div>

# Clair3 - Integrating pileup and full-alignment for high-performance long-read variant calling

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)  

Contact: Ruibang Luo  

Email: rbluo@cs.hku.hk  

---

## Introduction

Clair3 is a small variant caller for long-reads. Compare to PEPPER (r0.4), Clair3 (v0.1) shows a better SNP F1-score with ≤30-fold of ONT data (precisionFDA Truth Challenge V2), and a better Indel F1-score, while runs generally four times faster. Clair3 makes the best of both worlds of using pileup or full-alignment as input for deep-learning based long-read small variant calling. Clair3 is simple and modular for easy deployment and integration.

Clair3 is the 3<sup>rd</sup> generation of [Clair](https://github.com/HKU-BAL/Clair) (the 2<sup>nd</sup>) and [Clairvoyante](https://github.com/aquaskyline/Clairvoyante) (the 1<sup>st</sup>).

---

## Lastest Updates

We are actively fixing bugs and issues in Clair3 reported by users.

*v0.1-r2 (May 23)* : 1. Fixed BED file out of range error ([#12](https://github.com/HKU-BAL/Clair3/issues/12)). 2. Added support for both `.bam.bai` and `.bai` BAM index filename ([#10](https://github.com/HKU-BAL/Clair3/issues/10)). 3. Added some boundary checks on inputs. 4. Added version checks on required packages and utilities. 5. Increased pipeline robusity.

*v0.1-r1 (May 18)* : 1. Support relative path in Conda, but Docker and Singularity still require absolute path ([#5](https://github.com/HKU-BAL/Clair3/issues/5)). 2. Fix `taskset` CPU-core visibility and provide a Singularity image ([#6](https://github.com/HKU-BAL/Clair3/issues/6)).

*v0.1 (May 17)*: Initial release.

---

## We are working on ...

* A paper on detailed methods and benchmarks.
* A model trained with Guppy2 data. The available ONT models are tested and work well with Guppy3 and Guppy4 data, but perform even worse than Clair on Guppy2 data.
* A model trained with Guppy5 data. 

---

## Contents

* [Introduction](#introduction)
* [Latest Updates](#latest-updates)
* [What's New in Clair3](#whats-new-in-clair3)
* [Installation](#installation)
  + [Option 1. Docker pre-built image](#option-1--docker-pre-built-image)
  + [Option 2. Singularity](#option-2-singularity)
  + [Option 3. Build an anaconda virtual environment](#option-3-build-an-anaconda-virtual-environment)
  + [Option 4. Docker Dockerfile](#option-4-docker-dockerfile)
* [Quick Demo](#quick-demo)
* [Usage](#usage)
* [Folder Structure and Submodule Descriptions](#folder-structure-and-submodule-descriptions)
* [Training Data](#training-data)
* [VCF/GVCF Output Formats](#vcfgvcf-output-formats)
* [Pileup Model Training](docs/pileup_training.md)
* [Full-Alignment Model Training](docs/full_alignment_training.md)
* [Representation Unification](docs/representation_unification.md)
* [Visualization](docs)
  * [Model Input](docs/model_input_visualization.md)
  * [Representation Unification](docs/representation_unification_visualization.md)

---

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

**Run Clair3 ONT quick demo using pre-built docker image:**

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

A pre-built docker image is available [here](https://hub.docker.com/layers/hkubal/clair3/latest/images/sha256-769a241a9e1aab422d7309022ab14e8982d1e2af32c24ee7c16230c24b52cd74?context=explore). With it you can run Clair3 using a single command.

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR`. 

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
BIN_VERSION="v0.1-r2"

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:"${BIN_VERSION}" \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/ont" \     ## absolute model path prefix, change platform accordingly
  --output=${OUTPUT_DIR}               ## absolute output path prefix 
```

Check [Usage](#Usage) for more options.

### Option 2. Singularity

**Caution**: Absolute path is needed for both `INPUT_DIR` and `OUTPUT_DIR`. 

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
BIN_VERSION="v0.1-r2"

conda config --add channels defaults
conda create -n singularity-env -c conda-forge singularity -y
conda activate singularity-env

# singularity pull docker pre-built image
singularity pull docker://hkubal/clair3:v0.1-r2

# run clair3 like this afterward
singularity exec clair3_"${BIN_VERSION}".sif \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/ont" \     ## absolute model path prefix, change platform accordingly
  --output=${OUTPUT_DIR}               ## absolute output path prefix
```

### Option 3. Build an anaconda virtual environment

**Anaconda install**:

Please install anaconda using the official [guide](https://docs.anaconda.com/anaconda/install) or using the commands below:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./Miniconda3-latest-Linux-x86_64.sh 
./Miniconda3-latest-Linux-x86_64.sh
```

**Install Clair3 using anaconda step by step:**

*For using Clair3 on Illumina data, after the following steps, please also install the* [Boost Graph Library](https://www.boost.org/doc/libs/1_65_1/libs/graph/doc/index.html) *using this* [guidance](docs/quick_demo/illumina_quick_demo.md#step-2-install-boost-graph-library-for-illumina-realignment-process).

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
pypy3 -m pip install intervaltree==3.0.2
pypy3 -m pip install mpmath==1.2.1

# install python packages in environment
pip3 install tensorflow==2.2.0
pip3 install intervaltree==3.0.2  tensorflow-addons==0.11.2 tables==3.6.1
conda install -c anaconda pigz==2.4 -y
conda install -c conda-forge parallel=20191122 zstd=1.4.4 -y
conda install -c conda-forge -c bioconda samtools=1.10 -y
conda install -c conda-forge -c bioconda whatshap=1.0 -y

# clone Clair3
git clone https://github.com/HKU-BAL/Clair3.git
cd Clair3

# download pre-trained models
mkdir models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz 
tar -zxvf clair3_models.tar.gz -C ./models

# run clair3
./run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path=`pwd`"/models/ont" \    ## model path prefix, change platform accordingly
  --output=${OUTPUT_DIR}               ## output path prefix
```

### Option 4. Docker Dockerfile

This is the same as option 1 except that you are building a docker image yourself. Please refer to option 1 for usage. 

```bash
BIN_VERSION="v0.1-r2"

# clone Clair3
git clone https://github.com/hku-bal/Clair3.git
cd Clair3

# build a docker image named hkubal/clair3:v0.1-r2
# might require docker authentication to build docker image 
docker build -f ./Dockerfile -t hkubal/clair3:"${BIN_VERSION}" .

# run clair3 docker image like option 1
docker run -it hkubal/clair3:"${BIN_VERSION}" /opt/bin/run_clair3.sh --help
```

---

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
  -b, --bam_fn=FILE        BAM file input. The input file must be samtools indexed.
  -f, --ref_fn=FILE        FASTA reference file input. The input file must be samtools indexed.
  -m, --model_path=STR     The folder path containing a Clair3 model (requiring six files in the folder, including pileup.data-00000-of-00002, pileup.data-00001-of-00002 pileup.index, full_alignment.data-00000-of-00002, full_alignment.data-00001-of-00002  and full_alignment.index).
  -t, --threads=INT        Max threads to be used. The full genome will be divided into small chunks for parallel processing. Each chunk will use 4 threads. The chunks being processed simultaneously is ceil($threads/4)*3. 3 is the overloading factor.
  -p, --platform=STR       Select the sequencing platform of the input. Possible options: {ont,hifi,ilmn}.
  -o, --output=PATH        VCF/GVCF output directory.
```

**Other parameters:**

 **Caution**:  Use `=value` for optional parameters, e.g., `--bed_fn=fn.bed` instead of `--bed_fn fn.bed`

```bash
      --bed_fn=FILE        Call variants only in the provided bed regions.
      --vcf_fn=FILE        Candidate sites VCF file input, variants will only be called at the sites in the VCF file if provided.
      --ctg_name=STR       The name of the sequence to be processed.
      --sample_name=STR    Define the sample name to be shown in the VCF file.
      --qual=INT           If set, variants with >=$qual will be marked PASS, or LowQual otherwise.
      --samtools=STR       Path of samtools, samtools version >= 1.10 is required.
      --python=STR         Path of python, python3 >= 3.6 is required.
      --pypy=STR           Path of pypy3, pypy3 >= 3.6 is required.
      --parallel=STR       Path of parallel, parallel >= 20191122 is required.
      --whatshap=STR       Path of whatshap, whatshap >= 1.0 is required.
      --chunk_size=INT     The size of each chuck for parallel processing, default: 5Mbp.
      --pileup_only        Use the pileup model only when calling, default: disable.
      --print_ref_calls    Show reference calls (0/0) in vcf file, default: disable.
      --include_all_ctgs   Call variants on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}, default: disable.
      --gvcf               Enable GVCF output, default: disable.
      --snp_min_af=FLOAT   Minimum SNP AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.08,hifi:0.08,ilmn:0.08.
      --indel_min_af=FLOAT Minimum INDEL AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.15,hifi:0.08,ilmn:0.08.
      --var_pct_full=FLOAT EXPERIMENTAL: Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for full-alignment mode calling, default: 0.3.
      --ref_pct_full=FLOAT EXPERIMENTAL: Specify an expected percentage of low quality 0/0 variants called in the pileup mode for full-alignment mode calling, default: 0.3 for ilmn and hifi, 0.1 for ont.
      --fast_mode          EXPERIMENTAL: Skip variant candidates with AF <= 0.15, default: disable.
      --haploid_precise    EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant, default: disable.
      --haploid_sensitive  EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant, default: disable.
      --no_phasing_for_fa  EXPERIMENTAL: Call variants without whatshap phasing in full alignment calling, default: disable.
```

#### Call variants in a chromosome

```bash
CONTIGS_LIST="[YOUR_CONTIGS_LIST]"     # e.g "chr21" or "chr21,chr22"
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input  (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
BIN_VERSION="v0.1-r2"

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:"${BIN_VERSION}" \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/ont" \     ## absolute model path prefix, change platform accordingly
  --output=${OUTPUT_DIR} \             ## absolute output path prefix
  --ctg_name=${CONTIGS_LIST}
```

#### Call variants at known variant sites

```bash
KNOWN_VARIANTS_VCF="[YOUR_VCF_PATH]"   # e.g. /home/user1/known_variants.vcf.gz (absolute path needed)
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
BIN_VERSION="v0.1-r2"

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:"${BIN_VERSION}" \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/ont" \     ## absolute model path prefix, change platform accordingly
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
BIN_VERSION="v0.1-r2"

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:"${BIN_VERSION}" \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/ont" \     ## absolute model path prefix, change platform accordingly
  --output=${OUTPUT_DIR} \             ## absolute output path prefix
  --bed_fn=${BED_FILE_PATH}
```

#### Call variants in non-diploid organisms (Haploid calling)

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"            # e.g. 8
BIN_VERSION="v0.1-r2"

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:"${BIN_VERSION}" \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## change your reference name here
  --threads=${THREADS} \               ## maximum threads to be used
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="/opt/models/ont" \     ## absolute model path prefix, change platform accordingly
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

#### Pre-trained Models

Download models from [here](http://www.bio8.cs.hku.hk/clair3/clair3_models/) or click on the links below.

|      File       |  Platform   | Training samples | Included in the docker image | Release |   Date   |    Basecaller   |                             Link                             |
| :-------------: | :---------: | :--------------: | :----------------------------: | :-----: | :------: | :------: | :----------------------------------------------------------: |
|   ont.tar.gz    |     ONT     |   HG001,2,4,5    |              Yes               |    1    | 20210517 | Guppy3,4 | [Download](http://www.bio8.cs.hku.hk/clair3/clair3_models/ont.tar.gz) |
| ont_1235.tar.gz |     ONT     |   HG001,2,3,5    |                                |    1    | 20210517 | Guppy3,4 | [Download](http://www.bio8.cs.hku.hk/clair3/clair3_models/ont_1235.tar.gz) |
|   hifi.tar.gz   | PacBio HiFi |   HG001,2,4,5    |              Yes               |    1    | 20210517 | NA | [Download](http://www.bio8.cs.hku.hk/clair3/clair3_models/hifi.tar.gz) |
|   ilmn.tar.gz   |  Illumina   |   HG001,2,4,5    |              Yes               |    1    | 20210517 | NA | [Download](http://www.bio8.cs.hku.hk/clair3/clair3_models/ilmn.tar.gz) |

----

## VCF/GVCF Output Formats

Clair3 supports both VCF and GVCF output formats. Clair3 uses VCF version 4.2 specifications. Specifically, Clair3 adds a `P` INFO tag to the results called using a pileup model, and a `F` INFO tag to the results called using a full-alignment model.

Clair3 outputs a GATK-compatible GVCF format that passes GATK's `ValidateVariants` module. Different from DeepVariant that uses `<*>` to represent any possible alternative allele, Clair3 uses `<NON_REF>`, the same as GATK.

