
# Clair3 - Integrating pileup and full-alignment for high-performance long-read variant calling

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) \
Contact: Ruibang Luo \
Email: rbluo@cs.hku.hk

---

## Introduction



This is the formal release of Clair3, the successor of [Clair](https://github.com/HKU-BAL/Clair) and [Clairvoyante](https://github.com/aquaskyline/Clairvoyante).

---

## Contents

* [Introduction](#introduction)
* [What's New in Clair3](#whats-new-in-clair3)
* [Installation](#installation)
  + [Option 1. Docker pre-built image (recommended)](#option-1--docker-pre-built-image-recommended)
  + [Option 2. Docker Dockerfile](#option-2-docker-dockerfile)
  + [Option 3. Build an anaconda virtual environment](#option-3-build-an-anaconda-virtual-environment)
* [Quick Demo](#quick-demo)
* [Usage](#usage)
* [Folder Structure and Submodule Descriptions](#folder-structure-and-submodule-descriptions)
* [VCF Output Format](#vcf-output-format)
* [Pileup Model Training](docs/pileup_training.md)
* [Full-Alignment Model Training](docs/full_alignment_training.md)
* [Representation Unification](docs/representation_unification.md)
* [Visualization](docs)
  * [Model Input](docs/model_input_visualization.md)
  * [Representation Unification](docs/representation_unification_visualization.md)

---

## What's New in Clair3

* **Accuracy Improvements.** Clair3 outperforms Clair, reducing SNP errors by **~78%**,  and Indel errors by **~48%** in HG003 ONT case study. SNP F1-score reaches 99.69% and Indel F1-score reaches 80.58% in ONT HG003 ~85-fold coverage dataset.  
* **New Architecture.** Clair3 is an integration of pileup model and full-alignment model. Pileup model detects all candidate variants using summarized pileup input. Full-alignment model adopts more complete read-level representations with phased haplotype information to further decide the variant type of low-quality pileup candidates. 
* **High Efficiency.** 
  * Clair3 takes about ~8 hours for ONT ~50-fold coverage whole-genome-sequencing data using 36 CPUs, which is ~4.5x faster than PEPPER and ~14x faster than Medaka. Computational resource consumption using Clair3 is capped at 1 GB per CPU thread,  which is ~6 times lower than Clair and PEPPER. 
  * Clair3 takes about ~2 hours For PacBio HiFi ~35-fold coverage whole-genome-sequencing data using 36 CPUs, which is 13x faster than DeepVariant workflow.
* **New Base Caller Support.**  We support datasets base called using Guppy version 3.6.0~4.2.2  for ONT platform, check the [training data](docs/training_data.md) for specific datasets' link. We did not suggest using datasets base called using **Guppy version <= 3.6.0** for calling as we have discarded those datasets in model training.  

## Quick Demo

Clair3 supports germline variant calling in three sequencing platforms:

*   Oxford Nanopore (ONT) long-read data, see [ONT Quick Demo](docs/quick_demo/ont_quick_demo.md).
*   PacBio HiFi data, see [PaBio HiFi Quick Demo](docs/quick_demo/pacbio_hifi_quick_demo.md).
*   Illumina NGS data, see [Illumina Quick Demo](docs/quick_demo/illumina_quick_demo.md).

**Run Clair3 using pre-built docker image:**

```bash
cd ${HOME}
wget "http://www.bio8.cs.hku.hk/clair3/demo/clair3_ont_quick_demo.sh"
chmod +x clair3_ont_quick_demo.sh
./clair3_ont_quick_demo.sh
```

Check the results using `less ${HOME}/clair3_ont_quickDemo/output/merge_output.vcf.gz`

## Installation

### Option 1.  Docker pre-built image (recommended)

A pre-built docker image can be found [here](https://hub.docker.com/hku-bal/clair3). Then you can run Clair3 using one command:

```bash
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. input/
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. output/
THREADS="[MAXIMUM_THREADS]"            # e.g. 36
BIN_VERSION="v0.1"

docker run \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hku-bal/clair3:"${BIN_VERSION}" \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## Change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \          ## Change your reference name here
  --threads=${THREADS} \               ## Maximum threads to be used
  --platform="ont" \                   ## Options: {ont,hifi,ilmn}
  --model_path="/opt/models/ont" \     ## Options: {ont,hifi,ilmn}
  --output=${OUTPUT_DIR}
```

Check [Usage](#Usage)  for more options.

### Option 2. Docker Dockerfile

** **

```bash
# clone Clair3
git clone --depth 1 https://github.com/HKU-BAL/Clair3.git
cd Clair3

# build a docker image named clair3_docker
# might require docker authentication to build by docker
docker build -f ./Dockerfile -t clair3_docker .

# run clair docker image like this afterward
docker run -it clair3_docker /opt/bin/run_clair3.sh --help
```

### Option 3. Build an anaconda virtual environment

**Anaconda install**:

Please install anaconda using the installation guide at [Install](https://docs.anaconda.com/anaconda/install) or using the command below:

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod +x ./Miniconda3-latest-Linux-x86_64.sh 
./Miniconda3-latest-Linux-x86_64.sh
```

**Install Clair3 using anaconda step by step:**

```bash
# create and activate the environment named clair3
conda create -n clair3 python=3.6.10 -y
source activate clair3

# install pypy and packages on clair3 environemnt
conda install -c conda-forge pypy3.6 -y
pypy3 -m ensurepip
pypy3 -m pip install intervaltree==3.0.2
pypy3 -m pip install mpmath==1.2.1

# install python packages on clair3 environment
pip3 install tensorflow==2.2.0
pip3 install intervaltree==3.0.2  tensorflow-addons==0.11.2 tables==3.6.1
conda install -c anaconda pigz==2.4 -y
conda install -c conda-forge parallel=20191122 zstd=1.4.4 -y
conda install -c conda-forge -c bioconda samtools=1.10 -y
conda install -c conda-forge -c bioconda whatshap=1.0 -y

# clone Clair3
git clone --depth 1 https://github.com/HKU-BAL/Clair3.git
cd Clair3

# download pre-trained model
mkdir models
wget http://www.bio8.cs.hku.hk/clair3/clair3_models/clair3_models.tar.gz 
tar -zxvf clair3_models.tar.gz -C ./models

# run clair3 like this afterward
./run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## Change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## Change your reference name here
  --threads=${THREADS} \               ## Maximum threads to be used
  --platform="ont" \                   ## Options: {ont,hifi,ilmn}
  --model_path=`pwd`"/models/ont" \    ## Absolute model path prefix, change platform accordingly
  --output=${OUTPUT_DIR}
```

If using Clair3 Illumina model for variant calling, follows the [guidance](docs/quick_demo/illumina_quick_demo.md#step-2-install-boost-graph-library-for-illumina-realignment-process) to install [Boost Graph Library](https://www.boost.org/doc/libs/1_65_1/libs/graph/doc/index.html) additionally.

## Usage

### General Usage

```bash
# optional parameters should use "="
./run_clair3.sh \
  --bam_fn=${BAM} \
  --ref_fn=${REF} \
  --threads=${THREADS} \  		     
  --platform='ont' \               ## Options: {ont,hifi,ilmn}
  --model_path=${MODEL_PREFIX} \   ## Options: {ont,hifi,ilmn}
  --output=${OUTPUT_DIR}
  
## Pileup output file: ${OUTPUT_DIR}/pileup.vcf.gz
## Full alignment output file: ${OUTPUT_DIR}/full_alignment.vcf
## Final Clair3 output file: ${OUTPUT_DIR}/merge_output.vcf.gz
```

### Options

**Required parameters:**

```bash
  -b, --bam_fn FILE        BAM file input. The input file must be samtools indexed.
  -f, --ref_fn FILE        FASTA reference file input. The input file must be samtools indexed.
  -m, --model_path STR     The folder path containing a Clair3 model (requiring six files in the folder, including pileup.data-00000-of-00001, pileup.index, pileup.meta, full_alignment.data-00000-of-00001, full_alignment.index, and full_alignment.meta).
  -t, --threads INT        Max threads to be used. The full genome will be divided into small chunks for parallel processing. Each chunk will use 4 threads. The chunks being processed simultaneously is ceil($threads/4)*3. 3 is the overloading factor.
  -p, --platform STR       Select the sequencing platform of the input. Possible options: {ont,hifi,ilmn}.
  -o, --output PATH        VCF/GVCF output directory.
```

**Optional parameters:**

```bash
      --bed_fn FILE        Call variants only in the provided bed regions.
      --vcf_fn FILE        Candidate sites VCF file input, variants will only be called at the sites in the VCF file if provided.
      --ctg_name STR       The name of the sequence to be processed.
      --sample_name STR    Define the sample name to be shown in the VCF file.
      --qual INT           If set, variants with >=$qual will be marked PASS, or LowQual otherwise.
      --samtools STR       Path of samtools, samtools version >= 1.10 is required.
      --python STR         Path of python, python3 >= 3.6 is required.
      --pypy STR           Path of pypy3, pypy3 >= 3.6 is required.
      --parallel STR       Path of parallel, parallel >= 20191122 is required.
      --whatshap STR       Path of whatshap, whatshap >= 1.0 is required.
      --chunk_size INT     The size of each chuck for parallel processing, default: 5Mbp.
      --pileup_only        Use the pileup model only when calling, default: disable.
      --print_ref_calls    Show reference calls (0/0) in vcf file, default: disable.
      --include_all_ctgs   Call variants on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}, default: disable.
      --gvcf               Enable GVCF output, default: disable.
      --snp_min_af FLOAT   Minimum SNP AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.08,hifi:0.08,ilmn:0.08.
      --indel_min_af FLOAT Minimum INDEL AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.15,hifi:0.08,ilmn:0.08.
      --var_pct_full FLOAT EXPERIMENTAL: Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for full-alignment mode calling, default: 0.3.
      --ref_pct_full FLOAT EXPERIMENTAL: Specify an expected percentage of low quality 0/0 variants called in the pileup mode for full-alignment mode calling, default: 0.3 for ilmn and hifi, 0.1 for ont.
      --fast_mode          EXPERIMENTAL: Skip variant candidates with AF <= 0.15, default: disable.
      --haploid_precise    EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant, default: disable.
      --haploid_sensitive  EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant, default: disable.
      --no_phasing_for_fa  EXPERIMENTAL: Call variants without whatshap phasing in full alignment calling, default: disable.
```

#### Call variants in a chromosome

```bash
CONTIGS_LIST="[YOUR_CONTIGS_LIST]"     # e.g "chr21" or "chr21,chr22"
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. input/
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. output/
THREADS="[MAXIMUM_THREADS]"            # e.g. 36
BIN_VERSION="v0.1"

docker run \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hku-bal/clair3:"${BIN_VERSION}" \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## Change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## Change your reference name here
  --threads=${THREADS} \               ## Maximum threads to be used
  --platform="ont" \                   ## Options: {ont,hifi,ilmn}
  --model_path="/opt/models/ont" \     ## Options: {ont,hifi,ilmn}
  --output=${OUTPUT_DIR} \
  --ctg_name=${CONTIGS_LIST}
```

#### Call variants at known variant sites

```bash
KNOWN_VARIANTS_VCF="[YOUR_VCF_PATH]"   # e.g. known_variants.vcf
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. input/
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. output/
THREADS="[MAXIMUM_THREADS]"            # e.g. 36
BIN_VERSION="v0.1"

docker run \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hku-bal/clair3:"${BIN_VERSION}" \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## Change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## Change your reference name here
  --threads=${THREADS} \               ## Maximum threads to be used
  --platform="ont" \                   ## Options: {ont,hifi,ilmn}
  --model_path="/opt/models/ont" \     ## Options: {ont,hifi,ilmn}
  --output=${OUTPUT_DIR} \
  --vcf_fn=${KNOWN_VARIANTS_VCF}
```

#### Call variants at specific sites or bed regions

We highly recommended using bed format file to define single or multiple start like:

```shell
# define 0-based "ctg start end" if at specific sites
CONTIGS="[YOUR_CONTIGS_NAME]"          # e.g. chr22
START_POS="[YOUR_START_POS]"           # e.g. 0
END_POS="[YOUR_END_POS]"               # e.g 10000
echo -e "${CONTIGS}\t${START_POS}\t${END_POS}" > tmp.bed
```

Then run Clair3 like this:

```bash
BED_FILE_PATH="[YOUR_BED_FILE]"		   # e.g. tmp.bed
INPUT_DIR="[YOUR_INPUT_FOLDER]"        # e.g. input/
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"      # e.g. output/
THREADS="[MAXIMUM_THREADS]"            # e.g. 36
BIN_VERSION="v0.1"

docker run \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hku-bal/clair3:"${BIN_VERSION}" \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \    ## Change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \       ## Change your reference name here
  --threads=${THREADS} \               ## Maximum threads to be used
  --platform="ont" \                   ## Options: {ont,hifi,ilmn}
  --model_path="/opt/models/ont" \     ## Options: {ont,hifi,ilmn}
  --output=${OUTPUT_DIR} \
  --bed_fn=${BED_FILE_PATH}
```

## Folder Structure and Submodule Descriptions

Submodules in __`clair3/`__ are for variant calling and model training. Submodules in __`preprocess`__ are for data preparation.

*For all the submodules listed below, you can use the `-h` or `--help` option for available options.*

`clair3/` | Note: submodules under this folder are pypy incompatible, please run using python
---: | ---
`CallVariants` | Call variants using a trained model and tensors of candidate variants.
`CallVarBam` | Call variants using a trained model and a BAM file.
`Train` |  Training a model using warm-up learning rate startup with `RectifiedAdam` optimizer. We also use `Lookahead` optimizer to adjust `RectifiedAdam` parameters adaptively. The initial learning rate is `1e-3` with 0.1 learning rate warm-up. Input a binary tensors are created by `Tensor2Bin`. 

`preprocess/` | Note: submodules under this folder is Pypy compatible unless specified.
---: | ---
`CheckEnvs`| Check the environment and the validity of the input variables, preprocess the BED input if necessary, `--chunk_size` set the genome chuck size per job.
`CreateTensorPileup`| Generate variant candidate tensors using pileup for training or calling.
`CreateTensorFullAlignment`| Generate variant candidate tensors using phased full-alignment for training or calling.
`GetTruth`| Extract the variants from a truth VCF. Input: VCF; Reference FASTA if the VCF contains asterisks in ALT field.
`MergeVcf` | Merge pileup and full alignment VCF/GVCF.
`RealignReads` | Reads local realignment for Illumina platform.
`SelectCandidates`| Select pileup candidates for full alignment calling.
`SelectHetSnp` | Select heterozygous SNP candidates for whatshap phasing.
`SelectQual` | Select quality cut-off for phasing and full alignment calling globally from all candidates.
`SortVcf` | Sort VCF file according to variants start position and contig name.
`SplitExtendBed` | Split bed file regions according to the contig name and extend bed region.
`UnifyRepresentation` | Representation unification for candidate site and true variant.
`Tensor2Bin` | Combine the variant and non-variant tensors and convert them to a binary, using `blosc:lz4hc` meta-compressor, the overall training memory is 10~15G. (pypy incompatible)

## Training Data

More details about the training data and download links in [here](docs/training_data.md).

|  Platform   |   Reference   |      Aligner      | Training Samples |
| :---------: | :-----------: | :---------------: | :--------------: |
|     ONT     | GRCh38_no_alt |     minimap2      | HG001,2,(3/4),5  |
| PacBio HiFi | GRCh38_no_alt |       pbmm2       |   HG001,2,4,5    |
|  Illumina   |    GRCh38     | BWA-MEM/NovoAlign |   HG001,2,4,5    |

#### Pretained Model

Download models from [here](http://www.bio8.cs.hku.hk/clair3/clair3_models/) or click on the links below.

|   Folder    |  Platform   | Training Sample | Default in Docker |     Link     |
| :---------: | :---------: | :-------------: | :---------------: | :----------: |
|     ont     |     ONT     |   HG001,2,4,5   |        Yes        | [Download]() |
|  ont_hg004  |     ONT     |   HG001,2,3,5   |                   | [Download]() |
|  illumina   |  Illumina   |   HG001,2,4,5   |        Yes        | [Download]() |
| pacbio_hifi | PacBio HiFi |   HG001,2,4,5   |        Yes        | [Download]() |

## VCF Output Format

`clair3/CallVariants.py` outputs variants in VCF format with version 4.2 specifications.
Clair3 includes pileup calling and full-alignment calling submodule. Pileup calling results are denoted with a `P` INFO tag, while the full-alignment calling results are denoted with a `F` INFO tag.