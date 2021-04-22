
# Clair3 - Integrating pileup and full-alignment for high-performance long-read variant calling

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause) [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/clair/README.html) \
Contact: Ruibang Luo \
Email: rbluo@cs.hku.hk

---

## Introduction



This is the formal release of Clair3, the successor of Clair. Clair is published in Nature Machine Intelligence. A preprint is available in bioRxiv

---

## Contents

- [Installation](#installation)



---

## Installation

### Option 1. Build an anaconda virtual environment step by step

#### Please install anaconda using the installation guide at https://docs.anaconda.com/anaconda/install/

```bash
# create and activate the environment named clair3
conda create -n clair3 python=3.6
source activate clair3

# install pypy and packages on clair3 environemnt
conda install -c conda-forge pypy3.6
pypy3 -m ensurepip
pypy3 -m pip install intervaltree==3.0.2
pypy3 -m pip install mpmath==1.2.1 python-Levenshtein==0.12.0 

# install python packages on clair3 environment
pip3 install tensorflow==2.2.0
pip3 install intervaltree==3.0.2  tensorflow-addons==0.11.2 tables==3.6.1 python-Levenshtein==0.12.0
conda install -c anaconda pigz==2.4
conda install -c conda-forge parallel=20191122 zstd=1.4.4
conda install -c conda-forge -c bioconda samtools=1.10 vcflib=1.0.0 bcftools=1.10.2
conda install -c conda-forge -c bioconda whatshap=1.0

# clone Clair3
git clone --depth 1 https://github.com/HKU-BAL/Clair3.git
cd Clair3
chmod +x run_clair3.sh
chmod +x scripts/clair3.sh

# run clair3 like this afterwards
./run_clair3.sh --help
```

### Option 3. Docker

** **

```bash
# clone Clair3
git clone --depth 1 https://github.com/HKU-BAL/Clair3.git
cd Clair3

# build a docker image named clair3_docker
docker build -f ./Dockerfile -t clair3_docker . # 

# run docker image
# You might require docker authentication to build by docker
docker run -it clair3 # You might need root privilege

# run clair like this afterwards
./run_clair3.sh
```



## Usage

### General usage

```bash
#optional parameters should use "="
./run_clair3.sh \
    -b ${BAM} \
    -f ${REF} \
    -m ${MODEL_PREFIX} \
    -t ${THREADS} \
    -p ${PLATFORM} \  #{ont,hifi,ilmn}
    -o ${OUTPUT_DIR}

##Pileup output file: ${OUTPUT_DIR}/pileup.vcf.gz
##Full alignment output file: ${OUTPUT_DIR}/full_alignment.vcf
##Final merge output file: ${OUTPUT_DIR}/merge_output.vcf.gz
```

### Options

**Required parameters:**

```bash
  -b, --bam_fn FILE        BAM file input. The input file must be samtools indexed.
  -f, --ref_fn FILE        FASTA reference file input. The input file must be samtools indexed.
  -m, --model_path STR     The folder path containing a Clair3 model (requiring six files in the folder, including pileup.data-00000-of-00001, pileup.index, pileup.meta, full_alignment.data-00000-of-00001, full_alignment.index, and full_alignment.meta).
  -t, --threads INT        Max #threads to be used. The full genome will be divided into small chucks for parallel processing. Each chunk will use 4 threads. The #chucks being processed simaltaneously is ceil(#threads/4)*3. 3 is the overloading factor.
  -p, --platform STR       Selete the sequencing platform of the input. Possible options: {ont,hifi,ilmn}.
  -o, --output PATH        VCF/GVCF output directory.

```

**Optional parameters:**

```
Optional parameters:
      --bed_fn FILE        Call variants only in the provided bed regions.
      --vcf_fn FILE        Candidate sites VCF file input, variants will only be called at the sites in the VCF file if provided.
      --ctg_name STR       The name of the sequence to be processed.
      --sample_name STR    Define the sample name to be shown in the VCF file.
      --qual INT           If set, variants with >=$qual will be marked PASS, or LowQual otherwise.
      --samtools STR       Path of samtools, samtools verision >= 1.10 is required.
      --python STR         Path of python, python3 >= 3.6 is required.
      --pypy STR           Path of pypy3, pypy3 >= 3.6 is required.
      --parallel STR       Path of parallel, parallel >= 20191122 is required.
      --whatshap STR       Path of whatshap, whatshap >= 1.0 is required.
      --chunk_size INT     The size of each chuck for parallel processing, default: 5Mbp.
      --pileup_only        Use only the pileup mode for calling, default: disable.
      --print_ref_calls    Show reference calls (0/0) in vcf file, default: disable.
      --include_all_ctgs   Call variants on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}, default: disable.
      --gvcf               Enable GVCF output, default: disable.
      --snp_min_af FLOAT   Minimum SNP AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.08,hifi:0.08,ilmn:0.08.
      --indel_min_af FLOAT Minimum INDEL AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.15,hifi:0.08,ilmn:0.08.
      --var_pct_full FLOAT EXPERIMENTAL: Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for full-alignment mode calling, default: 0.3.
      --ref_pct_full FLOAT EXPERIMENTAL: Specify an expected percentage of low quality 0/0 variants called in the pileup mode for full-alignment mode calling, default: 0.3 for illumina and pb, 0.1 for ont.
      --fast_mode          EXPERIMENTAL: Skip variant candidates with AF <= 0.15, default: disable.
      --haploid_precise    EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant, default: disable.
      --haploid_sensitive  EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant, default: disable.
      --no_ra_for_illmn    EXPERIMENTAL: Call variants without reads realignment for illumina platform, default: disable.
      --no_phasing_for_fa  EXPERIMENTAL: Call variants without whatshap phasing in full alignment calling, default: disable.
```





## Folder Structure and Submodule Descriptions

Submodules in __`clair3/`__ are for variant calling and model training. Submodules in __`preprocess`__ are for data preparation.

*For all the submodules listed below, you can use the `-h` or `--help` option for available options.*

`clair3/` | Note: submodules under this folder are pypy incompatible, please run using python
---: | ---
`CallVariants` | Call variants using a trained model and tensors of candidate variants.
`CallVarBam` | Call variants using a trained model and a BAM file.
`Train` |  Training a model using warm-up learning rate startup with `RectifiedAdam` optimizer. We also use `Lookahead` optimizer to  adaptively adjust `RectifiedAdam` parameters. The initial learning rate is `1e-3` with 0.1 learning rate warm-up. Input a binary tensors are created by `Tensor2Bin`, using `blosc:lz4hc` meta-compressor, the overall training memory is only 10~15G.

`preprocess/` | Note: submodules under this folder is Pypy compatible unless specified.
---: | ---
`CheckEnvs`| Check the environment and the validity of the input variables, preprocess the BED input if necessary, `--chunk_size` set the genome chuck size per job.<br>
`CreateTensorPileup`| Generate variant candidate tensors using pileup for training or calling.
`CreateTensorFullAlignment`| Generate variant candidate tensors using phased full-alignment for training or calling.<br>
`GetTruth`| Extract the variants from a truth VCF. Input: VCF; Reference FASTA if the vcf contains asterisks in ALT field.<br>`RealignReads` | Reads local realignment for illumina platform.<br>`SelectCandidates` | Select pileup candidates for full alignment calling.<br>`SelectHetSnp` | Select heterozygous snp candidates for WhatsHap phasing.<br>`SelectQual` | Select quality cut-off for phasing and full alignment calling globally from all candidates.<br>`MergeVcf` | Merge pileup and full alignment VCF/GVCF.<br>`UnifyRepresentation` | Representation unification for candidate site and true variant.<br>

`Tensor2Bin` | Combine the variant and non-variant tensors and convert them to a binary.<br>(pypy incompatible)



## VCF Output Format

`clair3/CallVariants.py` outputs variants in VCF format with version 4.2 specifications.
Clair3 includes pileup calling and full-alignment calling submodule. For result of pileup calling are denoted with a `P` INFO tag, while the full-alignment calling result are denoted with a `F` INFO tag.