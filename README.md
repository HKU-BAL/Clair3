
# Clair3 - Integrate pileup and full-alignment for high-performance long-read variant calling

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



### Command to run Clair3

```bash
#optional parameters should use "="
./run_clair3.sh \
    -b "/mnt/bal36/zxzheng/testData/ont/data/input.bam" \
    -f "/mnt/bal36/zxzheng/testData/ont/data/ref.fa" \
    -m "/mnt/bal36/zxzheng/testData/ont/model"\
    -t 40 \
    -p "ont" \
    -o "/mnt/bal36/zxzheng/testData/ont/test/chr20_github" \
    --ctg_name='chr20'

##Pileup output file: ${OUTPUT_FOLDER}/pileup.vcf.gz
##Full alignment output file: ${OUTPUT_FOLDER}/full_alignment.vcf
##FInal merge output file: ${OUTPUT_FOLDER}/merge_output.vcf.gz
```

#### 