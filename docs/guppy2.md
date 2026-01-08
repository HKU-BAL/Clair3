# Guppy2 models

*Jun 28, 2021.*  The ONT model trained with Guppy3 and Guppy4 data performed even worse than Clair on Guppy2 data. We collected and trained a model specially for reads basecalled using Guppy2 or prior. The new model is available [here](https://github.com/HKU-BAL/Clair3#pre-trained-models) (first released with `v0.1-r4`). Below are some benchmarks (GIAB v4.2.1) showing a signficiant improvement when using the correct model for your data.

| Model/Caller    | Testing<br>Dataset      | Overall-<br>Precision | Overall-<br>Recall | Overall-<br>F1 | SNP-<br>Precision | SNP-<br>Recall | SNP-<br>F1 | Indel-<br>Precision | Indel-<br>Recall | Indel-<br>F1 |
| --------------- | ----------------------- | --------------------- | ------------------ | -------------- | ----------------- | -------------- | ---------- | ------------------- | ---------------- | ------------ |
| Guppy2 model    | Guppy2 HG002 chr20 ~64x |  98.365%               | 92.406%            | 95.292%        | 99.541%           | 99.080%        | 99.310%    | 85.847%             | 50.107%          | 63.279%      |
| Guppy3,4 model  | Guppy2 HG002 chr20 ~64x |  72.329%               | 90.896%            | 80.556%        | 99.337%           | 98.536%        | 98.935%    | 14.648%             | 42.475%          | 21.784%      |
| Guppy2 model    | Guppy4 HG002 chr20 ~64x |  86.897%               | 92.707%            | 89.708%        | 97.983%           | 99.579%        | 98.775%    | 35.773%             | 49.156%          | 41.410%      |
| Guppy3,4 model  | Guppy4 HG002 chr20 ~64x |  98.512%               | 94.601%            | 96.517%        | 99.759%           | 99.749%        | 99.754%    | 87.586%             | 61.976%          | 72.588%      |
| Clair ONT model | Guppy2 HG002 chr20 ~64x | 98.007%               | 87.739%            | 92.590%        | 99.290%           | 96.448%        | 97.848%    | 79.013%             | 32.552%          | 46.108%      |

----

## Use the Guppy2 model

```
INPUT_DIR="[YOUR_INPUT_FOLDER]"           # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"         # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"               # e.g. 8

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \       ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \          ## change your reference name here
  --threads=${THREADS} \                  ## maximum threads to be used
  --platform="ont" \                       
  --model_path="/opt/models/r941_prom_hac_g238" \  
  --output=${OUTPUT_DIR}                  ## absolute output path prefix 
```

Check [Usage](https://github.com/HKU-BAL/Clair3#Usage) for more options.

----

## Guppy2 training data

| Reference | Sample | Aligner  | Coverage |  Basecaller  |                             link                             |
| :-------: | :------: | :------: | :------: | :----------: | :----------------------------------------------------------: |
|  GRCh38   | HG001 | minimap2 |  ~124x   | Guppy v2.2.3 |                              -                               |
|  GRCh38   | HG002 | minimap2 |   ~64x   | Guppy v2.3.8 | [link](http://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/Ultralong_OxfordNanopore/guppy-V2.3.4_2019-06-26/ultra-long-ont_GRCh38_reheader.bam) |
|  GRCh38   | HG003 | minimap2 |   ~91x   | Guppy v2.3.8 |    [link](https://github.com/human-pangenomics/hpgp-data)    |
|  GRCh38   | HG004 | minimap2 |   ~90x   | Guppy v2.3.8 |    [link](https://github.com/human-pangenomics/hpgp-data)    |

