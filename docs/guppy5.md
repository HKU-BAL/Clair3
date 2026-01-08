# Guppy5 models

*Jun 9, 2021.* Guppy5 basecaller has been released recently. According to our observations, Guppy5 gives better base accuracy but has a substantially different error profile from the HAC mode in Guppy3 and Guppy4. While we are releasing Guppy5 models for Clair3 users to take advantage of the new Guppy5 basecaller, please aware that these models are not yet being heavily tested, and the Guppy5 data available for model training is still limited.

We **fine-tuned** the [HG001,2,4,5 (Guppy3,4)](http://www.bio8.cs.hku.hk/clair3/clair3_models/ont.tar.gz) model using ~70x HG002 data basecalled using the Guppy V5.0.7 using model `dna_r9.4.1_450bps_sup` from [ONT (visited May 25th, 2021)](https://labs.epi2me.io/gm24385_2021.05/). The new model is incorporated into `v0.1-r3` docker image with folder name `ont_guppy5`. The new model is also available [here](https://github.com/HKU-BAL/Clair3#pre-trained-models).

The benchmarks of the model and the data are as follows. With all precision increased slightly, the recall, espeicially the Indel recall, has improved significantly using Guppy5 data on a Guppy5 model (from 62.712% to 72.281%). The precision dropped when using Guppy5 data on the Guppy3,4 model, necessitating a Guppy5 model specially trained on Guppy5 data. The overall performance dropped signficiantly when using Guppy4 data on a Guppy5 model. We are not sure whether it's because the Guppy5 model was **fine-tuned** instead of trained from scratch, or it reflects the large differences between Guppy4 and Guppy5 data. A Guppy5 model trained from scratch will be provided when more Guppy5 data is available.

| Model          | Testing<br>Dataset | Overall<br>Precision | Overall<br>Recall | Overall<br>F1 | SNP<br>Precision | SNP<br>Recall | SNP<br>F1 | Indel<br>Precision | Indel<br>Recall | Indel<br>F1 |
| :--------------: | :---------------: | :----------: | :--------: | :---------------------: | :------------------: | :--------------: | :-----------------: | :--------------: | :----------: | :-------------------: |
| Guppy5 model   | Guppy5 HG002 chr20 ~70x  | 98.440%               | 96.100%            | 97.256%        | 99.812%           | 99.858%        | 99.835%    | 88.056%             | 72.281%          | 79.393%      |
| Guppy3,4 model | Guppy5 HG002 chr20 ~70x    | 97.294%               | 95.502%            | 96.390%        | 99.717%           | 99.858%        | 99.788%    | 79.616%             | 67.893%          | 73.289%      |
| Guppy3,4 model | Guppy4 HG002 chr20 70x    | 98.481%               | 94.709%            | 96.558%        | 99.780%           | 99.758%        | 99.769%    | 87.278%             | 62.713%          | 72.984%      |
| Guppy5 model   | Guppy4 HG002 chr20 70x    |94.184%               | 93.809%            | 93.996%        | 99.535%           | 99.634%        | 99.585%    | 59.359%             | 56.894%          | 58.100%      |

----

## Use the Guppy5 model

```
INPUT_DIR="[YOUR_INPUT_FOLDER]"           # e.g. /home/user1/input (absolute path needed)
OUTPUT_DIR="[YOUR_OUTPUT_FOLDER]"         # e.g. /home/user1/output (absolute path needed)
THREADS="[MAXIMUM_THREADS]"               # e.g. 8

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3::latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/input.bam \       ## change your bam file name here
  --ref_fn=${INPUT_DIR}/ref.fa \          ## change your reference name here
  --threads=${THREADS} \                  ## maximum threads to be used
  --platform="ont" \                       
  --model_path="/opt/models/r941_prom_sup_g506" \  
  --output=${OUTPUT_DIR}                  ## absolute output path prefix 
```

Check [Usage](https://github.com/HKU-BAL/Clair3#Usage) for more options.


