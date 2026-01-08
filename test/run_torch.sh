REF_GENOME="/autofs/bal19/zxzheng/testData/ont/data/GRCh38_no_alt_analysis_set.fasta"
THREADS="60"
MODEL_PATH="/autofs/nas5/xyu2/projects/Clair3_mv/hg002_hac/torch_model"
bash /autofs/nas5/xyu2/softwares/Clair3_no_change_v1.2_torch/run_clair3.sh \
    --bam_fn="/autofs/bal19/xyu/ont_open_data/HG001/hac_calls/downsample/HG001_30x.bam" \
    --ref_fn="${REF_GENOME}" \
    --threads=${THREADS} \
    --platform="ont" \
    --model_path="${MODEL_PATH}" \
    --output="/autofs/nas5/xyu2/softwares/Clair3_no_change_v1.2_torch/test/run_hac" \
    --enable_dwell_time \
    --ctg_name="chr20" 

# CLAIR3="${CLAIR3_PATH}/clair3.py"
# PYPY='/autofs/bal33/zxzheng/env/conda/envs/clair-somatic/bin/pypy3'
# WHATSHAP='/autofs/bal33/zxzheng/env/conda/envs/clair-somatic/bin/whatshap'
# PARALLEL='/autofs/bal33/zxzheng/env/conda/envs/clair-somatic/bin/parallel'
# TABIX='/autofs/bal33/zxzheng/env/conda/envs/clair-somatic/bin/tabix'
# SAMTOOLS='/autofs/bal33/zxzheng/env/conda/envs/clair-somatic/bin/samtools'

# PYTHON3='/autofs/bal19/zxzheng/env/conda/envs/mamba/envs/somatic/bin/python3'
# MODEL_PATH="/autofs/nas5/xyu2/projects/Clair3_mv/hg002_hac/model"
# bash /autofs/nas5/xyu2/softwares/Clair3_no_change_v1.2/run_clair3.sh \
#     --bam_fn="/autofs/bal19/xyu/ont_open_data/HG001/hac_calls/downsample/HG001_30x.bam" \
#     --ref_fn="${REF_GENOME}" \
#     --threads=${THREADS} \
#     --platform="ont" \
#     --model_path="${MODEL_PATH}" \
#     --output="/autofs/nas5/xyu2/softwares/Clair3_no_change_v1.2_torch/test/run_hac_tf" \
#     --enable_dwell_time \
#     --ctg_name="chr20" \
#     --samtools="${SAMTOOLS}" \
#     --pypy="${PYPY}" \
#     --parallel="${PARALLEL}" \
#     --whatshap="${WHATSHAP}" \
#     --python="${PYTHON3}" \


