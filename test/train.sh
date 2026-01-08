PYTHON3='python'
CLAIR3_PATH="/autofs/nas5/xyu2/softwares/Clair3_no_change_v1.2_torch"
CLAIR3="${CLAIR3_PATH}/clair3.py"
OUTPUT_DIR="/autofs/nas5/xyu2/softwares/Clair3_no_change_v1.2_torch/test/hac_train"
mkdir -p ${OUTPUT_DIR}

${PYTHON3} ${CLAIR3} Train \
    --bin_fn /autofs/nas5/xyu2/projects/Clair3_mv/hg002_hac/full_alignment/hg002/build/bins \
    --ochk_prefix ${OUTPUT_DIR} \
    --add_indel_length True \
    --random_validation \
    --platform ont \
    --learning_rate 1e-4 \
    --enable_dwell_time |& tee ${OUTPUT_DIR}/log.txt