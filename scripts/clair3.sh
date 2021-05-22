#!/bin/bash
SCRIPT_NAME=$(basename "$0")
Usage="Usage: ./${SCRIPT_NAME} --bam_fn=BAM --ref_fn=REF --output=OUTPUT_DIR --threads=THREADS --platform=PLATFORM --model_path=MODEL_PREFIX [--bed_fn=BED] [options]"
# INFO: whole calling workflow of clair3

set -e
ARGS=`getopt -o b:f:t:m:p:o:r::c::s::h::g \
-l bam_fn:,ref_fn:,threads:,model_path:,platform:,output:,\
bed_fn::,vcf_fn::,ctg_name::,sample_name::,help::,qual::,samtools::,python::,pypy::,parallel::,whatshap::,chunk_num::,chunk_size::,var_pct_full::,\
snp_min_af::,indel_min_af::,ref_pct_full::,pileup_only::,fast_mode::,gvcf::,print_ref_calls::,haploid_precise::,haploid_sensitive::,include_all_ctgs::,no_phasing_for_fa:: -n 'run_clair3.sh' -- "$@"`

if [ $? != 0 ] ; then echo"No input. Terminating...">&2 ; exit 1 ; fi
eval set -- "${ARGS}"

while true; do
   case "$1" in
    -b|--bam_fn ) BAM_FILE_PATH="$2"; shift 2 ;;
    -f|--ref_fn ) REFERENCE_FILE_PATH="$2"; shift 2 ;;
    -t|--threads ) THREADS="$2"; shift 2 ;;
    -m|--model_path ) MODEL_PATH="$2"; shift 2 ;;
    -p|--platform ) PLATFORM="$2"; shift 2 ;;
    -o|--output ) OUTPUT_FOLDER="$2"; shift 2 ;;
    --bed_fn ) BED_FILE_PATH="$2"; shift 2 ;;
    --vcf_fn ) VCF_FILE_PATH="$2"; shift 2 ;;
    --ctg_name ) CONTIGS="$2"; shift 2 ;;
    --sample_name ) SAMPLE="$2"; shift 2 ;;
    --chunk_num ) CHUNK_NUM="$2"; shift 2 ;;
    --chunk_size ) CHUNK_SIZE="$2"; shift 2 ;;
    --qual ) QUAL="$2"; shift 2 ;;
    --samtools ) SAMTOOLS="$2"; shift 2 ;;
    --python ) PYTHON="$2"; shift 2 ;;
    --pypy ) PYPY="$2"; shift 2 ;;
    --parallel ) PARALLEL="$2"; shift 2 ;;
    --whatshap ) WHATSHAP="$2"; shift 2 ;;
    --var_pct_full ) PRO="$2"; shift 2 ;;
    --ref_pct_full ) REF_PRO="$2"; shift 2 ;;
    --pileup_only ) PILEUP_ONLY="$2"; shift 2 ;;
    --fast_mode ) FAST_MODE="$2"; shift 2 ;;
    --print_ref_calls ) SHOW_REF="$2"; shift 2 ;;
    --gvcf ) GVCF="$2"; shift 2 ;;
    --snp_min_af ) SNP_AF="$2"; shift 2 ;;
    --indel_min_af ) INDEL_AF="$2"; shift 2 ;;
    --haploid_precise ) HAP_PRE="$2"; shift 2 ;;
    --haploid_sensitive ) HAP_SEN="$2"; shift 2 ;;
    --include_all_ctgs ) INCLUDE_ALL_CTGS="$2"; shift 2 ;;
    --no_phasing_for_fa ) NO_PHASING="$2"; shift 2 ;;

    -- ) shift; break; ;;
    -h|--help ) print_help_messages; break ;;
    * ) print_help_messages; exit 1 ;;
   esac
done


SHELL_FOLDER=$(cd "$(dirname "$0")";pwd)
CLAIR3="${SHELL_FOLDER}/../clair3.py"

if [ ${BED_FILE_PATH} = "EMPTY" ] ; then BED_FILE_PATH= ; fi

PILEUP_CHECKPOINT_PATH="${MODEL_PATH}/pileup"
FULL_ALIGNMENT_CHECKPOINT_PATH="${MODEL_PATH}/full_alignment"
LOG_PATH="${OUTPUT_FOLDER}/log"
TMP_FILE_PATH="${OUTPUT_FOLDER}/tmp"
SPLIT_BED_PATH="${TMP_FILE_PATH}/split_beds"
PILEUP_VCF_PATH="${TMP_FILE_PATH}/pileup_output"
GVCF_TMP_PATH="${TMP_FILE_PATH}/gvcf_tmp_output"
PHASE_OUTPUT_PATH="${TMP_FILE_PATH}/phase_output"
FULL_ALIGNMENT_OUTPUT_PATH="${TMP_FILE_PATH}/full_alignment_output"
PHASE_VCF_PATH="${PHASE_OUTPUT_PATH}/phase_vcf"
PHASE_BAM_PATH="${PHASE_OUTPUT_PATH}/phase_bam"
CANDIDATE_BED_PATH="${FULL_ALIGNMENT_OUTPUT_PATH}/candidate_bed"

echo "[INFO] Check envrionment variables"
${PYTHON} ${CLAIR3} CheckEnvs \
    --bam_fn ${BAM_FILE_PATH} \
    --bed_fn ${BED_FILE_PATH} \
    --output_fn ${OUTPUT_FOLDER} \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --vcf_fn ${VCF_FILE_PATH} \
    --ctg_name ${CONTIGS} \
    --chunk_num ${CHUNK_NUM} \
    --chunk_size ${CHUNK_SIZE} \
    --include_all_ctgs ${INCLUDE_ALL_CTGS} \
    --threads ${THREADS} \
    --python ${PYTHON} \
    --pypy ${PYPY} \
    --samtools ${SAMTOOLS} \
    --whatshap ${WHATSHAP} \
    --parallel ${PARALLEL} \
    --qual ${QUAL} \
    --var_pct_full ${PRO} \
    --ref_pct_full ${REF_PRO} \
    --snp_min_af ${SNP_AF} \
    --indel_min_af ${INDEL_AF}
readarray -t CHR < "${OUTPUT_FOLDER}/tmp/CONTIGS"
THREADS_LOW=$((${THREADS}*3/4))
if [[ ${THREADS_LOW} < 1 ]]; then THREADS_LOW=1; fi

cd ${OUTPUT_FOLDER}
# Pileup calling
#-----------------------------------------------------------------------------------------------------------------------
export CUDA_VISIBLE_DEVICES=""
echo "[INFO] 1/7 Calling variants using pileup model"
time ${PARALLEL} -C ' ' --joblog ${LOG_PATH}/parallel_1_call_var_bam_pileup.log -j ${THREADS_LOW} \
"${PYTHON} ${CLAIR3} CallVarBam \
    --chkpnt_fn ${PILEUP_CHECKPOINT_PATH} \
    --bam_fn ${BAM_FILE_PATH} \
    --call_fn ${PILEUP_VCF_PATH}/pileup_{1}_{2}.vcf \
    --sampleName ${SAMPLE} \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --extend_bed ${SPLIT_BED_PATH}/{1} \
    --bed_fn ${BED_FILE_PATH} \
    --vcf_fn ${VCF_FILE_PATH} \
    --ctgName {1} \
    --chunk_id {2} \
    --chunk_num {3} \
    --platform ${PLATFORM} \
    --fast_mode ${FAST_MODE} \
    --snp_min_af ${SNP_AF} \
    --indel_min_af ${INDEL_AF} \
    --gvcf ${GVCF} \
    --python ${PYTHON} \
    --pypy ${PYPY} \
    --samtools ${SAMTOOLS} \
    --temp_file_dir ${GVCF_TMP_PATH} \
    --pileup" :::: ${OUTPUT_FOLDER}/tmp/CHUNK_LIST |& tee ${LOG_PATH}/1_call_var_bam_pileup.log


echo "[INFO] Merge chunked contigs vcf files"
cat ${PILEUP_VCF_PATH}/pileup_*.vcf | ${PYPY} ${CLAIR3} SortVcf --output_fn ${OUTPUT_FOLDER}/pileup.vcf
bgzip -f ${OUTPUT_FOLDER}/pileup.vcf
tabix -f -p vcf ${OUTPUT_FOLDER}/pileup.vcf.gz

if [ ${PILEUP_ONLY} == True ]; then
    echo "[INFO] Only call pileup output with --pileup_only, output file: ${OUTPUT_FOLDER}/pileup.vcf.gz"
    echo "[INFO] Finish calling!"
    exit 1;
fi

# Whatshap phasing and haplotaging
#-----------------------------------------------------------------------------------------------------------------------
if [ ${NO_PHASING} == True ]
then
    echo "[INFO] 2/7 No phasing for full alignment calling"
    ${PARALLEL} -j${THREADS} ln -sf ${BAM_FILE_PATH} ${PHASE_BAM_PATH}/{1}.bam ::: ${CHR[@]}
    if [ -f ${BAM_FILE_PATH}.bai ]; then ${PARALLEL} -j${THREADS} ln -sf ${BAM_FILE_PATH}.bai ${PHASE_BAM_PATH}/{1}.bam.bai ::: ${CHR[@]}; fi
    if [ -f ${BAM_FILE_PATH%.*}.bai ]; then ${PARALLEL} -j${THREADS} ln -sf ${BAM_FILE_PATH%.*}.bai ${PHASE_BAM_PATH}/{1}.bam.bai ::: ${CHR[@]}; fi
else
    echo "[INFO] 2/7 Filter Hete SNP varaints for Whatshap phasing and haplotag"
    gzip -fdc ${OUTPUT_FOLDER}/pileup.vcf.gz | ${PYPY} ${CLAIR3} SelectQual --phase --output_fn ${PHASE_VCF_PATH}
    time ${PARALLEL} --joblog ${LOG_PATH}/parallel_2_select_hetero_snp.log -j${THREADS} \
    "${PYPY} ${CLAIR3} SelectHetSnp \
        --vcf_fn ${OUTPUT_FOLDER}/pileup.vcf.gz \
        --split_folder ${PHASE_VCF_PATH} \
        --ctgName {1}" ::: ${CHR[@]} ::: ${ALL_SAMPLE[@]} |& tee ${LOG_PATH}/2_select_hetero_snp.log

    echo "[INFO] 3/7 Whatshap phase vcf file"
    time ${PARALLEL} --joblog ${LOG_PATH}/parallel_3_phase.log -j${THREADS} \
    "${WHATSHAP} phase \
        --output ${PHASE_VCF_PATH}/phased_{1}.vcf.gz \
        --reference ${REFERENCE_FILE_PATH} \
        --chromosome {1} \
        --distrust-genotypes \
        --ignore-read-groups \
        ${PHASE_VCF_PATH}/{1}.vcf \
        ${BAM_FILE_PATH}" ::: ${CHR[@]} |& tee ${LOG_PATH}/3_phase.log
    ${PARALLEL} -j${THREADS} tabix -f -p vcf ${PHASE_VCF_PATH}/phased_{}.vcf.gz ::: ${CHR[@]}

    echo "[INFO] 4/7 Whatshap haplotag input bam file"
    time ${PARALLEL} --joblog ${LOG_PATH}/parallel_4_haplotag.log -j${THREADS} \
    "${WHATSHAP} haplotag \
        --output ${PHASE_BAM_PATH}/{1}.bam \
        --reference ${REFERENCE_FILE_PATH} \
        --ignore-read-groups \
        --regions {1} \
        ${PHASE_VCF_PATH}/phased_{1}.vcf.gz \
        ${BAM_FILE_PATH}" ::: ${CHR[@]} |& tee ${LOG_PATH}/4_haplotag.log
    time ${PARALLEL} -j${THREADS} ${SAMTOOLS} index -@12 ${PHASE_BAM_PATH}/{1}.bam ::: ${CHR[@]}
fi

# Full alignment calling
#-----------------------------------------------------------------------------------------------------------------------
echo "[INFO] 5/7 Select candidates for full alignment"
gzip -fdc ${OUTPUT_FOLDER}/pileup.vcf.gz | ${PYPY} ${CLAIR3} SelectQual --output_fn ${CANDIDATE_BED_PATH} \
--var_pct_full ${PRO} --ref_pct_full ${REF_PRO} --platform ${PLATFORM} --vcf_fn ${VCF_FILE_PATH}
time ${PARALLEL} --joblog ${LOG_PATH}/parallel_5_select_candidate.log -j${THREADS} \
"${PYPY} ${CLAIR3} SelectCandidates \
    --pileup_vcf_fn ${OUTPUT_FOLDER}/pileup.vcf.gz \
    --split_folder ${CANDIDATE_BED_PATH} \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --var_pct_full ${PRO} \
    --ref_pct_full ${REF_PRO} \
    --platform ${PLATFORM} \
    --ctgName {1}" ::: ${CHR[@]}  |& tee ${LOG_PATH}/5_select_candidate.log

echo "[INFO] 6/7 Calling variants using Full Alignment"
cat ${CANDIDATE_BED_PATH}/FULL_ALN_FILE_* > ${CANDIDATE_BED_PATH}/FULL_ALN_FILES
time ${PARALLEL} --joblog ${LOG_PATH}/parallel_6_call_var_bam_full_alignment.log -j ${THREADS_LOW} \
"${PYTHON} ${CLAIR3} CallVarBam \
    --chkpnt_fn ${FULL_ALIGNMENT_CHECKPOINT_PATH} \
    --bam_fn ${PHASE_BAM_PATH}/{1/.}.bam \
    --call_fn ${FULL_ALIGNMENT_OUTPUT_PATH}/full_alignment_{1/}.vcf \
    --sampleName ${SAMPLE} \
    --vcf_fn ${VCF_FILE_PATH} \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --full_aln_regions {1} \
    --ctgName {1/.} \
    --add_indel_length \
    --phasing_info_in_bam \
    --gvcf ${GVCF} \
    --python ${PYTHON} \
    --pypy ${PYPY} \
    --samtools ${SAMTOOLS} \
    --platform ${PLATFORM}" :::: ${CANDIDATE_BED_PATH}/FULL_ALN_FILES |& tee ${LOG_PATH}/6_call_var_bam_full_alignment.log

##Merge pileup and full alignment vcf
##-----------------------------------------------------------------------------------------------------------------------
cat ${FULL_ALIGNMENT_OUTPUT_PATH}/full_alignment_*.vcf | ${PYPY} ${CLAIR3} SortVcf --output_fn ${OUTPUT_FOLDER}/full_alignment.vcf
cat ${CANDIDATE_BED_PATH}/*.* > ${CANDIDATE_BED_PATH}/full_aln_regions
bgzip -f ${OUTPUT_FOLDER}/full_alignment.vcf
tabix -f -p vcf ${OUTPUT_FOLDER}/full_alignment.vcf.gz
if [ ${GVCF} == True ]; then cat ${GVCF_TMP_PATH}/*.tmp.g.vcf | ${PYPY} ${CLAIR3} SortVcf --output_fn ${GVCF_TMP_PATH}/non_var.gvcf; fi


echo "[INFO] 7/7 Merge pileup vcf and full alignment vcf"
time ${PARALLEL} --joblog ${LOG_PATH}/parallel_7_merge_vcf.log -j${THREADS} \
"${PYPY} ${CLAIR3} MergeVcf \
    --pileup_vcf_fn ${OUTPUT_FOLDER}/pileup.vcf.gz \
    --bed_fn ${CANDIDATE_BED_PATH}/full_aln_regions \
    --full_alignment_vcf_fn ${OUTPUT_FOLDER}/full_alignment.vcf.gz \
    --output_fn ${TMP_FILE_PATH}/merge_output/merge_{1}.vcf \
    --platform ${PLATFORM} \
    --print_ref_calls ${SHOW_REF} \
    --gvcf ${GVCF} \
    --haploid_precise ${HAP_PRE} \
    --haploid_sensitive ${HAP_SEN} \
    --gvcf_fn ${TMP_FILE_PATH}/merge_output/merge_{1}.gvcf \
    --non_var_gvcf_fn ${GVCF_TMP_PATH}/non_var.gvcf \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --ctgName {1}" ::: ${CHR[@]} |& tee ${LOG_PATH}/7_merge_vcf.log

cat ${TMP_FILE_PATH}/merge_output/merge_*.vcf | ${PYPY} ${CLAIR3} SortVcf --output_fn ${OUTPUT_FOLDER}/merge_output.vcf
if [ ${GVCF} == True ]; then cat ${TMP_FILE_PATH}/merge_output/merge_*.gvcf | ${PYPY} ${CLAIR3} SortVcf --output_fn ${OUTPUT_FOLDER}/merge_output.gvcf; fi
bgzip -f ${OUTPUT_FOLDER}/merge_output.vcf
tabix -f -p vcf ${OUTPUT_FOLDER}/merge_output.vcf.gz

echo "[INFO] Finish calling, output file: ${OUTPUT_FOLDER}/merge_output.vcf.gz"
