#!/bin/bash
SCRIPT_NAME=$(basename "$0")
Usage="\nUsage: ./${SCRIPT_NAME} -b BAM -f REF -o OUTPUT_DIR -t THREADS -p PLATFORM -m MODEL_PATH [--bed_fn=BED] [options]\n"

#./run_clair3.sh -b tmp.bam -f ref.fasta -t 32 -o tmp -p ont -m model_path
print_help_messages()
{
    echo -e ${Usage}
    echo $'Required input parameters:'
    echo $'  -b, --bam_fn FILE      Indexed BAM file input.\n'
    echo $'  -f, --ref_fn FILE      Indexed FASTA reference file input. \n'
    echo $'  -t, --threads INT      Number of additional threads to use.\n'
    echo $'  -m, --model_path STR   Model checkpoint path for calling.\n'
    echo $'  -p, --platform STR     Select which platform for variant calling, optional: [ont pb illumina].\n\n'

    echo $'Optional input parameters:'
    echo $'      --bed_fn FILE      Call variant only in provided bed regions, optional..\n'
    echo $'      --ctg_name STR     Checkpoint model path for calling, optional..\n'
    echo $'      --sample_name STR  Define the sample name to be shown in the VCF file, optional..\n'
    echo $'      --chunk_num INT    Total chunk of each number for parallel execution. Each chunk refer to a smaller reference regions, optional.\n'
    echo $'      --qual INT         If set, variant with equal or higher quality will be marked PASS, or LowQual otherwise, optional.\n'
    echo $'      --samtools STR     Path of samtools, samtools verision >= 1.10 is required.\n'
    echo $'      --python STR       Path of python, python3 >= 3.6 is required. \n'
    echo $'      --pypy STR         Path of pypy3, pypy3 >= 3.6 is required. \n'
    echo $'      --parallel STR     Path of parallel, parallel >= 20191122 is required. \n'
    echo $'      --whatshap STR     Path of whatshap, whatshap >= 1.0 is required. \n'
    echo $'      --chunk_size INT   Define the chunk size for each threads for processing, default 3000000. \n'
    echo $'      --chunk_num INT    Total chunk number for parallel execution. If set, will ignore the "chunk_size" , optional.\n'
    echo $'      --fast_mode        Ignore low allelic frequency <= 0.15 snp calling for ont platform,  optional.\n'
    echo $'      --proportion FLOAT Full alignment calling proportion. \n'
    echo $'      --ref_proportion FLOAT\n'
    echo $'                         Full alignment reference calling proportion. \n'
    echo $'      --threshold_for_snp_only FLOAT\n'
    echo $'                         SNP allele frequence below this threshold will be filter. \n'
    echo $'      --threshold_for_indel_only FLOAT\n'
    echo $'                         INDEL allele frequence below this threshold will be filter. \n'

    echo $'Output parameters:'
    echo $'  -o, --output PATH      Output (vcf/gvcf) output directory.\n'
    echo $'  -g, --gvcf             Whether to generate gvcf, default: False.\n'
    echo $'      --pileup_only      Only call pileup output, default: False.\n'
    exit 1
}

ARGS=`getopt -o b:f:t:m:p:o:r::c::s::h::g \
-l bam_fn:,ref_fn:,threads:,model_path:,platform:,output:,\
bed_fn::,ctg_name::,sample_name::,help::,qual::,samtools::,python::,pypy::,parallel::,whatshap::,chunk_num::,chunk_size::,proportion::,ref_proportion::,\
threshold_for_snp_only::,threshold_for_indel_only::,fast_mode,gvcf,pileup_only -n 'run_clair3.sh' -- "$@"`

if [ $? != 0 ] ; then echo"No input. Terminating...">&2 ; exit 1 ; fi
eval set -- "${ARGS}"

# default options
SAMPLE="EMPTY"
BED_FILE_PATH="EMPTY"
CONTIGS="EMPTY"
SAMTOOLS="samtools"
PYPY="pypy3"
PYTHON='python3'
PARALLEL='parallel'
WHATSHAP='whatshap'
CHUNK_NUM=0
CHUNK_SIZE=5000000
QUAL=0
PRO=0.3
REF_PRO=0.1
GVCF=False
PILEUP_ONLY=False
FAST_MODE=False
SNP_AF=0.0
INDEL_AF=0.0

while true; do
   case "$1" in
    -b|--bam_fn ) BAM_FILE_PATH="$2"; shift 2 ;;
    -f|--ref_fn ) REFERENCE_FILE_PATH="$2"; shift 2 ;;
    -t|--threads ) THREADS="$2"; shift 2 ;;
    -m|--model_path ) MODEL_PATH="$2"; shift 2 ;;
    -p|--platform ) PLATFORM="$2"; shift 2 ;;
    -o|--output ) OUTPUT_FOLDER="$2"; shift 2 ;;
    -r|--bed_fn ) BED_FILE_PATH="$2"; shift 2 ;;
    -c|--ctg_name ) CONTIGS="$2"; shift 2 ;;
    --sample_name ) SAMPLE="$2"; shift 2 ;;
    --chunk_num ) CHUNK_NUM="$2"; shift 2 ;;
    --chunk_size ) CHUNK_SIZE="$2"; shift 2 ;;
    --qual ) QUAL="$2"; shift 2 ;;
    --samtools ) SAMTOOLS="$2"; shift 2 ;;
    --python ) PYTHON="$2"; shift 2 ;;
    --pypy ) PYPY="$2"; shift 2 ;;
    --parallel ) PARALLEL="$2"; shift 2 ;;
    --whatshap ) WHATSHAP="$2"; shift 2 ;;
    --proportion ) PRO="$2"; shift 2 ;;
    --ref_proportion ) REF_PRO="$2"; shift 2 ;;
    --threshold_for_snp_only ) SNP_AF="$2"; shift 2 ;;
    --threshold_for_indel_only ) INDEL_AF="$2"; shift 2 ;;
    --gvcf ) GVCF=True; shift 1 ;;
    --pileup_only ) PILEUP_ONLY=True; shift 1 ;;
    --fast_mode ) FAST_MODE=True; shift 1 ;;

    -- ) shift; break; ;;
    -h|--help ) print_help_messages; break ;;
    * ) print_help_messages; exit 1 ;;
   esac
done

if [ -z ${BAM_FILE_PATH} ] || [ -z ${REFERENCE_FILE_PATH} ] || [ -z ${THREADS} ] || [ -z ${OUTPUT_FOLDER} ] || [ -z ${PLATFORM} ] || [ -z ${MODEL_PATH} ]; then
      print_help_messages
      echo "[ERROR] Required parameters missing";
fi

if [ ! -f ${BAM_FILE_PATH} ] || [ ! -f ${BAM_FILE_PATH}.bai ]; then echo "[ERROR] Bam file or Bam index bai file not found"; exit 1; fi
if [ ! -f ${REFERENCE_FILE_PATH} ] || [ ! -f ${REFERENCE_FILE_PATH}.fai ]; then echo "[ERROR] Reference file or Reference index fai file not found"; exit 1; fi
if [ ! ${BED_FILE_PATH} = "EMPTY" ] && [ ! -z ${BED_FILE_PATH} ] && [ ! -f ${BED_FILE_PATH} ] ; then echo "[ERROR] Bed file provides and not found"; exit 1; fi



mkdir -p ${OUTPUT_FOLDER}

#optional parameters should use "="
time (
echo "[INFO] BAM FILE PATH: ${BAM_FILE_PATH}"
echo "[INFO] REFERENCE FILE PATH: ${REFERENCE_FILE_PATH}"
echo "[INFO] MODEL PATH: ${MODEL_PATH}"
echo "[INFO] OUTPUT FOLDER: ${OUTPUT_FOLDER}"
echo "[INFO] PLATFORM: ${PLATFORM}"
echo "[INFO] THREADS: ${THREADS}"
echo "[INFO] BED FILE PATH: ${BED_FILE_PATH}"
echo "[INFO] CONTIGS: ${CONTIGS}"
echo "[INFO] SAMTOOLS PATH: ${SAMTOOLS}"
echo "[INFO] PYTHON PATH: ${PYTHON}"
echo "[INFO] PYPY PATH: ${PYPY}"
echo "[INFO] PARALLEL PATH: ${PARALLEL}"
echo "[INFO] WHATSHAP PATH: ${WHATSHAP}"
echo "[INFO] CHUNK SIZE: ${CHUNK_SIZE}"
echo "[INFO] CHUNK NUM: ${CHUNK_NUM}"
echo "[INFO] FULL ALIGN PROPORTION: ${PRO}"
echo "[INFO] FULL ALIGN RERFERENCE PROPORTION: ${REF_PRO}"
echo "[INFO] USER DEFINED SNP THRESHOLD: ${SNP_AF}"
echo "[INFO] USER DEFINED INDEL THRESHOLD: ${INDEL_AF}"
echo "[INFO] FILEUP ONLY CALLING: ${PILEUP_ONLY}"
echo "[INFO] FAST MODE CALLING: ${FAST_MODE}"
echo "[INFO] OUTPUT GVCF: ${GVCF}"
echo
scripts/clair3.sh \
    --bam_fn ${BAM_FILE_PATH} \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --threads ${THREADS} \
    --model_path ${MODEL_PATH} \
    --platform ${PLATFORM} \
    --output ${OUTPUT_FOLDER} \
    --bed_fn=${BED_FILE_PATH} \
    --ctg_name=${CONTIGS} \
    --sample_name=${SAMPLE} \
    --chunk_num=${CHUNK_NUM} \
    --chunk_size=${CHUNK_SIZE} \
    --samtools=${SAMTOOLS} \
    --python=${PYTHON} \
    --pypy=${PYPY} \
    --parallel=${PARALLEL} \
    --whatshap=${WHATSHAP} \
    --qual=${QUAL} \
    --proportion=${PRO} \
    --ref_proportion=${REF_PRO} \
    --threshold_for_snp_only=${SNP_AF} \
    --threshold_for_indel_only=${INDEL_AF} \
    --pileup_only=${PILEUP_ONLY} \
    --gvcf=${GVCF} \
    --fast_mode=${FAST_MODE}

) |& tee ${OUTPUT_FOLDER}/run_clair3.log






