#!/bin/bash
SCRIPT_NAME=$(basename "$0")
Usage="\nUsage: ./${SCRIPT_NAME} -b BAM -f REF -o OUTPUT_DIR -t THREADS -p PLATFORM -m MODEL_PREFIX [--bed_fn=BED] [options]\n"

set -e
#./run_clair3.sh -b tmp.bam -f ref.fasta -t 32 -o tmp -p ont -m model_path
print_help_messages()
{
    echo ${Usage}
    echo $'Required parameters:'
    echo $'  -b, --bam_fn FILE        BAM file input. The input file must be samtools indexed.'
    echo $'  -f, --ref_fn FILE        FASTA reference file input. The input file must be samtools indexed.'
    echo $'  -m, --model_path STR     The folder path containing a Clair3 model (requiring six files in the folder, including pileup.data-00000-of-00001, pileup.index, pileup.meta, full_alignment.data-00000-of-00001, full_alignment.index, and full_alignment.meta).'
    echo $'  -t, --threads INT        Max #threads to be used. The full genome will be divided into small chucks for parallel processing. Each chunk will use 4 threads. The #chucks being processed simaltaneously is ceil(#threads/4)*3. 3 is the overloading factor.'
    echo $'  -p, --platform STR       Selete the sequencing platform of the input. Possible options: {ont,hifi,ilmn}.'
    echo $'  -o, --output PATH        VCF/GVCF output directory.'
    echo $''
    echo $'Optional parameters:'
    echo $'      --bed_fn FILE        Call variants only in the provided bed regions.'
    echo $'      --vcf_fn FILE        Candidate sites VCF file input, variants will only be called at the sites in the VCF file if provided.'
    echo $'      --ctg_name STR       The name of the sequence to be processed.'
    echo $'      --sample_name STR    Define the sample name to be shown in the VCF file.'
    echo $'      --qual INT           If set, variants with >=$qual will be marked PASS, or LowQual otherwise.'
    echo $'      --samtools STR       Path of samtools, samtools verision >= 1.10 is required.'
    echo $'      --python STR         Path of python, python3 >= 3.6 is required.'
    echo $'      --pypy STR           Path of pypy3, pypy3 >= 3.6 is required.'
    echo $'      --parallel STR       Path of parallel, parallel >= 20191122 is required.'
    echo $'      --whatshap STR       Path of whatshap, whatshap >= 1.0 is required.'
    echo $'      --chunk_size INT     The size of each chuck for parallel processing, default: 5Mbp.'
    echo $'      --pileup_only        Use only the pileup mode for calling, default: disable.'
    echo $'      --print_ref_calls    Show reference calls (0/0) in vcf file, default: disable.'
    echo $'      --gvcf               Enable GVCF output, default: disable.'
    echo $'      --snp_min_af FLOAT   Minimum SNP AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.08,hifi:0.08,ilmn:0.08.'
    echo $'      --indel_min_af FLOAT Minimum INDEL AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.15,hifi:0.08,ilmn:0.08.'
    echo $'      --fast_mode          EXPERIMENTAL: Skip variant candidates with AF <= 0.15, default: disable.'
    echo $'      --haploid_precise    EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant, default: disable.'
    echo $'      --haploid_sensitive  EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant, default: disable.'
    echo $'      --var_pct_full FLOAT EXPERIMENTAL: Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for full-alignment mode calling, default: 0.3.'
    echo $'      --ref_pct_full FLOAT EXPERIMENTAL: Specify an expected percentage of low quality 0/0 variants called in the pileup mode for full-alignment mode calling, default: 0.3 for illumina and pb, 0.1 for ont.'
    echo $''

    exit 1
}

ARGS=`getopt -o b:f:t:m:p:o: \
-l bam_fn:,ref_fn:,threads:,model_path:,platform:,output:,\
bed_fn::,vcf_fn::,ctg_name::,sample_name::,help::,qual::,samtools::,python::,pypy::,parallel::,whatshap::,chunk_num::,chunk_size::,var_pct_full::,ref_pct_full::,\
snp_min_af::,indel_min_af::,fast_mode,gvcf,pileup_only,print_ref_calls,haploid_precise,haploid_sensitive -n 'run_clair3.sh' -- "$@"`

if [ $? != 0 ] ; then echo"No input. Terminating...">&2 ; exit 1 ; fi
eval set -- "${ARGS}"

# default options
SAMPLE="EMPTY"
BED_FILE_PATH="EMPTY"
VCF_FILE_PATH='EMPTY'
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
REF_PRO=0.3
GVCF=False
PILEUP_ONLY=False
FAST_MODE=False
SHOW_REF=False
SNP_AF=0.0
INDEL_AF=0.0
HAP_PRE=False
HAP_SEN=False

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
    --snp_min_af ) SNP_AF="$2"; shift 2 ;;
    --indel_min_af ) INDEL_AF="$2"; shift 2 ;;
    --gvcf ) GVCF=True; shift 1 ;;
    --pileup_only ) PILEUP_ONLY=True; shift 1 ;;
    --fast_mode ) FAST_MODE=True; shift 1 ;;
    --print_ref_calls ) SHOW_REF=True; shift 1 ;;
    --haploid_precise ) HAP_PRE=True; shift 1 ;;
    --haploid_sensitive ) HAP_SEN=True; shift 1 ;;

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
if [ ! ${VCF_FILE_PATH} = "EMPTY" ] && [ ! -z ${VCF_FILE_PATH} ] && [ ! -f ${VCF_FILE_PATH} ] ; then echo "[ERROR] Vcf file provides and not found"; exit 1; fi


mkdir -p ${OUTPUT_FOLDER}

#optional parameters should use "="
(time (
echo "[INFO] BAM FILE PATH: ${BAM_FILE_PATH}"
echo "[INFO] REFERENCE FILE PATH: ${REFERENCE_FILE_PATH}"
echo "[INFO] MODEL PATH: ${MODEL_PATH}"
echo "[INFO] OUTPUT FOLDER: ${OUTPUT_FOLDER}"
echo "[INFO] PLATFORM: ${PLATFORM}"
echo "[INFO] THREADS: ${THREADS}"
echo "[INFO] BED FILE PATH: ${BED_FILE_PATH}"
echo "[INFO] VCF FILE PATH: ${VCF_FILE_PATH}"
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
echo "[INFO] ENABLE FILEUP ONLY CALLING: ${PILEUP_ONLY}"
echo "[INFO] ENABLE FAST MODE CALLING: ${FAST_MODE}"
echo "[INFO] ENABLE PRINTING REFERENCE CALLS: ${SHOW_REF}"
echo "[INFO] ENABLE OUTPUT GVCF: ${GVCF}"
echo "[INFO] ENABLE HAPLOID PRECISE MODE: ${HAP_PRE}"
echo "[INFO] ENABLE HAPLOID SENSITIVE MODE: ${GVCF}"
echo $''
scripts/clair3.sh \
    --bam_fn ${BAM_FILE_PATH} \
    --ref_fn ${REFERENCE_FILE_PATH} \
    --threads ${THREADS} \
    --model_path ${MODEL_PATH} \
    --platform ${PLATFORM} \
    --output ${OUTPUT_FOLDER} \
    --bed_fn=${BED_FILE_PATH} \
    --vcf_fn=${VCF_FILE_PATH} \
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
    --var_pct_full=${PRO} \
    --ref_pct_full=${REF_PRO} \
    --snp_min_af=${SNP_AF} \
    --indel_min_af=${INDEL_AF} \
    --pileup_only=${PILEUP_ONLY} \
    --gvcf=${GVCF} \
    --fast_mode=${FAST_MODE} \
    --print_ref_calls=${SHOW_REF} \
    --haploid_precise=${HAP_PRE} \
    --haploid_sensitive=${HAP_SEN} \


)) |& tee ${OUTPUT_FOLDER}/run_clair3.log
