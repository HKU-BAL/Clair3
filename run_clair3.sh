#!/usr/bin/env bash
SCRIPT_NAME=$(basename "$0")
SCRIPT_PATH=$(dirname $(readlink -f "$0"))
VERSION='v1.1.2'
Usage="Usage: ${SCRIPT_NAME} --bam_fn=BAM --ref_fn=REF --output=OUTPUT_DIR --threads=THREADS --platform=PLATFORM --model_path=MODEL_PREFIX [--bed_fn=BED] [options]"

CMD="$0 $@"

set -e -o pipefail
#./run_clair3.sh -b tmp.bam -f ref.fasta -t 32 -o tmp -p ont -m model_path
print_help_messages()
{
    echo $''
    echo ${Usage}
    echo $''
    echo $'Required parameters:'
    echo $'  -b, --bam_fn=FILE             BAM file input. The input file must be samtools indexed.'
    echo $'  -f, --ref_fn=FILE             FASTA reference file input. The input file must be samtools indexed.'
    echo $'  -m, --model_path=STR          The folder path containing a Clair3 model (requiring six files in the folder, including pileup.data-00000-of-00002, pileup.data-00001-of-00002 pileup.index, full_alignment.data-00000-of-00002, full_alignment.data-00001-of-00002 and full_alignment.index).'
    echo $'  -t, --threads=INT             Max #threads to be used. The full genome will be divided into small chunks for parallel processing. Each chunk will use 4 threads. The #chunks being processed simultaneously is ceil(#threads/4)*3. 3 is the overloading factor.'
    echo $'  -p, --platform=STR            Select the sequencing platform of the input. Possible options: {ont,hifi,ilmn}.'
    echo $'  -o, --output=PATH             VCF/GVCF output directory.'
    echo $''
    echo $''
    echo $"Optional parameters (Use \"=value\" instead of \" value\". E.g., \"--bed_fn=fn.bed\" instead of \"--bed_fn fn.bed\".):"
    echo $'      --bed_fn=FILE             Call variants only in the provided bed regions.'
    echo $'      --vcf_fn=FILE             Candidate sites VCF file input, variants will only be called at the sites in the VCF file if provided.'
    echo $'      --ctg_name=STR            The name of the sequence to be processed.'
    echo $'      --sample_name=STR         Define the sample name to be shown in the VCF file.'
    echo $'      --qual=INT                If set, variants with >$qual will be marked PASS, or LowQual otherwise.'
    echo $'      --samtools=STR            Path of samtools, samtools version >= 1.10 is required.'
    echo $'      --python=STR              Path of python, python3 >= 3.6 is required.'
    echo $'      --pypy=STR                Path of pypy3, pypy3 >= 3.6 is required.'
    echo $'      --parallel=STR            Path of parallel, parallel >= 20191122 is required.'
    echo $'      --whatshap=STR            Path of whatshap, whatshap >= 1.0 is required.'
    echo $'      --longphase=STR           Path of longphase, longphase >= 1.0 is required.'
    echo $'      --chunk_size=INT          The size of each chuck for parallel processing, default: 5000000.'
    echo $'      --pileup_only             Use the pileup model only when calling, default: disable.'
    echo $'      --print_ref_calls         Show reference calls (0/0) in VCF file, default: disable.'
    echo $'      --include_all_ctgs        Call variants on all contigs, otherwise call in chr{1..22,X,Y} and {1..22,X,Y}, default: disable.'
    echo $'      --gvcf                    Enable GVCF output, default: disable.'
    echo $'      --use_whatshap_for_intermediate_phasing
                                           Phase high-quality heterozygous variants using whatshap for full-alignment model calling, default: enable.'
    echo $'      --use_longphase_for_intermediate_phasing
                                           Phase high-quality heterozygous variants using longphase for full-alignment model calling, default: disable.'
    echo $'      --use_whatshap_for_final_output_phasing
                                           Phase the output variants using whatshap, default: disable.'
    echo $'      --use_longphase_for_final_output_phasing
                                           Phase the output variants using longphase, default: disable.'
    echo $'      --use_whatshap_for_final_output_haplotagging
                                           Haplotag input BAM using output phased variants using whatshap, default: disable.'
    echo $'      --enable_phasing          It means `--use_whatshap_for_final_output_phasing`. The option is retained for backward compatibility.'
    echo $'      --longphase_for_phasing   It means `--use_longphase_for_intermediate_phasing`. The option is retained for backward compatibility.'
    echo $'      --disable_c_impl          Disable C implement with cffi for pileup and full-alignment create tensor, default: enable.'
    echo $'      --remove_intermediate_dir Remove intermediate directory, including intermediate phased BAM, pileup and full-alignment results. default: disable.'
    echo $'      --snp_min_af=FLOAT        Minimum SNP AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.08,hifi:0.08,ilmn:0.08.'
    echo $'      --indel_min_af=FLOAT      Minimum Indel AF required for a candidate variant. Lowering the value might increase a bit of sensitivity in trade of speed and accuracy, default: ont:0.15,hifi:0.08,ilmn:0.08.'
    echo $'      --var_pct_full=FLOAT      EXPERIMENTAL: Specify an expected percentage of low quality 0/1 and 1/1 variants called in the pileup mode for full-alignment mode calling, default: 0.3.'
    echo $'      --ref_pct_full=FLOAT      EXPERIMENTAL: Specify an expected percentage of low quality 0/0 variants called in the pileup mode for full-alignment mode calling, default: 0.3 for ilmn and hifi, 0.1 for ont.'
    echo $'      --var_pct_phasing=FLOAT   EXPERIMENTAL: Specify an expected percentage of high quality 0/1 variants used in WhatsHap phasing, default: 0.8 for ont guppy5 and 0.7 for other platforms.'
    echo $'      --pileup_model_prefix=STR EXPERIMENTAL: Model prefix in pileup calling, including $prefix.data-00000-of-00002, $prefix.data-00001-of-00002 $prefix.index. default: pileup.'
    echo $'      --fa_model_prefix=STR     EXPERIMENTAL: Model prefix in full-alignment calling, including $prefix.data-00000-of-00002, $prefix.data-00001-of-00002 $prefix.index, default: full_alignment.'
    echo $'      --min_mq=INT              EXPERIMENTAL: If set, reads with mapping quality with <$min_mq are filtered, default: 5.'
    echo $'      --min_coverage=INT        EXPERIMENTAL: Minimum coverage required to call a variant, default: 2.'
    echo $'      --min_contig_size=INT     EXPERIMENTAL: If set, contigs with contig size<$min_contig_size are filtered, default: 0.'
    echo $'      --fast_mode               EXPERIMENTAL: Skip variant candidates with AF <= 0.15, default: disable.'
    echo $'      --haploid_precise         EXPERIMENTAL: Enable haploid calling mode. Only 1/1 is considered as a variant, default: disable.'
    echo $'      --haploid_sensitive       EXPERIMENTAL: Enable haploid calling mode. 0/1 and 1/1 are considered as a variant, default: disable.'
    echo $'      --no_phasing_for_fa       EXPERIMENTAL: Call variants without whatshap phasing in full alignment calling, default: disable.'
    echo $'      --call_snp_only           EXPERIMENTAL: Call candidates pass SNP minimum AF only, ignore Indel candidates, default: disable.'
    echo $'      --enable_variant_calling_at_sequence_head_and_tail
                                           EXPERIMENTAL: Enable variant calling in sequence head and tail start or end regions that flanking 16bp windows having no read support. Default: disable.'
    echo $'      --output_all_contigs_in_gvcf_header
                                           EXPERIMENTAL: Enable output all contigs in gvcf header. Default: disable.'
    echo $'      --enable_long_indel       EXPERIMENTAL: Call long Indel variants(>50 bp), default: disable.'
    echo $'      --keep_iupac_bases        EXPERIMENTAL: Keep IUPAC reference and alternate bases, default: convert all IUPAC bases to N.'
    echo $'      --base_err=FLOAT          EXPERIMENTAL: Estimated base error rate when enabling gvcf option, default: 0.001.'
    echo $'      --gq_bin_size=INT         EXPERIMENTAL: Default gq bin size for merge non-variant block when enabling gvcf option, default: 5.'

    echo $''
}

print_version()
{
    echo "Clair3 ${VERSION}"
    exit 0
}

ERROR="\\033[31m[ERROR]"
WARNING="\\033[33m[WARNING]"
NC="\\033[0m"

ARGS=`getopt -o b:f:t:m:p:o:hv \
-l bam_fn:,ref_fn:,threads:,model_path:,platform:,output:,\
bed_fn::,vcf_fn::,ctg_name::,sample_name::,qual::,samtools::,python::,pypy::,parallel::,whatshap::,chunk_num::,chunk_size::,var_pct_full::,ref_pct_full::,var_pct_phasing::,longphase::,\
min_mq::,min_coverage::,min_contig_size::,snp_min_af::,indel_min_af::,pileup_model_prefix::,fa_model_prefix::,base_err::,gq_bin_size::,fast_mode,gvcf,pileup_only,print_ref_calls,haploid_precise,haploid_sensitive,include_all_ctgs,\
use_whatshap_for_intermediate_phasing,use_longphase_for_intermediate_phasing,use_whatshap_for_final_output_phasing,use_longphase_for_final_output_phasing,use_whatshap_for_final_output_haplotagging,keep_iupac_bases,\
remove_intermediate_dir,no_phasing_for_fa,call_snp_only,enable_variant_calling_at_sequence_head_and_tail,output_all_contigs_in_gvcf_header,enable_phasing,enable_long_indel,use_gpu,longphase_for_phasing,disable_c_impl,help,version -n 'run_clair3.sh' -- "$@"`

if [ $? != 0 ] ; then echo"No input. Terminating...">&2 ; exit 1 ; fi
eval set -- "${ARGS}"

# default options
SAMPLE="SAMPLE"
BED_FILE_PATH="EMPTY"
VCF_FILE_PATH='EMPTY'
CONTIGS="EMPTY"
SAMTOOLS="samtools"
PYPY="pypy3"
PYTHON='python3'
PARALLEL='parallel'
WHATSHAP='whatshap'
LONGPHASE='EMPTY'
CHUNK_NUM=0
CHUNK_SIZE=5000000
QUAL=2
MIN_MQ=5
MIN_COV=2
BASE_ERR=0.001
GQ_BIN_SIZE=5
MIN_CONTIG_SIZE=0
PHASING_PCT="0"
PRO="0"
REF_PRO="0"
GVCF=False
PILEUP_ONLY=False
FAST_MODE=False
SHOW_REF=False
SNP_AF="0"
INDEL_AF="0"
HAP_PRE=False
HAP_SEN=False
SNP_ONLY=False
INCLUDE_ALL_CTGS=False
NO_PHASING=False
RM_TMP_DIR=False
ENABLE_LONG_INDEL=False
TMP_LP_PHASING=False
TMP_WH_PHASING=True
FINAL_LP_PHASING=False
FINAL_WH_PHASING=False
FINAL_WH_HAPLOTAG=False
KEEP_IUPAC_BASES=False
USE_GPU=False
USE_LONGPHASE=False
ENABLE_C_IMPL=True
CALL_HT=False
OUTPUT_ALL_CONTIGS=False
PILEUP_PREFIX="pileup"
FA_PREFIX="full_alignment"

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
    --longphase ) LONGPHASE="$2"; shift 2 ;;
    --var_pct_full ) PRO="$2"; shift 2 ;;
    --ref_pct_full ) REF_PRO="$2"; shift 2 ;;
    --var_pct_phasing ) PHASING_PCT="$2"; shift 2 ;;
    --snp_min_af ) SNP_AF="$2"; shift 2 ;;
    --indel_min_af ) INDEL_AF="$2"; shift 2 ;;
    --min_mq ) MIN_MQ="$2"; shift 2 ;;
    --min_coverage ) MIN_COV="$2"; shift 2 ;;
    --min_contig_size ) MIN_CONTIG_SIZE="$2"; shift 2 ;;
    --pileup_model_prefix ) PILEUP_PREFIX="$2"; shift 2 ;;
    --fa_model_prefix ) FA_PREFIX="$2"; shift 2 ;;
    --base_err ) BASE_ERR="$2"; shift 2 ;;
    --gq_bin_size ) GQ_BIN_SIZE="$2"; shift 2 ;;
    --gvcf ) GVCF=True; shift 1 ;;
    --pileup_only ) PILEUP_ONLY=True; shift 1 ;;
    --fast_mode ) FAST_MODE=True; shift 1 ;;
    --call_snp_only ) SNP_ONLY=True; shift 1 ;;
    --enable_variant_calling_at_sequence_head_and_tail ) CALL_HT=True; shift 1 ;;
    --output_all_contigs_in_gvcf_header ) OUTPUT_ALL_CONTIGS=True; shift 1 ;;
    --print_ref_calls ) SHOW_REF=True; shift 1 ;;
    --haploid_precise ) HAP_PRE=True; shift 1 ;;
    --haploid_sensitive ) HAP_SEN=True; shift 1 ;;
    --include_all_ctgs ) INCLUDE_ALL_CTGS=True; shift 1 ;;
    --no_phasing_for_fa ) NO_PHASING=True; shift 1 ;;
    --remove_intermediate_dir ) RM_TMP_DIR=True; shift 1 ;;
    --use_whatshap_for_intermediate_phasing ) TMP_WH_PHASING=True; shift 1 ;;
    --use_longphase_for_intermediate_phasing ) USE_LONGPHASE=True; shift 1 ;;
    --use_whatshap_for_final_output_phasing ) FINAL_WH_PHASING=True; shift 1 ;;
    --use_longphase_for_final_output_phasing ) FINAL_LP_PHASING=True; shift 1 ;;
    --use_whatshap_for_final_output_haplotagging ) FINAL_WH_HAPLOTAG=True; shift 1 ;;
    --enable_phasing ) FINAL_WH_PHASING=True; shift 1 ;;
    --enable_long_indel ) ENABLE_LONG_INDEL=True; shift 1 ;;
    --keep_iupac_bases ) KEEP_IUPAC_BASES=True; shift 1 ;;
    --use_gpu ) USE_GPU=True; shift 1 ;;
    --longphase_for_phasing ) USE_LONGPHASE=True; shift 1 ;;
    --disable_c_impl ) ENABLE_C_IMPL=False; shift 1 ;;

    -- ) shift; break; ;;
    -h|--help ) print_help_messages; exit 0 ;;
    -v|--version ) print_version; exit 0 ;;
    * ) print_help_messages; break ;;
   esac
done

if [ -z ${BAM_FILE_PATH} ] || [ -z ${REFERENCE_FILE_PATH} ] || [ -z ${THREADS} ] || [ -z ${OUTPUT_FOLDER} ] || [ -z ${PLATFORM} ] || [ -z ${MODEL_PATH} ]; then
      if [ -z ${BAM_FILE_PATH} ] && [ -z ${REFERENCE_FILE_PATH} ] && [ -z ${THREADS} ] && [ -z ${OUTPUT_FOLDER} ] && [ -z ${PLATFORM} ] && [ -z ${MODEL_PATH} ]; then print_help_messages; exit 0; fi
      if [ -z ${BAM_FILE_PATH} ]; then echo -e "${ERROR} Require to define index BAM input by --bam_fn=BAM${NC}"; fi
      if [ -z ${REFERENCE_FILE_PATH} ]; then echo -e "${ERROR} Require to define FASTA reference file input by --ref_fn=REF${NC}"; fi
      if [ -z ${THREADS} ]; then echo -e "${ERROR} Require to define max threads to be used by --threads=THREADS${NC}"; fi
      if [ -z ${OUTPUT_FOLDER} ]; then echo -e "${ERROR} Require to define output folder by --output=OUTPUT_DIR${NC}"; fi
      if [ -z ${PLATFORM} ]; then echo -e "${ERROR} Require to define platform by --platform={ont,hifi,ilmn}${NC}"; fi
      if [ -z ${MODEL_PATH} ]; then echo -e "${ERROR} Require to define model path by --model_path=MODEL_PREFIX${NC}"; fi
      exit 1;
fi

# force to use absolute path when in docker or singularity environment
if [ `pwd` = "/opt/bin" ]; then
    if [[ ! "${BAM_FILE_PATH}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --bam_fn=FILE${NC}"; exit 1; fi
    if [[ ! "${REFERENCE_FILE_PATH}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --ref_fn=FILE${NC}"; exit 1; fi
    if [[ ! "${MODEL_PATH}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --model_path=PATH${NC}"; exit 1; fi
    if [[ ! "${OUTPUT_FOLDER}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --output=PATH${NC}"; exit 1; fi
    if [ "${BED_FILE_PATH}" != "EMPTY" ] &&  [ ! -z ${BED_FILE_PATH} ] && [[ ! "${BED_FILE_PATH}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --bef_fn=FILE${NC}"; exit 1; fi
    if [ "${VCF_FILE_PATH}" != "EMPTY" ] &&  [ ! -z ${VCF_FILE_PATH} ] && [[ ! "${VCF_FILE_PATH}" = /* ]]; then echo -e "${ERROR} Require to use absolute file path --vcf_fn=FILE${NC}"; exit 1; fi
fi

# relative path support
if [[ ! "${BAM_FILE_PATH}" = /* ]] && [ -f ${BAM_FILE_PATH} ]; then BAM_FILE_PATH=`pwd`/${BAM_FILE_PATH}; fi
if [[ ! "${REFERENCE_FILE_PATH}" = /* ]] && [ -f ${REFERENCE_FILE_PATH} ]; then REFERENCE_FILE_PATH=`pwd`/${REFERENCE_FILE_PATH}; fi
if [[ ! "${MODEL_PATH}" = /* ]] && [ -d ${MODEL_PATH} ]; then MODEL_PATH=`pwd`/${MODEL_PATH}; fi
if [ "${BED_FILE_PATH}" != "EMPTY" ] && [ ! -z ${BED_FILE_PATH} ] && [[ ! "${BED_FILE_PATH}" = /* ]] && [ -f ${BED_FILE_PATH} ]; then BED_FILE_PATH=`pwd`/${BED_FILE_PATH}; fi
if [ "${VCF_FILE_PATH}" != "EMPTY" ] && [ ! -z ${VCF_FILE_PATH} ] && [[ ! "${VCF_FILE_PATH}" = /* ]] && [ -f ${VCF_FILE_PATH} ]; then VCF_FILE_PATH=`pwd`/${VCF_FILE_PATH}; fi
if [[ ! "${OUTPUT_FOLDER}" = /* ]]; then echo -e "${WARNING} No absolute output path provided, using current directory as prefix${NC}"; OUTPUT_FOLDER=`pwd`/${OUTPUT_FOLDER}; fi

mkdir -p ${OUTPUT_FOLDER}
if [ ! -d ${OUTPUT_FOLDER} ]; then echo -e "${ERROR} Cannot create output folder ${OUTPUT_FOLDER}${NC}"; exit 1; fi

# show default reference proportion 0.3 for ilmn and hifi, 0.1 for ont
if [ "${PLATFORM}" = "ont" ] && [ "${REF_PRO}" = "0" ]; then REF_PRO=0.1; fi
if [ "${PLATFORM}" != "ont" ] && [ "${REF_PRO}" = "0" ]; then REF_PRO=0.3; fi

# show default variant proportion 0.3 for ilmn and hifi, 0.7 for ont
if [ "${PLATFORM}" = "ont" ] && [ "${PRO}" = "0" ]; then PRO=0.7; fi
if [ "${PLATFORM}" != "ont" ] && [ "${PRO}" = "0" ]; then PRO=0.3; fi

# set default af for ilmn and hifi and ont
if [ "${SNP_AF}" = "0" ]; then SNP_AF=0.08; fi
if [ "${PLATFORM}" = "ont" ] && [ "${INDEL_AF}" = "0" ]; then INDEL_AF=0.15; fi
if [ "${PLATFORM}" != "ont" ] && [ "${INDEL_AF}" = "0" ]; then INDEL_AF=0.08; fi

# show default high quality hete variant proportion for whatshap phasing, 0.8 for ont guppy5 and 0.7 for others
if [ "${PHASING_PCT}" = "0" ]; then PHASING_PCT=0.7; fi
BASE_MODEL=$(basename ${MODEL_PATH})
if [ "${BASE_MODEL}" = "r941_prom_sup_g5014" ] || [ "${BASE_MODEL}" = "r941_prom_hac_g5014" ] || [ "${BASE_MODEL}" = "ont_guppy5" ]; then PHASING_PCT=0.8; fi

# use the default longphase binary path
if ([ "${USE_LONGPHASE}" == True ] || [ "${FINAL_LP_PHASING}" == True ]) && [ "${LONGPHASE}" == "EMPTY" ]; then LONGPHASE="${SCRIPT_PATH}/longphase"; fi
if ([ "${USE_LONGPHASE}" == True ] || [ "${FINAL_LP_PHASING}" == True ]) && [ ! -f ${LONGPHASE} ]; then echo -e "${ERROR} Cannot find LongPhase path in ${LONGPHASE}, exit!${NC}"; exit 1; fi
if ([ "${USE_LONGPHASE}" == True ] || [ "${FINAL_LP_PHASING}" == True ]) && [ "${PLATFORM}" = "ilmn" ]; then echo -e "${WARNING} Illumina platform do not support longphase phasing, will enable whatshap phasing! ${NC}";  USE_LONGPHASE=False; fi

if [ "${FINAL_WH_HAPLOTAG}" == True ] && [ "${FINAL_WH_PHASING}" == False ] && [ "${FINAL_LP_PHASING}" == False ]; then FINAL_WH_PHASING=True; fi

# remove the last '/' character in directory input
OUTPUT_FOLDER=$(echo ${OUTPUT_FOLDER%*/})
MODEL_PATH=$(echo ${MODEL_PATH%*/})

# optional parameters should use "="
(time (
echo "[INFO] CLAIR3 VERSION: ${VERSION}"
echo "[INFO] BAM FILE PATH: ${BAM_FILE_PATH}"
echo "[INFO] REFERENCE FILE PATH: ${REFERENCE_FILE_PATH}"
echo "[INFO] MODEL PATH: ${MODEL_PATH}"
echo "[INFO] OUTPUT FOLDER: ${OUTPUT_FOLDER}"
echo "[INFO] PLATFORM: ${PLATFORM}"
echo "[INFO] THREADS: ${THREADS}"
echo "[INFO] BED FILE PATH: ${BED_FILE_PATH}"
echo "[INFO] VCF FILE PATH: ${VCF_FILE_PATH}"
echo "[INFO] CONTIGS: ${CONTIGS}"
echo "[INFO] CONDA PREFIX: ${CONDA_PREFIX}"
echo "[INFO] SAMTOOLS PATH: ${SAMTOOLS}"
echo "[INFO] PYTHON PATH: ${PYTHON}"
echo "[INFO] PYPY PATH: ${PYPY}"
echo "[INFO] PARALLEL PATH: ${PARALLEL}"
echo "[INFO] WHATSHAP PATH: ${WHATSHAP}"
echo "[INFO] LONGPHASE PATH: ${LONGPHASE}"
echo "[INFO] CHUNK SIZE: ${CHUNK_SIZE}"
if [ ${CHUNK_NUM} -gt 0 ]; then echo "[INFO] CHUNK NUM: ${CHUNK_NUM}"; fi
if [ ${MIN_CONTIG_SIZE} -gt 0 ]; then echo "[INFO] MIN CONTIG SIZE: ${CHUNK_NUM}"; fi
echo "[INFO] FULL ALIGN PROPORTION: ${PRO}"
echo "[INFO] FULL ALIGN REFERENCE PROPORTION: ${REF_PRO}"
echo "[INFO] PHASING PROPORTION: ${PHASING_PCT}"
echo "[INFO] MINIMUM MQ: ${MIN_MQ}"
echo "[INFO] MINIMUM COVERAGE: ${MIN_COV}"
echo "[INFO] SNP AF THRESHOLD: ${SNP_AF}"
echo "[INFO] INDEL AF THRESHOLD: ${INDEL_AF}"
echo "[INFO] BASE ERROR IN GVCF: ${BASE_ERR}"
echo "[INFO] GQ BIN SIZE IN GVCF: ${GQ_BIN_SIZE}"
echo "[INFO] ENABLE FILEUP ONLY CALLING: ${PILEUP_ONLY}"
echo "[INFO] ENABLE FAST MODE CALLING: ${FAST_MODE}"
echo "[INFO] ENABLE CALLING SNP CANDIDATES ONLY: ${SNP_ONLY}"
echo "[INFO] ENABLE PRINTING REFERENCE CALLS: ${SHOW_REF}"
echo "[INFO] ENABLE OUTPUT GVCF: ${GVCF}"
echo "[INFO] ENABLE HAPLOID PRECISE MODE: ${HAP_PRE}"
echo "[INFO] ENABLE HAPLOID SENSITIVE MODE: ${HAP_SEN}"
echo "[INFO] ENABLE INCLUDE ALL CTGS CALLING: ${INCLUDE_ALL_CTGS}"
echo "[INFO] ENABLE NO PHASING FOR FULL ALIGNMENT: ${NO_PHASING}"
echo "[INFO] ENABLE REMOVING INTERMEDIATE FILES: ${RM_TMP_DIR}"
echo "[INFO] ENABLE LONGPHASE FOR INTERMEDIATE VCF PHASING: ${USE_LONGPHASE}"
echo "[INFO] ENABLE PHASING FINAL VCF OUTPUT USING WHATSHAP: ${FINAL_WH_PHASING}"
echo "[INFO] ENABLE PHASING FINAL VCF OUTPUT USING LONGPHASE: ${FINAL_LP_PHASING}"
echo "[INFO] ENABLE HAPLOTAGGING FINAL BAM: ${FINAL_WH_HAPLOTAG}"
echo "[INFO] ENABLE LONG INDEL CALLING: ${ENABLE_LONG_INDEL}"
echo "[INFO] ENABLE C_IMPLEMENT: ${ENABLE_C_IMPL}"
echo $''

# file check
if [ ! -f ${BAM_FILE_PATH} ]; then echo -e "${ERROR} BAM file ${BAM_FILE_PATH} not found${NC}"; exit 1; fi
if [ ! -f ${BAM_FILE_PATH}.bai ] && [ ! -f ${BAM_FILE_PATH%.*}.bai ] && [ ! -f ${BAM_FILE_PATH}.csi ] && [ ! -f ${BAM_FILE_PATH%.*}.csi ] \
&& [ ! -f ${BAM_FILE_PATH}.crai ] && [ ! -f ${BAM_FILE_PATH%.*}.crai ]; then echo -e "${ERROR} BAM index bai file not found, please use 'samtools index \$BAM' first${NC}"; exit 1; fi
if [ ! -f ${REFERENCE_FILE_PATH} ]; then echo -e "${ERROR} Reference file ${REFERENCE_FILE_PATH} not found${NC}"; exit 1; fi
if [ ! -f ${REFERENCE_FILE_PATH}.fai ] && [ ! -f ${REFERENCE_FILE_PATH%.*}.fai ]; then echo -e "${ERROR} Reference index fai file not found, please use 'samtools faidx \$REF' first${NC}"; exit 1; fi

if [ "${BED_FILE_PATH}" != "EMPTY" ] && [ ! -z ${BED_FILE_PATH} ] && [ ! -f ${BED_FILE_PATH} ]; then echo -e "${ERROR} BED file ${BED_FILE_PATH} provides but not found${NC}"; exit 1; fi
if [ "${VCF_FILE_PATH}" != "EMPTY" ] && [ ! -z ${VCF_FILE_PATH} ] && [ ! -f ${VCF_FILE_PATH} ]; then echo -e "${ERROR} VCF file ${VCF_FILE_PATH} provides but not found${NC}"; exit 1; fi
if [ ! -d ${MODEL_PATH} ] && [ -z ${CONDA_PREFIX} ]; then echo -e "${ERROR} Conda prefix not found, please activate clair3 conda environment first, model path: ${MODEL_PATH}${NC}"; exit 1; fi
if [ ! -d ${MODEL_PATH} ]; then echo -e "${ERROR} Model path not found${NC}"; exit 1; fi

# max threads detection
if [ "$(uname)" = "Darwin" ]; then MAX_THREADS=$(sysctl -n hw.logicalcpu); else MAX_THREADS=$(nproc); fi
if [ "$(uname)" = "Darwin" ]; then SHELL_ENTRY=${SHELL}; else SHELL_ENTRY=""; fi
if [[ ! ${THREADS} =~ ^[\-0-9]+$ ]] || (( ${THREADS} <= 0)); then echo -e "${ERROR} Invalid threads input --threads=INT ${NC}"; exit 1; fi
if [[ ${THREADS} -gt ${MAX_THREADS} ]]; then echo -e "${WARNING} Threads setting exceeds maximum available threads ${MAX_THREADS}, set threads=${MAX_THREADS}${NC}"; THREADS=${MAX_THREADS}; fi

# max user ulimit threads detection
MAX_ULIMIT_THREADS=`ulimit -u`
if [ ! -z ${MAX_ULIMIT_THREADS} ]; then PER_ULIMIT_THREADS=$((${MAX_ULIMIT_THREADS}/30)); else MAX_ULIMIT_THREADS="unlimited"; PER_ULIMIT_THREADS=${THREADS}; fi
if [[ ${PER_ULIMIT_THREADS} < 1 ]]; then PER_ULIMIT_THREADS=1; fi
if [ "${MAX_ULIMIT_THREADS}" != "unlimited" ] && [[ ${THREADS} -gt ${PER_ULIMIT_THREADS} ]]; then echo -e "${WARNING} Threads setting exceeds maximum ulimit threads ${THREADS} * 30 > ${MAX_ULIMIT_THREADS} (ulimit -u), set threads=${PER_ULIMIT_THREADS}${NC}"; THREADS=${PER_ULIMIT_THREADS}; fi


# platform check
if [ ! ${PLATFORM} = "ont" ] && [ ! ${PLATFORM} = "hifi" ] && [ ! ${PLATFORM} = "ilmn" ]; then echo -e "${ERROR} Invalid platform input, optional: {ont, hifi, ilmn}${NC}"; exit 1; fi

# optional parameter detection
if [ -z ${BED_FILE_PATH} ]; then echo -e "${ERROR} Use '--bed_fn=FILE' instead of '--bed_fn FILE' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${VCF_FILE_PATH} ]; then echo -e "${ERROR} Use '--vcf_fn=FILE' instead of '--vcf_fn =FILE' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${CONTIGS} ]; then echo -e "${ERROR} Use '--ctg_name=STR' instead of '--ctg_name STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${SAMPLE} ]; then echo -e "${ERROR} Use '--sample_name=STR' instead of '--sample_name STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${QUAL} ]; then echo -e "${ERROR} Use '--qual=INT' instead of '--qual INT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${SAMTOOLS} ]; then echo -e "${ERROR} Use '--samtools=STR' instead of '--samtools STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PYTHON} ]; then echo -e "${ERROR} Use '--python=STR' instead of '--python STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PYPY} ]; then echo -e "${ERROR} Use '--pypy=STR' instead of '--pypy STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PARALLEL} ]; then echo -e "${ERROR} Use '--parallel=STR' instead of '--parallel STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${WHATSHAP} ]; then echo -e "${ERROR} Use '--whatshap=STR' instead of '--whatshap STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${CHUNK_SIZE} ]; then echo -e "${ERROR} Use '--chunk_size=INT' instead of '--chunk_size INT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${SNP_AF} ]; then echo -e "${ERROR} Use '--snp_min_af=FLOAT' instead of '--snp_min_af FLOAT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${INDEL_AF} ]; then echo -e "${ERROR} Use '--indel_min_af=FLOAT' instead of '--indel_min_af FLOAT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PRO} ]; then echo -e "${ERROR} Use '--var_pct_full=FLOAT' instead of '--var_pct_full FLOAT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${REF_PRO} ]; then echo -e "${ERROR} Use '--ref_pct_full=FLOAT' instead of '--ref_pct_full FLOAT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PHASING_PCT} ]; then echo -e "${ERROR} Use '--var_pct_phasing=FLOAT' instead of '--var_pct_phasing FLOAT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${PILEUP_PREFIX} ]; then echo -e "${ERROR} Use '--pileup_model_prefix=STR' instead of '--pileup_model_prefix STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${FA_PREFIX} ]; then echo -e "${ERROR} Use '--fa_model_prefix=STR' instead of '--fa_model_prefix STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${MIN_MQ} ]; then echo -e "${ERROR} Use '--min_mq=INT' instead of '--min_mq INT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${MIN_COV} ]; then echo -e "${ERROR} Use '--min_coverage=INT' instead of '--min_coverage INT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${MIN_CONTIG_SIZE} ]; then echo -e "${ERROR} Use '--min_contig_size=INT' instead of '--min_contig_size INT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${LONGPHASE} ]; then echo -e "${ERROR} Use '--longphase=STR' instead of '--longphase STR' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${BASE_ERR} ]; then echo -e "${ERROR} Use '--base_err=FLOAT' instead of '--base_err FLOAT' for optional parameters${NC}"; exit 1 ; fi
if [ -z ${GQ_BIN_SIZE} ]; then echo -e "${ERROR} Use '--gq_bin_size=INT' instead of '--gq_bin_size INT' for optional parameters${NC}"; exit 1 ; fi

# min mapping quality, min coverage and min contig size detection
if [[ ! ${THREADS} =~ ^[\-0-9]+$ ]] || (( ${THREADS} <= 0)); then echo -e "${ERROR} Invalid threads input --threads=INT ${NC}"; exit 1; fi
if [[ ! ${MIN_MQ} =~ ^[\-0-9]+$ ]] || (( ${MIN_MQ} < 5)); then echo -e "${WARNING} Invalid minimum mapping quality input --min_mq>=5 ${NC}"; MIN_MQ=5; fi
if [[ ! ${MIN_COV} =~ ^[\-0-9]+$ ]] || (( ${MIN_COV} < 2)); then echo -e "${WARNING} Invalid minimum coverage input --min_coverage>=2 ${NC}"; MIN_COV=2; fi
if [[ ! ${MIN_CONTIG_SIZE} =~ ^[\-0-9]+$ ]] || (( ${MIN_CONTIG_SIZE} < 0)); then echo -e "${WARNING} Invalid minimum contig size --min_contig_size>=0 ${NC}"; MIN_CONTIG_SIZE=0; fi
if [[ ! ${GQ_BIN_SIZE} =~ ^[\-0-9]+$ ]] || (( ${GQ_BIN_SIZE} < 0)); then echo -e "${WARNING} Invalid gq bin size --gq_bin_size>=0 ${NC}"; MIN_CONTIG_SIZE=0; fi

# in genotyping mode, set --snp_min_af and --indel_min_af to 0.0
if [ "${VCF_FILE_PATH}" != "EMPTY" ] && [ ! -z ${VCF_FILE_PATH} ] && [ -f ${VCF_FILE_PATH} ]; then SNP_AF=0.0; INDEL_AF=0.0; fi

# model prefix detection
if [ ! -f ${MODEL_PATH}/${PILEUP_PREFIX}.index ]; then echo -e "${ERROR} No pileup model found in provided model path and model prefix ${MODEL_PATH}/${PILEUP_PREFIX} ${NC}"; exit 1; fi
if [ ! -f ${MODEL_PATH}/${FA_PREFIX}.index ]; then echo -e "${ERROR} No full-alignment model found in provided model path and model prefix ${MODEL_PATH}/${FA_PREFIX} ${NC}"; exit 1; fi

CLAIR3_SCRIPT="clair3.sh"
if [ "${ENABLE_C_IMPL}" == True ] && [ "${PLATFORM}" = "ilmn" ]; then echo -e "${WARNING} Illumina platform will disable C implement to support short read realignment process! ${NC}";  ENABLE_C_IMPL=False; fi
if [ "${ENABLE_C_IMPL}" == True ]; then CLAIR3_SCRIPT="clair3_c_impl.sh"; fi

# keep command line info., for vcf header
mkdir -p ${OUTPUT_FOLDER}/tmp
echo "$CMD" > "${OUTPUT_FOLDER}/tmp/CMD"

set -x
${SHELL_ENTRY} ${SCRIPT_PATH}/scripts/${CLAIR3_SCRIPT} \
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
    --var_pct_phasing=${PHASING_PCT} \
    --snp_min_af=${SNP_AF} \
    --indel_min_af=${INDEL_AF} \
    --min_mq=${MIN_MQ} \
    --min_coverage=${MIN_COV} \
    --min_contig_size=${MIN_CONTIG_SIZE} \
    --pileup_only=${PILEUP_ONLY} \
    --gvcf=${GVCF} \
    --base_err=${BASE_ERR} \
    --gq_bin_size=${GQ_BIN_SIZE} \
    --fast_mode=${FAST_MODE} \
    --call_snp_only=${SNP_ONLY} \
    --enable_variant_calling_at_sequence_head_and_tail=${CALL_HT} \
    --output_all_contigs_in_gvcf_header=${OUTPUT_ALL_CONTIGS} \
    --print_ref_calls=${SHOW_REF} \
    --haploid_precise=${HAP_PRE} \
    --haploid_sensitive=${HAP_SEN} \
    --include_all_ctgs=${INCLUDE_ALL_CTGS} \
    --no_phasing_for_fa=${NO_PHASING} \
    --pileup_model_prefix=${PILEUP_PREFIX} \
    --fa_model_prefix=${FA_PREFIX} \
    --remove_intermediate_dir=${RM_TMP_DIR} \
    --enable_phasing=${FINAL_WH_PHASING} \
    --enable_long_indel=${ENABLE_LONG_INDEL} \
    --keep_iupac_bases=${KEEP_IUPAC_BASES} \
    --use_gpu=${USE_GPU} \
    --longphase_for_phasing=${USE_LONGPHASE} \
    --longphase=${LONGPHASE} \
    --use_whatshap_for_intermediate_phasing=${TMP_WH_PHASING} \
    --use_longphase_for_intermediate_phasing=${USE_LONGPHASE} \
    --use_whatshap_for_final_output_phasing=${FINAL_WH_PHASING} \
    --use_longphase_for_final_output_phasing=${FINAL_LP_PHASING} \
    --use_whatshap_for_final_output_haplotagging=${FINAL_WH_HAPLOTAG}

)) 2>&1 | tee ${OUTPUT_FOLDER}/run_clair3.log
