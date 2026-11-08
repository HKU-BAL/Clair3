#!/bin/bash
# =============================================================================
# Clair3 训练端到端参考脚本（chr21）
#
#
#
#
#   耗时（整条 chr21，2 个深度，1 epoch）: 约 40–60 分钟
#
# =============================================================================
#
# -----------------------------------------------------------------------------
# 关键路径：
#
# 【参考基因组 GRCh38】
#   REF          /autofs/bal19/zxzheng/testData/hifi/data/GRCh38_no_alt_analysis_set.fasta
#                （同名的 .fai 在同目录下）
#
# 【真值集 GIAB HG002 v5.0q】
#   TRUTH_VCF    /autofs/bal32/mghe/data/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz
#   TRUTH_VCF索引 /autofs/bal32/mghe/data/v5.0q/HG002_GRCh38_v5.0q_smvar.vcf.gz.tbi
#   TRUTH_BED    /autofs/bal32/mghe/data/v5.0q/HG002_GRCh38_v5.0q_smvar.benchmark.bed
#   CHM13 版同目录：HG002_CHM13v2.0_v5.0q_smvar.vcf.gz / .benchmark.bed
#
# 【BAM（PacBio Revio HiFi, HG002, GRCh38, phased）】
#   BAM 目录     /autofs/bal19/zxzheng/pb/revio_training/ru/hg002/phased_bam/
#   命名规则     hg002_revio_<depth>_<chr>.bam  (+ .bai)
#                深度depth ∈ {1000,778,667,556,444,333,222,184,147,111}；chr ∈ 1–19,21,22
#                （没有 chr20；1000 档是指向 hg002_revio_<chr>.bam 的符号链接）
#   本脚本示例使用   /autofs/bal19/zxzheng/pb/revio_training/ru/hg002/phased_bam/hg002_revio_1000_21.bam
#                即 chr21 全量 PacBio Revio HiFi，约 90x
#
# 【Clair3 代码（被测试的对象，PyTorch 版）】
#   REPO         /autofs/bal32/mghe/ClairS_issue/Clair3
#   入口         <REPO>/clair3.py          （预处理小工具 + Train）
#                <REPO>/run_clair3.py      （--pileup_only 产生候选）
#
# 【预训练权重（fine-tune 起点，这里用的是一个hifi_revio 模型）】
#   MODEL_SRC    /autofs/bal32/mghe/ClairS_issue/clair3_models_pytorch_all/hifi_revio/
#   Pileup       <MODEL_SRC>/pileup.pt
#   Full-aln     <MODEL_SRC>/full_alignment.pt
#
# 【Python / 工具环境】
#   pypy3        /autofs/bal33/zxzheng/env/conda/envs/clair3_1.1.1/bin/pypy3
#   python(pileup)/autofs/bal33/zxzheng/env/conda/envs/clair3_v2/bin/python3
#   python(train) /autofs/bal19/zxzheng/env/conda/envs/mamba/envs/somatic/bin/python3
#   samtools     /autofs/bal33/zxzheng/env/conda/envs/clair3_1.1.1/bin/samtools
#   bgzip        /autofs/bal33/zxzheng/env/conda/envs/clair3_1.1.1/bin/bgzip
#   tabix        /autofs/bal33/zxzheng/env/conda/envs/clair3_1.1.1/bin/tabix
#   parallel     /autofs/bal33/zxzheng/env/conda/envs/clair3_r6/bin/parallel
#   whatshap     /autofs/bal33/zxzheng/env/miniconda2/envs/clair2/bin/whatshap
#
# 【输出目录】
#   WORK         <本脚本所在目录>/chr21_quick
#                bam/ truth/ ru/ split/ var/ bins/{pileup,fa} candidates/ models/
# =============================================================================
set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

# -----------------------------------------------------------------------------
# 1. 配置
# -----------------------------------------------------------------------------
SAMPLE=${SAMPLE:-HG002}
CTG=chr21
CTG_NUM=${CTG#chr}          # SelectHetSnp 输出文件名会去掉 chr 前缀
REGION=${REGION:-$CTG}          # 默认整条染色体；也可设成 chr21:13000000-14000000

REF=${REF:?Set REF to an indexed reference FASTA}
TRUTH_VCF=${TRUTH_VCF:?Set TRUTH_VCF to an indexed truth VCF}
TRUTH_BED=${TRUTH_BED:?Set TRUTH_BED to a truth BED}
BAM_RAW=${BAM_RAW:?Set BAM_RAW to an indexed chr21 BAM}

# 深度档：1000=100%，222=22.2%。
# 加档位：例如 DEPTHS=(1000 444 222 111)， 脚本会从 1000 档降采样生成。
DEPTHS=(${DEPTHS_OVERRIDE:-1000 222})
VALIDATION_DEPTH=${VALIDATION_DEPTH:-222}

CHUNKS=${CHUNKS:-8}             # 每条染色体的分块数
THREADS=${THREADS:-16}
EPOCHS=${EPOCHS:-1}             # 冒烟测试 1 个 epoch 即可
BATCH_SIZE=${BATCH_SIZE:-1000}
CHUNK_SIZE=${CHUNK_SIZE:-200}
PLATFORM=${PLATFORM:-hifi}
DO_RU=${DO_RU:-1}               # 1 = 跑 Repun 表示统一；0 = 直接用真值 VCF

# -----------------------------------------------------------------------------
# 2. 环境与工具
# -----------------------------------------------------------------------------
CLAIR3_REPO=${CLAIR3_REPO:-$(cd "$SCRIPT_DIR/.." && pwd)}
CLAIR3=$CLAIR3_REPO/clair3.py
RUN_CLAIR3=$CLAIR3_REPO/run_clair3.py

# clair3.py accepts direct shell, Python, and PyPy invocation on this branch.
PYPY3=${PYPY3:-python3}
PY_BIN=${PY_BIN:-python3}
PY_TRAIN=${PY_TRAIN:-python3}
SAMTOOLS=${SAMTOOLS:-samtools}
BGZIP=${BGZIP:-bgzip}
TABIX=${TABIX:-tabix}
PARALLEL=${PARALLEL:-parallel}
WHATSHAP=${WHATSHAP:-whatshap}

# 预训练权重：fine-tune 的起点（公开 hifi_revio 模型）
MODEL_SRC=${MODEL_SRC:?Set MODEL_SRC to a directory containing pileup.pt and full_alignment.pt}
PRETRAIN_PILEUP=$MODEL_SRC/pileup.pt
PRETRAIN_FA=$MODEL_SRC/full_alignment.pt

WORK=${WORK:-$SCRIPT_DIR/chr21_quick}
CHUNK_LIST=$(seq 1 "$CHUNKS")

log() { echo "[$(date '+%F %T')] $*"; }
die() { echo "[ERROR] $*" >&2; exit 1; }
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}

QUICK_WORK="$WORK" \
QUICK_SCRIPT="$SCRIPT_DIR/$(basename "${BASH_SOURCE[0]}")" \
QUICK_CLAIR3="$CLAIR3" \
QUICK_RUN_CLAIR3="$RUN_CLAIR3" \
QUICK_SAMPLE="$SAMPLE" \
QUICK_CTG="$CTG" \
QUICK_REGION="$REGION" \
QUICK_DEPTHS="${DEPTHS[*]}" \
QUICK_VALIDATION_DEPTH="$VALIDATION_DEPTH" \
QUICK_CHUNKS="$CHUNKS" \
QUICK_PLATFORM="$PLATFORM" \
QUICK_DO_RU="$DO_RU" \
QUICK_REF="$REF" \
QUICK_TRUTH_VCF="$TRUTH_VCF" \
QUICK_TRUTH_BED="$TRUTH_BED" \
QUICK_BAM_RAW="$BAM_RAW" \
QUICK_PRETRAIN_PILEUP="$PRETRAIN_PILEUP" \
QUICK_PRETRAIN_FA="$PRETRAIN_FA" \
QUICK_PYPY3="$PYPY3" \
QUICK_PY_BIN="$PY_BIN" \
QUICK_PY_TRAIN="$PY_TRAIN" \
QUICK_SAMTOOLS="$SAMTOOLS" \
QUICK_BGZIP="$BGZIP" \
QUICK_TABIX="$TABIX" \
QUICK_PARALLEL="$PARALLEL" \
QUICK_WHATSHAP="$WHATSHAP" \
"$PY_TRAIN" - <<'PY' || die "缓存兼容性检查失败"
import hashlib
import json
import os
from pathlib import Path
import shutil


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(value):
    path = Path(value).expanduser().resolve(strict=True)
    return {"path": str(path), "sha256": sha256(path)}


def tool_identity(value):
    resolved = shutil.which(os.path.expanduser(value))
    identity = {"command": value, "path": resolved}
    if resolved:
        path = Path(resolved).resolve()
        identity["path"] = str(path)
        if path.is_file():
            identity["sha256"] = sha256(path)
    return identity


def changed_keys(before, after, prefix=""):
    if isinstance(before, dict) and isinstance(after, dict):
        return [
            key
            for name in sorted(set(before) | set(after))
            for key in changed_keys(
                before.get(name), after.get(name), f"{prefix}.{name}" if prefix else name
            )
        ]
    return [prefix] if before != after else []


work = Path(os.environ["QUICK_WORK"]).expanduser()
manifest_path = work / ".quick-training-manifest.json"
try:
    current = {
        "format": 1,
        "config": {
            "sample": os.environ["QUICK_SAMPLE"],
            "contig": os.environ["QUICK_CTG"],
            "region": os.environ["QUICK_REGION"],
            "depths": os.environ["QUICK_DEPTHS"].split(),
            "validation_depth": os.environ["QUICK_VALIDATION_DEPTH"],
            "chunks": os.environ["QUICK_CHUNKS"],
            "platform": os.environ["QUICK_PLATFORM"],
            "do_ru": os.environ["QUICK_DO_RU"],
        },
        "inputs": {
            name: file_identity(os.environ[f"QUICK_{name.upper()}"])
            for name in [
                "ref",
                "truth_vcf",
                "truth_bed",
                "bam_raw",
                "pretrain_pileup",
                "pretrain_fa",
            ]
        },
        "code": {
            name: file_identity(os.environ[f"QUICK_{name.upper()}"])
            for name in ["script", "clair3", "run_clair3"]
        },
        "tools": {
            name: tool_identity(os.environ[f"QUICK_{name.upper()}"])
            for name in [
                "pypy3",
                "py_bin",
                "py_train",
                "samtools",
                "bgzip",
                "tabix",
                "parallel",
                "whatshap",
            ]
        },
    }
except OSError as error:
    raise SystemExit(f"[ERROR] cannot fingerprint quick-workflow inputs/tools: {error}")

if manifest_path.exists():
    try:
        previous = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as error:
        raise SystemExit(
            f"[ERROR] invalid quick-workflow manifest {manifest_path}: {error}; "
            "use a new, empty WORK; existing data was not modified"
        )
    changed = changed_keys(previous, current)
    if changed:
        raise SystemExit(
            f"[ERROR] incompatible cached WORK {work}; changed identity: {', '.join(changed)}. "
            "Use a new, empty WORK; existing data was not modified"
        )

    depths = current["config"]["depths"]
    validation_depth = current["config"]["validation_depth"]
    sample = current["config"]["sample"]
    for kind in ["pileup", "fa"]:
        for split in ["train", "validation"]:
            expected = {
                f"bin_{sample}_{depth}"
                for depth in depths
                if (split == "validation") == (depth == validation_depth)
            }
            directory = work / "bins" / kind / split
            if not directory.is_dir():
                continue
            unexpected = sorted(
                path.name
                for path in directory.iterdir()
                if path.is_file() and path.name not in expected
            )
            if unexpected:
                raise SystemExit(
                    f"[ERROR] unexpected cached archive(s) in {directory}: "
                    f"{', '.join(unexpected)}; configured depths are {', '.join(depths)}. "
                    "Use a new, empty WORK; existing data was not modified"
                )
else:
    if work.exists() and (not work.is_dir() or any(work.iterdir())):
        raise SystemExit(
            f"[ERROR] nonempty unmanifested WORK {work} cannot be adopted as a compatible cache. "
            "Use a new, empty WORK; existing data was not modified"
        )
    work.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(current, indent=2, sort_keys=True) + "\n"
    temporary = manifest_path.with_name(f"{manifest_path.name}.{os.getpid()}.tmp")
    temporary.write_text(encoded)
    os.replace(temporary, manifest_path)
PY

mkdir -p "$WORK"/{bam,truth,ru,split,var,bins/pileup/train,bins/pileup/validation,bins/fa/train,bins/fa/validation,candidates,models}

log "=== Clair3 quick training: $SAMPLE $REGION ==="
log "  depths: ${DEPTHS[*]}   validation depth: $VALIDATION_DEPTH   chunks: $CHUNKS   epochs: $EPOCHS   DO_RU=$DO_RU"

# -----------------------------------------------------------------------------
# 3. 准备输入：BAM 切片 + 真值子集
# -----------------------------------------------------------------------------
log "[1/9] 准备 BAM / 参考 / 真值"

if [ "$REGION" = "$CTG" ]; then
    BAM=$BAM_RAW
else
    BAM=$WORK/bam/region.bam
    if [ ! -s "$BAM" ]; then
        "$SAMTOOLS" view -@8 -b -o "$BAM" "$BAM_RAW" "$REGION" || die "BAM 切片失败"
        "$SAMTOOLS" index -@8 "$BAM"
    fi
fi

# 真值 VCF（保留 header）
if [ ! -s "$WORK/truth/region.vcf.gz" ]; then
    "$TABIX" -h "$TRUTH_VCF" "$REGION" | "$BGZIP" -c > "$WORK/truth/region.vcf.gz"
    "$TABIX" -f -p vcf "$WORK/truth/region.vcf.gz"
fi

# 真值 BED（取与该区域重叠的区间；全染色体时等价于 grep）
if [ "$REGION" = "$CTG" ]; then
    awk -v c="$CTG" '$1==c' "$TRUTH_BED" > "$WORK/truth/region.bed"
else
    start=${REGION#*:}; start=${start%-*}; end=${REGION#*-}
    awk -v c="$CTG" -v s="$start" -v e="$end" \
        '$1==c && $3>s && $2<e' "$TRUTH_BED" > "$WORK/truth/region.bed"
fi

log "  真值变异: $(gzip -dc "$WORK/truth/region.vcf.gz" | grep -vc '^#')   BED 区间: $(wc -l < "$WORK/truth/region.bed")"

# -----------------------------------------------------------------------------
# 4. Haplotag + 降采样（用 v5.0q 自带的相位）
# -----------------------------------------------------------------------------
log "[2/9] whatshap haplotag + 降采样"
HT_BAM=$WORK/bam/$CTG.ht.bam
if [ ! -s "$HT_BAM" ]; then
    "$WHATSHAP" haplotag \
        --output "$HT_BAM" --reference "$REF" --regions "$CTG" \
        --ignore-read-groups "$WORK/truth/region.vcf.gz" "$BAM" \
        || die "haplotag 失败"
    "$SAMTOOLS" index -@8 "$HT_BAM"
fi

[ "$VALIDATION_DEPTH" != "1000" ] || die "VALIDATION_DEPTH 必须是非 1000 的独立深度"
case " ${DEPTHS[*]} " in
    *" $VALIDATION_DEPTH "*) ;;
    *) die "DEPTHS_OVERRIDE 必须包含 VALIDATION_DEPTH" ;;
esac

TRAIN_BAM=$WORK/bam/${CTG}_1000.training.bam
VALIDATION_BAM=$WORK/bam/${CTG}_${VALIDATION_DEPTH}.bam
if [ ! -s "$TRAIN_BAM" ] || [ ! -s "$VALIDATION_BAM" ]; then
    # -s 按 template/read pair 选择；-U 将未选择的互补 reads 写给训练集。
    "$SAMTOOLS" view -@8 -s "${VALIDATION_DEPTH}.${VALIDATION_DEPTH}" \
        -U "$TRAIN_BAM" -b -o "$VALIDATION_BAM" "$HT_BAM" "$CTG" \
        || die "独立训练/验证 reads 分割失败"
    "$SAMTOOLS" index -@8 "$TRAIN_BAM"
    "$SAMTOOLS" index -@8 "$VALIDATION_BAM"
fi
BAM_1000=$TRAIN_BAM
log "  depth 1000 训练使用验证集之外的互补 reads；depth $VALIDATION_DEPTH 为 held-out reads"

for d in "${DEPTHS[@]}"; do
    [ "$d" = "1000" ] || [ "$d" = "$VALIDATION_DEPTH" ] && continue
    out=$WORK/bam/${CTG}_${d}.bam
    if [ ! -s "$out" ]; then
        # 其他训练深度只能从互补训练 reads 中抽取，不能重新引入验证 reads。
        "$SAMTOOLS" view -@8 -s "${d}.${d}" -b -o "$out" "$BAM_1000" "$CTG" \
            || die "降采样 $d 失败"
        "$SAMTOOLS" index -@8 "$out"
    fi
done

# -----------------------------------------------------------------------------
# 5. Repun （可选）
#
#    跳过时直接用 region.vcf.gz 作为标签来源。
# -----------------------------------------------------------------------------
log "[3/9] Repun UnifyRepresentation (DO_RU=$DO_RU)"
UNIFIED_VCF=$WORK/truth/region.vcf.gz
if [ "$DO_RU" = "1" ]; then
    "$PYPY3" "$CLAIR3" SplitExtendBed --bed_fn "$WORK/truth/region.bed" \
        --output_fn "$WORK/split/$CTG" --ctgName "$CTG" || die "SplitExtendBed 失败"
    "$PYPY3" "$CLAIR3" GetTruth --vcf_fn "$WORK/truth/region.vcf.gz" \
        --ctgName "$CTG" --var_fn "$WORK/var/var_$CTG" || die "GetTruth 失败"
    $PARALLEL -j "$THREADS" \
        "\"$PYPY3\" \"$CLAIR3\" UnifyRepresentation \
            --bam_fn \"$BAM_1000\" --var_fn \"$WORK/var/var_$CTG\" --ref_fn \"$REF\" \
            --bed_fn \"$WORK/truth/region.bed\" --extend_bed \"$WORK/split/$CTG\" \
            --output_vcf_fn \"$WORK/ru/vcf_${CTG}_{1}\" \
            --min_af 0.08 --chunk_id {1} --chunk_num $CHUNKS \
            --platform $PLATFORM --ctgName $CTG" ::: $CHUNK_LIST \
        > "$WORK/ru/RU.log" 2>&1 || die "UnifyRepresentation 失败（见 $WORK/ru/RU.log）"

    grep "^#" "$WORK/ru/vcf_${CTG}_1" > "$WORK/ru/unified.vcf"
    for f in "$WORK"/ru/vcf_${CTG}_*; do grep -v "^#" "$f"; done \
        | sort -k1,1 -k2,2n >> "$WORK/ru/unified.vcf"
    "$BGZIP" -f "$WORK/ru/unified.vcf"
    "$TABIX" -f -p vcf "$WORK/ru/unified.vcf.gz"
    UNIFIED_VCF=$WORK/ru/unified.vcf.gz
    log "  unified 变异: $(gzip -dc "$UNIFIED_VCF" | grep -vc '^#')"
fi

# 后续 GetTruth 用的 var 文件
"$PYPY3" "$CLAIR3" GetTruth --vcf_fn "$UNIFIED_VCF" \
    --ctgName "$CTG" --var_fn "$WORK/var/truth_$CTG" || die "GetTruth(unified) 失败"

# =============================================================================
# 6. Pileup bins
# =============================================================================
log "[4/9] Pileup bins: CreateTrainingTensorDirect"
for d in "${DEPTHS[@]}"; do
    bam=$BAM_1000; [ "$d" != "1000" ] && bam=$WORK/bam/${CTG}_${d}.bam
    $PARALLEL --colsep ' ' -j "${THREADS_LOW:-12}" \
        "\"$PY_TRAIN\" \"$CLAIR3\" CreateTrainingTensorDirect \
            --bam_fn $bam --ref_fn $REF --var_fn $WORK/var/truth_$CTG \
            --bin_fn $WORK/bins/pileup/tensor_${d}_${CTG}_{1} --ctgName $CTG \
            --samtools $SAMTOOLS --snp_min_af 0.07 --indel_min_af 0.08 \
            --extend_bed $WORK/split/$CTG --bed_fn $WORK/truth/region.bed \
            --allow_duplicate_chr_pos --pileup --platform $PLATFORM --shuffle \
            --maximum_non_variant_ratio 5 \
            --chunk_id {1} --chunk_num $CHUNKS" ::: $CHUNK_LIST \
        > "$WORK/bins/pileup/tensor_${d}.log" 2>&1 \
        || die "pileup tensor depth=$d 失败"
done

log "[5/9] Pileup bins: MergeBin"
for d in "${DEPTHS[@]}"; do
    files=$(ls "$WORK"/bins/pileup/tensor_${d}_${CTG}_* 2>/dev/null)
    [ -n "$files" ] || die "depth=$d 没有 tensor"
    split=train
    [ "$d" = "$VALIDATION_DEPTH" ] && split=validation
    "$PY_TRAIN" "$CLAIR3" MergeBin $files \
        --out_fn "$WORK/bins/pileup/${split}/bin_${SAMPLE}_${d}" --pileup \
        > "$WORK/bins/pileup/merge_${d}.log" 2>&1 \
        || die "pileup MergeBin depth=$d 失败"
done

# =============================================================================
# 7. Full-alignment bins
# =============================================================================
log "[6/9] FA bins: run_clair3 --pileup_only 产生候选"
for d in "${DEPTHS[@]}"; do
    [ -s "$WORK/candidates/pileup_${d}/pileup.vcf.gz" ] && { log "  跳过 depth=${d}（已有候选）"; continue; }
    bam=$BAM_1000; [ "$d" != "1000" ] && bam=$WORK/bam/${CTG}_${d}.bam
    "$PY_BIN" "$RUN_CLAIR3" \
        --bam_fn "$bam" --ref_fn "$REF" --threads "$THREADS" \
        --platform "$PLATFORM" --model_path "$MODEL_SRC" \
        --output "$WORK/candidates/pileup_${d}" \
        --bed_fn "$WORK/truth/region.bed" --ctg_name "$CTG" --pileup_only \
        > "$WORK/candidates/pileup_${d}.log" 2>&1 \
        || die "pileup_only depth=$d 失败（见 $WORK/candidates/pileup_${d}.log）"
done

log "[7/9] FA bins: SelectHetSnp"
for d in "${DEPTHS[@]}"; do
    # 注意：SelectHetSnp 内部自己按 --chunk_num 分块，只需按 depth 调用一次，
    # 不要再套一层 parallel（否则 chunk 号会被当成多余的位置参数）。
    "$PYPY3" "$CLAIR3" SelectHetSnp \
        --alt_fn "$WORK/candidates/pileup_${d}/pileup.vcf.gz" \
        --split_folder "$WORK/candidates" --sampleName "$SAMPLE" --depth "$d" \
        --ref_pct_full 0.15 --var_pct_full 1.0 --chunk_num "$CHUNKS" \
        --phasing_info_in_bam --phase --ctgName "$CTG" \
        > "$WORK/candidates/select_${d}.log" 2>&1 \
        || die "SelectHetSnp depth=$d 失败（见 $WORK/candidates/select_${d}.log）"
done

log "[8/9] FA bins: CreateTrainingTensorDirect + MergeBin"
for d in "${DEPTHS[@]}"; do
    bam=$BAM_1000; [ "$d" != "1000" ] && bam=$WORK/bam/${CTG}_${d}.bam
    $PARALLEL --colsep ' ' -j "${THREADS_LOW:-12}" \
        "\"$PY_TRAIN\" \"$CLAIR3\" CreateTrainingTensorDirect \
            --bam_fn $bam --ref_fn $REF --var_fn $WORK/var/truth_$CTG \
            --bin_fn $WORK/bins/fa/tensor_${d}_${CTG}_{1} --ctgName $CTG \
            --extend_bed $WORK/split/$CTG \
            --full_aln_regions $WORK/candidates/${SAMPLE}_${d}_${CTG_NUM}_{1} \
            --bed_fn $WORK/truth/region.bed --phasing_info_in_bam \
            --add_no_phasing_data_training --allow_duplicate_chr_pos \
            --platform $PLATFORM --shuffle --maximum_non_variant_ratio 1 \
            --chunk_id {1} --chunk_num $CHUNKS" ::: $CHUNK_LIST \
        > "$WORK/bins/fa/tensor_${d}.log" 2>&1 \
        || die "FA tensor depth=$d 失败"

    files=$(ls "$WORK"/bins/fa/tensor_${d}_${CTG}_* 2>/dev/null)
    [ -n "$files" ] || die "FA depth=$d 没有 tensor"
    split=train
    [ "$d" = "$VALIDATION_DEPTH" ] && split=validation
    "$PY_TRAIN" "$CLAIR3" MergeBin $files \
        --out_fn "$WORK/bins/fa/${split}/bin_${SAMPLE}_${d}" --platform "$PLATFORM" \
        > "$WORK/bins/fa/merge_${d}.log" 2>&1 \
        || die "FA MergeBin depth=$d 失败"
done

# =============================================================================
# 8. 训练（fine-tune from 公开 hifi_revio 权重）
# =============================================================================
log "[9/9] 训练"

for kind in pileup fa; do
    for split in train validation; do
        files=$(ls "$WORK"/bins/${kind}/${split}/bin_* 2>/dev/null)
        [ -n "$files" ] || die "$kind $split 没有独立 bins；请包含 VALIDATION_DEPTH 对应深度"
    done
done

train_one() {
    local kind=$1 pretrain=$2 extra=$3
    local out=$WORK/models/${kind}
    mkdir -p "$out"
    local cuda_dev="${GPU:-0}"
    [ "$cuda_dev" = "cpu" ] && cuda_dev=""      # GPU=cpu 时会强制 CPU
    ( cd "$out" || exit 1
      # Train 把 training.log / best_val_loss.pt 写到当前目录，所以必须先 cd
      CUDA_VISIBLE_DEVICES="$cuda_dev" "$PY_TRAIN" "$CLAIR3" Train \
          --bin_fn "$WORK/bins/${kind}/train" \
          --validation_fn "$WORK/bins/${kind}/validation" \
          --ochk_prefix "$out/${kind}" \
          --chkpnt_fn "$pretrain" \
          --learning_rate 1e-5 --platform "$PLATFORM" \
          --batch_size "$BATCH_SIZE" --chunk_size "$CHUNK_SIZE" \
          --maxEpoch "$EPOCHS" $extra
    ) || die "$kind 训练失败"
}

train_one pileup "$PRETRAIN_PILEUP" "--pileup"
train_one fa     "$PRETRAIN_FA"     "--add_indel_length True"


echo
log "=== 完成 ==="
echo "  输出根目录: $WORK"
echo
echo "  Pileup bins:"
ls -la "$WORK/bins/pileup/train" 2>/dev/null | grep bin_ | awk '{printf "    %-45s %s\n", $9, $5}'
echo "  FA bins:"
ls -la "$WORK/bins/fa/train" 2>/dev/null | grep bin_ | awk '{printf "    %-45s %s\n", $9, $5}'
echo "  Models:"
for k in pileup fa; do
    echo "    $k:"
    ls -la "$WORK/models/$k"/*.pt 2>/dev/null | awk '{printf "      %-45s %s\n", $9, $5}'
done
echo
