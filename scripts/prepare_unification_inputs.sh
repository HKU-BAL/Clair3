#!/usr/bin/env bash
# Prepare the inputs for the differential representation-unification harness.
#
# The harness compares the original Python UnifyRepresentation with the native
# Rust implementation on identical inputs. Those inputs are normally produced
# by the original implementation itself, which used to require a machine-local
# checkout of the original Clair3 tree. This script materializes that
# implementation from this repository's own history at a pinned commit and
# generates the truth variants and the extended BED with it, so the harness
# can run from a clean clone.
#
# Usage:
#   scripts/prepare_unification_inputs.sh \
#       --bam_fn BAM --reference REF --truth_vcf VCF --bed BED \
#       --contig chr21 --output DIR [--original_ref COMMIT] [--python PYTHON]
#
# The pinned original implementation is never modified: it is extracted
# read-only from git and only executed.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

ORIGINAL_REF="1edf39d"
PYTHON_BIN="python3"
BAM=""
REFERENCE=""
TRUTH_VCF=""
BED=""
CONTIG=""
OUTPUT=""
ORIGINAL_ROOT=""

log() { printf '[prepare-unification] %s\n' "$*"; }
die() { printf '[prepare-unification] ERROR: %s\n' "$*" >&2; exit 1; }

while [ $# -gt 0 ]; do
  case "$1" in
    --bam_fn) BAM="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --truth_vcf) TRUTH_VCF="$2"; shift 2 ;;
    --bed) BED="$2"; shift 2 ;;
    --contig) CONTIG="$2"; shift 2 ;;
    --output) OUTPUT="$2"; shift 2 ;;
    --original_ref) ORIGINAL_REF="$2"; shift 2 ;;
    --original_root) ORIGINAL_ROOT="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

[ -n "$BAM" ] || die "--bam_fn is required"
[ -n "$REFERENCE" ] || die "--reference is required"
[ -n "$TRUTH_VCF" ] || die "--truth_vcf is required"
[ -n "$BED" ] || die "--bed is required"
[ -n "$CONTIG" ] || die "--contig is required"
[ -n "$OUTPUT" ] || die "--output is required"

mkdir -p "$OUTPUT"

if [ -z "$ORIGINAL_ROOT" ]; then
  ORIGINAL_ROOT="$OUTPUT/original-src"
  if [ ! -d "$ORIGINAL_ROOT" ]; then
    log "materializing the original Clair3 implementation at $ORIGINAL_REF"
    mkdir -p "$ORIGINAL_ROOT"
    git archive "$ORIGINAL_REF" | tar -x -C "$ORIGINAL_ROOT"
  fi
  ENTRY_POINT="$ORIGINAL_ROOT/clair3.py"
else
  ENTRY_POINT="$ORIGINAL_ROOT/clair3.py"
fi
[ -f "$ENTRY_POINT" ] || die "original Clair3 entry point not found: $ENTRY_POINT"

export PYTHONPATH="$ORIGINAL_ROOT${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONDONTWRITEBYTECODE=1

EXTENDED_BED="$OUTPUT/extended_bed_$CONTIG"
TRUTH_VARS="$OUTPUT/var_$CONTIG"

log "running the original SplitExtendBed"
"$PYTHON_BIN" "$ENTRY_POINT" SplitExtendBed \
  --bed_fn "$BED" \
  --output_fn "$EXTENDED_BED" \
  --ctgName "$CONTIG" \
  || die "original SplitExtendBed failed"

log "running the original GetTruth"
"$PYTHON_BIN" "$ENTRY_POINT" GetTruth \
  --vcf_fn "$TRUTH_VCF" \
  --ctgName "$CONTIG" \
  --var_fn "$TRUTH_VARS" \
  || die "original GetTruth failed"

[ -s "$EXTENDED_BED" ] || die "the extended BED is empty: $EXTENDED_BED"
[ -s "$TRUTH_VARS" ] || die "the truth variants file is empty: $TRUTH_VARS"

log "prepared unification inputs"
cat <<EOF
export CLAIR3_ORIGINAL_ROOT=$(cd "$(dirname "$ENTRY_POINT")" && pwd)
export CLAIR3_UNIFICATION_BAM=$BAM
export CLAIR3_UNIFICATION_TRUTH=$TRUTH_VARS
export CLAIR3_UNIFICATION_REFERENCE=$REFERENCE
export CLAIR3_UNIFICATION_BED=$BED
EOF
