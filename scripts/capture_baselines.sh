#!/usr/bin/env bash
# Capture the official Clair3 legacy baselines with the unmodified reference
# implementation and install them as committed test fixtures.
#
# The legacy baselines are the outputs of the original Python Clair3. They used
# to live only in a machine-local directory (target/legacy/baseline), which made
# every parity test that depends on them permanently skipped. This script
# regenerates those outputs with the official reference container and records a
# provenance manifest, so the expectations live in the repository and can be
# re-verified by anyone.
#
# Usage:
#   scripts/capture_baselines.sh [--platform ilmn|hifi] [--image IMAGE]
#                                 [--demo-dir DIR] [--work DIR] [--fixtures DIR]
#                                 [--verify-only]
#
# The reference implementation is never modified: the container runs the
# published image as-is, and the demo data is downloaded unmodified from the
# upstream Clair3 demo location.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

IMAGE="hkubal/clair3:v2.0.5"
PLATFORMS=("ilmn" "hifi")
DEMO_DIR="$REPO_ROOT/target/baselines/demo"
WORK_DIR="$REPO_ROOT/target/baselines/work"
FIXTURES_DIR="$REPO_ROOT/tests/fixtures/baselines"
VERIFY_ONLY=0
MAX_FIXTURE_BYTES=$((32 * 1024 * 1024))

DEMO_BASE_URL="https://www.bio8.cs.hku.hk/clair3/demo/quick_demo"
CONTIG="chr20"
REGION_START=100000
REGION_END=300000

log() { printf '[capture-baselines] %s\n' "$*"; }
die() { printf '[capture-baselines] ERROR: %s\n' "$*" >&2; exit 1; }

while [ $# -gt 0 ]; do
  case "$1" in
    --platform) PLATFORMS=("$2"); shift 2 ;;
    --image) IMAGE="$2"; shift 2 ;;
    --demo-dir) DEMO_DIR="$2"; shift 2 ;;
    --work) WORK_DIR="$2"; shift 2 ;;
    --fixtures) FIXTURES_DIR="$2"; shift 2 ;;
    --verify-only) VERIFY_ONLY=1; shift ;;
    -h|--help) sed -n '2,30p' "$0"; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

sha256() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | cut -d' ' -f1
  else
    shasum -a 256 "$1" | cut -d' ' -f1
  fi
}

# Upstream demo assets for the quick demo, unmodified.
platform_assets() {
  case "$1" in
    ilmn)
      printf '%s\n' \
        "GRCh38_chr20.fa" \
        "GRCh38_chr20.fa.fai" \
        "HG003_chr20_demo.bam" \
        "HG003_chr20_demo.bam.bai" \
        "HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz" \
        "HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz.tbi" \
        "HG003_GRCh38_chr20_v4.2.1_benchmark_noinconsistent.bed"
      ;;
    hifi)
      printf '%s\n' \
        "GRCh38_no_alt_chr20.fa" \
        "GRCh38_no_alt_chr20.fa.fai" \
        "HG003_chr20_demo.bam" \
        "HG003_chr20_demo.bam.bai" \
        "HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz" \
        "HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz.tbi" \
        "HG003_GRCh38_chr20_v4.2.1_benchmark_noinconsistent.bed"
      ;;
    *) die "unknown platform: $1" ;;
  esac
}

# Upstream directory name for the demo assets of a platform.
platform_upstream_dir() {
  case "$1" in
    ilmn) printf 'illumina' ;;
    hifi) printf 'pacbio_hifi' ;;
    *) die "unknown platform: $1" ;;
  esac
}

# Locally prepared copy of the same upstream assets.
platform_local_dir() {
  case "$1" in
    ilmn) printf 'illumina' ;;
    hifi) printf 'pacbio_hifi' ;;
    *) die "unknown platform: $1" ;;
  esac
}

platform_reference() {
  case "$1" in
    ilmn) printf 'GRCh38_chr20.fa' ;;
    hifi) printf 'GRCh38_no_alt_chr20.fa' ;;
    *) die "unknown platform: $1" ;;
  esac
}

# Prefer the already-prepared local demo data, otherwise download the upstream
# file unmodified. Every input is hashed so the manifest records exactly what
# produced the baselines.
prepare_demo() {
  local platform="$1" asset local_dir upstream_dir prepared_dir
  local_dir="$DEMO_DIR/$platform"
  upstream_dir="$(platform_upstream_dir "$platform")"
  prepared_dir="$REPO_ROOT/target/test-data/$(platform_local_dir "$platform")"
  mkdir -p "$local_dir"
  while IFS= read -r asset; do
    if [ -f "$local_dir/$asset" ]; then
      continue
    fi
    if [ -f "$prepared_dir/$asset" ]; then
      cp "$prepared_dir/$asset" "$local_dir/$asset"
      continue
    fi
    log "downloading $upstream_dir/$asset"
    curl -fsSL --retry 3 -o "$local_dir/$asset" "$DEMO_BASE_URL/$upstream_dir/$asset" \
      || die "failed to download $asset for $platform"
  done < <(platform_assets "$platform")
  printf '%s\t%s\t%s\t%s\n' \
    "$CONTIG" "$REGION_START" "$REGION_END" > "$local_dir/quick_demo.bed"
}

run_reference() {
  local platform="$1" demo_dir="$2" output_dir="$3"
  mkdir -p "$output_dir"
  docker run --rm --platform linux/amd64 \
    -v "$demo_dir:/input:ro" \
    -v "$output_dir:/output" \
    "$IMAGE" \
    /opt/bin/run_clair3.sh \
    --bam_fn="/input/HG003_chr20_demo.bam" \
    --ref_fn="/input/$(platform_reference "$platform")" \
    --threads=4 \
    --platform="$platform" \
    --model_path="/opt/models/$platform" \
    --output=/output \
    --bed_fn=/input/quick_demo.bed
}

# Files the parity tests consume, relative to the reference output directory.
# The Illumina run phases with WhatsHap and therefore also produces a phased
# BAM, while the HiFi run phases with LongPhase and keeps its phased VCF only.
expected_outputs() {
  cat <<'EOF'
pileup.vcf.gz
full_alignment.vcf.gz
merge_output.vcf.gz
tmp/split_beds/chr20
tmp/phase_output/phase_vcf/chr20.vcf
tmp/phase_output/phase_vcf/phase_qual
tmp/phase_output/phase_vcf/phased_chr20.vcf.gz
tmp/full_alignment_output/candidate_bed/chr20.0_1
tmp/full_alignment_output/candidate_bed/qual
tmp/full_alignment_output/candidate_bed/FULL_ALN_FILE_chr20
tmp/full_alignment_output/candidate_bed/FULL_ALN_FILES
tmp/full_alignment_output/full_alignment_chr20.0_1.vcf
tmp/merge_output/merge_chr20.vcf
EOF
}

# Extra reference outputs that only exist for one platform.
expected_platform_outputs() {
  case "$1" in
    ilmn)
      printf '%s\n' \
        "tmp/phase_output/phase_bam/chr20.bam" \
        "tmp/phase_output/phase_bam/chr20.bam.bai" \
        "tmp/full_alignment_output/candidate_bed/full_aln_regions_chr20"
      ;;
    hifi) : ;;
    *) die "unknown platform: $1" ;;
  esac
}

# The reference full-alignment tensor for the HiFi demo, captured with the
# reference C implementation, so the Rust tensor can be compared byte for byte.
capture_reference_tensor() {
  local platform="$1" demo_dir="$2" fixtures_dir="$3" output_dir="$4"
  [ "$platform" = "hifi" ] || return 0
  log "capturing the reference full-alignment tensor for $platform"
  cat > "$output_dir/capture_tensor.py" <<'PYTHON'
import sys
from argparse import Namespace

sys.path.insert(0, "/opt/bin")
from preprocess.CreateTensorFullAlignmentFromCffi import CreateTensorFullAlignment

CreateTensorFullAlignment(
    Namespace(
        ctgStart=1,
        ctgEnd=300000,
        full_aln_regions="/fixtures/tmp/full_alignment_output/candidate_bed/chr20.0_1",
        ref_fn="/input/GRCh38_no_alt_chr20.fa",
        ctgName="chr20",
        bam_fn="/input/HG003_chr20_demo.bam",
        platform="hifi",
        phased_vcf_fn="/fixtures/tmp/phase_output/phase_vcf/phased_chr20.vcf.gz",
        minMQ=5,
        minBQ=0,
        enable_long_indel=False,
        enable_dwell_time=False,
        no_phasing_for_fa=False,
        extend_bed=None,
        bed_fn=None,
        tensor_can_fn="/output/hifi-cffi.npy",
    )
)
PYTHON
  docker run --rm --platform linux/amd64 \
    -v "$demo_dir:/input:ro" \
    -v "$fixtures_dir:/fixtures:ro" \
    -v "$output_dir:/output" \
    "$IMAGE" python3 /output/capture_tensor.py \
    || die "reference tensor capture failed for $platform"
  cp "$output_dir/hifi-cffi.npy" "$fixtures_dir/hifi-cffi.npy"
  cp "$output_dir/hifi-cffi.npy.info" "$fixtures_dir/hifi-cffi.npy.info"
  log "installed $platform/hifi-cffi.npy and $platform/hifi-cffi.npy.info"
}

install_outputs() {
  local platform="$1" source_dir="$2" destination="$3" relative size
  mkdir -p "$destination"
  while IFS= read -r relative; do
    [ -n "$relative" ] || continue
    if [ ! -f "$source_dir/$relative" ]; then
      log "MISSING expected reference output: $relative"
      continue
    fi
    size=$(wc -c < "$source_dir/$relative" | tr -d ' ')
    if [ "$size" -gt "$MAX_FIXTURE_BYTES" ]; then
      log "SKIP oversized reference output ($size bytes): $relative"
      continue
    fi
    mkdir -p "$destination/$(dirname "$relative")"
    cp "$source_dir/$relative" "$destination/$relative"
    log "installed $platform/$relative ($size bytes)"
  done < <(expected_outputs; expected_platform_outputs "$platform")
}

vcf_record_count() {
  if [ ! -f "$1" ]; then
    printf '0'
    return
  fi
  case "$1" in
    *.gz) gzip -dc "$1" | grep -vc '^#' || true ;;
    *) grep -vc '^#' "$1" || true ;;
  esac
}

write_manifest() {
  local platform="$1" destination="$2" output_dir="$3" demo_dir="$4"
  local file relative size records
  {
    printf '{\n'
    printf '  "platform": "%s",\n' "$platform"
    printf '  "reference_image": "%s",\n' "$IMAGE"
    printf '  "reference_image_id": "%s",\n' "$(docker image inspect --format '{{.Id}}' "$IMAGE" 2>/dev/null || printf 'unknown')"
    printf '  "reference_entrypoint": "/opt/bin/run_clair3.sh",\n'
    printf '  "contig": "%s",\n' "$CONTIG"
    printf '  "region": "%s-%s",\n' "$REGION_START" "$REGION_END"
    printf '  "inputs": {\n'
    local first=1
    while IFS= read -r file; do
      [ -f "$demo_dir/$file" ] || continue
      [ $first -eq 1 ] || printf ',\n'
      first=0
      printf '    "%s": {"sha256": "%s", "bytes": %s}' \
        "$file" "$(sha256 "$demo_dir/$file")" "$(wc -c < "$demo_dir/$file" | tr -d ' ')"
    done < <(platform_assets "$platform")
    printf '\n  },\n'
    printf '  "outputs": {\n'
    first=1
    while IFS= read -r relative; do
      [ -n "$relative" ] || continue
      [ -f "$destination/$relative" ] || continue
      size=$(wc -c < "$destination/$relative" | tr -d ' ')
      records=$(vcf_record_count "$destination/$relative")
      [ $first -eq 1 ] || printf ',\n'
      first=0
      printf '    "%s": {"sha256": "%s", "bytes": %s, "records": %s}' \
        "$relative" "$(sha256 "$destination/$relative")" "$size" "$records"
    done < <(expected_outputs; expected_platform_outputs "$platform")
    printf '\n  }\n}\n'
  } > "$destination/manifest.json"
  log "wrote $destination/manifest.json"
}

verify_fixtures() {
  local platform="$1" destination="$2" relative
  while IFS= read -r relative; do
    [ -n "$relative" ] || continue
    if [ ! -f "$destination/$relative" ]; then
      log "MISSING fixture: $platform/$relative"
      return 1
    fi
  done < <(expected_outputs; expected_platform_outputs "$platform")
  return 0
}

for platform in "${PLATFORMS[@]}"; do
  demo_dir="$DEMO_DIR/$platform"
  output_dir="$WORK_DIR/$platform"
  destination="$FIXTURES_DIR/$platform"
  if [ "$VERIFY_ONLY" -eq 1 ]; then
    log "verifying installed fixtures for $platform"
    verify_fixtures "$platform" "$destination" || die "fixtures incomplete for $platform"
    continue
  fi
  log "preparing demo inputs for $platform"
  prepare_demo "$platform"
  if [ ! -d "$output_dir/tmp" ]; then
    log "running the reference implementation for $platform (image $IMAGE)"
    run_reference "$platform" "$demo_dir" "$output_dir" \
      || die "reference run failed for $platform"
  else
    log "reusing existing reference output in $output_dir"
  fi
  install_outputs "$platform" "$output_dir" "$destination"
  capture_reference_tensor "$platform" "$demo_dir" "$destination" "$output_dir"
  write_manifest "$platform" "$destination" "$output_dir" "$demo_dir"
done

log "done"
