# Quick-testing Clair3-Rust

Use these checks after following [Installing Clair3-Rust](rust_installation.md). Run the commands from the repository root unless otherwise noted. The small fixture exercises native model inference and VCF generation without downloading a large BAM file.

## 1. Check the installed tools

```bash
clair3 --version
whatshap --version
samtools --version
```

`clair3` should be on `PATH` when installed with `make clair3 PREFIX=...`. To test a binary built in the repository instead, replace `clair3` below with `target/release/clair3`.

## 2. Run the small end-to-end fixture

Choose a HiFi model directory containing `pileup.pt` and `full_alignment.pt`. The golden fixture uses `hifi_sequel2`; a different checkpoint, including `hifi_revio`, can produce different calls:

```bash
# After make models PREFIX="${CONDA_PREFIX}":
MODEL_DIR="${CONDA_PREFIX}/bin/models/hifi_sequel2"

# Contributors with the prepared test assets can instead use:
# MODEL_DIR=target/test-data/models/hifi

test -f "${MODEL_DIR}/pileup.pt"
test -f "${MODEL_DIR}/full_alignment.pt"
```

Run without phasing for a minimal calling test:

```bash
clair3 run \
  --bam_fn tests/fixtures/data/toy.bam \
  --ref_fn tests/fixtures/data/toy.fa \
  --threads 2 \
  --platform hifi \
  --model_path "${MODEL_DIR}" \
  --output target/quick-test/hifi \
  --ctg_name chrToy \
  --no_phasing_for_fa
```

The command uses the native Rust pipeline and should produce `target/quick-test/hifi/merge_output.vcf.gz` and its tabix index. Display the output:

```bash
gzip -dc target/quick-test/hifi/merge_output.vcf.gz
gzip -dc target/quick-test/hifi/merge_output.vcf.gz | grep -vc '^#'
```

The tested fixture produces **eight** records. The first two are:

```text
chrToy  500  .  A     C     25.31  PASS  P  GT:GQ:DP:AD:AF  0/1:25:12:6,6:0.5000
chrToy  700  .  T     TGGG  19.22  PASS  F  GT:GQ:DP:AD:AF  1/1:19:12:0,12:1.0000
```

Different model checkpoints can produce different calls, so use the fixture model when checking exact values.

## 3. Test phasing with stock WhatsHap

Run the same fixture without `--no_phasing_for_fa` and choose a new output directory:

```bash
clair3 run \
  --bam_fn tests/fixtures/data/toy.bam \
  --ref_fn tests/fixtures/data/toy.fa \
  --threads 2 \
  --platform hifi \
  --model_path "${MODEL_DIR}" \
  --output target/quick-test/phasing \
  --ctg_name chrToy
```

With WhatsHap 2.8, the log includes independent BAM regions, native allele loading, and read selection. This tiny fixture has too few informative reads to form phase blocks; the full benchmark below exercises independent phase blocks. In both cases, the installed upstream DP solver is used and no replacement core is required. The final fixture VCF again contains eight records for the fixture model.

`CLAIR3_BIN` does not need to be set when invoking `clair3` directly: native phasing workers use the running executable. When using `run_clair3.sh`, set `CLAIR3_BIN` only to select a different binary explicitly.

## 4. Run build and fixture checks

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --locked -- -D warnings
cargo test --release --locked native_phase_command_uses_current_executable
cargo test --test pileup_parity --release --locked
cargo test --test model_parity --release --locked
```

The model-parity test uses official checkpoints when they are available under `target/test-data/models`. It reports a skip for unavailable checkpoints.

High-depth full-alignment fixtures have Linux-specific `.linux.npy` goldens because the original C implementation's `rand()` downsampling differs between macOS and Linux. The Linux goldens were generated with the unmodified original CFFI implementation and are compared exactly with Rust; the original macOS full-alignment goldens are unchanged. CPU model parity selects original-Python/PyTorch goldens by platform: `.macos.npy` on macOS, `.linux-amd64.npy` on Linux x86_64, and `.linux.npy` for Linux ARM64 pileup cases. The nine macOS model goldens were regenerated with original Python and official PyTorch 2.13 in the validated macOS environment. The official model-parity tolerance remains `1e-4`; no model-parity assertion was removed or weakened. Model-parity inference runs inside `tch::no_grad`, matching production and Python's `inference_mode`. CUDA continues to use its separate `.cuda.npy` goldens. The wrapped-checkpoint test runs when `target/wrapped-state-dict.pt` is available and reports a skip otherwise.

Run the unfiltered test suite:

```bash
cargo test --all-targets --release --locked
```

The unfiltered suite passed in the final native macOS validation, including `official_lstm_layers_match`. Tests with external-data ignore markers remain skipped by default.

## 5. Run the supplied HG002 chr20 benchmark

This is a full benchmark, not a quick smoke test. If the supplied files are available, use a separate output directory:

```bash
python3 run_clair3.py \
  --bam_fn ~/Temp/clairBenchmark/hg002_revio_20.bam \
  --ref_fn ~/Temp/clairBenchmark/GRCh38_no_alt_analysis_set.fasta \
  --threads 8 \
  --platform hifi \
  --model_path ~/Temp/clairBenchmark/v0.1_chm13_mix_torch \
  --output ~/Temp/clairBenchmark/output-rust-quick-test \
  --bed_fn ~/Temp/clairBenchmark/chr20.benchmark.bed \
  --ctg_name chr20
```

For stage timings, add `CLAIR3_PROFILE=1` to the environment before running it.

The validated runs contain:

| Output | Current native macOS Rust | Linux/ARM64 original and Rust | Linux/x86_64 original and Rust |
| --- | ---: | ---: | ---: |
| `pileup.vcf.gz` | 135,180 | 135,180 | 135,180 |
| `full_alignment.vcf.gz` | 40,494 | 40,494 | 40,493 |
| `merge_output.vcf.gz` | 94,551 | 94,550 | 94,559 |
| `tmp/phase_output/phase_vcf/phased_chr20.vcf.gz` | 33,824 | 33,824 | 33,824 |

Record counts alone are not a parity check. Compare contig, position, reference, alternate, genotype, QUAL, GQ, and ordering against the chosen reference. The original Python output is not completely identical to the optimized-Rust output, and the original implementation itself has platform-dependent numerical differences. Use a same-platform original run as the Linux reference; see the measured differences in the [optimization and verification report](rust_optimization.md).

## 6. Optional Docker test

Build the CPU image and run the phased fixture inside it. This checks Linux installation, model loading, VCF generation, and native phasing without relying on the host's libraries:

```bash
docker build -f Dockerfile -t clair3-rust:local .
docker run --rm -v "${PWD}:/workspace:ro" clair3-rust:local \
  --bam_fn /workspace/tests/fixtures/data/toy.bam \
  --ref_fn /workspace/tests/fixtures/data/toy.fa \
  --threads 2 \
  --platform hifi \
  --model_path /opt/models/hifi_sequel2 \
  --output /tmp/clair3-quick-test \
  --ctg_name chrToy
```

The fixture should contain eight records at `/tmp/clair3-quick-test/merge_output.vcf.gz` inside the container. The Dockerfile supports Linux x86_64 and ARM64; ARM64 uses an official PyTorch wheel and an upstream LongPhase source build. Both architectures passed this fixture and the selector serialization check. x86_64 was tested under `linux/amd64` emulation on an ARM64 host; add `--platform linux/amd64` to the Docker build and run commands to reproduce that configuration.

With the prepared model assets from section 2, run the complete local Docker checks as follows:

```bash
docker build --target build -t clair3-rust:build .
docker run --rm --entrypoint /bin/bash \
  -v "${PWD}/tests:/opt/clair3-src/tests:ro" \
  -v "${PWD}/target/test-data/models:/opt/clair3-src/target/test-data/models:ro" \
  -w /opt/clair3-src clair3-rust:build -c '
    rustup component add rustfmt clippy &&
    cargo fmt --all -- --check &&
    cargo clippy --all-targets --locked -- -D warnings &&
    cargo test --all-targets --release --locked
  '
```

The Linux/ARM64 and emulated Linux/x86_64 runs passed these checks and the phased fixture, including native selector initialization with the unmodified WhatsHap 2.8 core. Tests requiring additional external data remain marked as ignored, and the optional wrapped-checkpoint test reports a skip if its input is unavailable. With official WhatsHap 2.8 installed, explicitly run the selector regression as follows; install `whatshap==2.8` in the test container first if necessary:

```bash
cargo test --test pipeline_parity readset_serialization_preserves_variant_fields \
  --release --locked -- --ignored --exact
```

The test verifies position, allele, and quality through the installed core and the first-party selector. It does not modify WhatsHap or its DP solver.

## 7. Quick model-training workflow

The repository includes `scripts/run_quick_train_chr21.sh` for an end-to-end HiFi training smoke test. It exercises BAM/truth preparation, haplotagging, downsampling, representation unification, tensor generation, bin merging, candidate calling, and one-epoch fine-tuning of both pileup and full-alignment models. It requires external indexed input data and is not part of the small fixture test.

Prepare indexed reference/BAM/truth inputs, a directory of pretrained `pileup.pt` and `full_alignment.pt` checkpoints, and the `samtools`, `bgzip`, `tabix`, `parallel`, and `whatshap` tools. The `clair3.py` launcher accepts direct shell, Python, or PyPy invocation. Set `CLAIR3_BIN` if you need to select a particular native executable.

From the repository root, configure paths and run the bounded example:

```bash
REGION=chr21:13000000-14000000 \
DEPTHS_OVERRIDE="1000 222" \
VALIDATION_DEPTH=222 \
BATCH_SIZE=1000 CHUNK_SIZE=200 \
EPOCHS=1 GPU=cpu DO_RU=1 \
REF=/path/to/GRCh38_no_alt_analysis_set.fasta \
TRUTH_VCF=/path/to/HG002_GRCh38_v5.0q_smvar.vcf.gz \
TRUTH_BED=/path/to/HG002_GRCh38_v5.0q_smvar.benchmark.bed \
BAM_RAW=/path/to/hg002_revio_1000_21.bam \
MODEL_SRC=/path/to/hifi_revio \
WORK="${PWD}/target/quick-test/training" \
bash scripts/run_quick_train_chr21.sh
```

The supplied data's first chr21 interval is not covered by the truth BED; use the populated interval shown above or choose another interval with reads, truth variants, and BED regions. The default whole-chromosome run is much longer. The workflow uses `samtools view -s 222.222 -U` to select depth-222 validation templates/read pairs and writes their complementary, non-overlapping QNAME set to `chr21_1000.training.bam`. Thus the depth-1000 label retains the nominal depth metadata but contains approximately 77.8% of the haplotagged reads, not the full BAM; additional training depths are drawn only from this complement. `VALIDATION_DEPTH` selects a different non-1000 depth that must be present in `DEPTHS_OVERRIDE`. The bounded defaults `BATCH_SIZE=1000` and `CHUNK_SIZE=200` ensure that each read-disjoint split has complete batches; larger datasets can use the standard 2000/200 settings. Use `GPU=0` for CUDA training on a configured NVIDIA host; `GPU=cpu` forces CPU training. MPS training is not currently validated. `CHUNKS`, `THREADS`, `THREADS_LOW`, tool paths, and `CLAIR3_REPO` can be overridden when needed.

Successful one-epoch training produces:

```text
WORK/models/pileup/pileup.01.pt
WORK/models/pileup/best_val_loss.pt
WORK/models/pileup/training.log
WORK/models/fa/fa.01.pt
WORK/models/fa/best_val_loss.pt
WORK/models/fa/training.log
```

Check that both training and validation losses are finite and greater than zero, the disjoint merged bins are nonempty, and both checkpoints load successfully. Training-side `_20_` and sample exclusions do not remove independent validation archives. An explicit validation directory that yields no complete batches is an error rather than a zero-loss success. Resume with `--chkpnt_fn` is a weight-only restart with a fresh optimizer, epoch schedule, and RNG state, not exact continuation. See the detailed [pileup training](pileup_training.md) and [full-alignment training](full_alignment_training_r1.md) guides for production workflows.

The bounded workflow passed on native macOS ARM64 CPU, local Linux ARM64 CPU, emulated Linux x86_64 CPU, and native Linux x86_64 CUDA on an RTX 4090. Original Python/PyTorch successfully read the resulting LZF archives and loaded the generated checkpoints. Training remains stochastic; these checks establish working artifacts, not byte-identical model weights. See the [training validation results](rust_optimization.md#end-to-end-model-training-validation) for timings and known differences.

## 8. Optional GPU smoke test

CPU is the default. Use an explicit device to test an optional GPU backend:

```bash
# Apple Silicon:
clair3 run ... --use_gpu --device mps

# NVIDIA:
clair3 run ... --use_gpu --device cuda:0
```

Pileup remains on CPU for GPU runs to preserve candidate selection. Full-alignment inference uses dedicated GPU workers. See [GPU quick start](gpu_quick_start.md) for setup. There is currently no AMD-specific backend.

## Expected results and next checks

- `clair3 --version` runs and the native executable can be located.
- The fixture produces a compressed, indexed VCF with the expected records.
- Phasing completes using upstream WhatsHap 2.8 and the native helpers.
- Build, formatting, and applicable fixture tests pass.
- The full benchmark, when run, produces the platform-appropriate counts above; use the same-platform original-Python output to assess numerical and genotype parity.

For optimization history, benchmark results, and known limitations, read [Rust optimization and verification](rust_optimization.md).
