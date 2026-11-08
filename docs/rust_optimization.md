# Rust rewrite, optimization, and verification report

## Scope and final policy

This report covers the `clair3-rust` branch. The `main` branch was not changed.

The implementation retains the native Rust caller and the validated CPU tensor and inference optimizations described below. CPU remains the default. MPS and CUDA are optional; there is no AMD-specific backend. No quantization or candidate-filter relaxation is used.

The final implementation uses pristine third-party sources and the installed upstream WhatsHap 2.8 DP solver. It restores four pipeline-level phasing optimizations: independent phase-block scheduling, parallel BAM preparation, native allele loading, and native read selection. DP-specific solver changes, including a replacement phase core and altered cost calculations, are excluded. Gloo is downloaded from a pinned upstream commit, and `torch-sys` is resolved from crates.io without a patch. Other WhatsHap versions use the stock phasing path.

The experimental 30-second CPU target was canceled before completion. The best measured runtime with the former custom DP changes was 42.62 seconds. The current implementation retains only the four optimizations above; its benchmark is reported separately below.

## Benchmark

The following command was used for the original Python implementation and the Rust implementation on the supplied HG002 Revio chr20 data with eight requested threads:

```bash
python3 run_clair3.py \
  --bam_fn ~/Temp/clairBenchmark/hg002_revio_20.bam \
  --ref_fn ~/Temp/clairBenchmark/GRCh38_no_alt_analysis_set.fasta \
  --threads 8 \
  --platform hifi \
  --model_path ~/Temp/clairBenchmark/v0.1_chm13_mix_torch \
  --output ~/Temp/clairBenchmark/output \
  --bed_fn ~/Temp/clairBenchmark/chr20.benchmark.bed \
  --ctg_name chr20
```

The output directory was changed between runs to preserve the results. Timing is wall time for the complete command, not an isolated model call. Measurements on this host vary with background load; large changes were accepted only after output comparison and a clean confirmation.

| Version | Wall time | Notes |
| --- | ---: | --- |
| Original Python implementation | 309.10 s | 626.80 s user CPU; 41.70 s system CPU; 1.97 GB peak RSS |
| Initial native Rust implementation | 301.89 s | Native end-to-end port before the optimization milestones |
| Best optimized CPU run with former custom DP changes | 42.62 s | Historical result at `02311b1`; independent confirmation 43.19 s |
| Pristine dependency baseline without restored phasing optimizations | 191.09 s | Stock WhatsHap 2.8; phasing 160.7 s, including 113.8 s in its DP solver |
| Final four phasing optimizations with the upstream DP solver | 72.01 s | Phasing 29.2 s, including 20.6 s in the original DP solver |

The final implementation is the last configuration. Restoring the four pipeline-level optimizations reduced the pristine dependency baseline by 62.3% without changing the DP implementation. The final run used 384.01 seconds of user CPU, 20.42 seconds of system CPU, and 3.93 GB peak resident memory. The 42.62-second historical result included DP solver changes and is not the deliverable.

## Final output verification

Comparisons are made by `(contig, position, reference, alternate)` rather than by row number. The comparison checks record keys, genotypes, QUAL, and GQ. The phased comparison also checks GT/PS fields and record ordering.

### Preservation of the validated optimized-Rust reference

The exact-source reference is `output-rust-baseline-rebuilt`. This is the appropriate control for isolating optimization and cleanup changes from earlier Python-to-Rust numerical differences.

| Output | Reference records | Final records | Result |
| --- | ---: | ---: | --- |
| Pileup VCF | 135,180 | 135,180 | Exact: no missing/extra keys, genotype, QUAL, or GQ differences |
| Phased intermediate VCF | 33,824 | 33,824 | Exact GT/PS fields and record ordering |
| Full-alignment VCF | 40,493 | 40,493 | Exact: no missing/extra keys, genotype, QUAL, or GQ differences |
| Merged VCF | 94,572 | 94,572 | Exact: no missing/extra keys, genotype, QUAL, or GQ differences |

### Comparison with the original Python run

The original output is not byte-identical to the optimized-Rust reference. It would be incorrect to claim otherwise.

| Comparison | Original Python | Optimized-Rust reference | Observed difference |
| --- | ---: | ---: | --- |
| Pileup records | 135,180 | 135,180 | Same keys and genotypes; 78 QUAL values differ, maximum absolute difference 0.01; GQ unchanged |
| Phased GT/PS fields | 33,824 | 33,824 | Identical, including record ordering |
| Full-alignment records | 40,494 | 40,493 | 293 Python-only and 292 Rust-only keys; 22 genotype differences among 40,201 common keys; 29,794 QUAL differences (maximum 26.85) and 15,822 GQ differences |
| Merged records | 94,551 | 94,572 | 206 Python-only and 227 Rust-only keys; 22 genotype differences among 94,345 common keys; 20,649 QUAL differences (maximum 15.42) and 8,206 GQ differences |

The residual differences are concentrated in full-alignment model decoding and indel representation. They are not introduced by the final dependency cleanup, which preserves the optimized-Rust reference exactly. The exact source of every per-record difference was not isolated. Consequently, the final version is **not** claimed to be byte-identical to the original Python output.

### Current native macOS comparison

The current branch was revalidated natively on macOS after the CUDA work. Both original Python and Rust used the current environment with official PyTorch 2.13 and WhatsHap 2.8. The complete supplied HG002 Rust run finished successfully in 206.62 seconds. This was a functional validation under concurrent test load, not a new controlled optimization benchmark.

| Output | Original and Rust records | Positions, alleles, genotypes, and ordering | QUAL differences | GQ differences |
| --- | ---: | --- | ---: | ---: |
| Pileup | 135,180 | Exact match | 83; maximum absolute difference 0.01 | 1 |
| Full alignment | 40,494 | Exact match | 50; maximum absolute difference 0.07 | 0 |
| Merged | 94,551 | Exact match | 58; maximum absolute difference 0.07 | 0 |

All 33,824 phased GT/PS fields and their record ordering match exactly. This comparison uses a different original-Python/PyTorch environment from the historical optimized-Rust reference above. Production code is unchanged in this validation, so the difference is not evidence of a new production-code regression; the source of every residual per-record numerical difference was not isolated. The full unfiltered test suite passed, followed by the explicit WhatsHap selector serialization regression and the phased HiFi fixture. The previously reported native macOS `official_lstm_layers_match` assertion now passes without being removed or skipped.

Two test-only corrections were needed. Nine `.macos.npy` model goldens were generated using original Python/PyTorch in the current macOS environment. In particular, original Python reproduced the former ONT golden difference of `0.0012764633`, so this was not a Rust model regression; the strict `1e-4` model tolerance was retained. The macOS stride-two NNPACK comparison encountered one batch-35 difference of `5.9604645e-7`; its test tolerance is now `1e-6`, instead of requiring bit-for-bit equality. No production model, calling, phasing, or third-party algorithm was changed.

## Retained implementation optimizations

### Native implementation and data handling

- Replaced Python, PyPy, and CFFI calling with Rust, direct HTSlib BAM/CRAM access, and LibTorch model execution.
- Implemented pileup and full-alignment calling, candidate selection, decoding, VCF/GVCF merging, realignment, representation unification, post-processing, metrics, checkpoint conversion, and training helpers natively.
- Preserved shell-compatible `run_clair3.sh` and `run_clair3.py` entry points and legacy model formats.
- Used bounded tensor generation, independent job scheduling, per-job model ownership, and job-indexed result collection to preserve read, candidate, and output order.
- Used deterministic per-job random sampling and native macOS/Linux sampling bridges.

### Phasing with the upstream DP solver

- Schedule independent connected phase blocks in forked workers, largest first, and merge haplotypes in genomic order. Each block calls the installed upstream `PedigreeDPTable`; no solver source or cost calculation is changed.
- Capture alignment spans while pileup calling proceeds, partition BAM reads across workers, and use native Rust/HTSlib allele loading with exact-match shortcuts, bounded edit-distance workspaces, and compatible supplementary-read merging.
- Use the first-party C++ read selector with unmodified upstream headers and the installed WhatsHap 2.8 core. Preserve the original 15× selection and priority-queue behavior.
- Load the official core to resolve selector symbols; do not substitute a bundled `whatshap.core`. Accelerated behavior is limited to the tested 2.8 API/ABI, with stock behavior for other versions.

### Parallel CPU inference

- Scheduled the requested number of independent calling jobs with one LibTorch intra-op thread per model worker.
- Split pileup tensor generation and inference for non-GVCF CPU runs, then distributed unchanged 200-record inference batches over independently owned models.
- Split full-alignment tensor generation from inference, using four tensor workers and eight inference workers with unchanged 100-record batches.
- Added scheduler tests covering job identity, batch boundaries, and record order.
- Used four full-alignment BAM decoder threads and ordered read preparation.

### Native model kernels

- Cached supported NNPACK convolution kernels and spatially grouped independent samples while retaining per-sample normalization and remainder handling.
- Added a direct im2col/CBLAS backend for full-alignment 3x3 convolutions, with correct stride, padding, NCHW layout, and per-channel bias. The original backend remains a fallback.
- Added a native bidirectional LSTM path with batched input projections, recurrent matrix multiplication, zero initial states, separate input/recurrent biases, and batch-first views over time-first storage.
- Retained LibTorch-compatible sigmoid and tanh evaluation and the original LibTorch path for training, unsupported systems, and opt-out use.
- Kept CPU as the default and retained optional, dedicated MPS/CUDA workers. Pileup remains on CPU for GPU runs to preserve candidate selection.

### Tensor construction and BAM access

- Enabled libdeflate through `rust-htslib` to reduce BGZF decompression and CRC overhead without changing BAM records.
- Used binary variant lookup, bounded CIGAR scans, decoded sequence slicing, exact single-edit shortcuts, and precomputed sampling in full-alignment tensor construction.
- Replaced repeated HTSlib pileup work with direct, cache-windowed CIGAR walking in 4096-base windows.
- Used range-delta counts for reference-matching runs, then reconstructed depth, strand, and reference counts before applying the original candidate rules.
- Preserved mismatch, deletion, insertion, reference-skip, boundary, high-depth, and Illumina tie-breaking behavior.
- Reused indexed BAM readers where appropriate and reduced unnecessary record copies and container overhead.
- Retained `CLAIR3_DIRECT_PILEUP=0` as the HTSlib fallback and `CLAIR3_CONV_GEMM=0` / `CLAIR3_NATIVE_LSTM=0` as backend opt-outs.

## Optimization milestones

These are historical measurements from the optimization process, not claims that every intermediate revision used the final dependency policy.

| Commit | Main change | Reported measurement |
| --- | --- | --- |
| `4754c23` | Native Rust rewrite | 301.89 s |
| `7bf7de9` | Parallel independent WhatsHap phase blocks | Phasing 192.4 to 117.66 s; restored with the upstream DP solver |
| `182f317` | Parallel BAM allele loading for phasing | 162.61 s with six workers; restored |
| `3c94f82` | Requested worker count for independent calling | 145.01 s with eight workers |
| `e745991` | Optional dedicated MPS/CUDA inference workers | Optional MPS run 105.94 s; CPU remains default |
| `dee6460` | Parallel CPU inference and streamlined tensors | 112.76 s |
| `928aa07` | Decoupled pileup inference and native phasing | 94.21 s; DP changes excluded, other optimizations retained |
| `a0526f5` | Cached NNPACK, tensor scheduling, native read selection | 74.55 s; native read selection restored |
| `d616e2c` | Native LSTM and native allele loading | 70.20 s; native allele loading restored |
| `bbe5a8b` | libdeflate BAM decompression | 57.69 s |
| `f1d2e6b` | Direct full-alignment GEMM convolution | 51.72 s |
| `fdcac30` | Cache-windowed CIGAR walking | 44.62 s |
| `02311b1` | Range-counted reference matches | 42.62 s; clean confirmation 43.19 s |

Each accepted end-to-end improvement greater than 10% was committed with implementation, benchmark, and validation details. Later commits below 10% were also documented individually. The 30-second target was not reached and was explicitly canceled; no result is presented as meeting it.

## Third-party dependency cleanup

The final revision makes these installation changes:

- Removed `third_party/whatshap` and the replacement DP phase core. The first-party read selector now lives in `csrc/fast_readselect.cpp`; it is built against the official WhatsHap 2.8 headers, downloaded from PyPI and verified with SHA-256 `073ee7946d563e0125dbce58d9fc86a7a3fae377f27859dacb398309374b2cbe`.
- Removed `third_party/torch-sys` and the Cargo source patch. `Cargo.lock` now identifies the official crates.io `torch-sys 0.26.0` package and checksum.
- Removed `third_party/gloo`. `build.rs` downloads upstream commit `bcd1672ee07538123ea8f4fac76832efc58fb8ef`; its tree was compared with the previous vendored tree before pinning.
- Restored background BAM span capture, native allele reading, and independent phase-block scheduling. The Python wrapper calls the installed `PedigreeDPTable` and never substitutes a replacement DP core. Accelerated paths are restricted to WhatsHap 2.8; other versions retain stock behavior.
- Added `make python_tools`, which installs unmodified `whatshap==2.8`. On macOS it also installs the official `torch==2.13.0` wheel for the upstream `LIBTORCH_USE_PYTORCH=1` build path.
- Updated CPU and GPU Docker images to install upstream WhatsHap 2.8 and removed their vendored-source copy steps.
- Retained the existing upstream LibTorch, LongPhase, model, samtools, bgzip, and tabix download/install paths.

No downloaded third-party source file is edited by the final build. CMake configuration, compilation, and linking are performed on the downloaded sources without source patches.

## Rejected or superseded experiments

These changes were not included because they were slower, unproven, exceeded numerical tolerance, or conflicted with the final dependency policy.

| Experiment | Outcome |
| --- | --- |
| Larger inference batches (1,000 pileup / 400 full alignment) | Exact output comparison, but 53.75 s and about 7.9 GB peak RSS; reverted |
| Profile-guided Rust optimization | Exact output comparison; 53.07 s versus a 52.15 s fresh generic control; not enabled |
| Architecture-specific native compiler flags | Exact output comparison but slower, 52.77 s; reverted |
| Parallel forward/reverse LSTM directions | Focused tests and full VCF comparisons passed; 18.00 s direct model time versus 17.00 s for ordinary workers; rejected |
| Four full-alignment inference workers | Exact output comparison; 45.95 s total and 13.805 s inference versus about 10.4 s with eight; rejected |
| Native read-selection bookkeeping rewrite | Exact output comparison but no measurable end-to-end gain; reverted |
| Single-sample generic-phasing specialization | Identical GT/PS output but no measurable improvement; reverted before the phasing cleanup |
| Tiled convolution packing | Identical tensors in the microbenchmark but no speed improvement; not applied |
| Custom stride-one Winograd convolution | 17.15 s full inference versus 10.43 s; rejected |
| Approximate `vvtanhf`, `vvexpf`, and SLEEF `u35` math | Faster microbenchmarks but approximately 0.0013–0.0017 model error, exceeding the 0.0001 tolerance; not enabled |
| Column-major CBLAS | Correct matrices but slightly slower in the microbenchmark; rejected |
| BNNS convolution | Filter creation returned null for supported test layouts; no replacement enabled |
| Raw-record caching and full/shared-stage overlap | No validated end-to-end gain; reverted |
| NNPACK/CBLAS pool tuning and four workers with two BLAS threads | No validated gain; reverted |
| Pre-transposed weights and reusable GEMM scratch | No validated end-to-end gain; reverted |
| Uninitialized LSTM scratch | No validated end-to-end gain; reverted |
| NHWC and per-sample direct NCHW GEMM | No validated gain; reverted |
| Range-tree read selection and generic/HTSlib buffer pools | No validated gain; reverted |
| SNP realignment shortcuts and compact read-selection indexing | Correct output after fixes, but no meaningful end-to-end gain; not retained |
| Native phasing result serialization and VCF sorting/serialization changes | Correct output but no measurable gain; reverted |
| Full read-preparation overlap | No validated gain; reverted |
| LongPhase source build | Build failed; not used to replace verified phasing |
| Quantization and candidate-filter relaxation | Not permitted because they change model precision or biological inputs |

## Verification commands and status

The final build used a new target directory so the build could not reuse the removed local `torch-sys` or Gloo copies:

```bash
LIBTORCH_USE_PYTORCH=1 cargo build --release --locked
```

Additional checks:

- `git diff --check`: passed during cleanup and checked again before commit.
- Fresh release build with pristine upstream dependencies: passed.
- `make python_tools` in an isolated virtual environment: passed; installed official PyTorch 2.13.0 and WhatsHap 2.8.
- `make clair3` installation and `clair3 --version` smoke test: passed with the installed libraries and stock WhatsHap 2.8.
- Full supplied HG002 chr20 command and VCF comparisons: see above.
- `cargo fmt --all -- --check`: passed.
- `cargo clippy --all-targets --locked -- -D warnings`: passed.
- `cargo test --all-targets --release --locked`: passed in the final native macOS revalidation with 59 unit tests and one ignored, followed by all integration targets with their documented external-data ignores. No test filter or LSTM skip was used. The explicit selector regression was then run separately and passed.
- Direct installed `clair3 run` without a `CLAIR3_BIN` override: passed with native phasing. The small HiFi fixture produced eight expected records both with and without phasing.
- Native macOS does not require a Docker or LongPhase build for the validated `clair3 run` path. Linux/ARM64 and emulated Linux/x86_64 Docker validation are recorded below, followed by remote NVIDIA RTX 4090 validation. AMD GPU execution remains unverified.

## Local Docker Linux validation

The GitHub-based Linux workflow was removed in favor of testing with the user's local Docker daemon. The Linux/ARM64 CPU image builds successfully with the official PyTorch CPU wheel and an unmodified upstream LongPhase 1.7.3 source build. It installs official WhatsHap 2.8 and Debian samtools 1.16.1; no third-party source is patched. The build stage also contains samtools for pipeline tests. The earlier LongPhase source-build failure listed under rejected experiments was on macOS, not a failure of this validated Linux build.

Local Docker formatting, Clippy with warnings treated as errors, and `cargo test --all-targets --release --locked` passed in the final ARM64 revalidation. This includes 60 unit tests with one ignored, all full-alignment tensor fixtures, model parity, pipeline tests, and training tests. External-data tests retain their documented ignore markers. The wrapped-checkpoint test was also exercised with its generated input and skips when that input is unavailable. The explicit WhatsHap-dependent selector serialization regression passed, the phased HiFi fixture produced the expected eight variants, and the native selector initialized against the unmodified upstream WhatsHap core.

High-depth full-alignment golden tensors are platform-specific because the original C implementation uses the platform C library's `rand()` for downsampling. Six Linux `.linux.npy` fixtures were generated with the unmodified original CFFI implementation; Rust matches them exactly. Existing macOS full-alignment fixtures were not modified. Pileup model Linux goldens were generated with the original Python/PyTorch implementation, preserving the strict `1e-4` tolerance rather than loosening it. Every official model case also has an original-Python x86_64 `.linux-amd64.npy` golden; the new macOS model goldens are described above. Model-parity inference uses `tch::no_grad`, matching production and Python's `inference_mode`; gradient-enabled test execution was not treated as a production numerical failure.

### Original-versus-Rust HG002 comparison on Linux/ARM64

Both implementations ran the supplied HG002 chr20 command with eight threads and the supplied `v0.1_chm13_mix_torch` checkpoint in local Docker. The original Python run used unmodified source from `main` with official dependencies. Results below compare Rust with that same-platform original, not with the macOS baseline:

| Output | Records in each version | Positions, alleles, genotypes, and ordering | QUAL differences | GQ differences |
| --- | ---: | --- | --- | --- |
| Pileup | 135,180 | Exact match | 3; maximum absolute difference 0.01 | 0 |
| Full alignment | 40,494 | Exact match | 65; maximum absolute difference 0.11 | 3 |
| Merged | 94,550 | Exact match | 50; maximum absolute difference 0.11 | 3 |

All 33,824 phased GT/PS fields and their record ordering match the original exactly. A batch-100 comparison with the supplied full-alignment checkpoint found a maximum absolute probability difference of `1.66893e-6`, accounting for the small remaining numerical differences. The larger macOS-versus-Linux differences also occur between two original Python runs and must not be attributed to the Rust port. No new optimization or production-code change was introduced during this Linux validation.

The complete original run took 375.24 seconds and the earlier Rust validation run took 209.22 seconds, approximately 1.79 times faster. The Rust run overlapped other validation work, so these timings document the test environment rather than a controlled performance benchmark; they do not revise the macOS optimization measurements above. The post-CUDA ARM64 rerun completed in 967.04 seconds while other Docker validation jobs were active and reproduced every comparison above. This heavily loaded timing is not a performance regression measurement.

### Original-versus-Rust HG002 comparison on Linux/x86_64

The x86_64 validation used `linux/amd64` containers on the ARM64 host. It verifies the x86_64 build and behavior, not native x86_64 hardware performance. Formatting, Clippy, and the complete test suite passed with 60 unit tests and one ignored, followed by an explicitly executed WhatsHap-dependent serialization regression. The phased fixture and native selector check passed on both x86_64 and ARM64.

In the post-CUDA x86_64 rerun, formatting, Clippy, all test targets, and the explicit selector regression passed again. The subsequent repeated fixture/selector smoke command was not run because the user stopped the validation job after the regression test. This interruption was not a test failure; the earlier x86_64 fixture and selector smoke results remain documented above.

The official x86_64 WhatsHap wheel exports the old C++ string ABI, while the tested ARM64 source-built core uses the new ABI. Linux now builds the same first-party selector with `_GLIBCXX_USE_CXX11_ABI=0` and `_GLIBCXX_USE_CXX11_ABI=1` and selects the appropriate module from symbols exported by the installed WhatsHap 2.8 core. The upstream DP solver and all third-party sources remain unmodified.

The first full x86_64 run also exposed an evaluation-order bug in first-party deserialization: three `reader.integer()` calls in one `addVariant` expression could be evaluated right-to-left, assigning quality to position. The fields are now read sequentially. A saved-job check retained 4,083 positions across 26,825 reads, selected 3,507 reads, and produced no unsorted reads. The permanent regression test verifies position, allele, and quality through serialization and deserialization on the installed core.

The original Python and final Rust x86_64 outputs compare as follows:

| Output | Records in each version | Positions, alleles, genotypes, and ordering | QUAL differences | GQ differences |
| --- | ---: | --- | --- | --- |
| Pileup | 135,180 | Exact match | 1; maximum absolute difference 0.01 | 0 |
| Full alignment | 40,493 | Exact match | 29; maximum absolute difference 0.11 | 0 |
| Merged | 94,559 | Exact match | 28; maximum absolute difference 0.11 | 0 |

All 33,824 phased GT/PS fields and their record ordering match the original exactly. Optimized phasing completed successfully rather than falling back to the installed selector. The complete original run took 878.63 seconds and the earlier Rust validation run took 386.89 seconds, approximately 2.27 times faster. Both ran under emulation with other validation work active, so these are validation timings rather than a controlled performance benchmark or a new optimization result. The post-CUDA x86_64 rerun completed in 1,222.20 seconds while the emulated test build ran concurrently and reproduced every comparison above. It validates behavior under emulation, not native x86_64 performance.

### CUDA validation on an NVIDIA RTX 4090

The CUDA version was built and run on a remote native Linux/x86_64 server with an RTX 4090, driver 580.82.09, Docker, and the NVIDIA Container Toolkit. Validation used an isolated checkout and SHA-256-verified copies of the supplied HG002 inputs and checkpoints. The server had sufficient free space after the user cleaned it; no unrelated server images were removed.

The CUDA packaging now installs CMake 3.31.6, preserves the compiled binary outside the BuildKit cache, discovers CUDA headers through `CUDA_HOME`/`CUDA_PATH`, and adds the required library search and link paths. Official unmodified NVIDIA cuSPARSE-Lt 0.8.1, NCCL 2.29.7, NVSHMEM 3.4.5, and CUPTI packages provide the companion libraries required by LibTorch 2.13. The runtime also builds unmodified upstream samtools 1.24 because Ubuntu 22.04's samtools 1.13 is below Clair3's minimum. Explicit CUDA linkage keeps CUDA operator registration available in model-test binaries.

Formatting, Clippy, and the full CUDA-enabled test suite passed with 60 unit tests and one ignored. The opt-in CUDA model-parity test passed for every official model case against original-Python CUDA goldens generated on the same GPU, with the existing `0.001` tolerance unchanged. The explicit WhatsHap selector serialization regression passed. The phased HiFi fixture completed with eight expected variants and working selector serialization. Direct original-Python-versus-Rust ONT inference on CUDA was exactly equal at every probability.

The final supplied HG002 workflow completed in **60.89 seconds** with CUDA. This is a functional validation timing, not a controlled CPU-versus-GPU benchmark. The following comparison uses the original Python **x86_64 CPU** output, not an original-Python GPU run:

| Output | Records | Numerical and genotype differences from original x86_64 CPU |
| --- | --- | --- |
| Pileup | 135,180 in each | Exact keys, genotypes, and ordering; 111 QUAL differences, maximum 0.01; no GQ differences |
| Full alignment | 40,493 in each | 10,215 QUAL differences, maximum 0.11; 115 GQ differences; one genotype difference and one ref/alt representation difference |
| Merged | 94,559 original; 94,560 CUDA | 4,764 QUAL differences, maximum 0.11; 60 GQ differences; one genotype difference and one extra CUDA record |

All 33,824 phased GT/PS fields match the original exactly. The residual full-alignment differences are at `chr20:1352240` and `chr20:16576770`, where genotypes or ref/alt representation cross decision thresholds while QUAL/GQ remain very close. Original Python also differs at these sites across CPU platforms. This is not exact original-Python parity, and no claim is made that every residual per-record difference has been eliminated. No production model algorithm or third-party DP solver was changed for this validation.

Historically, in an earlier macOS test environment, the strict `official_lstm_layers_match` assertion checked the official ONT pileup model and found an LSTM2 difference of `0.053098708391189575`, with first-timestep difference `0.0008534938097000122`. It also reproduced at the earlier `f1d2e6b` milestone and was not a result of the range-counting or dependency-cleanup changes. This historical failure is not the final macOS result: the assertion passes in the current official PyTorch 2.13 environment, and the complete unfiltered suite passed. It was not removed or skipped.

A direct comparison with `CLAIR3_NATIVE_LSTM=0` on the supplied HiFi benchmark showed that the LSTM backend affects results: 51 pileup QUAL values changed (maximum absolute difference 0.01), one GQ value changed, and downstream full-alignment/merged outputs had 22 genotype differences and altered representations. Although the historical failing fixture was ONT-specific, this comparison confirms that LSTM numerical differences can be consequential. The assertion remains in place rather than being deleted to hide numerical differences.

## Documentation and final cleanup verification

The documentation pass corrects current installation, Docker, GPU, quick-demo, dwell-time, and CLI guidance while retaining historical benchmark tables as historical references. Documentation checks validate local link targets and balanced code fences. The [project prompt history](clair3_rust_prompts.md) preserves all 38 user-authored prompts verbatim and in chronological order; automated messages and assistant or tool output are excluded.

Code cleanup removed the obsolete `unused_parens` suppression and 14 redundant parenthesized `break` expressions in `src/call/decode.rs`, plus the unreferenced `contains_index` helper in `csrc/fast_readselect.cpp`. No functional expression, model algorithm, compatibility option, fallback, or third-party source changed.

Native macOS formatting, Clippy with warnings treated as errors, and the unfiltered release test suite passed: 59 unit tests passed and one was ignored, with all integration targets passing their applicable tests. The explicit WhatsHap selector serialization regression passed. The phased HiFi fixture produced eight `chrToy` records with exact alleles and genotypes; its comparison correctly excludes the separate `chr:colon` fixture record because the command selects only `chrToy`.

A fresh supplied HG002 chr20 run completed in 62.23 seconds. All 135,180 pileup, 40,494 full-alignment, 94,551 merged, and 33,824 phased records exactly match the current native macOS Rust reference in full record content and ordering. This is a functional cleanup check, not a new optimization milestone or controlled benchmark. Earlier original-Python numerical differences and platform limitations documented above remain unchanged.

## End-to-end model-training validation

### Current correctness-first results

The current evidence comes from the Task 4 review-fix state at `b06121c`. It uses complementary, read-disjoint training and validation splits. Earlier one-epoch figures at the end of this section are retained only as historical initial smoke results; their overlapping-read validation metrics are superseded and must not be used as current evidence.

The supplied chr21 workflow uses `chr21:13000000-14000000`, nominal depths 1000 and 222, representation unification, and bounded batches. The originally suggested 1–2 Mb interval has no overlapping truth BED intervals and no reads; the populated 13–14 Mb interval was selected from the supplied data instead. The workflow covers haplotagging, complementary-read splitting, representation unification, tensor generation, bin merging, candidate calling, and fine-tuning of both models. `scripts/run_quick_train_chr21.sh` provides the configurable native workflow without modifying the original Clair3 checkout.

#### Fixed mismatches and actual root causes

| Area | Demonstrated root cause | Corrected behavior and covering evidence |
| --- | --- | --- |
| Candidate CIGAR decoding | `create_candidate_details` did not advance `query_pos` for `Cigar::Ins`, so later SNP bases and inserted sequences were read at stale offsets. | Advance by the insertion length. A real indexed `5M2I5M` BAM regression checks the insertion and both following bases exactly; all four real-data anchors pass. |
| Phased SNP/insertion/deletion matching | `reference_base` was a `u8`; `to_string()` produced decimal ASCII such as `"67"` instead of `"C"`, causing allele lookup and phased genotype matching to fail. | Convert the byte to a character in all three arms. Separate SNP, insertion, and deletion regressions preserve haplotype 0 versus haplotypes 1/2 behavior. |
| Multi-allelic ALT order | Native collected selected alleles in a lexical `BTreeSet`, changing the deterministic count-ranked candidate order. | Emit selected alleles in `candidate.alternates` order. Exact allele and GT/GQ/DP/AF/PS assertions pass at the four anchors. Original Python hash-set order remains a separate residual, not a stable reference ordering. |
| Deterministic dropout stream | Rust seeded before model construction, so initialization consumed random values that original training reserved for the dropout stream. | Seed after model construction and checkpoint loading. This is existing trainer ordering, not a normal-training seed override; explicit seeds remain test-only. |
| AdamW weight decay | Rust hardcoded `1e-7`, but original pileup uses `1e-4` and full alignment uses `1e-7`. | Select decay by model type. Exact checkpoint comparison exposed the pileup-only residual after the seed fix. |
| Independent validation archives | Training `_20_` and `exclude_training_samples` filename exclusions were also applied to explicit validation directories, discarding independent data; an empty batch set could return zero loss. | Apply exclusions only to training and fail explicitly when validation has no complete batch. The regression retains `_20_` validation data, checks both training exclusions, and rejects an empty directory. |
| Quick-workflow read independence | The initial workflow used the full depth-1000 BAM for training while selecting depth 222 from those same reads, yielding 1,340 shared native QNAMEs. | Select validation with samtools `-s` and write the `-U` complement for training; additional training depths use only that complement. Per-QNAME conservation and zero intersection are asserted on all four platforms. |
| Validation diagnostics and bounded smoke settings | A missing validation path was described as a training archive directory; the standard 2000-row batch yielded empty bounded splits, and an existing-candidate log hit an unbound-variable error. | Select the correct directory kind in the error, use smoke defaults `BATCH_SIZE=1000` and `CHUNK_SIZE=200`, and brace `${d}`. CLI defaults for larger workflows are unchanged. |

The original `SplitExtendBed` final-interval omission was **not** reproduced. Native retains the final merged interval, guarded by a fixture with both merged and final disjoint intervals. Controlled unification uses the same original extended BED on both sides so that this original boundary bug cannot change chunking and conceal candidate differences.

#### Deterministic original-versus-Rust training

Both fixtures create 40 labeled examples with `BinWriter`, read them through `BinArchive`, and supply identical arrays and initial checkpoints to unmodified original model/loss/epoch code and native Rust. Seed 0, one fixed 20-example batch, `batch_size=chunk_size=20`, one optimizer step, and no validation pass make the comparison controlled without changing normal RNG, dropout, or public training interfaces.

- Pileup compares both task losses; full alignment compares all four task losses, including indel-length tasks. Final per-task and actual total losses match exactly.
- Every updated checkpoint tensor matches original with **zero tolerance**. Both implementations reload the checkpoints and produce finite `(2,24)` pileup and `(2,90)` full-alignment outputs.
- The explicit two-process native C10d test uses nonzero, rank-specific focal-loss gradients and compares every all-reduced tensor on both ranks with independently reproduced rank-local arithmetic means. Maximum gradient difference is exactly 0. The prior test's discarded reductions and zero-derivative softmax-sum were corrected in the test, not hidden by a merged-batch comparison.
- The original deterministic core consumes arrays exported from the native archive because the initial test environment lacked h5py. Original h5py/DataLoader execution was subsequently performed separately as described below; it is not inferred from the deterministic fixture.

#### Complementary reads and checkpoint restart semantics

The quick script runs `samtools view -s 222.222 -U chr21_1000.training.bam ... -o chr21_222.bam`. Selection is by template/read pair, and `-U` writes reads not selected for validation. The depth-1000 filename remains nominal metadata: training now contains approximately 77.8% of the haplotagged reads rather than the full BAM. All additional training depths are subsampled only from this complement.

| Platform | Full QNAMEs | Training QNAMEs | Validation QNAMEs | Shared | Per-QNAME record conservation |
| --- | ---: | ---: | ---: | ---: | --- |
| Native macOS ARM64 | 6,147 | 4,807 | 1,340 | 0 | Holds for every QNAME |
| Linux ARM64 CPU | 6,147 | 4,756 | 1,391 | 0 | Holds for every QNAME |
| Emulated Linux x86_64 CPU | 6,147 | 4,756 | 1,391 | 0 | Holds for every QNAME |
| Linux x86_64 CUDA, RTX 4090 | 6,147 | 4,756 | 1,391 | 0 | Holds for every QNAME |

Native and Linux tool/runtime stacks select different template counts from identical raw inputs; each split is internally complementary, but identical cross-platform datasets are not claimed. Two repeated native samtools probes were identical, and a paired-SAM control kept both mates on the same side.

Checkpoint resume is explicitly a **weight-only restart** with a fresh AdamW optimizer, epoch schedule, and RNG schedule, not exact optimizer/epoch/RNG continuation. Controlled tests through the actual `train()` interface used seed 29 and exactly reproduced fresh one-epoch loss, validation loss, complete log, and every checkpoint tensor: pileup loss `3.3610267639160156`, full-alignment loss `10.915680885314941`. A different-starting-checkpoint control changed tensors, proving that the selected checkpoint is used. Real-data restarts also created a fresh `.01.pt` and one log row. No production restart change was needed.

#### Current bounded platform results

All four runs used the same mounted input sources with disjoint output directories and the corrected complementary-read split. Times are functional validation timings, not controlled performance benchmarks. Losses are `training / validation` after one epoch.

| Rust platform | Device | Review-fix image build | Bounded workflow | Pileup loss | Full-alignment loss |
| --- | --- | ---: | ---: | ---: | ---: |
| Native macOS ARM64 | CPU | Native release build | 38.594 s | 0.768210 / 0.441530 | 9.129483 / 0.651172 |
| Local Docker Linux ARM64 | CPU | 219.335 s | 75.597 s | 0.785272 / 0.416018 | 8.996529 / 0.503815 |
| Local Docker Linux x86_64, emulated on ARM64 | CPU | 222.079 s | 82.491 s | 0.776526 / 0.416016 | 8.823021 / 0.503768 |
| Docker Linux x86_64 | RTX 4090 CUDA | 101.704 s, including source synchronization | 36.766 s, including SSH | 0.784594 / 0.416052 | 8.933035 / 0.503819 |

The current archive matrices and nonzero-label totals are:

| Platform | Pileup train / validation | Full-alignment train / validation |
| --- | --- | --- |
| Native macOS ARM64 | `(3610,33,18)` / `(5396,33,18)`; 14,440 / 21,584 | `(3397,55,33,8)` / `(3492,55,33,8)`; 13,588 / 13,968 |
| Linux ARM64 and emulated x86_64 | `(3649,33,18)` / `(4310,33,18)`; 14,596 / 17,240 | `(3398,55,33,8)` / `(3543,55,33,8)`; 13,592 / 14,172 |
| Linux x86_64 CUDA | Same as Linux CPU | Same as Linux CPU |

Every label matrix has shape `(rows,90)` and exactly four nonzero labels per row. Complete 90-class distributions and SHA-256 values remain in the Task 4 review-fix logs. These shapes and totals do not imply identical matrices, example order, or weights.

#### Convergence, early stopping, loaders, and checkpoint compatibility

On the current native disjoint data, three-epoch `(training, validation)` losses were:

| Model | Epoch 1 | Epoch 2 | Epoch 3 |
| --- | --- | --- | --- |
| Pileup | `(0.696778, 0.441970)` | `(0.574981, 0.440305)` | `(0.654248, 0.438492)` |
| Full alignment | `(9.195758, 0.650995)` | `(7.702460, 0.986795)` | `(8.630635, 1.689801)` |

The two explicit real-data tests passed together in 422.14 seconds. Both models stopped at epoch 11 of 12 in the early-stopping runs and passed the restart checks. Pileup training improves overall and validation decreases in this bounded run. **Full-alignment validation worsens during the bounded native and CUDA tests even though training improves overall; no full-alignment validation-convergence claim is made.** A zero-learning-rate probe changed full-alignment validation from `0.658174` to `1.181838` with unchanged weights, demonstrating training-mode state effects; constant validation was therefore not used as a substitute for the actual early-stopping check.

CUDA three-epoch results on the corrected bins were:

| Model | Training losses | Validation losses |
| --- | --- | --- |
| Pileup | `0.737139, 0.599273, 0.720201` | `0.416417, 0.414775, 0.412990` |
| Full alignment | `9.363139, 7.571423, 8.508289` | `0.504872, 0.900081, 1.585159` |

Original h5py/DataLoader training read the generated LZF archives on each platform's own disjoint bins. On native macOS its training/validation losses were `0.312487 / 0.442511` for pileup and `8.875472 / 0.651737` for full alignment. On CUDA they were `0.719845 / 0.416259` and `9.286971 / 0.507250`. Native macOS stock multiprocessing spawn cannot pickle the original module's open h5py handles; the unmodified original loader was exercised with its Linux-style `fork` start method. This spawn limitation remains.

Both native and original-trained checkpoints reload in Rust and original Python on the tested platforms, with 26 pileup and 81 full-alignment tensors and finite real-input outputs of `(2,24)` and `(2,90)`. CUDA reload used `cuda:0`. **Original CUDA DDP was not executed**: the host has one RTX 4090, while original rank-to-`cuda:<rank>` mapping requires distinct devices. Only original single-process CUDA DataLoader training and native two-process gradient averaging are claimed.

#### Current held-out chr20 F1 evidence

Upstream GIAB v5.0q truth/BED files were used without modification. Every record in the chr21 13–14 Mb training interval exactly matches the supplied training truth. Calling used the supplied HG002 chr20 BAM/reference/BED. Rust and original models were each fine-tuned for one epoch on the same current native disjoint bins; the pretrained rows remain available as the common starting-checkpoint comparison. Upstream `jmcdani20/hap.py:v0.3.12` with vcfeval and `--pass-only` produced these PASS metrics (identical to ALL here):

| Model | INDEL recall / precision / F1 | SNP recall / precision / F1 |
| --- | --- | --- |
| Pretrained, unchanged reference | 0.959621 / 0.965542 / **0.962572** | 0.997248 / 0.997287 / **0.997268** |
| Rust disjoint-read one-epoch | 0.945211 / 0.962045 / **0.953554** | 0.996734 / 0.997168 / **0.996951** |
| Original disjoint-read one-epoch | 0.946526 / 0.963052 / **0.954718** | 0.997129 / 0.997050 / **0.997090** |

There is **no invented F1 acceptance threshold** and no improvement claim. These one-epoch runs do not improve on pretrained INDEL F1. Current Rust and original calling times were 82.811 and 70.636 seconds; the unchanged pretrained calling run took 71.916 seconds. These are evidence timings, not a new optimization campaign.

#### Executed tests, named skips, and residual differences

The Task 4 review-fix native `cargo test --all-targets --release --locked` suite passed **113 tests with 0 failed and 15 ignored**. Task 5 final verification on the documentation and test-alias state again passed formatting, Clippy with warnings denied, the same 113/0/15 suite, the explicit selector regression, and the eight-record phased HiFi golden comparison; shell syntax, local-link/fence, and whitespace checks also passed. Ignored tests were not counted as executed passes. Explicit additional evidence includes:

- The four real anchor assertions, independent-validation regression, two controlled restart regressions, two real-data convergence/stop/restart tests, and four complementary-read platform checks passed.
- The explicit WhatsHap 2.8 selector regression and colon-contig phasing test passed with installed upstream tools; no WhatsHap source or DP solver changed.
- `official_realigned_bam_matches` passed with all 28,730 BAM records matching the native golden. The original candidate BED had been generated before the original workflow later lacked its realigner binary.
- Ordinary fixture tests cover ONT, HiFi, Illumina, dwell time, CRAM, GVCF, and long indels. Native, Linux ARM64, emulated Linux x86_64, and CUDA completed the phased HiFi fixture. The native fixture's eight `chrToy` records exactly match its golden.
- The supplied native HG002 run took 64.082 seconds. Against the unchanged `output-rust-macos-training` reference, pileup 135,180, full-alignment 40,494, merged 94,551, and phased 33,824 records each had **0 differing records**, including complete record content and ordering. This is regression preservation, not a claim of identical original-versus-native unification.

Four official-baseline comparisons remain **named skips, not passes**: `candidate_and_merge_outputs_match`, `official_illumina_calls_match`, `production_full_call_includes_flanking_candidates`, and `official_hifi_phased_context_matches`. Their messages enumerate missing official legacy baseline outputs/CFFI capture. An attempted unmodified original Illumina baseline regeneration reached full alignment but could not complete because `preprocess/realign/realigner` is absent; partial output was not promoted to a baseline. Python 3.14 also exposed an original argparse `None`/`%(default)f` incompatibility; the successful partial generation used Python 3.11.

**The exact eight-chunk original/native unification differential remains RED.** The latest explicitly executed Task 4 run used all eight chunks, produced 1,547 original and 1,552 native records, and failed. Its two diagnosed residual classes are:

1. Five native-only records at `chr21:13600450`, `13600728`, `13626244`, `13698384`, and `13735785`. Their indel AFs exceed the harness's explicit `--min_af 0.08` but fall below original CTFA's separately configured effective default `--indel_min_af 0.15`; original generated no candidate rows there. Giving original the native candidate rows reproduced the five records, isolating candidate selection rather than `unify_regions`. These are original-side indel-threshold policy differences, not fixed or hidden native defects.
2. Run-dependent original Python hash-order ALT differences. Original collects selected alternates in a `set`, so exact ALT sequence order changes between runs. REF, the ALT allele set, and compared metadata agree at these sites, but the harness intentionally continues to compare exact ordered records. Task 2's final run had five such differences; Task 4's run had eight, demonstrating run dependence. Native deterministic candidate order is not a claim of identical original/native unification.

The differential was not rerun in Task 5 merely to recreate this expected failure, and no threshold, assertion, ALT normalization, or original source was changed to make it pass. Other numerical original-versus-Rust calling differences documented earlier in this report also remain subject to their platform-specific evidence.

Known warnings are retained: Cargo reports the pre-existing `proc-macro-error2 v2.0.1` future-incompatibility notice, and emulated x86_64 emits NNPACK `Unsupported hardware` warnings after otherwise passing checks. Neither warning is presented as a test failure or suppressed. Training order, dropout randomness, and numerical backends still prevent a general claim of byte-identical datasets or weights.

### Historical initial smoke results (superseded)

The following table records the earlier initial workflow before the complementary-read review fix. Its training and validation reads overlapped, so these losses are **not current independent-validation evidence**. They are retained only as historical functional timings. Losses are `training / validation` after one epoch.

| Implementation and platform | Device | Time | Pileup loss | Full-alignment loss |
| --- | --- | ---: | ---: | ---: |
| Original Python, Linux x86_64 | RTX 4090 | 97.61 s | 0.601942 / 0.804108 | 8.894713 / 0.540924 |
| Rust, native macOS ARM64 | CPU | 35.91 s | 0.552349 / 0.624374 | 10.663239 / 1.117319 |
| Rust, local Linux ARM64 | CPU | 77.23 s | 0.412903 / 1.062350 | 10.244903 / 0.875078 |
| Rust, emulated Linux x86_64 | CPU | 93.62 s | 0.412619 / 1.062332 | 10.170233 / 0.875076 |
| Rust, native Linux x86_64 | RTX 4090 | 31.70 s | 0.394398 / 1.062423 | 10.109503 / 0.875232 |

At that initial milestone, original archives contained 3,182/3,984 pileup rows and 3,203/3,501 full-alignment rows at depths 1000/222. Rust macOS contained 3,411/5,278 and 3,299/3,486; Linux and CUDA contained 3,411/4,235 and 3,299/3,539. Original Python/h5py read the LZF archives, original PyTorch loaded the checkpoints, and both trained models completed the native calling fixture. These artifacts demonstrated cross-compatibility, not identical data or weights.

Earlier training-specific fixes retained by the current implementation include:

- Added the `CreateTrainingTensorDirect` alias and Python/PyPy-compatible `clair3.py` launcher, preserving direct shell execution and `CLAIR3_BIN` selection.
- Separated the internal Clap IDs for `call_fn` and `bin_fn`, which previously shared `output` and incorrectly required an unrelated output argument.
- Implemented native training-mode `SelectHetSnp` options, quality selection, contig-prefix handling, chunk filenames, BED intervals, and phased metadata, while retaining the existing phasing-only VCF filter.
- Detect gzip and Zstandard text input by magic bytes rather than filename suffix, fixing gzip-compressed `GetTruth` outputs without `.gz` extensions.
- Accept legacy case-insensitive boolean spellings such as `--add_indel_length True` for training.
- Replaced Blosc-only HDF5 output with original-compatible LZF-plus-shuffle filters so original Python/h5py can read the archives without an additional plugin.
- Added the official GNU Parallel package to both Docker runtime images for the optional training workflow.

The historical native macOS suite then passed 63 unit tests and one ignored plus applicable integration tests. A fresh HG002 calling run at that milestone took 71.42 seconds and preserved all 135,180 pileup, 40,494 full-alignment, 94,551 merged, and 33,824 phased records exactly relative to the then-current native macOS Rust reference. Neither this historical timing nor the current validation resumes or revises the canceled optimization goal.

## Reproduction notes

- Keep output directories distinct when comparing implementations.
- Compare records by contig, position, reference, and alternate; do not rely only on row positions or record counts.
- Check genotype, QUAL, GQ, and phased GT/PS fields explicitly.
- Benchmark one configuration at a time on an otherwise quiet host. Several rejected changes looked beneficial only under different background load.
- The supplied benchmark is a HiFi CPU case on Apple Silicon. It is not a performance claim for other platforms, sequencing chemistries, or GPU backends.
