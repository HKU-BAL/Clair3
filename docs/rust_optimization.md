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

The final documentation pass corrects current installation, Docker, GPU, quick-demo, dwell-time, and CLI guidance while retaining historical benchmark tables as historical references. All 25 Markdown files have valid local link targets and balanced code fences. The [project prompt history](clair3_rust_prompts.md) preserves all 34 user-authored prompts verbatim and in chronological order; automated messages and assistant or tool output are excluded.

Code cleanup removed the obsolete `unused_parens` suppression and 14 redundant parenthesized `break` expressions in `src/call/decode.rs`, plus the unreferenced `contains_index` helper in `csrc/fast_readselect.cpp`. No functional expression, model algorithm, compatibility option, fallback, or third-party source changed.

Native macOS formatting, Clippy with warnings treated as errors, and the unfiltered release test suite passed: 59 unit tests passed and one was ignored, with all integration targets passing their applicable tests. The explicit WhatsHap selector serialization regression passed. The phased HiFi fixture produced eight `chrToy` records with exact alleles and genotypes; its comparison correctly excludes the separate `chr:colon` fixture record because the command selects only `chrToy`.

A fresh supplied HG002 chr20 run completed in 62.23 seconds. All 135,180 pileup, 40,494 full-alignment, 94,551 merged, and 33,824 phased records exactly match the current native macOS Rust reference in full record content and ordering. This is a functional cleanup check, not a new optimization milestone or controlled benchmark. Earlier original-Python numerical differences and platform limitations documented above remain unchanged.

## Reproduction notes

- Keep output directories distinct when comparing implementations.
- Compare records by contig, position, reference, and alternate; do not rely only on row positions or record counts.
- Check genotype, QUAL, GQ, and phased GT/PS fields explicitly.
- Benchmark one configuration at a time on an otherwise quiet host. Several rejected changes looked beneficial only under different background load.
- The supplied benchmark is a HiFi CPU case on Apple Silicon. It is not a performance claim for other platforms, sequencing chemistries, or GPU backends.
