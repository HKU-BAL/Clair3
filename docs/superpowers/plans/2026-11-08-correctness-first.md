# Correctness-first differential validation implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Identify and fix real unification and training defects with differential tests, then prove real model quality without conflating working artifacts with identical training data or weights.

**Architecture:** Compare original Python and Rust per chunk, move from candidate generation through partition matching and final records, and carry validated examples into deterministic trainer comparisons. Preserve the existing native architecture and use small, test-driven corrections at the responsible boundaries.

**Tech Stack:** Rust 2024 with Rust 1.88, LibTorch/tch 0.26 and PyTorch 2.13, HDF5, HTSlib, the original Python Clair3 reference, local Docker on macOS, and the existing NVIDIA CUDA server.

**Spec:** `docs/superpowers/specs/2026-11-08-correctness-first-design.md`

## Global Constraints

- Work only in `.worktrees/clair3-rust-correctness` on `clair3-rust-correctness`; keep `main` unchanged and preserve `b22b066`.
- Every production change begins with a failing test and an observed RED result. Do not write implementation before a test fails for the intended reason.
- Use identical BAM, truth, reference, BED, and extended-BED inputs for original and Rust comparisons; never modify the original implementation or third-party source.
- Keep the original final-interval omission separate from genuine defects. Rust must retain the final merged interval.
- Do not change normal dropout or random behavior for production. Any seed override or fixed batch ordering is limited to the test harness.
- Do not relax assertions, hide errors, substitute synthetic success for real workflows, or resume the canceled optimization goal.
- Preserve `TrainingOptions`, `TrainingSelectOptions`, existing CLI options, and all valid calling paths unless a demonstrated defect requires a change.
- Commit each completed, independently tested task before starting the next.
- Final verification includes the complete native suite, applicable external-data tests, exact HG002 record comparison, and documentation link/fence checks.

## Review Focus

- Final BED interval: it must not disappear or be clipped. Task 1 writes a fixture with adjacent intervals and asserts the final merged interval is emitted.
- Candidate representation and ordering: multi-allelic and nearby positions must preserve ref/alt alleles, GT/PS metadata, and record ordering. Task 1 compares complete chunk records, not counts alone.
- Haplotype-specific matching: reads with `0`, `1`, `2`, or missing phase information must not be assigned to the wrong haplotype. Task 2 inspects candidate/read mappings for a mismatching candidate.
- High-depth partition boundaries: candidate/truth partitions must not duplicate or drop positions when distance or combination limits change. Task 2 covers these partitions before fixing the responsible logic.
- Deterministic training and distributed reduction: identical batches and initial weights must not silently change through seed handling, batch normalization, task selection, or gradient reduction. Task 3 compares all tasks and all checkpoint tensors.
- Missing or malformed validation data: the harness must skip optional real-data tests explicitly or fail with a useful input error, not pass vacuously. Tasks 1, 3, and 4 check this behavior.

---

### Task 1: Differential unification harness

**Files:**
- Create: `tests/unification_parity.rs`
- Create: `tests/fixtures/data/final_interval.bed`

**Interfaces:**
- Consumes: `UnifyRepresentationArgs`, `split_extend_bed`, and the original `UnifyRepresentation.py` implementation used only as a validation reference.
- Produces: a reusable per-chunk comparison harness and `split_extend_bed` regression for the final interval.

- [ ] **Step 1: Write the failing tests**

Add an ignored real-data `chunked_unification_matches_original` test that runs all chunks with identical inputs and compares `(contig, position, ref, alt, genotype, GQ, DP, AF)` and PS when present in order, plus a non-ignored `split_extend_bed_keeps_final_interval` unit test. Compare actual metadata; do not add PS or other fields merely to make the harness fit. The unit fixture must contain adjacent intervals whose expansion merges them and a final disjoint interval; the expected output contains both. Optional real-data tests must report a named missing-input skip, never an empty pass.

- [ ] **Step 2: Verify RED and the baseline guard**

Run: `cargo test --test unification_parity split_extend_bed_keeps_final_interval --release --locked` and the ignored real-data test with the prepared original inputs.

Expected: the final-interval guard passes because the current implementation is already correct. The real-data test exposes specific chunk/candidate differences rather than only a different record count; that RED result drives Task 2. Do not change correct production code merely to force a failure.

- [ ] **Step 3: Implement the harness and final-interval guard**

Use `env!("CARGO_BIN_EXE_clair3")` to execute the native `UnifyRepresentation` command with the same arguments for every chunk. Run the original command only in the optional differential test. Keep original and Rust extended-BED generation separate when measuring the original final-interval defect. Retain the existing correct `split_extend_bed` behavior and assert that it writes the last current region after all input rows are consumed.

- [ ] **Step 4: Verify GREEN for the harness and baseline guard**

Run: `cargo test --test unification_parity split_extend_bed_keeps_final_interval --release --locked` and `cargo test --lib pipeline::training::tests::truth_and_bed_commands_work --release --locked`.

Expected: both PASS; the original final interval is retained and existing truth/BED behavior is unchanged. Task 2 owns the real-data mismatch fixes.

- [ ] **Step 5: Commit**

```bash
git add tests/unification_parity.rs tests/fixtures/data/final_interval.bed
git commit -m "Add differential unification harness and interval guard"
```

### Task 2: Resolve demonstrated unification defects

**Files:**
- Modify: `src/pipeline/unification.rs`
- Modify: `src/tensor/full_alignment.rs` (candidate-detail CIGAR decoding, after its own failing regression)
- Modify: `tests/unification_parity.rs`

**Interfaces:**
- Consumes: the per-chunk harness from Task 1 and the observed mismatching positions/candidates.
- Produces: focused regressions for the affected candidate, partition, allele, or haplotype paths.

- [ ] **Step 1: Add regressions for concrete mismatches**

Start with the observed original-only and Rust-only records, including `chr21:13251076` and `chr21:13048299`, plus the multi-allelic records at `chr21:13093614` and `chr21:13245144`. Assert an evidence-backed expected representation for the responsible candidate/partition, and assert relevant read haplotype and metadata mappings. Do not treat every original result as automatically correct, replace the real mismatch with a record-count assertion, or accept an explanation as a substitute for fixing a demonstrated first-party defect.

- [ ] **Step 2: Verify RED**

Run: `cargo test --test unification_parity <selected_mismatch_test> --release --locked`.

Expected: each new regression fails for a reproducible implementation reason and names the affected candidate or metadata.

- [ ] **Step 3: Isolate and fix the responsible behavior**

Trace the failing candidate through `candidate_rows`, allele decoding, `lock_variant`, `phased_genotype`, `partitions`, `find_match_pairs`, and `unify_regions`, as applicable. Compare the original reference implementation at that boundary. Make the smallest first-party change that corrects the demonstrated behavior without changing partition limits, hiding the candidate, or changing unrelated output.

- [ ] **Step 4: Verify GREEN**

Run: `cargo test --test unification_parity --release --locked -- --ignored` and `cargo test --lib --release --locked`.

Expected: the selected mismatch tests and existing library tests PASS; all four concrete records are fixed or have an evidence-backed explanation documented. The ignored real-data harness must be executed explicitly rather than assumed from the ordinary target run.

- [ ] **Step 5: Commit**

```bash
git add src/pipeline/unification.rs tests/unification_parity.rs
git commit -m "Fix unification mismatches found by differential testing"
```

### Task 3: Deterministic trainer parity

**Files:**
- Modify: `tests/training.rs`
- Modify: `src/train/trainer.rs`
- Modify: `src/train/data.rs`

**Interfaces:**
- Consumes: `TrainingOptions`, `BinArchive`, `BinWriter`, and the validated examples/checkpoints from Tasks 1 and 2.
- Produces: deterministic pileup/full-alignment loss and checkpoint comparisons, including distributed gradient reduction where available.

- [ ] **Step 1: Add failing deterministic comparisons**

For each model type, create the same initial checkpoint and labeled examples in both harnesses. Use `seed=Some(0)`, a fixed batch sequence, and `random_offset=0` in the test harness. Compare per-task losses and all checkpoint tensors after one optimizer step. For distributed mode, compare the distributed all-reduce result with the arithmetic mean of two independently reproduced rank-local gradients, not generally with a merged single-process batch whose batch normalization and dropout differ. Missing original/DDP support must produce an explicit skip, not a vacuous pass.

- [ ] **Step 2: Verify RED**

Run: `cargo test --test training deterministic_pileup_matches_original --release --locked` and `cargo test --test training deterministic_full_alignment_matches_original --release --locked`.

Expected: at least one comparison exposes a mismatch in task losses, tensor changes, or RNG/batch handling; if all comparisons already pass, inspect the test for an untested assumption rather than creating a failure.

- [ ] **Step 3: Correct the responsible trainer behavior**

Change only behavior shown to differ from the original: archive ordering, chunk offsets, task slicing, dropout/batch normalization, AdamW, or distributed averaging. Preserve the public `TrainingOptions` interface and ordinary random training behavior. Do not disable dropout or force a match by weakening tolerances.

- [ ] **Step 4: Verify GREEN**

Run: `cargo test --test training --release --locked` and `cargo test --test distributed --release --locked`.

Expected: all three tests pass. Checkpoints reload with both Rust and original Python, with finite outputs and the expected 24-value pileup tasks and 90-value full-alignment tasks.

- [ ] **Step 5: Commit**

```bash
git add tests/training.rs src/train/trainer.rs src/train/data.rs
git commit -m "Validate deterministic trainer parity"
```

### Task 4: Real training and external-data validation

**Files:**
- Modify: `src/train/trainer.rs`
- Modify: `tests/pipeline_parity.rs`
- Modify: `tests/calling_parity.rs`
- Modify: `tests/full_alignment_parity.rs`
- Modify: `scripts/run_quick_train_chr21.sh`
- Modify: `docs/rust_quick_test.md`

**Interfaces:**
- Consumes: the corrected training workflow and deterministic comparison helpers from Tasks 1–3, plus the `archives` validation-directory loader.
- Produces: platform-specific training/validation results, non-vacuous validation, and held-out model-quality comparisons.

- [ ] **Step 1: Add non-vacuous validation and real-data checks**

Add `validation_directory_preserves_independent_inputs`: a validation directory containing `_20_` must not lose its archives through the training-contig exclusion, and empty validation batches must produce a useful error rather than zero loss. Add real-data convergence, early-stopping, and weight-only checkpoint-restart checks for both models. Resume means loading weights and starting a fresh optimizer/epoch schedule, not exact optimizer/epoch/RNG continuation. Exercise applicable ONT, HiFi, Illumina, dwell-time, CRAM, GVCF, long-indel, and realignment workflows with explicit input-dependent skips. Add held-out chr20 F1 comparisons as evidence, not as an invented numerical acceptance threshold.

- [ ] **Step 2: Verify RED and execute optional evidence explicitly**

Run: `cargo test --test training validation_directory_preserves_independent_inputs --release --locked`, the ignored per-chunk differential test with `--ignored`, the two-process gradient comparison with `--ignored`, and the named external-data tests.

Expected: the validation-directory test fails because `_20_` is filtered incorrectly or validation is empty. The optional differential and distributed tests execute and report their actual result. No optional test is inferred from an ordinary target run.

- [ ] **Step 3: Correct validation behavior and run native macOS validation**

Keep the training `_20_` exclusion for training archives but do not apply it to independent validation archives. If an explicit validation directory yields no batches, return a clear error instead of reporting a zero validation loss. Then run `cargo test --all-targets --release --locked` and the supplied quick-training workflow with independent validation bins.

Expected: the new validation test passes; all applicable tests pass, losses are finite, and checkpoint/log/best-validation artifacts exist. No external-data test is silently converted to a pass without its named inputs.

- [ ] **Step 4: Run Linux and CUDA validation**

Run: the same bounded workflow in freshly rebuilt ARM64 and x86_64 Docker CPU images, then the CUDA image on the RTX 4090. Extend convergence evaluation on the CUDA host as resources permit. Use identical staged inputs and record native versus emulated timing separately.

Expected: all four platform workflows complete and both model checkpoints reload in Rust and original Python. Record archive shapes, label distributions, losses, epochs, and held-out F1 metrics rather than claiming identical weights. Record executed external-data tests and named skips explicitly.

- [ ] **Step 5: Verify calling and documentation behavior**

Run: `cargo test --all-targets --release --locked`, the explicit WhatsHap selector regression, the phased HiFi fixture, and the supplied HG002 command. Check all four VCF record sets against the current native macOS Rust reference, and check all Markdown links and fences.

Expected: applicable tests pass; all 135,180 pileup, 40,494 full-alignment, 94,551 merged, and 33,824 phased records match exactly. Do not authorize out-of-scope calling changes or replace the reference to hide a difference.

- [ ] **Step 6: Commit**

```bash
git add src/train/trainer.rs tests/pipeline_parity.rs tests/calling_parity.rs tests/full_alignment_parity.rs scripts/run_quick_train_chr21.sh docs/rust_quick_test.md
git commit -m "Validate real training and external-data workflows"
```

### Task 5: Final report, review, and integration preparation

**Files:**
- Modify: `docs/rust_optimization.md`
- Modify: `docs/clair3_rust_prompts.md`
- Modify: `docs/pileup_training.md`
- Modify: `docs/full_alignment_training_r1.md`
- Modify: `tests/unification_parity.rs` (Clippy-only type alias; preserve all assertions and behavior)

**Interfaces:**
- Consumes: completed task logs and measured results from Tasks 1–4.
- Produces: a reviewable final diff and an evidence-backed verification report.

- [ ] **Step 1: Update the prompt record and report**

Append `Clean up docker images and builds`, `Suggest what else tasks can be done to improve the comprehensiveness and quality of the project`, and `proceed with the work package` as prompts 36–38. Record every fixed mismatch, residual difference, test result, platform timing, and held-out metric. Do not call any unresolved difference resolved.

- [ ] **Step 2: Run final verification**

Run: `cargo fmt --all -- --check`, `cargo clippy --all-targets --locked -- -D warnings`, `cargo test --all-targets --release --locked`, the selector regression, `bash -n scripts/run_quick_train_chr21.sh`, documentation checks, and `git diff --check`.

Expected: all commands exit 0; no test failures are hidden and no generated artifacts enter the commit.

- [ ] **Step 3: Commit**

```bash
git add docs/rust_optimization.md docs/clair3_rust_prompts.md docs/pileup_training.md docs/full_alignment_training_r1.md
git commit -m "Document correctness validation and residual differences"
```

- [ ] **Step 4: Review the complete branch**

Run the prescribed whole-branch review against the plan and spec. Address Critical and Important findings with new RED→GREEN tests, then rerun the complete suite before asking the user to choose the integration method. Do not merge or push without that choice.
