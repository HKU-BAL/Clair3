# Correctness-first work package design

## Intent

Close the remaining correctness gaps in Clair3-Rust representation unification and model training. The first milestone is a per-chunk differential harness that identifies real implementation differences, followed by deterministic trainer parity and held-out model-quality validation. Broader validation must not certify the wrong labels or confuse successful file creation with scientific parity.

## Current state

- Native calling and the bounded chr21 training workflow run on macOS ARM64, Linux ARM64, Linux x86_64, and CUDA.
- The isolated baseline passes 63 unit tests and one ignored test, plus all applicable integration targets.
- Original `GetTruth` output and native training artifacts are compatible with original Python/PyTorch after the earlier fixes.
- Original `SplitExtendBed` omits the final merged interval. Rust includes it. That is an original bug and must not be reintroduced.
- Some representation-unification differences remain with identical haplotagged BAM and truth inputs. Their sources have not all been isolated.
- Training order, dropout randomness, and numerical backends differ. Existing checks prove working artifacts, not byte-identical datasets or weights.
- The unrelated publication-figures commit `b22b066` is present and must be preserved.

## Goals

1. Compare original and Rust unification per chunk with identical BAM, truth, reference, BED, and extended-BED inputs. Record candidate rows, alleles, haplotypes, matched pairs, and final records.
2. Classify and fix demonstrated first-party defects with focused regressions. Do not use changed thresholds to hide mismatches.
3. Compare pileup and full-alignment training using identical initial checkpoints, examples, batch order, and explicit test-only seeds. Check every full-alignment loss task, batch normalization, AdamW, and distributed gradient averaging.
4. Validate real training with independent validation data, convergence, early stopping, and checkpoint resume. Evaluate held-out chr20 F1 metrics against original training and the pretrained baseline.
5. Run applicable external-data tests for ONT, HiFi, Illumina, dwell-time, CRAM, GVCF, long-indel, and realignment workflows without restoring GitHub-based Linux testing.
6. Document every fixed difference and every residual limitation.

## Non-goals and anti-goals

- Do not reproduce the original final-interval omission merely to obtain matching counts.
- Do not resume the canceled 30-second optimization goal or start new performance experiments.
- Do not claim byte-identical weights or datasets across different random streams and numerical backends.
- Do not relax tolerances, delete failing assertions, or replace real validation with a synthetic success.
- Do not modify third-party source, WhatsHap's DP solver, or the original Clair3 checkout.
- Do not rewrite working unification, training, or model code merely to enlarge the change.

## Approaches

### Selected: unification first, then trainer parity

Start with the per-chunk differential harness. Fix defects at the responsible function, then use the validated examples and initial checkpoints to compare trainer behavior. This finds incorrect labels before optimization behavior is judged and keeps each change testable.

### Rejected: trainer-first comparison

Comparing losses and weights before unification is validated would mix label differences with optimizer and numerical differences. That makes root-cause isolation slower and can certify the wrong training data.

### Rejected: broad rewrite of unification and training

A wholesale rewrite would be larger than the observed defects justify and would make regression attribution difficult. Preserve existing architecture and make only changes required by differential failures.

## Constraints

- Work only in `.worktrees/clair3-rust-correctness` on `clair3-rust-correctness`. Keep `main` unchanged.
- Every production change starts with a failing test whose failure is observed before implementation. Run the relevant suite after each fix.
- Use identical input files for original and Rust comparisons. Do not change original code to manufacture agreement.
- Keep original final-interval behavior separate from genuine unification defects. The native implementation retains the interval.
- Use explicit seeds only inside the differential harness; normal training semantics, including dropout, remain unchanged.
- Tolerances must be justified by measured floating-point differences and must not hide biological or representation changes.
- Preserve validated calling behavior and the exact native macOS Rust HG002 reference.
- Use local Docker and the existing CUDA server for platform validation; no GitHub Linux workflow.
- Keep training logs, validation outputs, and investigation artifacts outside the commit unless a small fixture is required for a permanent regression.

## Deliverables

- A reusable per-chunk differential unification harness and focused regression tests.
- Minimal first-party fixes for demonstrated unification defects, each with a recorded root cause.
- Deterministic pileup and full-alignment trainer comparisons, including checkpoint reload and cross-Python compatibility checks.
- Real training and held-out evaluation results, including convergence, early stopping, and platform limitations.
- Updated training and verification documentation and a maintained chronological prompt history.

## Acceptance evidence

- Each targeted mismatch is either fixed with a regression or explained explicitly; no blanket "unisolated" conclusion is used for the targeted fixtures.
- Deterministic fixtures have expected record content, labels, archive schemas, and output ordering. Non-deterministic ordering is compared separately from biological content.
- Both pileup and full-alignment checkpoints load with Rust and original Python and run real calling workflows.
- Final native formatting, Clippy, all applicable tests, and HG002 record comparison pass. Any remaining difference is reported rather than presented as resolved.

## Remaining builder decisions

- Exact failing candidate and partition examples will be selected from the differential harness, including the already observed original-only and Rust-only records.
- The responsible function for each defect will be chosen from evidence at the corresponding pipeline boundary, not guessed from the symptom.
- Any numerical tolerance will be chosen from measured comparisons and explained in the verification report.
- Checkpoint resume means a weight-only restart with a fresh optimizer and epoch schedule, not exact optimizer, epoch, or RNG continuation.
- Held-out F1 is reported as evidence without an invented numerical acceptance threshold. Independent validation directories must not inherit the training contig exclusion, and empty validation batches must fail explicitly.
