# Reference baselines for Clair3-Rust parity testing

The parity tests compare the native Rust implementation against the original
Python Clair3. Those comparisons used to depend on a machine-local directory
(`target/legacy/baseline/...`) that only existed on the machine where the
baselines were first generated, so the tests were permanently skipped
everywhere else.

This document describes how the reference outputs are captured with the
unmodified reference implementation and committed to the repository, and how to
reproduce them.

## What is committed

`tests/fixtures/baselines/<platform>/` holds the outputs of the original Clair3
for the official quick-demo datasets:

| Path | Contents |
| --- | --- |
| `pileup.vcf.gz` | reference pileup calling result |
| `full_alignment.vcf.gz` | reference full-alignment calling result |
| `merge_output.vcf.gz` | reference merged result |
| `tmp/split_beds/chr20` | reference extended BED |
| `tmp/phase_output/phase_bam/chr20.bam` | reference phased BAM (Illumina only) |
| `tmp/phase_output/phase_vcf/chr20.vcf` | reference phase VCF |
| `tmp/phase_output/phase_vcf/phased_chr20.vcf.gz` | reference phased VCF |
| `tmp/full_alignment_output/candidate_bed/chr20.0_1` | reference candidate regions |
| `tmp/full_alignment_output/candidate_bed/qual` | reference quality thresholds |
| `tmp/full_alignment_output/full_alignment_chr20.0_1.vcf` | reference per-chunk full-alignment VCF |
| `tmp/merge_output/merge_chr20.vcf` | reference per-chunk merge VCF |
| `hifi-cffi.npy`, `hifi-cffi.npy.info` | reference full-alignment tensor from the reference C implementation (HiFi only) |

`manifest.json` records the reference image and its digest, the SHA-256 and
size of every input and every captured output, and the record count of each
VCF. It is the provenance record for the fixtures.

## Regenerating the baselines

```bash
scripts/capture_baselines.sh --platform ilmn
scripts/capture_baselines.sh --platform hifi
```

The script

1. prepares the demo inputs, reusing `target/test-data/<platform>` when it
   exists and otherwise downloading the upstream demo files unmodified from
   `https://www.bio8.cs.hku.hk/clair3/demo/quick_demo/`,
2. runs the published reference container (`hkubal/clair3:v2.0.0` by default)
   over the quick-demo region `chr20:100000-300000`,
3. installs the small reference outputs as committed fixtures, and
4. captures the reference full-alignment tensor for HiFi.

Useful options: `--image` selects another reference image, `--demo-dir`,
`--work`, and `--fixtures` relocate the scratch and fixture directories, and
`--verify-only` checks that the committed fixtures are complete.

The reference implementation is never modified: the container runs the
published image as-is and the demo data is downloaded unmodified.

### Reference image version

Use `hkubal/clair3:v2.0.0`. The `v2.0.5` image cannot complete the pipeline
because upstream marks the realigner's `free_memory` helper `static` in
`preprocess/realign/realigner.cpp`, so the symbol is not exported and the
Python wrapper fails to load it. This is an upstream packaging defect, not a
Clair3-Rust issue, and it is why the realignment step of the reference
pipeline fails in that image.

## Tests that use the committed baselines

| Test | Comparison |
| --- | --- |
| `tests/calling_parity.rs::official_illumina_calls_match` | pileup and full-alignment VCFs against the reference, including GT, GQ, DP, AF, and PS |
| `tests/calling_parity.rs::production_full_call_includes_flanking_candidates` | production calling against the reference full-alignment VCF |
| `tests/pipeline_parity.rs::candidate_and_merge_outputs_match` | candidate selection, heterozygous-SNP selection, and merging for Illumina and HiFi |
| `tests/full_alignment_parity.rs::official_hifi_phased_context_matches` | full-alignment tensor against the reference C implementation, byte for byte |
| `tests/realign_workflow_parity.rs::official_realigned_bam_matches` | realigned BAM against the committed reference BAM |

The tests skip with a named message when the external demo inputs (BAM,
reference FASTA, and the pretrained models) are not present. The expectations
themselves now live in the repository.

## Preparing unification inputs

The differential representation-unification harness compares the original
`UnifyRepresentation` with the native implementation on identical inputs. It
needs a truth variant file and an extended BED, which the reference
implementation produces. `scripts/prepare_unification_inputs.sh` materializes
the reference implementation from this repository's own history at a pinned
commit and generates those inputs, so no separate checkout is required:

```bash
scripts/prepare_unification_inputs.sh \
  --bam_fn target/test-data/illumina/HG003_chr20_demo.bam \
  --reference target/test-data/illumina/GRCh38_chr20.fa \
  --truth_vcf target/test-data/illumina/HG003_GRCh38_chr20_v4.2.1_benchmark.vcf.gz \
  --bed target/test-data/illumina/HG003_GRCh38_chr20_v4.2.1_benchmark_noinconsistent.bed \
  --contig chr20 \
  --output target/baselines/unification-demo
```

The script prints the environment variables for the harness:

```bash
export CLAIR3_ORIGINAL_ROOT=<output>/original-src
export CLAIR3_UNIFICATION_BAM=...
export CLAIR3_UNIFICATION_TRUTH=...
export CLAIR3_UNIFICATION_REFERENCE=...
export CLAIR3_UNIFICATION_BED=...
export CLAIR3_UNIFICATION_CONTIG=chr20
```

With those variables set, the differential runs anywhere the demo data is
available:

```bash
cargo test --test unification_parity chunked_unification_matches_original \
  --release --locked -- --ignored --exact --nocapture
```

The differential is expected to report the residual differences documented in
the [optimization and verification report](rust_optimization.md): the
reference applies its own indel allele-frequency threshold separately from
`--min_af`, and it emits alternate alleles in hash order, while the native
implementation uses the deterministic candidate order.
