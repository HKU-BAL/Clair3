# Bioinformatics submission package

**Working title:** Modernising a long-read variant caller in Rust: preserving scientific behaviour and attributing performance in the Clair3 rewrite

**Target journal:** *Bioinformatics* (Oxford University Press)
**Recommended article type:** Original Paper
**Manuscript status:** Complete draft for author completion and submission
**Authors:** ⟨Author list and ORCIDs⟩
**Corresponding author:** ⟨name, e-mail⟩

This file contains the complete submission package in three parts:

- **Part A — Author notes (not part of the manuscript).** Journal fit, formatting decisions applied, and residual risks.
- **Part B — Manuscript.** Submit this portion (Parts B and C of the manuscript body), after resolving the ⟨placeholders⟩.
- **Part C — Submission package.** Cover letter, pre-submission checklist, the recommended accuracy experiment, and figure/table production notes.

**Estimated lengths:** main text ≈4,900 words (excluding figures, tables and references); structured abstract ≈148 words; 19 references; 7 tables; 2 figures.

**Resolving placeholders:** every ⟨⟨angle-bracket⟩⟩ marker denotes information the authors must supply (hardware models, accessions, DOIs, URLs, funding, contributions) or text that depends on the outcome of the recommended accuracy experiment in Part C.

---

# Part A — Author notes

## A.1 Journal fit and article-type decision

**Recommended article type: Original Paper.** *Bioinformatics* caps Original Papers at 7 journal-template pages (≈5,000 words excluding figures) and Application Notes at 4 pages; manuscripts exceeding the limit by ≥20% are routinely returned without review ([formatting summary](https://manusights.com/blog/bioinformatics-formatting-requirements)). This paper is a retrospective methodology study with a large evidence base, not a tool announcement, so an Application Note would be too short for the content; an Original Paper fits if the main text is held to ≈5,000 words with detail moved to supplementary material.

**What was changed in this rewrite relative to the earlier long-form draft**

| Aspect | Previous draft | This version |
| --- | --- | --- |
| Abstract | Unstructured, ~250 words | Structured with the journal's own headings: **Motivation / Results / Availability and Implementation / Contact / Supplementary information**, ≤150 words ([OUP instructions](https://academic.oup.com/bioinformatics/pages/instructions_for_authors)) |
| Framing | "Rustification" as the headline | Biological stakes first: variant calls and training data underpin population-scale and clinical interpretation; the systems rewrite is the vehicle |
| Length | ~6,300 words | ~4,900 words main text; tables carry the data, detail moved to Supplementary Material |
| Structure | Generic IMRaD | Journal order: Introduction → Materials and methods → Results → Discussion → Conclusions, then Availability, Data availability, CRediT, Competing interests, Funding, Supplementary data |
| References | Numbered, journal-agnostic | Numeric style as used by recent published papers, e.g. "Zheng Z, Li S, Su J, et al. … Nat Comput Sci 2022;2:797–803" (see a [2025 Bioinformatics article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12806030/)) |
| Validation | Software-engineering framing | Five-layer differential-validation framework, benchmark table including runtime, memory and hardware, as the journal expects |
| Availability | Internal repository mentions | Explicit Availability and implementation + Data availability statements, as the journal enforces at the desk stage ([submission guidance](https://manusights.com/blog/bioinformatics-submission-guide)) |

## A.2 Formatting facts that must be confirmed before submission

The journal's own author pages returned HTTP 403 to automated retrieval during preparation of this draft. The following were verified from indexed copies of the official instructions and from recently published articles; each should be re-checked on the OUP site before submission.

1. **Structured abstract headings.** The official instructions specify, for original papers, five headings: *Motivation*, *Results*, *Availability and Implementation*, *Contact*, *Supplementary information*, with some flexibility permitted. Many recent papers instead use *Background / Methods / Results / Conclusions* plus an *Availability* line. This draft uses the five official headings.
2. **Section order.** Title page → Structured abstract → Introduction → System and methods (Algorithm, Implementation) → Results → Discussion → References, with Funding stated separately and Supplementary data supplied as a separate file. This draft uses "Materials and methods" in place of "System and methods" and states explicitly that no new algorithm is proposed.
3. **Reference style.** Recent published *Bioinformatics* articles use numbered references in the form "Author AB, Author CD, et al. Title. Journal Year;Volume:pages". The 2026 OUP LaTeX template advertises an `author-date` bibliography option ([template notes](https://github.com/RagnarGrootKoerkamp/oxford-bioinformatics-template)), which conflicts with published practice. This draft follows the numbered style of recent articles; if the current submission sheet requires author–date, the conversion is mechanical in the `.bib` file.
4. **Length caps.** Original Papers: 7 journal-template pages (≈5,000 words, excluding figures). Application Notes: 4 pages (≈2,600 words). Manuscript length is enforced administratively, so compile the PDF and count pages before submission.
5. **Figure and table production.** Colour/halftone at 350 dpi, line art at 1,200 dpi; single-column width 86 mm, double-column 178 mm; tables embedded in the manuscript in editable form, numbered consecutively with descriptive titles.

## A.3 Editorial-risk items and the single most valuable addition

The journal requires that new methods be compared with the state of the art on real biological data, that machine-learning papers report performance on an independent test set, that software be freely available at a stable URL with documentation and test data, and that runtime and memory be reported with hardware specifications. It also requires a cover letter.

This manuscript satisfies the software-availability, benchmark, real-data and reporting requirements, but its *accuracy* evidence is limited to genotype-level equivalence of outputs. The single most valuable addition is therefore an independent accuracy evaluation, specified in Part C.4. It would allow the manuscript to replace a limitation with a positive, evidence-backed claim about biological equivalence.

A second, smaller risk is fit framing. *Bioinformatics* editors look for biological relevance to be visible early. Part B therefore states the biological stakes of variant-call stability in the first paragraph of the Introduction and in the Motivation section of the abstract, and Figure 1 presents the methodological contribution as the paper's central object rather than an implementation report.

---

# Part B — Manuscript

**Title:** Modernising a long-read variant caller in Rust: preserving scientific behaviour and attributing performance in the Clair3 rewrite

**Authors:** ⟨Author 1⟩^1, ⟨Author 2⟩^1, ⟨Author 3⟩^2
**Affiliations:** ^1⟨Department, Institution, City, Country⟩; ^2⟨Department, Institution, City, Country⟩
**Corresponding author:** ⟨Name⟩, ⟨email⟩, ORCID ⟨0000-0000-0000-0000⟩
**Category suggested on submission:** Sequence analysis

---

### Structured abstract

**Motivation.** Long-read variant callers underpin population-scale and clinical variant interpretation, so their outputs must remain scientifically trustworthy while performance and hardware support evolve. We report how a systems-language rewrite of Clair3 was managed to preserve behaviour.
**Results.** Differential validation uncovered nine compatibility defect classes, seven of which component tests could not detect. The first native implementation matched the original runtime (301.89 s versus 309.10 s), whereas profiling and pipeline restructuring reduced a pristine-dependency baseline from 191.09 s to 72.01 s with the upstream phasing solver unmodified. Variant keys, genotypes and all 33,824 phased fields were preserved across CPU and CUDA runs.
**Availability and Implementation.** Rust implementation, training workflow, containers and validation logs: ⟨URL⟩. Clair3 v2 (Python/PyTorch) remains supported.
**Contact.** ⟨email⟩
**Supplementary information.** ⟨URL⟩

**Keywords:** long-read variant calling; software maintenance; reproducibility; differential testing; performance engineering; phasing

---

## 1 Introduction

Long-read sequencing has moved from niche to routine, and variant callers built on it now sit inside pipelines that genotype tens of thousands of human genomes, feed genome-wide association studies, and support clinical interpretation of rare variants [1]. The consequence is that a caller's implementation is scientific infrastructure: a change in its numerical behaviour propagates into association statistics, into benchmark results reported by other groups, and potentially into individual diagnoses.

This creates a durable tension. Variant-calling software must be modernised—its dependencies evolve, its hardware targets multiply from Apple silicon to GPU servers, and its users expect faster turnaround—yet its scientific outputs must remain stable. A modernisation that quietly changes which variants are called, or how they are phased, is a regression even if every test passes.

Language-level rewrite is one response now available to bioinformatics, because Rust offers memory safety, explicit data layout, native interoperability and reproducible builds. It is also an expensive response: the caller Clair3 combines deep-learning inference, haplotype phasing, dynamic-programming optimisation, alignment I/O, reference handling and model training [1,2]. Every one of these layers can perturb the result.

We report a retrospective case study of rewriting Clair3 in Rust and, more generally, of what it takes to modernise scientific software without losing scientific control. We contribute: a five-layer differential-validation framework; a taxonomy of nine compatibility defect classes with the layer that detected each; a performance attribution analysis showing that the language change itself contributed almost nothing and that measurable gains came from profiling and pipeline restructuring; an explicit accounting of unresolved numerical differences; and ten practices that we argue transfer to other scientific rewrites.

We use the informal term *rustification* for the replacement of an implementation with a Rust implementation; we make no claim that the term is established in the literature.

---

## 2 Materials and methods

### 2.1 Case study, scope and the biological stakes

Clair3 combines a fast pileup model with a haplotype-aware full-alignment model applied to unresolved candidates [1]. Phasing is provided by WhatsHap [2] and, for high-coverage HiFi data, optionally by LongPhase [3]. Alignment data are handled through samtools/htslib [4,5]. The original implementation comprised Python pipeline code, a C tensor-generation library accessed through CFFI, and external command-line tools.

Two properties make Clair3 informative. Its outputs are consumed as measurements, so behavioural drift has scientific consequences rather than merely operational ones. And it is numerically heterogeneous: inference, phasing, optimisation, I/O and reference handling all contribute, so a rewrite must preserve the interaction of layers rather than any single layer.

The scope of the rewrite was deliberately bounded. Calling, tensor generation, representation unification, realignment, post-processing, checkpoint handling and training were reimplemented natively [6–12]. The WhatsHap DP solver, LibTorch inference kernels and alignment file formats were retained unmodified. No new variant-calling algorithm is proposed in this paper; the object of study is the engineering process and its evidence.

### 2.2 Modernisation strategy

The original implementation was preserved unmodified on the `main` branch throughout, with the rewrite on a dedicated branch. This provided an always-available differential oracle, allowed outputs to be compared under identical inputs and environments, and kept the question "did the rewrite change the science?" separate from "did the optimisation change the science?".

Work proceeded as vertical slices rather than layer-by-layer translation: calling, tensor generation, unification, realignment, post-processing and training were each ported, validated and merged independently. External processes were replaced only where the interface was well defined; scientific cores owned by other groups were retained.

### 2.3 Dependency and boundary policy

Third-party code was never patched in place. LibTorch, HTSlib [5], WhatsHap [2], LongPhase [3], Gloo, samtools [4] and bgzip/tabix were obtained from official distributions, pinned by version or commit and checksum-verified where practical. The Rust build resolved dependencies from a lockfile with locked versions, and the upstream DP solver was used without modification, with acceleration restricted to the tested dependency version and a stock path retained for others.

Where acceleration could not be reproduced exactly, an explicit opt-out to stock behaviour was retained and documented rather than emulated silently.

### 2.4 A layered differential-validation framework

Five layers were used, and Table 1 lists which defect classes each layer detects.

**Table 1.** Differential-validation layers and the defect classes they detect.

| Layer | Artifacts compared | Defect classes detected |
| --- | --- | --- |
| L1 Component | Tensors, decoded records, selectors, LSTM layers, checkpoints | Localised numerical and logic errors |
| L2 Backend golden | Original-Python outputs regenerated per numerical backend | Backend-specific arithmetic, tolerance calibration |
| L3 Full workflow | Complete calling workflow on real GIAB data | Argument models, launcher assumptions, value spellings, compression detection, ABI assumptions |
| L4 Artifact | Training archives and checkpoints consumed by original Python/PyTorch | Filter and format incompatibility |
| L5 Consumption | Artifacts used in downstream calling and inference | Unusable or silently wrong outputs |

Because the original implementation's high-depth sampling used the platform C library's `rand()`, a single cross-platform golden is not scientifically meaningful. We therefore generated original-Python goldens per platform and per backend. Tolerances were set from measurement rather than convenience: official-model parity remained at 1e-4 on CPU and 0.001 on CUDA; one macOS convolution comparison was widened to 1e-6 after observing a 5.96 × 10⁻⁷ difference, with 2 × 10⁻⁴ retained for other backends; and approximate transcendental functions producing 1.3–1.7 × 10⁻³ model error were rejected outright.

Comparisons used record keys (contig, position, reference, alternate) and separately assessed genotype, QUAL, GQ and phased GT/PS fields. Record counts alone were never treated as a parity criterion.

### 2.5 Performance methodology

Wall-clock time was measured for the complete command on the GIAB HG002 Revio chr20 dataset with eight requested threads, within a single measurement campaign on the hardware given in Table 2. Configurations were compared with one another rather than across campaigns. A configuration that changed upstream scientific behaviour was excluded from the deliverable and reported separately. Platform timings obtained under emulation or on differing hardware are labelled functional timings.

**Table 2.** Benchmark platform.

| Platform | CPU | RAM | Accelerators |
| --- | --- | --- | --- |
| macOS (native) | ⟨Apple silicon model⟩ | ⟨GB⟩ | — |
| Linux ARM64 (Docker, native host) | ⟨host CPU model⟩ | ⟨GB⟩ | — |
| Linux x86_64 (Docker, emulated) | ⟨host CPU model⟩ | ⟨GB⟩ | — |
| Linux x86_64 (remote, native) | ⟨server CPU model⟩ | ⟨GB⟩ | 1 × NVIDIA RTX 4090 (24 GB), driver 580.82.09 |

### 2.6 Training validation workflow

Model training was validated separately with a complete chr21 workflow covering haplotagging, downsampling, representation unification, tensor generation, bin merging, candidate calling and one-epoch fine-tuning of both models at two depth levels. The same region-specific BAM, truth set, reference and pretrained weights were used for the original reference and every Rust platform. The region chr21:13–14 Mb was selected because the truth BED contains no intervals in the initially suggested 1–2 Mb interval. Archives were inspected with the original Python/h5py stack and checkpoints with original PyTorch, and trained models were then used for downstream calling.

### 2.7 Threats to validity

This is a single-project retrospective and cannot isolate the effect of the language, since implementation, dependencies, concurrency structure and hardware changed together. Timings varied with concurrent validation work; configurations, not headline ratios, are reported. Several residual numerical differences have plausible but unconfirmed causes. Training outputs depend on sampling, dropout and ordering, so weight-level equality was neither expected nor claimed. MPS was validated for calling only; no AMD backend was implemented; native Linux x86_64 CPU performance was not measured on x86_64 hardware. Linux validation was executed locally in containers rather than in hosted CI on native runners, which we do not recommend.

---

## 3 Results

### 3.1 The expensive defects were compatibility defects, not algorithmic ones

Nine distinct compatibility defect classes were found; seven were detected only at L3 or later, and none by L1 component tests. Table 3 summarises them.

**Table 3.** Compatibility defect classes, resolutions and detection layer.

| Class | Symptom | Cause | Resolution | Layer |
| --- | --- | --- | --- | --- |
| Argument model | `CreateTrainingTensorDirect` rejected | Alias absent | Alias to the existing command | L3 |
| Argument model | `--bin_fn` demanded an unrelated output | `call_fn` and `bin_fn` shared the internal ID `output` | Distinct internal IDs | L3 |
| Launcher model | `python clair3.py …` failed | Shell-only launcher | Polyglot shell/Python launcher | L3 |
| Value model | `--add_indel_length True` rejected | Case-sensitive boolean parsing | Legacy case-insensitive parser | L3 |
| Data format | Gzip input without `.gz` rejected | Compression chosen by extension | Magic-byte detection | L3 |
| Data format | Original h5py could not read native archives | Blosc filter absent in target environment | Original LZF plus shuffle filters [15] | L4 |
| Evaluation order | Quality written into position on x86_64 | Unspecified order of three reads in one expression | Explicit sequential reads | L2/L3 |
| ABI | Selector failed to load | Old versus new C++ string ABI in dependency wheels | Build both; select by exported symbol | L3 |
| Numerical | ONT LSTM golden mismatch | Backend difference, not a rewrite error | Golden regenerated with original Python | L2 |

Two findings merit emphasis. First, the argument-identifier defect could not have been found by inspection of the algorithm: it arose from a data model shared between two command-line options. Second, the evaluation-order defect was compiler-dependent, so a single-platform test suite would have certified it as correct.

### 3.2 Performance attribution: the language was not the lever

**Table 4.** Measured runtime on the HG002 chr20 benchmark.

| Configuration | Wall time | Status |
| --- | ---: | --- |
| Original Python/CFFI implementation | 309.10 s | Baseline |
| First native Rust implementation | 301.89 s | Port only; no meaningful speedup |
| Best configuration with a custom DP solver | 42.62 s (43.19 s confirmation) | Excluded: altered upstream solver behaviour |
| Pristine dependencies, stock DP solver, no restored pipeline optimisations | 191.09 s | Reference for the effect of pipeline work |
| Final: stock DP solver with four restored pipeline optimisations | 72.01 s | Deliverable |

The first native implementation was within 2.3% of the original, which is the paper's central performance result: the language change itself contributed almost nothing. The final configuration reduced the pristine-dependency baseline by 62.3% (191.09 s → 72.01 s) while leaving the upstream DP solver untouched, through four changes confined to the pipeline: independent phase-block scheduling with concurrent workers calling the unmodified solver; background capture of alignment spans; native allele loading with exact-match shortcuts and bounded edit-distance workspaces; and a first-party C++ read selector reproducing the original selection behaviour against upstream headers.

Two classes of candidate speedup were rejected on scientific grounds: approximate transcendental functions (model error 1.3–1.7 × 10⁻³ against a 1e-4 tolerance) and any change to model precision or candidate selection. Others were rejected on measurement despite being plausible: larger inference batches (7.9 GB peak), profile-guided optimisation, architecture-specific flags, parallel forward/reverse LSTM evaluation, four rather than eight full-alignment workers, and custom Winograd convolution. Peak memory rose from 1.97 GB to 3.93 GB relative to the original — an explicit cost of the parallelism that delivered the speedup.

### 3.3 Scientific outputs were preserved at the granularity the evidence supports

**Table 5.** Comparison with the original implementation on the same platform.

| Comparison | Records (each) | Keys and order | Genotypes | QUAL differences | GQ differences | Phased GT/PS |
| --- | ---: | --- | --- | --- | --- | --- |
| macOS | 135,180 pileup; 40,494 full; 94,551 merged | Exact | Exact | 83 / 50 / 58 (max 0.07) | 1 / 0 / 0 | 33,824 exact |
| Linux ARM64 | 135,180; 40,494; 94,550 | Exact | Exact | 3 / 65 / 50 (max 0.11) | 0 / 3 / 3 | 33,824 exact |
| Linux x86_64 | 135,180; 40,493; 94,559 | Exact | Exact | 1 / 29 / 28 (max 0.11) | 0 / 0 / 0 | 33,824 exact |
| CUDA versus original x86_64 CPU | 135,180; 40,493; 94,559 versus 94,560 | One extra CUDA record | One difference | 111 / 10,215 / 4,764 | 0 / 115 / 60 | 33,824 exact |

A controlled comparison of 100 full-alignment records against the supplied checkpoint found a maximum absolute probability difference of 1.66893 × 10⁻⁶, consistent with the small QUAL and GQ differences on CPU. The CUDA comparison used the original CPU output, because no original-Python GPU full-pipeline run was available; it is not a GPU-versus-GPU equivalence claim.

Three distinct claims must therefore be kept apart, and only the first is supportable as "identical": variant sites, alleles, genotypes and record order were preserved; QUAL and GQ values were not identical; and on CUDA one genotype and one record representation differed. For downstream use, genotype-level stability and exact preservation of the 33,824 phased fields are the properties that matter most, since phasing errors propagate into haplotype-aware calling and into association analyses that treat haplotypes as units.

### 3.4 Cross-platform validation

**Table 6.** Platform validation.

| Platform | Mode | Outcome |
| --- | --- | --- |
| macOS ARM64 | Native | Full suite (63 unit tests, 1 ignored), full HG002 workflow, exact match to the current native reference |
| Linux ARM64 | Native container | Formatting, Clippy, full suite, phased fixture, selector ABI test, full HG002 workflow |
| Linux x86_64 | Emulated container | Full suite, phased fixture, full HG002 workflow; functional only, not native performance |
| Linux x86_64 + RTX 4090 | Native container, CUDA | CUDA-enabled suite, CUDA model parity, full HG002 workflow |

Two platform-specific behaviours required engineering rather than documentation. The x86_64 WhatsHap wheel exports the pre-C++11 string ABI while the ARM64 source build exports the new ABI, so the selector is built in both variants and the appropriate one is selected by inspecting symbols exported by the installed core. And per-platform goldens were required because original sampling depends on the platform C library.

### 3.5 The training workflow is functional and interoperable, but not identical

**Table 7.** Complete bounded training workflow.

| Implementation and platform | Device | Time | Pileup train / validation loss | Full-alignment train / validation loss |
| --- | --- | ---: | ---: | ---: |
| Original Python, Linux x86_64 | RTX 4090 | 97.61 s | 0.601942 / 0.804108 | 8.894713 / 0.540924 |
| Rust, native macOS ARM64 | CPU | 35.91 s | 0.552349 / 0.624374 | 10.663239 / 1.117319 |
| Rust, Linux ARM64 | CPU | 77.23 s | 0.412903 / 1.062350 | 10.244903 / 0.875078 |
| Rust, emulated Linux x86_64 | CPU | 93.62 s | 0.412619 / 1.062332 | 10.170233 / 0.875076 |
| Rust, Linux x86_64 | RTX 4090 (CUDA) | 31.70 s | 0.394398 / 1.062423 | 10.109503 / 0.875232 |

These are functional timings across different devices, not a controlled training benchmark. They establish that the workflow completes, produces finite losses, and yields usable artifacts.

Training exposed defects that calling validation could not. The original `GetTruth` writes gzip content without a `.gz` suffix, and the original `SplitExtendBed` does not write the final merged interval of a region; both became visible only when the complete training script ran end to end. Original Python/h5py subsequently read every native archive, original PyTorch loaded every checkpoint with finite tensors, and on each platform both trained models completed a downstream calling fixture.

Archive row counts nevertheless differed: the original produced 3,182/3,984 pileup and 3,203/3,501 full-alignment rows at the two depths; Rust macOS produced 3,411/5,278 and 3,299/3,486; Rust Linux and CUDA produced 3,411/4,235 and 3,299/3,539. One contributing cause is identified — the omitted final merged interval — but other representation-unification differences remained even with identical haplotagged input. We therefore claim functional equivalence and artifact interoperability, not identical training datasets or weights.

### 3.6 Maintenance artifacts that outlived the rewrite

Four artifacts were retained deliberately. A negative-results register records every rejected optimisation with its measurement and reason, including experiments faster in isolation but slower end to end. A requirements-provenance record preserves, verbatim and chronologically, the 35 user requests that defined the scope, providing an executable trace of why each capability exists. A layered documentation set separates installation, quick testing, GPU guidance, per-model training guidance and a verification report that records what was measured and what remains different. Finally, compatibility shims are narrow and documented: legacy options are accepted as no-ops, and accelerated paths have opt-outs so an unsupported configuration can be reverted to stock behaviour without a rebuild.

---

## 4 Discussion

### 4.1 What this case study shows

**Preserving behaviour is a verification programme, not a translation effort.** The dominant defect classes were interface-level: argument models, launcher invocation, value spellings, compression detection, serialization order and ABI assumptions. These are exercised by callers rather than by functions, which is precisely why component testing missed seven of the nine.

**Performance is attributable to structure, not to language.** The first native implementation was within 2.3% of the original. The defensible claims are configuration-specific: 72.01 s against a 191.09 s pristine-dependency baseline in the same campaign, and a 42.62 s configuration excluded from the deliverable because it modified the upstream solver. Claims such as "Rust is *n*× faster" are not supportable from these data.

**Dependency boundaries determine upgradeability.** Keeping the upstream DP solver unmodified meant inheriting its semantics and future upgrades for free, at the price of the excluded 42.62 s configuration. Every upstream modification we removed bought maintainability. Because the expensive interactions were in the orchestration layer, we still gained performance without touching the solver.

**Training must be validated as a production workflow.** The training script alone exposed five compatibility defects, and it is where the remaining dataset differences became visible. A rewrite validated only on the flagship calling command would have been certified incomplete.

### 4.2 Ten practices

Each practice below is stated as an actionable rule with the evidence motivating it.

1. **Write the behavioural contract before the first line of new code.** Enumerate interfaces, artifact formats, ordering guarantees, tolerances and supported platforms, and keep the original implementation executable. *Evidence: retaining the original made every comparison possible and prevented "the reference stopped working" from becoming an argument against change.*
2. **Modernise the boundaries, not the scientific cores.** Replace orchestration, data movement and integration; retain trusted algorithm implementations unless a measured and scientifically acceptable reason exists. *Evidence: excluding DP-solver modifications kept the 42.62 s result out of the deliverable while preserving upgradeability.*
3. **Validate in layers, and make the highest layer mandatory.** Component tests localise defects; complete workflows reveal them. *Evidence: seven of nine defect classes surfaced only at L3 or beyond.*
4. **Generate golden data per numerical backend from the original implementation.** Do not impose a cross-platform golden on a pipeline whose original behaviour is platform-dependent. *Evidence: original `rand()`-based sampling made a single golden scientifically meaningless.*
5. **Set tolerances from evidence and record the evidence.** *Evidence: a macOS tolerance was widened to 1e-6 only after observing 5.96 × 10⁻⁷; model parity stayed at 1e-4; approximate arithmetic was rejected at 1.3–1.7 × 10⁻³ error.*
6. **Investigate failing tests before touching them.** A failing assertion is a measurement. *Evidence: the ONT LSTM discrepancy was reproduced with original Python before attribution, and the assertion was never deleted.*
7. **Prefer structural equivalence to numerical approximation.** Remove process boundaries, parallelise independent work and eliminate redundant passes before approximating arithmetic. *Evidence: a 62.3% reduction was achieved without altering any numerical kernel.*
8. **Test against the installed dependency, not an idealisation of it.** *Evidence: the C++ string ABI differed between the x86_64 wheel and the ARM64 source build, and unspecified evaluation order produced a quality-in-position defect on one compiler only.*
9. **Treat training as a first-class workflow.** Validate data preparation, tensor generation, archive schema, loss finiteness, checkpoint reload and downstream use. *Evidence: the training script alone exposed five compatibility defects.*
10. **Publish negative results and residual differences at honest granularity.** *Evidence: exact keys and genotypes, non-identical QUAL, and one CUDA genotype difference are three different claims, and only the first is "identical".*

Two process-level recommendations follow from the case but were not fully implemented in it. First, institutionalise platform validation in hosted CI on native runners: local container execution does not substitute for scheduled CI, and emulation must never be reported as native performance. Second, separate scientific-accuracy evaluation from software-equivalence evaluation. Our equivalence evidence is strong; our accuracy evidence is limited to the equivalence of genotype-level outputs, and the accuracy of Rust-trained relative to original-trained checkpoints has not been measured.

### 4.3 When modernisation is appropriate

The evidence supports modernisation when at least one condition holds: a structural bottleneck dominates user experience; concurrency safety is a recurring defect source; the project must support multiple architectures with one codebase; native integration with C libraries dominates the design; or the existing implementation cannot be maintained because its behavioural dependencies are implicit. The evidence argues against it when the original performance is adequate, when the scientific core is large, subtle and owned elsewhere, when the output is already reproducible, or when no trustworthy reference exists for differential testing. A useful decision rule is to ask what can be proved afterwards that could not be proved before. For Clair3 the answer was a per-backend behavioural contract, an executable cross-platform workflow, and quantified performance attribution; without such an answer a rewrite is an act of faith.

### 4.4 Limitations

The study is retrospective, single-site and lacks a controlled language comparison. Several QUAL/GQ and representation-unification differences remain incompletely explained, and no original-Python GPU full-pipeline baseline exists. Training was validated for functional correctness and artifact interoperability rather than statistical equivalence of learned models, and the accuracy of Rust-trained checkpoints relative to original-trained checkpoints has not been evaluated on an independent truth set. The most valuable next experiments are an accuracy evaluation of both trained checkpoints against the GIAB truth set, a controlled benchmark of the final configuration on native Linux x86_64 hardware, and replication of the validation framework on a second scientific codebase to test whether the ten practices transfer.

---

## 5 Conclusions

Modernising scientific software in a systems language is usually described as a performance intervention. This case study suggests it should be described as a verification and maintenance intervention. The Rust rewrite of Clair3 produced its largest gains not from the language but from profiling and restructuring, and its most instructive failures were contractual rather than computational: mismatched argument models, launcher assumptions, value spellings, compression conventions, serialization order, ABI assumptions and filter formats. The practices that emerged are transferable — write the behavioural contract first, keep a trustworthy reference, validate up to complete real workflows including training, generate per-backend goldens from the original implementation, justify every tolerance from measurement, prefer structural to numerical optimisation, and publish negative results and residual differences at honest granularity. Scientific software earns trust by making its behaviour legible and its maintenance decisions reviewable, and a rewrite that does both is a contribution to the longevity of the science it serves.

---

### Availability and implementation

The Rust implementation of Clair3, the configurable end-to-end training workflow, container definitions and validation logs are available at:

- Source: ⟨https://github.com/HKU-BAL/Clair3⟩, branch `clair3-rust`, commit ⟨⟨commit⟩⟩
- Archived release (software, test data, documentation): ⟨https://doi.org/10.5281/zenodo.⟨⟨zenodo-doi⟩⟩⟩
- Containers: ⟨registry/image references with digests⟩ for CPU Linux x86_64 and ARM64, and CUDA on Linux x86_64

The implementation requires Rust 1.88 or newer, CMake 3.26 or newer, a C++20-capable compiler, samtools 1.15.1 or newer, bgzip, tabix, Python 3 with pip, and GNU Parallel for the training workflow. Core calling and training do not require Python at run time beyond the WhatsHap phasing dependency. The upstream Python/PyTorch implementation (Clair3 v2, archived at ⟨https://doi.org/10.5281/zenodo.15064816⟩) remains available and is the reference implementation used for all comparisons in this paper [17].

### Data availability

All datasets used are publicly available. The GIAB HG002 v5.0q small-variant truth VCF and benchmark BED are distributed by the NIST Genome in a Bottle reference materials programme [16]; the PacBio Revio HiFi HG002 alignments are publicly available under accession ⟨SRA accession⟩. The exact command lines, checksums of the reference and truth files, the region used for training validation, and the derived training archives are provided in Supplementary Material and in the archived artifact.

### Supplementary data

- **Supplementary Figure S1.** Full chronology of the validation ladder, with each of the nine defect classes placed at the layer that detected it, and the port/optimisation timeline with configuration identifiers.
- **Supplementary Table S1.** Per-platform comparison of variant records, genotypes, QUAL and GQ differences, with the maximum absolute probability difference over the controlled 100-record comparison.
- **Supplementary Table S2.** Complete register of rejected optimisations, with measurements and reasons.
- **Supplementary Table S3.** Training-archive schemas, row counts, label distributions and checkpoint tensor inventories for all platforms.
- **Supplementary Material S4.** End-to-end training workflow, environment specifications, and the 35 user requests that defined scope, preserved verbatim as a requirements-provenance record.

### Acknowledgements

⟨Acknowledge colleagues, infrastructure, and the reference implementation maintainers.⟩

### Funding

⟨Grant numbers and funders; or: "This work received no specific grant from any funding agency in the public, commercial, or not-for-profit sectors."⟩

### Competing interests

⟨None declared.⟩

### Author contributions

⟨CRediT: Conceptualisation—⟨⟩; Methodology—⟨⟩; Software—⟨⟩; Validation—⟨⟩; Investigation—⟨⟩; Data curation—⟨⟩; Writing (original draft)—⟨⟩; Writing (review and editing)—⟨⟩; Visualization—⟨⟩; Supervision—⟨⟩; Project administration—⟨⟩.⟩

### Ethics approval

Not applicable; no human participants or animal subjects were involved.

---

### Figures

**Figure 1.** The differential-validation ladder and the defect class first detected at each layer. Five layers (L1 component, L2 backend golden, L3 full workflow, L4 artifact, L5 consumption) run from unit-level checks to scientific-output-level validation; each of the nine compatibility defect classes is placed at the layer that first detected it and is colour-coded by defect category. Seven classes were detected only at L3 or beyond, and none at L1. Files: `figures/figure1_validation_ladder.pdf` and `.svg` (vector masters), with PNG at 600 and 1,200 dpi; 178 mm double column, 117 mm high.

**Figure 2.** Performance attribution for the HG002 chr20 benchmark. Bars show the original implementation (309.10 s), the first native port (301.89 s), the pristine-dependency configuration with the stock solver (191.09 s), and the final configuration with the stock solver and four restored pipeline optimisations (72.01 s); the excluded 42.62 s configuration with a modified DP solver is drawn separately and annotated as not shipped. The annotation between the first two bars isolates the effect of the language change (−2.3%), and the annotation between the third and fourth bars the effect of profiling and pipeline restructuring (−62.3%), with the upstream solver unmodified. Files: `figures/figure2_performance_attribution.pdf` and `.svg` (vector masters), with PNG at 600 and 1,200 dpi; 178 mm double column, 111 mm high.

### Tables

Tables 1–7 above are embedded in the manuscript, numbered consecutively, with descriptive titles; the journal requests tables in an editable format at the end of the main text. Table 2 (benchmark platform) and the hardware entries in Tables 5–7 require the ⟨placeholders⟩ to be filled from the archived run logs before submission.

### References

1. Zheng Z, Li S, Su J, Leung AWS, Lam TW, Luo R. Symphonizing pileup and full-alignment for deep learning-based long-read variant calling. Nat Comput Sci 2022;2:797–803. https://doi.org/10.1038/s43588-022-00448-7
2. Patterson M, Marschall T, Pisanti N, van Iersel L, Stougie J, Klau GW, Schönthuth A. WhatsHap: weighted haplotype assembly for future-generation sequencing reads. J Comput Biol 2015;22:498–509. https://doi.org/10.1089/cmb.2014.0157
3. Lin J-H, Chen L-C, Yu S-C, Huang Y-T. LongPhase: an ultra-fast chromosome-scale phasing algorithm for small and large variants. Bioinformatics 2022;38:1816–22. https://doi.org/10.1093/bioinformatics/btac058
4. Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics 2009;25:2078–9. https://doi.org/10.1093/bioinformatics/btp616
5. Bonfield JK, Marshall J, Danecek P, et al. HTSlib: C library for reading/writing high-throughput sequencing data. GigaScience 2021;10:giab007. https://doi.org/10.1093/gigascience/giab007
6. rust-htslib: Rust bindings for HTSlib. https://github.com/rust-bio/rust-htslib (accessed ⟨date⟩).
7. tch: Rust bindings for PyTorch/LibTorch. https://github.com/LaurentMazare/tch (accessed ⟨date⟩).
8. hdf5-metno: HDF5 bindings for Rust. https://crates.io/crates/hdf5-metno (accessed ⟨date⟩).
9. Sandve GK, Nekrutenko A, Taylor J, Hovig E. Ten simple rules for reproducible computational research. PLoS Comput Biol 2013;9:e1003285. https://doi.org/10.1371/journal.pcbi.1003285
10. Prlić A, Procter JB. Ten simple rules for the open development of scientific software. PLoS Comput Biol 2012;8:e1002802. https://doi.org/10.1371/journal.pcbi.1002802
11. Wilkinson MD, Dumontier M, Aalbersberg IJ, et al. The FAIR Guiding Principles for scientific data management and stewardship. Sci Data 2016;3:160018. https://doi.org/10.1038/sdata.2016.18
12. Wilson G, Bryan JD, Cranston K, et al. Good enough practices in scientific computing. PLoS Comput Biol 2017;13:e1005510. https://doi.org/10.1371/journal.pcbi.1005510
13. The Rust Programming Language. https://www.rust-lang.org (accessed ⟨date⟩).
14. NNPACK: neural network package. https://github.com/goupsiod/nnpack (accessed ⟨date⟩).
15. The HDF Group. HDF5 File Format Specification. https://support.hdfgroup.org/HDF5 (accessed ⟨date⟩).
16. NIST Genome in a Bottle. Genome in a Bottle reference materials. https://www.nist.gov/genomics-research-project/genome-assembly-1000g-reference (accessed ⟨date⟩).
17. Zheng Z, Li S, Su J, et al. Clair3 v1.0.11. Zenodo 2025. https://doi.org/10.5281/zenodo.15064816
18. Clair3-Rust: native Rust implementation of Clair3. ⟨repository URL⟩ ⟨commit/tag⟩ ⟨DOI⟩.
19. Zheng Z, Yu X, Lei L, et al. Clair3-RNA: a deep learning-based small variant caller for long-read RNA sequencing data. Nat Commun 2025;16:11553. https://doi.org/10.1038/s41467-025-67237-y

*(Ordering above follows first citation; renumber once the text is final. Add the local-Docker and CUDA validation reports, the upstream WhatsHap package checksum, and the pinned Gloo commit as Supplementary references rather than in-text citations.)*

---

# Part C — Submission package

## C.1 Cover letter (required by the journal)

> Dear Editor,
>
> We submit an Original Paper entitled "Modernising a long-read variant caller in Rust: preserving scientific behaviour and attributing performance in the Clair3 rewrite" for consideration in *Bioinformatics*.
>
> Long-read variant callers are used at population scale and in clinical interpretation, so the scientific consequence of an implementation change is rarely trivial. Yet modernising such software is unavoidable as dependencies, hardware and performance expectations evolve. Our paper asks a question the community faces repeatedly and answers it with a complete, auditable case: we rewrote the Clair3 calling, phasing and training pipelines in Rust and report what was required to preserve scientific behaviour.
>
> The contribution is methodological. We formalise a five-layer differential-validation framework and show that the expensive defects in such a rewrite are compatibility defects — argument models, launcher assumptions, compression conventions, serialization order and dependency ABIs — rather than algorithmic ones, with seven of nine classes invisible to component testing. We further show that the language change itself contributed almost no speedup (301.89 s versus 309.10 s), that measurable gains came from profiling and pipeline restructuring (191.09 s to 72.01 s with the upstream phasing solver unmodified), and that the fastest configuration had to be excluded because it altered upstream solver behaviour. We publish the negative results, the residual numerical differences and the training-data differences we could not fully explain.
>
> All code, containers, training workflow, validation logs and the underlying evidence are freely available at ⟨URL⟩ with an archival release at ⟨DOI⟩, and the software and test data are usable from the stated instructions on the platforms reported.
>
> We believe the paper suits *Bioinformatics* because it addresses a shared bottleneck in the field — the maintenance of widely used computational methods whose outputs are treated as measurements — and because it offers a reusable validation framework rather than an implementation report. ⟨Related preprint disclosure, if any.⟩
>
> Sincerely,
> ⟨Corresponding author⟩

## C.2 Pre-submission checklist for this manuscript

**Editorial-risk items (fix first)**

- [ ] Software available and installable at the stated URL at submission time, with a working example and documented outputs; archived version and test data at Zenodo or equivalent.
- [ ] Add the missing independent accuracy evaluation (Section C.4). This is the most likely reviewer objection.
- [ ] Fill every ⟨placeholder⟩: hardware models, RAM, SRA/benchmark accessions, commit/tag/DOI, contact e-mail, ORCIDs, funding, CRediT.
- [ ] Ensure the abstract's Availability line contains the working URL, as the journal checks it.

**Formatting and structure**

- [ ] Structured abstract with Motivation / Results / Availability and Implementation / Contact / Supplementary information, ≤150 words (current: ~148).
- [ ] Main text within 7 journal-template pages; convert to the OUP template and check the compiled page count (Figures 1–2 and Tables 1–7 are excluded from the word budget but consume page space — consider moving Tables 5 and 7 to Supplementary Material if over).
- [ ] Section order as drafted; tables embedded in the manuscript in editable form; figures at 350 dpi (colour) / 1,200 dpi (line art); single column 86 mm, double column 178 mm.
- [ ] Reference style: confirm numbered versus author–date against the current submission sheet; the list above follows recent published articles.
- [ ] Supplementary material submitted as a separate file, cross-referenced from the main text.

**Reproducibility and data**

- [ ] Data availability statement with accession numbers for the GIAB v5.0q truth set and the Revio HG002 alignments; note that the internal training inputs used here are equivalent to public accessions and state the mapping.
- [ ] Include the exact command lines, container digests, `Cargo.lock`, and pinned dependency versions (LibTorch 2.13, WhatsHap 2.8 with published checksum, pinned Gloo commit) in the archive.
- [ ] State the licence (e.g. BSD-3-Clause) explicitly in the Availability section.
- [ ] If the journal's AI-disclosure policy applies, add the standard statement for AI-assisted writing/analysis.

## C.3 The one experiment that converts equivalence into a biological claim

The paper currently proves that the rewrite preserves software behaviour, but not that it preserves *accuracy*. A single additional experiment would close that gap and answer the reviewer's first question.

1. **Calling accuracy on independent truth.** On the benchmark region and, ideally, a second independent GIAB sample, compare original-Python and Rust VCFs against the GIAB v5.0q truth set with hap.py (as used in the project's quick demos) and report SNP/INDEL precision, recall and F1 for each implementation, plus the direct concordance between them. Expected outcome: identical F1 to within rounding, with the observed 1–3 GQ differences having no effect.
2. **Trained-checkpoint accuracy.** Train a second checkpoint from the Rust pipeline and one from the original pipeline on the same chr21 data, call the same region with each, and compare both against the truth set. This directly addresses the requirement that machine-learning results be reported on an independent test set and replaces our current, weaker statement that the training workflows are "functionally equivalent and cross-compatible".
3. **Report the residual differences as biological quantities.** For the two documented residual differences (one CUDA genotype; the training representation-unification differences), report the allele-supporting read counts and the genotype likelihoods, and state explicitly whether either would change an interpretation at a locus of interest.

The resulting sentences can replace the current limitation wording with a positive, evidence-backed claim — for example: "Rust-called variants achieve SNP F1 *X* and INDEL F1 *Y* on GIAB HG002, indistinguishable from the original implementation, and checkpoints trained by either implementation yield equivalent accuracy on an independent region".

## C.4 Figure and table production notes

- **Figure 1** (validation ladder with defect classes) is the paper's argument in one image and should be the first figure; editors look for the contribution to be visible early. **Both figures are now rendered** at 178 mm double column and are provided as vector masters (PDF and SVG) with PNG at 600 and 1,200 dpi: `figures/figure1_validation_ladder.*` and `figures/figure2_performance_attribution.*`. They are generated by `scripts/make_paper_figures.py`, which reads the same numbers as `docs/rust_optimization.md`, so the figures can be regenerated from the archived artifact rather than edited by hand.
- **Figure 2** (performance attribution) should isolate the language effect from the restructuring effect; the 42.62 s bar is visually separated and annotated as not shipped, with the reason stated in the caption rather than only in the text.
- **Tables 3, 5, 7** are the evidentiary core. Consider moving Table 7 (training results) and the per-platform detail of Table 5 into Supplementary Material, retaining summary rows in the main text, if page count becomes tight.
- All numbers already exist in the repository's verification report (`docs/rust_optimization.md`), quick-test guide and the archived run logs, so figures and tables can be regenerated from the artifact rather than transcribed by hand — worth stating in the Supplementary Material for reproducibility.

## C.5 Source of every quantitative claim

| Claim | Repository source |
| --- | --- |
| Runtime milestones (309.10 s, 301.89 s, 42.62 s, 191.09 s, 72.01 s) | `docs/rust_optimization.md` → Benchmark table |
| Compatibility defect classes and detection layers | `docs/rust_optimization.md` → compatibility findings; commit history for `fb11a82` and `586a827` |
| Per-platform record counts, QUAL/GQ differences, phased-field equality | `docs/rust_optimization.md` → platform comparison tables |
| Batch-100 maximum probability difference (1.66893 × 10⁻⁶) | `docs/rust_optimization.md` → Linux/ARM64 comparison section |
| Test totals (63 unit tests, 1 ignored) | `target/tmp/training-full-native-tests.log` (archived run log) |
| Training timings, losses, archive rows, checkpoint validation | `docs/rust_optimization.md` → End-to-end model-training validation |
| Requirements provenance (35 requests) | `docs/clair3_rust_prompts.md` |
| Tooling, dependency versions and checksums | `docs/rust_installation.md`, `docs/gpu_quick_start.md`, `docs/rust_optimization.md` → dependency policy |

---

*Prepared with the *Bioinformatics* house style in view: structured abstract with Motivation/Results/Availability and Implementation/Contact/Supplementary information, ~5,000-word main text, numeric reference style, explicit availability and data statements, runtime and memory reporting with hardware, and real-data validation. Confirm the current guidelines on the OUP site before submission, since the journal's own author pages were not machine-readable during preparation.*
