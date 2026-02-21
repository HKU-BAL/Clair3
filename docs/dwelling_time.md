# Dwelling Time Feature

_New in Clair3 v2.0_

## Overview

Clair3 v2.0 supports incorporating per-base dwell time from ONT sequencing to improve variant calling accuracy. When enabled, dwell time is extracted from the BAM `mv` (move table) tag and added as an additional channel to the full-alignment input tensor.

---

## The `mv` Tag

The `mv` tag is a move table produced by Dorado basecalling. It encodes a coarse sequence-to-signal mapping with the following format:

```
mv:B:c,[block_stride],[signal_block_move_list]
```

- `block_stride`: the number of raw signal samples per move table element (set to the basecalling model's input stride).
- `signal_block_move_list`: a list of `int8_t` values. `1` indicates a base boundary (a new base was emitted), `0` indicates a continuation of the current base.

For example, given `mv:B:c,5,1,0,1,0,1`:
- The block stride is `5`
- The move list `1,0,1,0,1` indicates three bases were emitted at the 1st, 3rd, and 5th strided blocks

Clair3 parses the `mv` tag to compute the number of strided signal blocks (dwell time) for each aligned base.

To include `mv` tags in your BAM, run Dorado with the `--emit-moves` flag. See the [Dorado move table documentation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/move_table/) for details.

---

## Prerequisites

- **BAM files with `mv` tags** produced by Dorado (`--emit-moves`).
- **Platform**: ONT only. Not supported for PacBio HiFi or Illumina.
- **C implementation**: must be enabled (the default). `--disable_c_impl` is incompatible with dwell time.

---

## Usage

Add the `--enable_dwell_time` flag to your Clair3 command:

```bash
python run_clair3.py \
  --bam_fn=input.bam \
  --ref_fn=ref.fa \
  --threads=8 \
  --platform="ont" \
  --model_path="${MODEL_PATH}" \
  --output=${OUTPUT_DIR} \
  --enable_dwell_time
```

GPU acceleration is also supported by adding `--use_gpu`.

---

## Verifying mv Tags

Before enabling dwell time, verify that your BAM file contains `mv` tags:

```bash
samtools view input.bam | head -5 | tr '\t' '\n' | grep "^mv:"
```

Expected output format:
```
mv:B:c,5,1,0,1,0,1,...
```

If no `mv` tags are found, re-run basecalling with Dorado using the `--emit-moves` flag.

---

## Quick Demo

See [ONT Dwelling Time Quick Demo](quick_demo/ont_mv_quick_demo.md).
