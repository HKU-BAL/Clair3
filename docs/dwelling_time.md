# Dwelling Time Feature

_New in Clair3 v2.0_

## Overview

Clair3 v2.0 introduces **signal-aware variant calling** for Oxford Nanopore (ONT) data. During nanopore sequencing, each nucleotide produces an electrical signal whose duration (the "dwell time") varies depending on the base identity and its local sequence context. By incorporating per-base dwell time information into the variant calling model, Clair3 can leverage signal-level evidence to improve calling accuracy.

When enabled, the dwell time is extracted from the BAM `mv` (move table) tag and added as a 9th channel to the full-alignment input tensor (in addition to the existing 8 channels: reference base, alternative base, mapping quality, base quality, strand info, variant type, insert base, and phasing info).

---

## Prerequisites

1. **BAM files with `mv` tags** — The input BAM must contain move table tags produced by Dorado. The move table is a record of the model's base emissions in strided signal space and gives a coarse sequence-to-signal mapping.
   - **Dorado**: Use the `--emit-moves` flag to include move table tags in the output. See the [Dorado move table documentation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/move_table/) for details.

2. **Platform** — Dwelling time is supported for the `ont` platform only. It is **not** supported for PacBio HiFi (`hifi`) or Illumina (`ilmn`).

3. **C implementation** — The dwelling time feature requires the C implementation to be enabled (this is the default). If you have explicitly set `--disable_c_impl`, Clair3 will automatically re-enable it when `--enable_dwell_time` is used.

4. **Dwell-time-aware model** — For best results, use a pre-trained model that was trained with dwell time enabled. Standard models will still work but may not fully benefit from the additional signal information.

---

## Usage

### Variant Calling

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

### Variant Calling with GPU

Dwelling time is fully compatible with GPU-accelerated calling:

```bash
python run_clair3.py \
  --bam_fn=input.bam \
  --ref_fn=ref.fa \
  --threads=8 \
  --platform="ont" \
  --model_path="${MODEL_PATH}" \
  --output=${OUTPUT_DIR} \
  --enable_dwell_time \
  --use_gpu
```

### Training with Dwelling Time

To train a full-alignment model with dwell time support, add `--enable_dwell_time` to the training command. This adjusts the model input shape from `[89, 33, 8]` to `[89, 33, 9]` to accommodate the additional dwell time channel.

```bash
# Single-GPU training
python clair3.py Train \
  --bin_fn ${BINS_FOLDER_PATH} \
  --ochk_prefix ${OUTPUT_PREFIX} \
  --add_indel_length True \
  --random_validation \
  --platform ont \
  --learning_rate 1e-4 \
  --enable_dwell_time

# Multi-GPU training with DDP
torchrun --nproc_per_node=4 clair3.py Train \
  --bin_fn ${BINS_FOLDER_PATH} \
  --ochk_prefix ${OUTPUT_PREFIX} \
  --add_indel_length True \
  --random_validation \
  --platform ont \
  --learning_rate 1e-4 \
  --scale_lr \
  --enable_dwell_time
```

**Note**: Training data (bin files) must have been generated with dwell time tensors. When creating training tensors, the dwell time channel is automatically included if the source BAM contains `mv` tags and the `--enable_dwell_time` flag is passed to the tensor generation step.

---

## How It Works

### Signal Extraction

During nanopore sequencing, the raw electrical signal is segmented into per-base events by the basecaller. The `mv` tag in the BAM file encodes this segmentation as a **move table** — a record of the model's base emissions in strided signal space that gives a coarse sequence-to-signal mapping.

#### Move Table Format

The move table SAM/BAM tag has the following format:

```
mv:B:c,[block_stride],[signal_block_move_list]
```

- **`block_stride`**: An `int8_t` containing the number of source signal samples that each element in the move list corresponds to. This is set to the input stride of the basecalling model.
- **`signal_block_move_list`**: A list of `int8_t` values, each containing a single move table element. Each element corresponds to `block_stride` samples of the raw source signal. A value of `1` indicates a base boundary (a new base was emitted), and `0` indicates a continuation of the current base.

For example, given `mv:B:c,5,0,0,1,0,1`:
- The block stride is `5` (the first value)
- The remaining values `0,0,1,0,1` indicate that bases were emitted in the 3rd and 5th strided blocks
- Converting to signal space (`[0-4, 5-9, 10-14, 15-19, 20-24]`), the bases were emitted from the 10th-14th and 20th-24th signal samples respectively

For full details, see the [Dorado move table documentation](https://software-docs.nanoporetech.com/dorado/latest/basecaller/move_table/).

#### How Clair3 Uses the Move Table

Clair3 parses the `mv` tag to compute the number of strided signal blocks (dwell time) for each aligned base. The implementation:

1. Reads the `mv` tag from each alignment record
2. Skips the first value (block stride) and iterates over the move list
3. Counts the number of signal blocks between consecutive base boundaries to derive per-base dwell time
4. Handles strand orientation (reverses the array for reverse-strand reads)
5. Accumulates dwell time for inserted bases into the preceding aligned position
6. Stores the dwell time value in the full-alignment tensor as channel 9

### Tensor Representation

The standard full-alignment tensor has shape `[89, 33, 8]` where:
- `89` is the matrix depth (number of reads)
- `33` is the flanking window size (16bp on each side + 1 center position)
- `8` channels encode alignment features

With dwell time enabled, the tensor shape becomes `[89, 33, 9]`, where channel 9 contains the per-base dwell time values clipped to the range `[0, 127]` for int8 representation.

### Model Architecture

The `Clair3_F` (full-alignment) model automatically adjusts its input layer to accept either 8 or 9 channels based on the `input_channels` parameter. The dwell time channel is processed through the same ResNet-based convolutional architecture as the other channels, allowing the model to learn signal-level patterns that correlate with variant presence.

The pileup model (`Clair3_P`) is not affected by the dwell time feature, as dwell time information is only incorporated at the full-alignment stage.

---

## Verifying mv Tags

Before enabling dwell time, verify that your BAM file contains `mv` tags:

```bash
# Check for mv tags in the first few reads
samtools view input.bam | head -5 | tr '\t' '\n' | grep "^mv:"
```

Expected output format:
```
mv:B:c,5,0,0,1,0,1,...
```

The tag follows the format `mv:B:c,[block_stride],[signal_block_move_list]`:
- `mv:B:c` indicates a signed byte array
- The first value (`5` in this example) is the block stride — the number of raw signal samples per move table element
- The remaining values are the move table entries (`1` = base boundary, `0` = continuation)

If no `mv` tags are found, your BAM was produced without move table output. Re-run basecalling with Dorado using the `--emit-moves` flag.

---

## Limitations

- **Platform support**: Dwelling time is only supported for the `ont` platform with Dorado-basecalled reads. It is not supported for PacBio HiFi or Illumina data.
- **C implementation required**: The dwell time feature relies on the C implementation for efficient `mv` tag parsing. It cannot be used with `--disable_c_impl`.
- **BAM compatibility**: If the BAM file does not contain `mv` tags, the dwell time channel will contain zero values. This will not cause an error but may reduce accuracy compared to running without the flag.
- **Pileup model**: Dwell time is only used in the full-alignment stage. The pileup model is unaffected.
- **Model compatibility**: For optimal results, use a model that was trained with `--enable_dwell_time`. Using a standard model with dwell time input, or a dwell-time model without dwell time input, may produce suboptimal results.

---

## Quick Demo

For a step-by-step quick demo of signal-aware variant calling, see [ONT Dwelling Time Quick Demo](quick_demo/ont_mv_quick_demo.md).
