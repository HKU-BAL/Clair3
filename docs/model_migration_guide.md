# Clair3 Model Migration Guide: TensorFlow to PyTorch

## Overview

Convert Clair3 variant-calling model checkpoints from TensorFlow to PyTorch `.pt` format.

## Environment Dependencies

The `clair3-rust` branch performs conversion natively. Install or build the Rust executable as described in [Installing Clair3-Rust](rust_installation.md); TensorFlow, Python, and PyTorch packages are not required at runtime. The build must retain the default `tf-checkpoint` Cargo feature.

Both Keras 2 (TensorFlow ≤ 2.15) and Keras 3 (TensorFlow ≥ 2.16) checkpoint layouts are supported. The converter reads the checkpoint tensors directly instead of relying on Keras `load_weights()`, so conversion is independent of the installed Keras version. If any weight cannot be matched, it stops with an error rather than writing an untrained `.pt` (see [#488](https://github.com/HKU-BAL/Clair3/issues/488)).

Run the command from the Clair3 project root, or use an installed `clair3` executable.

## CLI Reference

**Native command:** `clair3 ConvertCheckpoint` (also available as `clair3 convert-tf-checkpoint`)

### Single-Model Mode

| Argument | Required | Description |
|---|---|---|
| `--checkpoint` | Yes | TF checkpoint prefix or directory |
| `--output` | Yes | Output `.pt` file path |
| `--model-type` | Yes | `pileup` or `full_alignment` |
| `--platform` | No | `ont` (default) / `hifi` / `ilmn` |
| `--add-indel-length` | No | Include indel-length heads |
| `--enable-dwell-time` | No | Add dwell-time channel (full-alignment only) |

### Batch Mode

| Argument | Required | Description |
|---|---|---|
| `--checkpoint-dir` | Yes | Dir containing `pileup.*` and `full_alignment.*` checkpoint files |
| `--output-dir` | Yes | Output directory for `.pt` files |
| `--platform` | No | `ont` (default) / `hifi` / `ilmn` |
| `--pileup-prefix` | No | Checkpoint prefix for pileup (default: `pileup`) |
| `--fa-prefix` | No | Checkpoint prefix for full-alignment (default: `full_alignment`) |
| `--pileup-add-indel-length` | No | Indel-length heads for pileup |
| `--fa-add-indel-length` | No | Indel-length heads for full-alignment |
| `--enable-dwell-time` | No | Dwell-time channel (full-alignment only) |

## Usage Examples

Expected input directory structure — each model directory contains TF checkpoint files with `pileup` and `full_alignment` prefixes:

```
<model_dir>/
├── pileup.index
├── pileup.data-00000-of-NNNNN
├── full_alignment.index
└── full_alignment.data-00000-of-NNNNN
```

Each checkpoint consists of an `.index` file and one or more `.data-*` shard files. The `--checkpoint` argument should point to the checkpoint prefix (e.g. `<model_dir>/pileup`), and the native command resolves it automatically.

### 1. Single pileup model

```bash
clair3 ConvertCheckpoint \
    --checkpoint /path/to/tf_models/r1041_e82_400bps_sup_v500/pileup \
    --output /path/to/output/pileup.pt \
    --model-type pileup \
    --platform ont
```

### 2. Single full-alignment model (with indel length)

```bash
clair3 ConvertCheckpoint \
    --checkpoint /path/to/tf_models/r1041_e82_400bps_sup_v500/full_alignment \
    --output /path/to/output/full_alignment.pt \
    --model-type full_alignment \
    --platform ont \
    --add-indel-length
```

### 3. Batch-convert both models from one directory

```bash
clair3 ConvertCheckpoint \
    --checkpoint-dir /path/to/tf_models/r1041_e82_400bps_sup_v500 \
    --output-dir /path/to/output/r1041_e82_400bps_sup_v500 \
    --platform ont \
    --fa-add-indel-length
```

Output:
```
r1041_e82_400bps_sup_v500/
├── pileup.pt
└── full_alignment.pt
```

## Key Notes

- Pileup models do **not** use indel-length heads; full-alignment models **do**.
- `tests/fixtures/tensorflow/dense` is a tiny real TensorFlow 2.20 checkpoint generated with upstream `tensorflow-cpu` in an x86_64 Python 3.11 container. The non-ignored `restores_dense_checkpoint` test checks its exact dense-weight transposition, so checkpoint parsing no longer depends on a machine-local TensorFlow fixture.
