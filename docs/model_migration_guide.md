# Clair3 Model Migration Guide: TensorFlow to PyTorch

## Overview

Convert Clair3 variant-calling model checkpoints from TensorFlow to PyTorch `.pt` format.

## Environment Dependencies

```bash
pip install tensorflow-cpu torch numpy
```

The script also requires Clair3's internal modules (`shared.param_p`, `shared.param_f`, `clair3.model`). Run from the Clair3 project root is recommended.

## CLI Reference

**Script:** `convert_tf_checkpoint_to_torch.py`

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

Each checkpoint consists of an `.index` file and one or more `.data-*` shard files. The `--checkpoint` argument should point to the checkpoint prefix (e.g. `<model_dir>/pileup`), and the script resolves it automatically.

### 1. Single pileup model

```bash
cd Clair3_no_change_v2.0

python convert_tf_checkpoint_to_torch.py \
    --checkpoint /path/to/tf_models/r1041_e82_400bps_sup_v500/pileup \
    --output /path/to/output/pileup.pt \
    --model-type pileup \
    --platform ont
```

### 2. Single full-alignment model (with indel length)

```bash
python convert_tf_checkpoint_to_torch.py \
    --checkpoint /path/to/tf_models/r1041_e82_400bps_sup_v500/full_alignment \
    --output /path/to/output/full_alignment.pt \
    --model-type full_alignment \
    --platform ont \
    --add-indel-length
```

### 3. Batch-convert both models from one directory

```bash
python convert_tf_checkpoint_to_torch.py \
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
