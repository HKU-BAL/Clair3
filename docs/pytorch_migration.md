# Clair3 TensorFlow → PyTorch Migration Notes

_Last updated: 2026-01-07_

Status: core training and inference now run on PyTorch, and HDF5 (h5py + hdf5plugin) replaces PyTables for bin I/O. The sections below reflect the original TensorFlow audit and should be refreshed as we finish cleanup.

This document captures the current TensorFlow (TF) footprint inside the Clair3
codebase and enumerates the concrete steps required to complete the PyTorch
migration outlined in the roadmap. It should remain the canonical reference
while we touch code, environments, and documentation.

---

## 1. Current TensorFlow Surface Area

| Area | File(s) | TF Usage | Notes |
| --- | --- | --- | --- |
| Model definitions | `clair3/model.py` | `tf.keras` layers, `tf.keras.Model`, `tf.keras.regularizers`, `tf.keras.layers.Softmax` | Implements both pileup (`Clair3_P`) and full-alignment (`Clair3_F`) networks plus helper blocks (ResNet-style convs, pyramid pooling). |
| Training entrypoint | `clair3/Train.py` | `tensorflow` + `tensorflow_addons` (Lookahead, RectifiedAdam, F1 metrics), `tf.keras.callbacks`, `tf.config.threading` | Entire training loop uses `tf.keras.Model.fit`. Also relies on TF graph warm-up and checkpoint formats (`.index`, `.data-0000x-of-00002`). |
| Direct inference | `clair3/CallVariants.py` | `tensorflow` session-level API (`tf.config.threading`, `tf.nn`, etc.) | Executes forward passes, softmaxes, probability post-processing, and GPU thread configuration. |
| CFFI GPU wrappers | `clair3/CallVariantsFromCffi.py`, `clair3/CallVariantsFromCffiGPU.py` | Imports TF to load checkpoints, configure GPUs, and run `model.predict_on_batch`. | Also toggles `TF_ENABLE_ONEDNN_OPTS`, `CUDA_VISIBLE_DEVICES`, and uses TF for multi-GPU virtualization. |
| Triton interface | `CallVariantsFromCffi.py` (Triton branch) | Feeds TF logits into Triton when `--use_triton_gpu` is enabled. | Needs TorchScript/ONNX equivalent once migration happens. |
| TFRecord tooling | `clair3/bin2tf.py` | Uses `tf.io.TFRecordWriter`, `tf.train.Example` to emit sharded datasets. | Entire script is TF-specific; Dataloader rewrite must replace this. |
| Parameter knobs | `shared/param_p.py`, `shared/param_f.py` | `tensorflow_threads`, training defaults sized for TF graph | Need torch-friendly knobs (e.g., `torch_num_threads`, `default_amp`). |
| CLI wrappers | `run_clair3.sh`, `scripts/clair3.sh`, `preprocess/CheckEnvs.py` | Expect TF checkpoints/env vars, warn about TF versions, pass `--tensorflow_threads` around | Bash + Python wrappers must surface PyTorch env validation, new CLI flags, and no longer assume TF artifacts. |
| Utilities/docs | `clair3/utils.py`, `docs/gpu_quick_start.md`, `README.md`, `colab/*.ipynb` | Explanations reference TensorFlow, TensorFlow-Metal, TF install commands, and GPU behavior. | All user-facing docs must shift to PyTorch instructions. |
| Docker/conda envs | `Dockerfile`, `Dockerfile.gpu`, README install guides | Base images (`tensorflow/tensorflow:2.15.0-gpu`), packages (`tensorflow==2.2.0`, `tensorflow-addons`, `tensorflow-cpu`), and Apple Silicon guidance rely on TF. | Need PyTorch 2.2 + CUDA 12.1 environment plus CPU/MPS support. |

### Other secondary references
- `preprocess/CreateTrainingTensor.py` still mentions `--tensorflow_threads` in comments.
- `private_script/unified_training_multisample.sh` depends on TF checkpoints (`*.index`/`.data` pairs) and env exports tuned for TF.
- `shared/utils.py` and `preprocess/CheckEnvs.py` emit warnings that explicitly name TensorFlow when validating user environments.
- Multiple unit/benchmark docs and release notes explicitly cite TensorFlow performance numbers.

---

## 2. Dependency Inventory

| Layer | Location | Current Packages | Migration Impact |
| --- | --- | --- | --- |
| Conda/mamba env | README Option 4, docs/gpu_quick_start.md | `tensorflow=2.15.0`, `tensorflow-addons`, `pypy3.6`, `python=3.9` | Replace with `python=3.11`, `pytorch=2.2.*`, `pytorch-cuda=12.1`, `torchvision`, `torchaudio`, `torchmetrics`, `onnx`, `tritonclient`, plus existing CLI deps (samtools, whatshap, parallel, pypy3.10). |
| Docker CPU | `Dockerfile` | `continuumio/miniconda3` base + conda TF install | Swap to mamba-based build that installs PyTorch CPU wheels and refreshed CLI deps. |
| Docker GPU | `Dockerfile.gpu` | `tensorflow/tensorflow:2.15.0-gpu`, pip TF extras | Change base to `pytorch/pytorch:2.2.2-cuda12.1-cudnn8-runtime`, install whatshap/cffi/longphase as before. |
| Apple Silicon | `docs/gpu_quick_start.md` | `tensorflow-macos`, `tensorflow-metal` instructions | Need equivalent `pip install torch==...` with `--index-url https://download.pytorch.org/whl/cpu` or `mps` guidance. |
| Colab notebooks | `colab/clair3_*_quick_demo.ipynb` | `pip install tensorflow==2.2.0`, `tensorflow-addons` | Update cells to install PyTorch + CLI deps, adjust narrative text, and revalidate outputs. |

---

## 3. Runtime & Checkpoint Expectations

1. **Checkpoint format** — TF stores pileup/full-alignment weights as `.index` plus `.data-0000x-of-00002`. Migration requires a converter that loads TF checkpoints and writes PyTorch `.pt` (state_dict) plus TorchScript/ONNX exports for inference. Existing users must be able to re-use downloaded models.
2. **Threading** — `CallVariants*` currently calls `tf.config.threading.set_intra_op_parallelism_threads(...)` and `tf.config.threading.set_inter_op_parallelism_threads(...)` while also exporting `OMP_NUM_THREADS`. Equivalent PyTorch knobs (`torch.set_num_threads`, `torch.set_num_interop_threads`) need to be wired through CLI flags.
3. **GPU binding** — GPU wrappers rely on TF to enumerate devices and split memory. We must replicate this logic with `torch.cuda.mem_get_info`, `torch.cuda.device_count`, or `nvidia-smi` parsing.
4. **Data serialization** — Training optionally converts `.bin` tensors into TFRecords to unblock tf.data streaming. PyTorch Datasets can read PyTables directly, so TFRecord generation either becomes obsolete or should emit `.pt`/NumPy shards.
5. **Metrics/optimizers** — `tensorflow_addons` provides Lookahead + RectifiedAdam + F1 metrics. PyTorch equivalents can be built using `torch.optim` + `torchmetrics` (F1) + community Lookahead wrappers or by implementing custom classes.
6. **CFFI path** — C implementations in `preprocess/` feed numpy arrays directly to TF. Torch inference should accept the same numpy buffers (converted to tensors) without rewriting the C layer.

---

## 4. Migration Tasks (Detailed)

### 4.1 Model & Checkpoint Port
- Re-implement `Clair3_P` and `Clair3_F` as `torch.nn.Module` classes in `clair3/model.py`.
- Factor out shared building blocks (`BasicConv2D`, `BasicBlock`, `PyramidPooling`, LSTM stacks) using Torch layers.
- Provide utility scripts:
  - `clair3/model_export.py` — loads `.pt` weights and exports TorchScript (`torch.jit.trace`) and optional ONNX for Triton.
  - `scripts/convert_tf_checkpoint_to_torch.py` — loads TF checkpoints (via `tensorflow` or `tf2onnx` in a temporary env) and emits `.pt` + TorchScript. This is needed during transition so users do not lose pre-trained models.

### 4.2 Training Stack
- Replace `tf.keras` training loop with a PyTorch loop that uses:
  - `torch.utils.data.Dataset` to stream from `.bin` (PyTables) directly without TFRecords.
  - Mixed precision via `torch.cuda.amp.autocast` when `--use_gpu` or `--enable_dwell_time`.
  - Optimizers: `torch.optim.AdamW` + optional Lookahead wrapper; scheduler replicating warm-up/warm-down.
  - Metrics via `torchmetrics.F1Score`.
- Maintain CLI parity (`Train.py` flags) so automation and scripts keep working.
- Update `private_script/unified_training_multisample.sh` to activate the new mamba env, pass PyTorch-specific knobs (AMP, DDP world size), and consume `.pt` checkpoints.

### 4.3 Inference & Pipeline
- Update `CallVariants.py` and `CallVariantsFromCffi*.py`:
  - Load TorchScript models once per process, pin them to CPU/GPU/Apple MPS according to CLI flags.
  - Replace `model.predict_on_batch(X)` with `with torch.inference_mode(): model(torch.from_numpy(X))`.
  - Mirror the current scaling logic (depth clipping, per-chromosome batching).
  - Expose new CLI switches: `--torch_num_threads`, `--precision {fp32,fp16,bf16}`, `--amp`.
- Rework GPU wrappers to use PyTorch for device enumeration and memory-awareness.
- Ensure Triton flows ingest ONNX exported from PyTorch rather than TensorFlow SavedModels.

### 4.4 Environment & Tooling
- Add `env/clair3-torch.yml` describing the mamba environment (Python 3.11, PyTorch 2.2, CUDA 12.1 toolkit, CLI deps).
- Update both Dockerfiles to install from that environment file (or mimic it) and clean up TF packages.
- Update `README.md`, `docs/gpu_quick_start.md`, and quick demos to show new install commands (`mamba env create -f env/clair3-torch.yml && conda activate clair3-torch`).
- Refresh Colab notebooks so they install PyTorch wheels and demonstrate Torch-based inference.

### 4.5 Documentation & Communication
- Document the migration story (this file) and summarize user-facing changes in `README.md` and release notes.
- Add a “Migration FAQ” section covering checkpoint conversion, reproducibility, and how to keep using TF if needed (short term).
- Update `docs/training_data.md`, `docs/full_alignment_training*.md`, etc., with new CLI examples referencing PyTorch.

### 4.6 Validation & CI
- Create `tests/test_model_parity.py` that loads a frozen TF checkpoint (via converter) and a PyTorch checkpoint to compare logits on deterministic tensors (tolerances per task).
- Build smoke tests for calling/training (subset of HG002) to run inside CI (CPU-only plus GPU job).
- Extend CI workflows to install PyTorch env, run linting, execute critical smoke tests, and, if GPUs unavailable, at least exercise CPU inference.

---

## 5. Risks & Mitigations

| Risk | Impact | Mitigation |
| --- | --- | --- |
| Numerical drift between TF and PyTorch models | Variant quality/regression | Implement parity tests, lock seeds, and use double precision comparisons during bring-up. |
| Checkpoint conversion failures | Users cannot reuse pre-trained models | Provide verified conversion scripts + ship newly trained PyTorch checkpoints. Allow fallback to TF for one release. |
| Performance regressions (GPU utilization, multi-threading) | Longer runtimes, unhappy users | Benchmark each platform (ONT, HiFi, Illumina) post-migration; tune `torch.set_num_threads`, use `torch.compile` or `channels_last` as needed. |
| Environment conflicts (CUDA vs. PyTorch versions) | Install failures | Pin versions in `env/clair3-torch.yml`, document `nvidia-smi` matrix, and add env validation logic in `CheckEnvs`. |
| Triton/ONNX integration gaps | GPU calling inside wrappers may break | Add ONNX export tests and ensure Triton model repos accept Torch outputs. |

---

## 6. Immediate Next Steps
1. **Finalize audit (this doc)** – keep it updated as we touch modules.
2. **Branch management** – create `feature/pytorch-migration`, tag current TF baseline, and enforce PR checkpoints per workstream.
3. **Kick off `torch-models` task** – start porting `Clair3_P/F` with parity unit tests before touching training/inference flows.
4. **Plan user communication** – draft migration announcement for release notes once PyTorch builds are ready.

Please keep this document in sync with the actual migration progress so new contributors can ramp quickly.

