use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
    process::{Command, Output},
    sync::Mutex,
};

use clair3::{
    config::Platform,
    model::{Model, ModelKind, load_python_state_dict},
    tensor::labels::{labels_from_reference, labels_from_vcf_columns},
    train::{
        checkpoint::save_checkpoint,
        data::{BinArchive, BinWriter, TensorExample},
        trainer::{TrainingOptions, train},
    },
};
use serde::Deserialize;
use tch::{Device, Kind, Tensor};

const DETERMINISTIC_BATCH_SIZE: usize = 20;
const DETERMINISTIC_LOSS_TOLERANCE: f64 = 0.0;
const DETERMINISTIC_CHECKPOINT_TOLERANCE: f64 = 0.0;
static TORCH_TEST_LOCK: Mutex<()> = Mutex::new(());

fn directory() -> PathBuf {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/verify/training/smoke");
    fs::create_dir_all(&path).unwrap();
    path
}

#[test]
fn checkpoint_roundtrip_preserves_every_tensor() {
    let _guard = TORCH_TEST_LOCK.lock().unwrap();
    let path = directory().join("checkpoint.pt");
    tch::manual_seed(11);
    let model = Model::new(ModelKind::Pileup, Device::Cpu, Some(18), false, None).unwrap();
    save_checkpoint(&model.var_store, &path).unwrap();
    let loaded = Model::new(ModelKind::Pileup, Device::Cpu, Some(18), false, Some(&path)).unwrap();
    let expected = model.var_store.variables();
    let actual = loaded.var_store.variables();
    assert_eq!(expected.len(), actual.len());
    for (name, tensor) in expected {
        assert_eq!(
            tensor.eq_tensor(&actual[&name]).all().int64_value(&[]),
            1,
            "tensor {name}"
        );
    }
}

#[test]
fn optimizer_step_writes_a_loadable_checkpoint() {
    let _guard = TORCH_TEST_LOCK.lock().unwrap();
    let root = directory();
    let data = root.join("data");
    fs::create_dir_all(&data).unwrap();
    let mut writer = BinWriter::create(&data.join("toy.bin"), &[33, 18], true).unwrap();
    let examples: Vec<_> = (0..120)
        .map(|index| TensorExample {
            position: format!("chrToy:{}:{}", index + 1, "A".repeat(33)),
            matrix: (0..33 * 18)
                .map(|value| ((value + index) % 101) - 50)
                .collect(),
            alt_info: "20-RA 20".to_string(),
            labels: if index % 3 == 0 {
                labels_from_vcf_columns(&["chrToy", "1", "A", "C", "0", "1"]).unwrap()
            } else {
                labels_from_reference("A").unwrap()
            },
        })
        .collect();
    writer.append(&examples).unwrap();
    drop(writer);
    let options = TrainingOptions {
        pileup: true,
        add_indel_length: false,
        max_epochs: 1,
        mini_epochs: 1,
        batch_size: 10,
        chunk_size: 5,
        random_validation: true,
        learning_rate: 1e-4,
        seed: Some(7),
        log_path: root.join("training.log"),
        best_checkpoint: root.join("best.pt"),
        ..TrainingOptions::new(data, root.join("model"), Platform::Ont)
    };
    let report = train(&options).unwrap();
    assert_eq!(report.epochs, 1);
    assert!(report.final_loss.is_finite());
    assert_eq!(report.checkpoints.len(), 1);
    assert!(report.best_validation_loss.is_some());
    let model = Model::new(
        ModelKind::Pileup,
        Device::Cpu,
        Some(18),
        false,
        Some(&report.checkpoints[0]),
    )
    .unwrap();
    let input = Tensor::from_slice(&vec![10_i32; 2 * 33 * 18]).reshape([2, 33, 18]);
    let output = model.forward_tasks(&input, false);
    assert_eq!(output.len(), 2);
    assert_eq!(output[0].size(), vec![2, 21]);
    assert_eq!(output[1].size(), vec![2, 3]);
    assert_eq!(output[0].isfinite().all().int64_value(&[]), 1);
}

#[test]
fn full_alignment_fine_tuning_writes_a_loadable_checkpoint() {
    let _guard = TORCH_TEST_LOCK.lock().unwrap();
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/verify/training/full_smoke");
    fs::create_dir_all(&root).unwrap();
    let data = root.join("data");
    fs::create_dir_all(&data).unwrap();
    let mut writer = BinWriter::create(&data.join("toy.bin"), &[55, 33, 8], false).unwrap();
    let examples: Vec<_> = (0..120)
        .map(|index| TensorExample {
            position: format!("chrToy:{}:{}", index + 1, "A".repeat(55)),
            matrix: (0..55 * 33 * 8)
                .map(|value| ((value + index) % 101) - 50)
                .collect(),
            alt_info: "20-RA 20".to_string(),
            labels: if index % 3 == 0 {
                labels_from_vcf_columns(&["chrToy", "1", "A", "C", "0", "1"]).unwrap()
            } else {
                labels_from_reference("A").unwrap()
            },
        })
        .collect();
    writer.append(&examples).unwrap();
    drop(writer);

    tch::manual_seed(13);
    let initial = Model::new(ModelKind::FullAlignment, Device::Cpu, Some(8), true, None).unwrap();
    let pretrained = root.join("pretrained.pt");
    save_checkpoint(&initial.var_store, &pretrained).unwrap();
    let options = TrainingOptions {
        checkpoint: Some(pretrained),
        pileup: false,
        add_indel_length: true,
        max_epochs: 1,
        mini_epochs: 1,
        batch_size: 10,
        chunk_size: 5,
        random_validation: true,
        learning_rate: 1e-4,
        seed: Some(7),
        log_path: root.join("training.log"),
        best_checkpoint: root.join("best.pt"),
        ..TrainingOptions::new(data, root.join("model"), Platform::Hifi)
    };
    let report = train(&options).unwrap();
    assert_eq!(report.epochs, 1);
    assert!(report.final_loss.is_finite());
    assert!(report.final_loss > 0.0);
    assert_eq!(report.checkpoints.len(), 1);
    assert!(report.best_validation_loss.is_some());

    let model = Model::new(
        ModelKind::FullAlignment,
        Device::Cpu,
        Some(8),
        true,
        Some(&report.checkpoints[0]),
    )
    .unwrap();
    let input = Tensor::from_slice(&vec![10_i32; 2 * 55 * 33 * 8]).reshape([2, 55, 33, 8]);
    let output = model.forward_tasks(&input, false);
    assert_eq!(output.len(), 4);
    for (tensor, classes) in output.iter().zip([21, 3, 33, 33]) {
        assert_eq!(tensor.size(), vec![2, classes]);
        assert_eq!(tensor.isfinite().all().int64_value(&[]), 1);
    }
}

#[test]
fn validation_directory_preserves_independent_inputs() {
    let _guard = TORCH_TEST_LOCK.lock().unwrap();
    let root = tempfile::tempdir().unwrap();
    let data = write_deterministic_archive(root.path(), &[33, 18], true);
    fs::write(data.join("ignored_20_.bin"), "not an archive").unwrap();
    fs::write(data.join("excluded-sample.bin"), "not an archive").unwrap();
    let validation = root.path().join("validation");
    fs::create_dir_all(&validation).unwrap();
    fs::copy(data.join("toy.bin"), validation.join("bin_HG002_20_chr20")).unwrap();
    let options = TrainingOptions {
        pileup: true,
        add_indel_length: false,
        max_epochs: 1,
        mini_epochs: 1,
        batch_size: DETERMINISTIC_BATCH_SIZE,
        chunk_size: DETERMINISTIC_BATCH_SIZE,
        validation_dir: Some(validation),
        exclude_training_samples: ["excluded-sample".to_string()].into_iter().collect(),
        seed: Some(7),
        log_path: root.path().join("training.log"),
        best_checkpoint: root.path().join("best.pt"),
        ..TrainingOptions::new(data.clone(), root.path().join("model"), Platform::Ont)
    };
    let report = train(&options).unwrap();
    let validation_loss = report.best_validation_loss.unwrap();
    assert!(
        validation_loss.is_finite() && validation_loss > 0.0,
        "independent _20_ validation archive was lost: validation loss {validation_loss}"
    );
    assert!(options.best_checkpoint.is_file());
    assert!(options.log_path.is_file());

    let empty = root.path().join("empty-validation");
    fs::create_dir_all(&empty).unwrap();
    let error = train(&TrainingOptions {
        validation_dir: Some(empty.clone()),
        output_prefix: root.path().join("empty-model"),
        log_path: root.path().join("empty-training.log"),
        best_checkpoint: root.path().join("empty-best.pt"),
        ..options.clone()
    })
    .unwrap_err();
    assert!(
        error.to_string().contains("no validation batches")
            && error.to_string().contains(&empty.display().to_string()),
        "unexpected empty-validation error: {error}"
    );

    let missing = root.path().join("missing-validation");
    let error = train(&TrainingOptions {
        validation_dir: Some(missing),
        ..options
    })
    .unwrap_err();
    assert!(
        error
            .to_string()
            .contains("failed to list validation archive directory"),
        "validation path error used the wrong directory kind: {error}"
    );
}

#[derive(Debug, Deserialize)]
struct OriginalMetrics {
    task_losses: Vec<f64>,
    total_loss: f64,
}

const ORIGINAL_PARITY_HARNESS: &str = r#"
import json
from itertools import accumulate
import sys
import types

import numpy as np
import torch

original_root = sys.argv[1]
sys.path.insert(0, original_root)
sys.modules.setdefault("h5py", types.ModuleType("h5py"))
import clair3.Train as trainer
from clair3.model import Clair3_F, Clair3_P

input_path, label_path, initial_path, updated_path, metrics_path, reload_path = sys.argv[2:8]
kind = sys.argv[8]
weight_decay = float(sys.argv[9])
pileup = kind == "pileup"
model_class = Clair3_P if pileup else Clair3_F
add_indel_length = not pileup
inputs = np.load(input_path)
labels = np.load(label_path)
assert inputs.shape[0] == 20 and labels.shape == (20, 90)

device = torch.device("cpu")
torch.manual_seed(0)
model = model_class(add_indel_length=add_indel_length, input_channels=inputs.shape[-1])
torch.save(model.state_dict(), initial_path)
trainer._load_checkpoint(model, initial_path, device)
torch.manual_seed(0)
optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=weight_decay)
label_shapes = [21, 3] if pileup else [21, 3, 33, 33]
label_shape_cum = list(accumulate(label_shapes))
recorded_losses = []

class RecordingFocalLoss(trainer.FocalLoss):
    def forward(self, y_true, y_pred):
        losses = super().forward(y_true, y_pred)
        recorded_losses.append(float(losses.mean().detach().cpu()))
        return losses

loss_funcs = [
    RecordingFocalLoss(label_shape_cum, task, None)
    for task in range(len(label_shapes))
]
total_loss = trainer._run_epoch(
    model,
    [(inputs, labels)],
    optimizer,
    loss_funcs,
    label_shapes,
    device,
    train=True,
    distributed=False,
)
assert len(recorded_losses) == len(label_shapes)
assert np.isfinite(recorded_losses).all() and np.isfinite(total_loss)
torch.save(model.state_dict(), updated_path)

reloaded = model_class(add_indel_length=add_indel_length, input_channels=inputs.shape[-1])
trainer._load_checkpoint(reloaded, updated_path, device)
reloaded.eval()
with torch.no_grad():
    outputs = reloaded(torch.from_numpy(inputs[:2]))
output = torch.cat(outputs, dim=1).cpu().numpy()
expected_width = 24 if pileup else 90
assert output.shape == (2, expected_width)
assert np.isfinite(output).all()
np.save(reload_path, output)
metrics = {
    "task_losses": recorded_losses,
    "total_loss": total_loss,
    "random_offset": 0,
    "batch_count": 1,
}
with open(metrics_path, "w") as output_file:
    json.dump(metrics, output_file)
print(json.dumps(metrics))
"#;

fn original_support(test_name: &str) -> Option<(PathBuf, PathBuf)> {
    let python =
        PathBuf::from("/Users/Aqua/Temp/Clair3/target/tmp/macos-validation-venv/bin/python");
    let original = PathBuf::from("/Users/Aqua/Temp/Clair3/target/tmp/original-main");
    let missing: Vec<_> = [
        ("validation Python", python.clone()),
        ("original Clair3", original.join("clair3.py")),
        ("original trainer", original.join("clair3/Train.py")),
    ]
    .into_iter()
    .filter(|(_, path)| !path.is_file())
    .map(|(name, path)| format!("{name} ({})", path.display()))
    .collect();
    if !missing.is_empty() {
        eprintln!("SKIP {test_name}: missing {}", missing.join(", "));
        return None;
    }
    let probe = Command::new(&python)
        .args(["-c", "import numpy, torch"])
        .output()
        .unwrap();
    if !probe.status.success() {
        eprintln!(
            "SKIP {test_name}: original NumPy/PyTorch support unavailable: {}",
            String::from_utf8_lossy(&probe.stderr)
        );
        return None;
    }
    Some((python, original))
}

fn run_original(command: &mut Command) -> Output {
    let display = format!("{command:?}");
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "{display} failed\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn write_deterministic_archive(root: &Path, shape: &[usize], pileup: bool) -> PathBuf {
    let data = root.join("data");
    fs::create_dir_all(&data).unwrap();
    let mut writer = BinWriter::create(&data.join("toy.bin"), shape, pileup).unwrap();
    let matrix_size: usize = shape.iter().product();
    let examples: Vec<_> = (0..2 * DETERMINISTIC_BATCH_SIZE)
        .map(|index| {
            let columns: Vec<&str> = match index % 4 {
                0 => vec!["chrToy", "1", "A", "C", "0", "1"],
                1 => vec!["chrToy", "1", "A", "AG", "0", "1"],
                2 => vec!["chrToy", "1", "A", "ATG", "0", "1"],
                _ => vec!["chrToy", "1", "A", "A", "0", "0"],
            };
            TensorExample {
                position: format!(
                    "chrToy:{}:{}",
                    index + 1,
                    "A".repeat(shape[shape.len() - 2])
                ),
                matrix: (0..matrix_size)
                    .map(|value| ((value * 17 + index * 31) % 101) as i32 - 50)
                    .collect(),
                alt_info: "20-RA 20".to_string(),
                labels: if index % 4 == 3 {
                    labels_from_reference("A").unwrap()
                } else {
                    labels_from_vcf_columns(&columns).unwrap()
                },
            }
        })
        .collect();
    writer.append(&examples).unwrap();
    drop(writer);
    data
}

fn tensor_map(tensors: Vec<(String, Tensor)>) -> BTreeMap<String, Tensor> {
    tensors.into_iter().collect()
}

fn compare_tensor_maps(
    actual: &BTreeMap<String, Tensor>,
    expected: &BTreeMap<String, Tensor>,
    tolerance: f64,
    context: &str,
) -> f64 {
    assert_eq!(actual.len(), expected.len(), "{context} tensor count");
    let mut maximum = 0.0_f64;
    for (name, expected) in expected {
        let actual = actual
            .get(name)
            .unwrap_or_else(|| panic!("{context} missing {name}"));
        assert_eq!(actual.size(), expected.size(), "{context} shape for {name}");
        if actual.kind() == Kind::Int64 {
            assert_eq!(
                actual.eq_tensor(expected).all().int64_value(&[]),
                1,
                "{context} integer tensor {name}"
            );
            continue;
        }
        let difference = (actual - expected).abs().max().double_value(&[]);
        maximum = maximum.max(difference);
        assert!(
            difference <= tolerance,
            "{context} tensor {name} maximum absolute difference {difference} exceeds {tolerance}"
        );
    }
    maximum
}

fn tensors_changed(before: &BTreeMap<String, Tensor>, after: &BTreeMap<String, Tensor>) -> bool {
    after
        .iter()
        .any(|(name, after)| after.ne_tensor(&before[name]).any().int64_value(&[]) != 0)
}

fn focal_loss(labels: &Tensor, predictions: &Tensor) -> Tensor {
    let predictions = predictions.clamp(1e-9, 1.0 - 1e-9);
    let cross_entropy = -labels * predictions.log();
    let weights = (1.0_f64 - &predictions).square() * labels;
    (cross_entropy * weights)
        .sum_dim_intlist(&[-1i64][..], false, Kind::Float)
        .mean(Kind::Float)
}

#[allow(clippy::too_many_arguments)]
fn deterministic_matches_original(
    test_name: &str,
    kind: ModelKind,
    platform: Platform,
    pileup: bool,
    shape: &[usize],
    weight_decay: f64,
) {
    let _guard = TORCH_TEST_LOCK.lock().unwrap();
    let Some((python, original)) = original_support(test_name) else {
        return;
    };
    let temporary_root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/verify/training");
    fs::create_dir_all(&temporary_root).unwrap();
    let root = tempfile::Builder::new()
        .prefix(&format!("{test_name}-"))
        .tempdir_in(temporary_root)
        .unwrap();
    let root = root.path();
    let data = write_deterministic_archive(root, shape, pileup);
    let archive = BinArchive::open(&data.join("toy.bin")).unwrap();
    assert_eq!(
        archive.chunk_count(DETERMINISTIC_BATCH_SIZE, DETERMINISTIC_BATCH_SIZE),
        1
    );
    let input_channels = *shape.last().unwrap();
    let (input, labels) = archive
        .read_batch(0, DETERMINISTIC_BATCH_SIZE, input_channels)
        .unwrap();
    let input_path = root.join("input.npy");
    let label_path = root.join("labels.npy");
    input.write_npy(&input_path).unwrap();
    labels.write_npy(&label_path).unwrap();

    let initial = root.join("initial.pt");
    let original_updated = root.join("original.pt");
    let metrics = root.join("metrics.json");
    let original_reload = root.join("original-reload.npy");
    run_original(
        Command::new(python)
            .arg("-c")
            .arg(ORIGINAL_PARITY_HARNESS)
            .arg(original)
            .arg(&input_path)
            .arg(&label_path)
            .arg(&initial)
            .arg(&original_updated)
            .arg(&metrics)
            .arg(&original_reload)
            .arg(if pileup { "pileup" } else { "full_alignment" })
            .arg(weight_decay.to_string()),
    );
    let metrics: OriginalMetrics = serde_json::from_slice(&fs::read(&metrics).unwrap()).unwrap();
    let add_indel_length = !pileup;
    let model = Model::new(
        kind,
        Device::Cpu,
        Some(input_channels as i64),
        add_indel_length,
        Some(&initial),
    )
    .unwrap();
    let initial_tensors = tensor_map(load_python_state_dict(&initial).unwrap());
    compare_tensor_maps(
        &tensor_map(model.var_store.variables().into_iter().collect()),
        &initial_tensors,
        0.0,
        "initial checkpoint",
    );

    tch::manual_seed(0);
    let outputs = model.forward_tasks(&input, true);
    let offsets = [0, 21, 24, 57];
    let widths = [21, 3, 33, 33];
    let task_losses: Vec<_> = outputs
        .iter()
        .enumerate()
        .map(|(task, output)| {
            focal_loss(&labels.narrow(1, offsets[task], widths[task]), output).double_value(&[])
        })
        .collect();
    assert_eq!(task_losses.len(), metrics.task_losses.len());
    for (task, (actual, expected)) in task_losses.iter().zip(&metrics.task_losses).enumerate() {
        let difference = (actual - expected).abs();
        assert!(
            difference <= DETERMINISTIC_LOSS_TOLERANCE,
            "task {task} loss difference {difference} exceeds {DETERMINISTIC_LOSS_TOLERANCE}: rust={actual}, original={expected}"
        );
    }

    let options = TrainingOptions {
        checkpoint: Some(initial.clone()),
        pileup,
        add_indel_length,
        max_epochs: 1,
        mini_epochs: 1,
        batch_size: DETERMINISTIC_BATCH_SIZE,
        chunk_size: DETERMINISTIC_BATCH_SIZE,
        random_validation: false,
        learning_rate: 1e-3,
        seed: Some(0),
        log_path: root.join("training.log"),
        best_checkpoint: root.join("best.pt"),
        ..TrainingOptions::new(data, root.join("rust-model"), platform)
    };
    let report = train(&options).unwrap();
    assert_eq!(report.epochs, 1);
    assert_eq!(report.checkpoints.len(), 1);
    let total_difference = (report.final_loss - metrics.total_loss).abs();
    assert!(
        total_difference <= DETERMINISTIC_LOSS_TOLERANCE,
        "total loss difference {total_difference} exceeds {DETERMINISTIC_LOSS_TOLERANCE}: rust={}, original={}",
        report.final_loss,
        metrics.total_loss
    );

    let rust_updated = &report.checkpoints[0];
    let trained = Model::new(
        kind,
        Device::Cpu,
        Some(input_channels as i64),
        add_indel_length,
        Some(rust_updated),
    )
    .unwrap();
    let actual_tensors = tensor_map(trained.var_store.variables().into_iter().collect());
    let expected_tensors = tensor_map(load_python_state_dict(&original_updated).unwrap());
    assert!(tensors_changed(&initial_tensors, &expected_tensors));
    assert!(tensors_changed(&initial_tensors, &actual_tensors));
    let checkpoint_difference = compare_tensor_maps(
        &actual_tensors,
        &expected_tensors,
        DETERMINISTIC_CHECKPOINT_TOLERANCE,
        "updated checkpoint",
    );

    let (reload_input, _) = archive.read_batch(0, 2, input_channels).unwrap();
    let rust_reload = trained.forward(&reload_input, false);
    let expected_width = if pileup { 24 } else { 90 };
    assert_eq!(rust_reload.size(), vec![2, expected_width]);
    assert_eq!(rust_reload.isfinite().all().int64_value(&[]), 1);
    let original_reload = Tensor::read_npy(&original_reload).unwrap();
    assert_eq!(original_reload.size(), vec![2, expected_width]);
    assert_eq!(original_reload.isfinite().all().int64_value(&[]), 1);
    eprintln!(
        "{test_name}: task_losses={task_losses:?}, total_loss_difference={total_difference:.3e}, maximum_checkpoint_difference={checkpoint_difference:.3e}"
    );
}

#[test]
fn deterministic_pileup_matches_original() {
    deterministic_matches_original(
        "deterministic_pileup_matches_original",
        ModelKind::Pileup,
        Platform::Ont,
        true,
        &[33, 18],
        1e-4,
    );
}

#[test]
fn deterministic_full_alignment_matches_original() {
    deterministic_matches_original(
        "deterministic_full_alignment_matches_original",
        ModelKind::FullAlignment,
        Platform::Hifi,
        false,
        &[55, 33, 8],
        1e-7,
    );
}

fn real_training_inputs(test_name: &str, kind: &str) -> Option<(PathBuf, PathBuf, PathBuf)> {
    let Some(root) = std::env::var_os("CLAIR3_REAL_TRAINING_ROOT").map(PathBuf::from) else {
        eprintln!("SKIP {test_name}: CLAIR3_REAL_TRAINING_ROOT is not set");
        return None;
    };
    let train = root.join(kind).join("train");
    let validation = root.join(kind).join("validation");
    let checkpoint = root.join("models").join(if kind == "pileup" {
        "pileup.pt"
    } else {
        "full_alignment.pt"
    });
    let missing: Vec<_> = [train.clone(), validation.clone(), checkpoint.clone()]
        .into_iter()
        .filter(|path| !path.exists())
        .collect();
    if !missing.is_empty() {
        eprintln!(
            "SKIP {test_name}: missing {}",
            missing
                .iter()
                .map(|path| path.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
        return None;
    }
    Some((train, validation, checkpoint))
}

fn real_archive_metrics(test_name: &str, path: &Path, input_channels: usize) {
    let mut archives: Vec<_> = fs::read_dir(path)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.is_file())
        .collect();
    archives.sort();
    assert!(!archives.is_empty(), "{test_name}: no archives in {path:?}");
    for path in archives {
        let archive = BinArchive::open(&path).unwrap();
        let rows = archive.len();
        assert!(rows > 0, "{test_name}: empty archive {path:?}");
        let mut label_counts = Tensor::zeros([90], (Kind::Float, Device::Cpu));
        for start in (0..rows).step_by(200) {
            let (_, labels) = archive
                .read_batch(start, (start + 200).min(rows), input_channels)
                .unwrap();
            label_counts += labels.sum_dim_intlist(&[0i64][..], false, Kind::Float);
        }
        let label_counts = Vec::<f32>::try_from(&label_counts).unwrap();
        let expected_nonzero = rows * 4;
        assert_eq!(
            label_counts.iter().sum::<f32>() as usize,
            expected_nonzero,
            "{test_name}: labels in {path:?}"
        );
        assert!(
            label_counts.iter().filter(|count| **count > 0.0).count() > 2,
            "{test_name}: labels in {path:?} are not representative"
        );
        eprintln!(
            "{test_name}: archive={} shape={:?} rows={} label_counts={label_counts:?}",
            path.display(),
            archive.tensor_shape(),
            rows
        );
    }
}

fn training_log(path: &Path) -> Vec<(f64, f64)> {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some("epoch\tloss\tval_loss"));
    lines
        .enumerate()
        .map(|(index, line)| {
            let fields: Vec<_> = line.split('\t').collect();
            assert_eq!(fields.len(), 3, "malformed training log line: {line}");
            let epoch = fields[0].parse::<usize>().unwrap();
            assert_eq!(epoch, index + 1);
            let loss = fields[1].parse::<f64>().unwrap();
            let validation = fields[2].parse::<f64>().unwrap();
            assert!(loss.is_finite() && loss > 0.0);
            assert!(validation.is_finite() && validation > 0.0);
            (loss, validation)
        })
        .collect()
}

fn real_training_workflow(test_name: &str, kind: ModelKind, directory: &str, pileup: bool) {
    let _guard = TORCH_TEST_LOCK
        .lock()
        .unwrap_or_else(|poison| poison.into_inner());
    let Some((training, validation, checkpoint)) = real_training_inputs(test_name, directory)
    else {
        return;
    };
    let input_channels = if pileup { 18 } else { 8 };
    real_archive_metrics(test_name, &training, input_channels);
    real_archive_metrics(test_name, &validation, input_channels);
    let output = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("target/verify/training/task4-real")
        .join(format!("{directory}-{}", std::process::id()));
    fs::create_dir_all(&output).unwrap();
    let options = TrainingOptions {
        checkpoint: Some(checkpoint),
        pileup,
        add_indel_length: !pileup,
        max_epochs: 3,
        mini_epochs: 1,
        batch_size: 1000,
        chunk_size: 200,
        learning_rate: 1e-5,
        validation_dir: Some(validation.clone()),
        seed: Some(19),
        log_path: output.join("convergence.log"),
        best_checkpoint: output.join("convergence-best.pt"),
        ..TrainingOptions::new(training.clone(), output.join("convergence"), Platform::Hifi)
    };
    let report = train(&options).unwrap();
    assert_eq!(report.epochs, 3);
    assert_eq!(report.checkpoints.len(), 3);
    assert!(options.best_checkpoint.is_file());
    let losses = training_log(&options.log_path);
    assert_eq!(losses.len(), 3);
    assert!(
        losses.iter().skip(1).any(|(loss, _)| *loss < losses[0].0),
        "{test_name}: training loss never improved: {losses:?}"
    );
    eprintln!("{test_name}: convergence losses={losses:?}");

    let restart_checkpoint = report.checkpoints.last().unwrap().clone();
    let restart = TrainingOptions {
        checkpoint: Some(restart_checkpoint.clone()),
        max_epochs: 1,
        seed: Some(23),
        output_prefix: output.join("restart"),
        log_path: output.join("restart.log"),
        best_checkpoint: output.join("restart-best.pt"),
        ..options.clone()
    };
    let restarted = train(&restart).unwrap();
    assert_eq!(restarted.epochs, 1);
    assert_eq!(restarted.checkpoints.len(), 1);
    assert!(
        restarted.checkpoints[0]
            .file_name()
            .unwrap()
            .to_string_lossy()
            .ends_with(".01.pt")
    );
    assert_eq!(training_log(&restart.log_path).len(), 1);
    let before = tensor_map(load_python_state_dict(&restart_checkpoint).unwrap());
    let after = tensor_map(load_python_state_dict(&restarted.checkpoints[0]).unwrap());
    assert!(tensors_changed(&before, &after));
    let model = Model::new(
        kind,
        Device::Cpu,
        Some(input_channels as i64),
        !pileup,
        Some(&restarted.checkpoints[0]),
    )
    .unwrap();
    let archive = BinArchive::open(
        &fs::read_dir(&training)
            .unwrap()
            .map(|entry| entry.unwrap().path())
            .find(|path| path.is_file())
            .unwrap(),
    )
    .unwrap();
    let (input, _) = archive.read_batch(0, 2, input_channels).unwrap();
    let output_tensor = model.forward(&input, false);
    assert_eq!(output_tensor.size()[0], 2);
    assert_eq!(output_tensor.size()[1], if pileup { 24 } else { 90 });
    assert_eq!(output_tensor.isfinite().all().int64_value(&[]), 1);

    let stopping = TrainingOptions {
        checkpoint: restart.checkpoint.clone(),
        max_epochs: 12,
        learning_rate: 0.0,
        output_prefix: output.join("stopping"),
        log_path: output.join("stopping.log"),
        best_checkpoint: output.join("stopping-best.pt"),
        ..restart.clone()
    };
    let stopped = train(&stopping).unwrap();
    assert_eq!(
        stopped.epochs, 11,
        "{test_name}: early stopping did not trigger"
    );
    assert_eq!(stopped.checkpoints.len(), 11);
    assert!(stopping.best_checkpoint.is_file());
    let stopping_losses = training_log(&stopping.log_path);
    assert_eq!(stopping_losses.len(), 11);
    eprintln!("{test_name}: early-stopping losses={stopping_losses:?}");
    eprintln!(
        "{test_name}: weight-only restart and early stopping passed; outputs={}",
        output.display()
    );
}

#[test]
#[ignore = "requires independent real training archives and pretrained checkpoints"]
fn real_data_pileup_converges_stops_and_restarts() {
    real_training_workflow(
        "real_data_pileup_converges_stops_and_restarts",
        ModelKind::Pileup,
        "pileup",
        true,
    );
}

#[test]
#[ignore = "requires independent real training archives and pretrained checkpoints"]
fn real_data_full_alignment_converges_stops_and_restarts() {
    real_training_workflow(
        "real_data_full_alignment_converges_stops_and_restarts",
        ModelKind::FullAlignment,
        "fa",
        false,
    );
}

fn controlled_weight_only_restart(
    test_name: &str,
    kind: ModelKind,
    platform: Platform,
    pileup: bool,
    shape: &[usize],
) {
    let _guard = TORCH_TEST_LOCK.lock().unwrap();
    let root = tempfile::tempdir().unwrap();
    let data = write_deterministic_archive(root.path(), shape, pileup);
    let validation = root.path().join("validation");
    fs::create_dir_all(&validation).unwrap();
    fs::copy(data.join("toy.bin"), validation.join("validation_20.bin")).unwrap();
    let input_channels = *shape.last().unwrap();
    let selected = root.path().join("selected.pt");
    let different = root.path().join("different.pt");
    tch::manual_seed(101);
    let model = Model::new(
        kind,
        Device::Cpu,
        Some(input_channels as i64),
        !pileup,
        None,
    )
    .unwrap();
    save_checkpoint(&model.var_store, &selected).unwrap();
    tch::manual_seed(202);
    let model = Model::new(
        kind,
        Device::Cpu,
        Some(input_channels as i64),
        !pileup,
        None,
    )
    .unwrap();
    save_checkpoint(&model.var_store, &different).unwrap();
    assert!(tensors_changed(
        &tensor_map(load_python_state_dict(&selected).unwrap()),
        &tensor_map(load_python_state_dict(&different).unwrap()),
    ));

    let options = TrainingOptions {
        checkpoint: Some(selected.clone()),
        pileup,
        add_indel_length: !pileup,
        max_epochs: 1,
        mini_epochs: 1,
        batch_size: DETERMINISTIC_BATCH_SIZE,
        chunk_size: DETERMINISTIC_BATCH_SIZE,
        learning_rate: 1e-3,
        validation_dir: Some(validation),
        seed: Some(29),
        log_path: root.path().join("restart.log"),
        best_checkpoint: root.path().join("restart-best.pt"),
        ..TrainingOptions::new(data, root.path().join("restart"), platform)
    };
    let restarted = train(&options).unwrap();
    assert_eq!(restarted.epochs, 1);
    assert_eq!(restarted.checkpoints.len(), 1);
    assert_eq!(training_log(&options.log_path).len(), 1);
    assert!(
        restarted.checkpoints[0]
            .file_name()
            .unwrap()
            .to_string_lossy()
            .ends_with(".01.pt")
    );

    let reproduced = train(&TrainingOptions {
        output_prefix: root.path().join("reproduced"),
        log_path: root.path().join("reproduced.log"),
        best_checkpoint: root.path().join("reproduced-best.pt"),
        ..options.clone()
    })
    .unwrap();
    assert_eq!(restarted.final_loss, reproduced.final_loss);
    assert_eq!(
        restarted.best_validation_loss,
        reproduced.best_validation_loss
    );
    assert_eq!(
        fs::read_to_string(&options.log_path).unwrap(),
        fs::read_to_string(root.path().join("reproduced.log")).unwrap()
    );
    let restarted_tensors = tensor_map(load_python_state_dict(&restarted.checkpoints[0]).unwrap());
    compare_tensor_maps(
        &restarted_tensors,
        &tensor_map(load_python_state_dict(&reproduced.checkpoints[0]).unwrap()),
        0.0,
        "fresh weight-only restart reproduction",
    );

    let control = train(&TrainingOptions {
        checkpoint: Some(different),
        output_prefix: root.path().join("different-control"),
        log_path: root.path().join("different-control.log"),
        best_checkpoint: root.path().join("different-control-best.pt"),
        ..options
    })
    .unwrap();
    assert!(tensors_changed(
        &restarted_tensors,
        &tensor_map(load_python_state_dict(&control.checkpoints[0]).unwrap()),
    ));
    eprintln!(
        "{test_name}: exact fresh-run loss={}, different-checkpoint control changed tensors",
        restarted.final_loss
    );
}

#[test]
fn pileup_weight_only_restart_reproduces_fresh_training() {
    controlled_weight_only_restart(
        "pileup_weight_only_restart_reproduces_fresh_training",
        ModelKind::Pileup,
        Platform::Ont,
        true,
        &[33, 18],
    );
}

#[test]
fn full_alignment_weight_only_restart_reproduces_fresh_training() {
    controlled_weight_only_restart(
        "full_alignment_weight_only_restart_reproduces_fresh_training",
        ModelKind::FullAlignment,
        Platform::Hifi,
        false,
        &[55, 33, 8],
    );
}

fn bam_qname_counts(path: &Path) -> BTreeMap<Vec<u8>, usize> {
    use rust_htslib::bam::{self, Read as _};

    let mut reader = bam::Reader::from_path(path).unwrap();
    let mut counts = BTreeMap::new();
    for record in reader.records() {
        *counts.entry(record.unwrap().qname().to_vec()).or_default() += 1;
    }
    counts
}

#[test]
#[ignore = "requires a completed bounded quick-training workflow"]
fn quick_workflow_uses_complementary_training_and_validation_reads() {
    let Some(work) = std::env::var_os("CLAIR3_QUICK_WORK").map(PathBuf::from) else {
        eprintln!(
            "SKIP quick_workflow_uses_complementary_training_and_validation_reads: CLAIR3_QUICK_WORK is not set"
        );
        return;
    };
    let full = work.join("bam/chr21.ht.bam");
    let validation = work.join("bam/chr21_222.bam");
    let fixed_training = work.join("bam/chr21_1000.training.bam");
    let training = if fixed_training.is_file() {
        fixed_training
    } else {
        full.clone()
    };
    let missing: Vec<_> = [&full, &validation, &training]
        .into_iter()
        .filter(|path| !path.is_file())
        .collect();
    if !missing.is_empty() {
        eprintln!(
            "SKIP quick_workflow_uses_complementary_training_and_validation_reads: missing {}",
            missing
                .iter()
                .map(|path| path.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
        return;
    }
    let full_counts = bam_qname_counts(&full);
    let training_counts = bam_qname_counts(&training);
    let validation_counts = bam_qname_counts(&validation);
    assert!(!full_counts.is_empty());
    assert!(!training_counts.is_empty());
    assert!(!validation_counts.is_empty());
    let shared: Vec<_> = training_counts
        .keys()
        .filter(|name| validation_counts.contains_key(*name))
        .map(|name| String::from_utf8_lossy(name).into_owned())
        .collect();
    assert!(
        shared.is_empty(),
        "training and validation share {} QNAMEs; first {:?}; training path {}",
        shared.len(),
        &shared[..shared.len().min(5)],
        training.display()
    );
    assert_eq!(
        training_counts
            .keys()
            .chain(validation_counts.keys())
            .collect::<BTreeSet<_>>()
            .len(),
        full_counts.len()
    );
    for (name, count) in &full_counts {
        assert_eq!(
            training_counts.get(name).copied().unwrap_or(0)
                + validation_counts.get(name).copied().unwrap_or(0),
            *count,
            "QNAME {}",
            String::from_utf8_lossy(name)
        );
    }
    eprintln!(
        "quick workflow disjoint QNAMEs: full={} training={} validation={} shared=0; per-QNAME record complement holds",
        full_counts.len(),
        training_counts.len(),
        validation_counts.len()
    );
}
