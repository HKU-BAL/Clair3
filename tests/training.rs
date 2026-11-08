use std::{
    collections::BTreeMap,
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
