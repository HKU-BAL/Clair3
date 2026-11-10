use std::{fs, path::PathBuf};

use clair3::{
    model::{Model, ModelKind},
    train::data::BinArchive,
};
use tch::{Device, Kind};

/// Per-task metrics of the full-alignment model on the independent validation
/// archives, for every checkpoint of a training run.
///
/// The aggregate validation loss of the full-alignment model can rise while
/// the training loss falls. This diagnostic reports the error of every one of
/// the 90 model outputs separately, for each epoch, so that the tasks that
/// degrade are visible instead of being hidden in the aggregate.
fn env_path(name: &str, default: &str) -> PathBuf {
    std::env::var_os(name)
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(default))
}

fn archives(directory: &PathBuf) -> Vec<PathBuf> {
    let mut paths: Vec<_> = fs::read_dir(directory)
        .unwrap_or_else(|error| panic!("failed to read {directory:?}: {error}"))
        .map(|entry| entry.unwrap().path())
        .filter(|path| path.is_file())
        .collect();
    paths.sort();
    paths
}

/// The mean absolute error of every model output, over the archives.
fn per_task_error(model: &Model, directory: &PathBuf, input_channels: i64) -> Vec<f64> {
    let mut errors = vec![0.0_f64; 90];
    let mut rows = 0_usize;
    for path in archives(directory) {
        let archive = BinArchive::open(&path).unwrap();
        let length = archive.len();
        for start in (0..length).step_by(200) {
            let (input, labels) = archive
                .read_batch(start, (start + 200).min(length), input_channels as usize)
                .unwrap();
            let output = model.forward(&input, false);
            let difference = (output - labels).abs().to_kind(Kind::Float);
            let per_task = difference.sum_dim_intlist(&[0i64][..], false, Kind::Float);
            let values = Vec::<f64>::try_from(&per_task).unwrap();
            for (index, value) in values.iter().enumerate() {
                errors[index] += value;
            }
            rows += input.size()[0] as usize;
        }
    }
    assert!(rows > 0, "no rows in {directory:?}");
    errors.iter().map(|value| value / rows as f64).collect()
}

#[test]
fn full_alignment_validation_tasks_are_reported_per_epoch() {
    let root = env_path(
        "CLAIR3_VALIDATION_RUN",
        "/Users/Aqua/Temp/Clair3/target/correctness-validation-20261108/verify/training/task4-real/fa-403",
    );
    let validation = env_path(
        "CLAIR3_VALIDATION_ARCHIVES",
        "/Users/Aqua/Temp/Clair3/target/correctness-validation-20261108/verify/training/task4-real-inputs/fa/validation/bin_HG002_20_222",
    );
    let checkpoints: Vec<_> = ["convergence.01.pt", "convergence.02.pt", "convergence.03.pt"]
        .iter()
        .map(|name| root.join(name))
        .collect();
    let missing: Vec<_> = checkpoints
        .iter()
        .chain(std::iter::once(&validation))
        .filter(|path| !path.exists())
        .collect();
    if !missing.is_empty() {
        eprintln!(
            "SKIP full_alignment_validation_tasks_are_reported_per_epoch: missing {}",
            missing
                .iter()
                .map(|path| path.display().to_string())
                .collect::<Vec<_>>()
                .join(", ")
        );
        return;
    }

    let mut per_epoch = Vec::new();
    for checkpoint in &checkpoints {
        let model = Model::new(
            ModelKind::FullAlignment,
            Device::Cpu,
            Some(8),
            true,
            Some(checkpoint),
        )
        .unwrap();
        let errors = per_task_error(&model, &validation, 8);
        assert!(
            errors.iter().all(|value| value.is_finite()),
            "non-finite per-task error for {}",
            checkpoint.display()
        );
        let mean: f64 = errors.iter().sum::<f64>() / errors.len() as f64;
        eprintln!(
            "{}: mean per-task error {mean:.6}, worst {:?}",
            checkpoint.display(),
            errors
                .iter()
                .enumerate()
                .max_by(|left, right| left.1.total_cmp(right.1))
                .map(|(index, value)| (index, *value))
        );
        per_epoch.push(errors);
    }

    // Report the tasks whose error grows the most between the first and the
    // last epoch, so that a degrading task class is visible by name.
    let first = &per_epoch[0];
    let last = per_epoch.last().unwrap();
    let mut changes: Vec<_> = first
        .iter()
        .enumerate()
        .map(|(index, value)| (index, last[index] - value))
        .collect();
    changes.sort_by(|left, right| right.1.total_cmp(&left.1));
    eprintln!(
        "tasks with the largest error increase between the first and the last epoch: {:?}",
        &changes[..changes.len().min(10)]
    );
    let decreasing = changes
        .iter()
        .filter(|(_, change)| *change < 0.0)
        .count();
    eprintln!(
        "tasks that improved: {decreasing}, tasks that degraded: {}",
        changes.len() - decreasing
    );
}
