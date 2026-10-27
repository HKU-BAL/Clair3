use std::{fs, path::PathBuf};

use clair3::{
    config::Platform,
    model::{Model, ModelKind},
    tensor::labels::{labels_from_reference, labels_from_vcf_columns},
    train::{
        checkpoint::save_checkpoint,
        data::{BinWriter, TensorExample},
        trainer::{TrainingOptions, train},
    },
};
use tch::{Device, Tensor};

fn directory() -> PathBuf {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/verify/training/smoke");
    fs::create_dir_all(&path).unwrap();
    path
}

#[test]
fn checkpoint_roundtrip_preserves_every_tensor() {
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
