use std::path::{Path, PathBuf};

use clair3::model::{Model, ModelKind};
use tch::{Device, Tensor};

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn cases() -> Vec<(&'static str, ModelKind, i64, bool, &'static str)> {
    vec![
        ("pileup_ont", ModelKind::Pileup, 18, false, "ont/pileup.pt"),
        (
            "pileup_hifi",
            ModelKind::Pileup,
            18,
            false,
            "hifi/pileup.pt",
        ),
        (
            "pileup_ilmn",
            ModelKind::Pileup,
            18,
            false,
            "ilmn/pileup.pt",
        ),
        (
            "full_ont",
            ModelKind::FullAlignment,
            8,
            true,
            "ont/full_alignment.pt",
        ),
        (
            "full_hifi",
            ModelKind::FullAlignment,
            8,
            true,
            "hifi/full_alignment.pt",
        ),
        (
            "full_ilmn",
            ModelKind::FullAlignment,
            8,
            true,
            "ilmn/full_alignment.pt",
        ),
        (
            "full_dwell",
            ModelKind::FullAlignment,
            9,
            true,
            "r1041_e82_400bps_hac_with_mv/full_alignment.pt",
        ),
        (
            "full_v500",
            ModelKind::FullAlignment,
            8,
            true,
            "r1041_e82_400bps_hac_v500/full_alignment.pt",
        ),
        (
            "pileup_indel",
            ModelKind::Pileup,
            18,
            true,
            "tests/fixtures/models/pileup_indel.pt",
        ),
    ]
}

fn checkpoint(reference: &str) -> PathBuf {
    if reference.starts_with("tests/") {
        root().join(reference)
    } else {
        root().join("target/test-data/models").join(reference)
    }
}

fn tensor(path: PathBuf) -> Tensor {
    Tensor::read_npy(path).unwrap()
}

fn check(case: &(&str, ModelKind, i64, bool, &str), device: Device, tolerance: f64) {
    let (name, kind, channels, add_indel, reference) = case;
    let checkpoint = checkpoint(reference);
    if !checkpoint.exists() {
        eprintln!("skipping {name}: {} is not available", checkpoint.display());
        return;
    }
    let model = Model::new(
        *kind,
        device,
        Some(*channels),
        *add_indel,
        Some(&checkpoint),
    )
    .unwrap();
    let input = tensor(root().join(format!("tests/fixtures/models/npy/{name}_input.npy")));
    let expected = tensor(root().join(format!("tests/fixtures/models/npy/{name}_cpu.npy")));
    let actual = model
        .forward(&input.to_device(device), false)
        .to_device(Device::Cpu);
    assert_eq!(actual.size(), expected.size(), "output shape for {name}");
    let difference = (actual - expected).abs().max().double_value(&[]);
    assert!(
        difference <= tolerance,
        "{name} maximum absolute difference {difference} exceeds {tolerance}"
    );
}

#[test]
fn official_models_match_cpu_golden_outputs() {
    for case in cases() {
        check(&case, Device::Cpu, 1e-4);
    }
}

#[test]
fn official_models_match_mps_golden_outputs() {
    if std::env::var("CLAIR3_TEST_MPS").is_err() {
        return;
    }
    for case in cases().into_iter().filter(|case| {
        matches!(
            case.0,
            "pileup_ont" | "full_ont" | "full_dwell" | "pileup_indel"
        )
    }) {
        check(&case, Device::Mps, 0.001);
    }
}

#[test]
fn bloom_checkpoint_loads_and_preserves_legacy_shape_failure() {
    let checkpoint = root().join("tests/fixtures/models/fb.pt");
    let model = Model::new(
        ModelKind::FullAlignmentBloom,
        Device::Cpu,
        Some(8),
        false,
        Some(&checkpoint),
    )
    .unwrap();
    let input = tensor(root().join("tests/fixtures/models/npy/fb_input.npy"));
    let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        model.forward(&input, false)
    }));
    assert!(
        result.is_err(),
        "Clair3_FB unexpectedly accepted its incompatible dense input"
    );
}

#[test]
fn checkpoint_loader_rejects_wrong_architecture() {
    let checkpoint = root().join("tests/fixtures/models/pileup_indel.pt");
    let error = Model::new(
        ModelKind::Pileup,
        Device::Cpu,
        Some(18),
        false,
        Some(&checkpoint),
    );
    assert!(error.is_err());
}

#[test]
fn checkpoint_loader_supports_ddp_prefixes() {
    let source = root().join("tests/fixtures/models/pileup_indel.pt");
    if !source.exists() {
        return;
    }
    let tensors = clair3::model::load_python_state_dict(&source).unwrap();
    let prefixed: Vec<(String, Tensor)> = tensors
        .into_iter()
        .map(|(name, tensor)| (format!("module.{name}"), tensor))
        .collect();
    let path = root().join("target/ddp-prefixed.safetensors");
    Tensor::write_safetensors(&prefixed, &path).unwrap();
    Model::new(ModelKind::Pileup, Device::Cpu, Some(18), true, Some(&path)).unwrap();
}

#[test]
fn checkpoint_loader_supports_wrapped_state_dict() {
    let checkpoint = root().join("target/wrapped-state-dict.pt");
    Model::new(
        ModelKind::Pileup,
        Device::Cpu,
        Some(18),
        true,
        Some(&checkpoint),
    )
    .unwrap();
}

#[test]
fn model_checkpoint_path_is_used() {
    let path = root().join("target/test-data/models/ont/pileup.pt");
    if Path::new(&path).exists() {
        let model = Model::new(ModelKind::Pileup, Device::Cpu, None, false, Some(&path)).unwrap();
        assert!(!model.var_store.is_empty());
    }
}

#[test]
#[ignore = "requires captured production full-alignment tensors"]
fn full_production_model_matches_original() {
    let model = Model::new(
        ModelKind::FullAlignment,
        Device::Cpu,
        Some(8),
        true,
        Some(&root().join("target/test-data/models/ont/full_alignment.pt")),
    )
    .unwrap();
    let input = Tensor::read_npy(root().join("target/verify/full-python.npy"))
        .unwrap()
        .reshape([5, 89, 33, 8]);
    let expected = Tensor::read_npy(root().join("target/verify/full-python-output.npy")).unwrap();
    let difference = (model.forward(&input, false) - expected)
        .abs()
        .max()
        .double_value(&[]);
    assert!(
        difference <= 1e-4,
        "full-alignment probabilities differ by {difference}"
    );
}
