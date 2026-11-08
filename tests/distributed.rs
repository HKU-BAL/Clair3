#![cfg(feature = "distributed")]

use std::{
    collections::BTreeMap,
    fs,
    net::TcpListener,
    path::{Path, PathBuf},
    process::{Command, Stdio},
    time::{Duration, Instant},
};

use clair3::{
    model::{Model, ModelKind},
    train::{collective::C10dCollective, distributed::Collective},
};
use tch::{Device, Kind, Tensor};

fn focal_loss(labels: &Tensor, predictions: &Tensor) -> Tensor {
    let predictions = predictions.clamp(1e-9, 1.0 - 1e-9);
    let cross_entropy = -labels * predictions.log();
    let weights = (1.0_f64 - &predictions).square() * labels;
    (cross_entropy * weights)
        .sum_dim_intlist(&[-1i64][..], false, Kind::Float)
        .mean(Kind::Float)
}

fn rank_local_gradients(rank: usize) -> Vec<(String, Tensor)> {
    tch::manual_seed(0);
    let model = Model::new(ModelKind::Pileup, Device::Cpu, Some(18), false, None).unwrap();
    let values: Vec<i32> = (0..2 * 33 * 18)
        .map(|index| ((index + rank * 37) % 41) as i32 - 20)
        .collect();
    let input = Tensor::from_slice(&values).reshape([2, 33, 18]);
    let mut labels = vec![0.0_f32; 2 * 90];
    for sample in 0..2 {
        labels[sample * 90 + (sample + rank) % 21] = 1.0;
        labels[sample * 90 + 21 + (sample + rank) % 3] = 1.0;
    }
    let labels = Tensor::from_slice(&labels).reshape([2, 90]);
    let outputs = model.forward_tasks(&input, true);
    let losses: Vec<_> = outputs
        .iter()
        .enumerate()
        .map(|(task, output)| {
            let offset = if task == 0 { 0 } else { 21 };
            let width = if task == 0 { 21 } else { 3 };
            focal_loss(&labels.narrow(1, offset, width), output)
        })
        .collect();
    Tensor::stack(&losses, 0).sum(Kind::Float).backward();
    let mut variables: Vec<_> = model.var_store.variables().into_iter().collect();
    variables.sort_by(|left, right| left.0.cmp(&right.0));
    variables
        .into_iter()
        .map(|(name, variable)| {
            let gradient = variable.grad();
            let gradient = if gradient.defined() {
                gradient
            } else {
                variable.zeros_like()
            };
            (name, gradient)
        })
        .collect()
}

fn tensor_map(tensors: Vec<(String, Tensor)>) -> BTreeMap<String, Tensor> {
    tensors.into_iter().collect()
}

fn assert_reduced_gradients(directory: &Path) {
    let rank0 = rank_local_gradients(0);
    let rank1 = rank_local_gradients(1);
    assert_eq!(rank0.len(), rank1.len());
    assert!(rank0.iter().any(|(name, left)| {
        let right = &rank1.iter().find(|(other, _)| other == name).unwrap().1;
        left.ne_tensor(right).any().int64_value(&[]) != 0
    }));
    assert!(
        rank0
            .iter()
            .any(|(_, gradient)| gradient.abs().max().double_value(&[]) > 0.0)
    );
    let reduced0 = tensor_map(
        Tensor::read_safetensors(directory.join("rank0-gradients.safetensors")).unwrap(),
    );
    let reduced1 = tensor_map(
        Tensor::read_safetensors(directory.join("rank1-gradients.safetensors")).unwrap(),
    );
    let mut maximum = 0.0_f64;
    for ((name, gradient0), (_, gradient1)) in rank0.iter().zip(&rank1) {
        let expected = (gradient0 + gradient1) / 2.0;
        for (rank, reduced) in [&reduced0, &reduced1].into_iter().enumerate() {
            let actual = &reduced[name];
            let difference = (actual - &expected).abs().max().double_value(&[]);
            maximum = maximum.max(difference);
            assert_eq!(
                actual.eq_tensor(&expected).all().int64_value(&[]),
                1,
                "rank {rank} reduced gradient {name} differs from the independent arithmetic mean by {difference}"
            );
        }
    }
    eprintln!("maximum distributed gradient difference: {maximum:.3e}");
}

#[test]
fn distributed_worker() {
    let Ok(rank) = std::env::var("CLAIR3_TEST_RANK") else {
        return;
    };
    let rank: usize = rank.parse().unwrap();
    eprintln!("rank {rank}: initializing");
    let collective = C10dCollective::new(rank, 2, rank, false).unwrap();
    eprintln!("rank {rank}: broadcasting");
    assert_eq!(collective.broadcast_i64(1234, 0).unwrap(), 1234);
    let input = Tensor::from_slice(&[rank as f32 + 1.0, 10.0]);
    eprintln!("rank {rank}: reducing");
    let output = collective.all_reduce_mean(&input).unwrap();
    assert_eq!(Vec::<f32>::try_from(&output).unwrap(), vec![1.5, 10.0]);
    let input = Tensor::from_slice(&vec![rank as f32 + 1.0; 640 * 256]).reshape([640, 256]);
    let output = collective.all_reduce_mean(&input).unwrap();
    assert_eq!(output.double_value(&[0, 0]), 1.5);
    assert_eq!(output.double_value(&[639, 255]), 1.5);
    for _ in 0..32 {
        let output = collective.all_reduce_mean(&input).unwrap();
        assert_eq!(output.double_value(&[0, 0]), 1.5);
    }
    if std::env::var_os("CLAIR3_TEST_AUTOGRAD").is_some() {
        let gradients = rank_local_gradients(rank);
        let mut reduced = Vec::new();
        for (name, gradient) in gradients {
            eprintln!("rank {rank}: gradient {name} {:?}", gradient.size());
            reduced.push((name, collective.all_reduce_mean(&gradient).unwrap()));
        }
        if let Some(directory) = std::env::var_os("CLAIR3_TEST_GRADIENT_DIR") {
            Tensor::write_safetensors(
                &reduced,
                Path::new(&directory).join(format!("rank{rank}-gradients.safetensors")),
            )
            .unwrap();
        }
    }
    collective.barrier().unwrap();
    eprintln!("rank {rank}: finished");
}

#[test]
#[ignore = "requires two local LibTorch processes"]
fn two_process_collectives_work() {
    let port = TcpListener::bind("127.0.0.1:0")
        .unwrap()
        .local_addr()
        .unwrap()
        .port();
    let executable = std::env::current_exe().unwrap();
    let directory = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/verify/distributed");
    fs::create_dir_all(&directory).unwrap();
    for rank in 0..2 {
        let path = directory.join(format!("rank{rank}-gradients.safetensors"));
        if path.exists() {
            fs::remove_file(path).unwrap();
        }
    }
    let mut children = Vec::new();
    for rank in 0..2 {
        let log = fs::File::create(directory.join(format!("rank{rank}.log"))).unwrap();
        children.push((
            Command::new(&executable)
                .arg("distributed_worker")
                .arg("--exact")
                .arg("--nocapture")
                .env("CLAIR3_TEST_RANK", rank.to_string())
                .env("MASTER_ADDR", "127.0.0.1")
                .env("MASTER_PORT", port.to_string())
                .env("LIBTORCH_USE_PYTORCH", "1")
                .env("RUST_BACKTRACE", "1")
                .env("CLAIR3_TEST_AUTOGRAD", "1")
                .env("CLAIR3_TEST_GRADIENT_DIR", &directory)
                .stdout(Stdio::from(log.try_clone().unwrap()))
                .stderr(Stdio::from(log))
                .spawn()
                .unwrap(),
            rank,
        ));
    }
    for (child, rank) in &mut children {
        let started = Instant::now();
        loop {
            if child.try_wait().unwrap().is_some() {
                break;
            }
            if started.elapsed() > Duration::from_secs(30) {
                child.kill().unwrap();
                child.wait().unwrap();
                break;
            }
            std::thread::sleep(Duration::from_millis(100));
        }
        let status = child.wait().unwrap();
        assert!(
            status.success(),
            "worker {rank} failed: {}",
            fs::read_to_string(directory.join(format!("rank{rank}.log"))).unwrap()
        );
    }
    assert_reduced_gradients(&directory);
}
