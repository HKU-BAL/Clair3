#![cfg(feature = "distributed")]

use std::{
    fs,
    net::TcpListener,
    path::PathBuf,
    process::{Command, Stdio},
    time::{Duration, Instant},
};

use clair3::{
    model::{Model, ModelKind},
    train::{collective::C10dCollective, distributed::Collective},
};
use tch::{Device, Tensor};

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
        let model = Model::new(ModelKind::Pileup, Device::Cpu, Some(18), false, None).unwrap();
        let input = Tensor::from_slice(&vec![10_i32; 10 * 33 * 18]).reshape([10, 33, 18]);
        let loss = model.forward(&input, true).sum(tch::Kind::Float);
        loss.backward();
        let mut variables: Vec<_> = model.var_store.variables().into_iter().collect();
        variables.sort_by(|left, right| left.0.cmp(&right.0));
        for (name, variable) in variables {
            let gradient = variable.grad();
            let gradient = if gradient.defined() {
                gradient
            } else {
                variable.zeros_like()
            };
            eprintln!("rank {rank}: gradient {name} {:?}", gradient.size());
            let _ = collective.all_reduce_mean(&gradient).unwrap();
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
}
