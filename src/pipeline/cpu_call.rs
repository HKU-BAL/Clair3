use std::{
    collections::BTreeMap,
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
        mpsc,
    },
};

use anyhow::{Context, Result, bail};
use tch::Device;

use crate::{
    call::{
        caller::{TensorInput, call_tensor_batch},
        decode::DecodeOptions,
    },
    config::FULL_ALIGNMENT_PREDICT_BATCH_SIZE,
    io::{reference::Reference, variant::VcfRecord},
    model::{Model, ModelKind},
    pipeline::call::{CallConfig, CallOutput, full_alignment_batch, full_decode_options},
    tensor::full_alignment::FullAlignmentBatch,
};

struct GeneratedBatch {
    values: Arc<FullAlignmentBatch>,
    decode: DecodeOptions,
}

struct InferenceBatch {
    job: usize,
    values: Arc<FullAlignmentBatch>,
    decode: DecodeOptions,
    start: usize,
    end: usize,
}

fn generate_batches(jobs: &[CallConfig], threads: usize) -> Result<Vec<GeneratedBatch>> {
    let jobs = Arc::new(jobs);
    let next = Arc::new(AtomicUsize::new(0));
    let (sender, receiver) = mpsc::channel();
    let workers = (threads / 2).max(1).min(jobs.len());
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..workers {
            let jobs = Arc::clone(&jobs);
            let next = Arc::clone(&next);
            let sender = sender.clone();
            handles.push(scope.spawn(move || {
                let worker = (|| -> Result<()> {
                    let reference = Reference::open(&jobs[0].reference)?;
                    loop {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        if index >= jobs.len() {
                            break;
                        }
                        let result = full_alignment_batch(&jobs[index], &reference).map(|batch| {
                            GeneratedBatch {
                                values: Arc::new(batch),
                                decode: full_decode_options(&jobs[index]),
                            }
                        });
                        if sender.send((index, result)).is_err() {
                            break;
                        }
                    }
                    Ok(())
                })();
                if let Err(error) = worker {
                    let _ = sender.send((usize::MAX, Err(error)));
                }
            }));
        }
        drop(sender);
        let mut outputs = BTreeMap::new();
        let mut errors = Vec::new();
        for (index, result) in receiver {
            match result {
                Ok(output) => {
                    outputs.insert(index, output);
                }
                Err(error) => errors.push(error),
            }
        }
        for handle in handles {
            if handle.join().is_err() {
                bail!("a full-alignment tensor worker panicked");
            }
        }
        if let Some(error) = errors.into_iter().next() {
            return Err(error);
        }
        if outputs.len() != jobs.len() {
            bail!(
                "tensor workers returned {} of {} jobs",
                outputs.len(),
                jobs.len()
            );
        }
        Ok(outputs.into_values().collect())
    })
}

fn configured_batch_size() -> usize {
    std::env::var("CLAIR3_FULL_ALIGNMENT_BATCH")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .filter(|size| *size > 0)
        .unwrap_or(FULL_ALIGNMENT_PREDICT_BATCH_SIZE)
}

fn infer_batch(batch: &InferenceBatch, model: &Model) -> Result<Vec<VcfRecord>> {
    let values = &batch.values;
    let count = batch.end - batch.start;
    let shape = [count, values.shape[1], values.shape[2], values.shape[3]];
    let item_size = values.shape[1..].iter().product::<usize>();
    let batch_size = configured_batch_size();
    tch::no_grad(|| {
        call_tensor_batch(
            model,
            TensorInput {
                shape: &shape,
                values: &values.values[batch.start * item_size..batch.end * item_size],
                positions: &values.positions[batch.start..batch.end],
                alt_info: &values.alt_info[batch.start..batch.end],
            },
            &batch.decode,
            batch_size,
        )
    })
}

fn execute_inference(
    jobs: &[CallConfig],
    batches: Vec<GeneratedBatch>,
    threads: usize,
) -> Result<Vec<CallOutput>> {
    let batch_size = configured_batch_size();
    let mut tasks = Vec::new();
    for (job, batch) in batches.into_iter().enumerate() {
        for start in (0..batch.values.positions.len()).step_by(batch_size) {
            tasks.push(InferenceBatch {
                job,
                values: Arc::clone(&batch.values),
                decode: batch.decode.clone(),
                start,
                end: (start + batch_size).min(batch.values.positions.len()),
            });
        }
    }
    if tasks.is_empty() {
        return Ok((0..jobs.len()).map(|_| CallOutput::default()).collect());
    }
    let tasks = Arc::new(tasks);
    let next = Arc::new(AtomicUsize::new(0));
    let (sender, receiver) = mpsc::channel();
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..threads.min(tasks.len()) {
            let tasks = Arc::clone(&tasks);
            let next = Arc::clone(&next);
            let sender = sender.clone();
            handles.push(scope.spawn(move || {
                let worker = (|| -> Result<()> {
                    let config = jobs.first().context("no full-alignment jobs")?;
                    let model = Model::new(
                        ModelKind::FullAlignment,
                        Device::Cpu,
                        Some(if config.enable_dwell_time { 9 } else { 8 }),
                        config.add_indel_length,
                        Some(&config.checkpoint),
                    )?;
                    loop {
                        let index = next.fetch_add(1, Ordering::Relaxed);
                        if index >= tasks.len() {
                            break;
                        }
                        let result = infer_batch(&tasks[index], &model);
                        if sender.send((index, result)).is_err() {
                            break;
                        }
                    }
                    Ok(())
                })();
                if let Err(error) = worker {
                    let _ = sender.send((usize::MAX, Err(error)));
                }
            }));
        }
        drop(sender);
        let mut outputs = BTreeMap::new();
        let mut errors = Vec::new();
        for (index, result) in receiver {
            match result {
                Ok(output) => {
                    outputs.insert(index, output);
                }
                Err(error) => errors.push(error),
            }
        }
        for handle in handles {
            if handle.join().is_err() {
                bail!("a full-alignment inference worker panicked");
            }
        }
        if let Some(error) = errors.into_iter().next() {
            return Err(error);
        }
        if outputs.len() != tasks.len() {
            bail!(
                "inference workers returned {} of {} batches",
                outputs.len(),
                tasks.len()
            );
        }
        let mut results: Vec<Vec<VcfRecord>> = (0..jobs.len()).map(|_| Vec::new()).collect();
        for (task, batch) in outputs {
            results[tasks[task].job].extend(batch);
        }
        Ok(results
            .into_iter()
            .map(|records| CallOutput {
                records,
                nonvariants: Vec::new(),
            })
            .collect())
    })
}

pub fn execute_full_cpu_jobs(jobs: &[CallConfig], threads: usize) -> Result<Vec<CallOutput>> {
    if jobs.is_empty() {
        return Ok(Vec::new());
    }
    let started = std::time::Instant::now();
    let batches = generate_batches(jobs, threads)?;
    let generated = started.elapsed();
    let outputs = execute_inference(jobs, batches, threads)?;
    if std::env::var_os("CLAIR3_PROFILE").is_some() {
        eprintln!(
            "full tensor stage={:.3}s inference stage={:.3}s",
            generated.as_secs_f64(),
            started.elapsed().as_secs_f64() - generated.as_secs_f64(),
        );
    }
    Ok(outputs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{config::Platform, model::FullAlignmentModel};

    #[test]
    fn parallel_inference_preserves_job_and_record_order() {
        tch::set_num_threads(1);
        let directory = tempfile::tempdir().unwrap();
        let checkpoint = directory.path().join("model.safetensors");
        let store = tch::nn::VarStore::new(Device::Cpu);
        let _network = FullAlignmentModel::new(&store.root(), false, 8);
        store.save(&checkpoint).unwrap();
        let root = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
        let mut config = CallConfig::new(
            root.join("tests/fixtures/data/toy.bam"),
            root.join("tests/fixtures/data/toy.fa"),
            checkpoint.clone(),
            "chrToy".to_string(),
            Platform::Hifi,
        );
        config.add_indel_length = false;
        let jobs = vec![config.clone(), config.clone(), config];
        let count = 450;
        let batch = FullAlignmentBatch {
            shape: vec![count, 55, 33, 8],
            values: (0..count * 55 * 33 * 8)
                .map(|index| (index % 100) as i8)
                .collect(),
            positions: (0..count)
                .map(|index| format!("chrToy:{}:A", index + 1))
                .collect(),
            alt_info: (0..count).map(|_| "10-XC 5 RA 5".to_string()).collect(),
        };
        let batch = Arc::new(batch);
        let decode = full_decode_options(&jobs[0]);
        let model = Model::new(
            ModelKind::FullAlignment,
            Device::Cpu,
            Some(8),
            false,
            Some(&checkpoint),
        )
        .unwrap();
        let expected = tch::no_grad(|| {
            call_tensor_batch(
                &model,
                TensorInput {
                    shape: &batch.shape,
                    values: &batch.values,
                    positions: &batch.positions,
                    alt_info: &batch.alt_info,
                },
                &decode,
                configured_batch_size(),
            )
            .unwrap()
        });
        let expected: Vec<String> = expected.iter().map(|record| record.to_line()).collect();
        let batches = (0..jobs.len())
            .map(|_| GeneratedBatch {
                values: Arc::clone(&batch),
                decode: decode.clone(),
            })
            .collect();
        let outputs = execute_inference(&jobs, batches, 4).unwrap();
        assert_eq!(outputs.len(), jobs.len());
        for output in outputs {
            let actual: Vec<String> = output.records.iter().map(|record| record.to_line()).collect();
            assert_eq!(actual, expected);
        }
    }
}
