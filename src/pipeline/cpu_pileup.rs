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
    config::PREDICT_BATCH_SIZE,
    io::{interval::Region, reference::Reference, variant::VcfRecord},
    model::{Model, ModelKind},
    pipeline::call::{CallConfig, CallOutput, pileup_batch, pileup_decode_options},
    tensor::pileup::PileupBatch,
};

struct GeneratedBatch {
    values: Arc<PileupBatch>,
    region: Region,
    decode: DecodeOptions,
}

struct InferenceBatch {
    job: usize,
    values: Arc<PileupBatch>,
    decode: DecodeOptions,
    start: usize,
    end: usize,
}

fn generate_batches(jobs: &[CallConfig], threads: usize) -> Result<Vec<GeneratedBatch>> {
    let jobs = Arc::new(jobs);
    let next = Arc::new(AtomicUsize::new(0));
    let (sender, receiver) = mpsc::channel();
    let workers = if threads >= 8 {
        threads
    } else {
        (threads / 2).max(1)
    }
    .min(jobs.len());
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
                        let mut config = jobs[index].clone();
                        config.pileup_parallelism = if workers >= 8 { 1 } else { 2 };
                        let result = pileup_batch(&config, &reference).map(|(batch, region)| {
                            GeneratedBatch {
                                values: Arc::new(batch),
                                region,
                                decode: pileup_decode_options(&config),
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
                bail!("a pileup tensor worker panicked");
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

fn infer_batch(batch: &InferenceBatch, model: &Model) -> Result<Vec<VcfRecord>> {
    let values = &batch.values;
    let count = batch.end - batch.start;
    let shape = [count, values.shape[1], values.shape[2]];
    let item_size = values.shape[1..].iter().product::<usize>();
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
            PREDICT_BATCH_SIZE,
        )
    })
}

fn execute_inference(
    jobs: &[CallConfig],
    batches: Vec<GeneratedBatch>,
    threads: usize,
) -> Result<Vec<CallOutput>> {
    let mut tasks = Vec::new();
    for (job, batch) in batches.iter().enumerate() {
        for start in (0..batch.values.positions.len()).step_by(PREDICT_BATCH_SIZE) {
            tasks.push(InferenceBatch {
                job,
                values: Arc::clone(&batch.values),
                decode: batch.decode.clone(),
                start,
                end: (start + PREDICT_BATCH_SIZE).min(batch.values.positions.len()),
            });
        }
    }
    if tasks.is_empty() {
        return Ok((0..jobs.len()).map(|_| CallOutput::default()).collect());
    }
    let tasks = Arc::new(tasks);
    let next = Arc::new(AtomicUsize::new(0));
    let (sender, receiver) = mpsc::channel();
    let outputs = std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..threads.min(tasks.len()) {
            let tasks = Arc::clone(&tasks);
            let next = Arc::clone(&next);
            let sender = sender.clone();
            handles.push(scope.spawn(move || {
                let worker = (|| -> Result<()> {
                    let config = jobs.first().context("no pileup jobs")?;
                    let model = Model::new(
                        ModelKind::Pileup,
                        Device::Cpu,
                        Some(18),
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
                bail!("a pileup inference worker panicked");
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
        Ok(outputs)
    })?;
    let mut results: Vec<Vec<VcfRecord>> = (0..jobs.len()).map(|_| Vec::new()).collect();
    for (task, batch) in outputs {
        results[tasks[task].job].extend(batch);
    }
    Ok(batches
        .into_iter()
        .zip(results)
        .map(|(batch, records)| CallOutput {
            records: records
                .into_iter()
                .filter(|record| {
                    record.contig == batch.region.contig
                        && record.position > batch.region.start
                        && record.position <= batch.region.end
                })
                .collect(),
            nonvariants: Vec::new(),
        })
        .collect())
}

pub fn execute_pileup_cpu_jobs(jobs: &[CallConfig], threads: usize) -> Result<Vec<CallOutput>> {
    if jobs.is_empty() {
        return Ok(Vec::new());
    }
    let started = std::time::Instant::now();
    let batches = generate_batches(jobs, threads)?;
    let generated = started.elapsed();
    let outputs = execute_inference(jobs, batches, threads)?;
    if std::env::var_os("CLAIR3_PROFILE").is_some() {
        eprintln!(
            "pileup tensor stage={:.3}s inference stage={:.3}s",
            generated.as_secs_f64(),
            started.elapsed().as_secs_f64() - generated.as_secs_f64(),
        );
    }
    Ok(outputs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{config::Platform, model::PileupModel};

    #[test]
    fn parallel_pileup_inference_preserves_job_and_record_order() {
        tch::set_num_threads(1);
        let directory = tempfile::tempdir().unwrap();
        let checkpoint = directory.path().join("model.safetensors");
        let store = tch::nn::VarStore::new(Device::Cpu);
        let _network = PileupModel::new(&store.root(), false, 18);
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
        let batch = PileupBatch {
            shape: vec![count, 33, 18],
            values: (0..count * 33 * 18)
                .map(|index| (index % 100) as i32)
                .collect(),
            positions: (0..count)
                .map(|index| format!("chrToy:{}:A", index + 1))
                .collect(),
            alt_info: (0..count).map(|_| "10-XC 5 RA 5".to_string()).collect(),
            ..Default::default()
        };
        let batch = Arc::new(batch);
        let decode = pileup_decode_options(&jobs[0]);
        let model = Model::new(
            ModelKind::Pileup,
            Device::Cpu,
            Some(18),
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
                PREDICT_BATCH_SIZE,
            )
            .unwrap()
        });
        let expected: Vec<String> = expected.iter().map(|record| record.to_line()).collect();
        let batches = (0..jobs.len())
            .map(|_| GeneratedBatch {
                values: Arc::clone(&batch),
                region: Region::new("chrToy", 0, count as i64).unwrap(),
                decode: decode.clone(),
            })
            .collect();
        let outputs = execute_inference(&jobs, batches, 4).unwrap();
        assert_eq!(outputs.len(), jobs.len());
        for output in outputs {
            let actual: Vec<String> = output
                .records
                .iter()
                .map(|record| record.to_line())
                .collect();
            assert_eq!(actual, expected);
        }
    }
}
