use std::{
    collections::BTreeSet,
    fs,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use tch::{Device, Kind, Tensor, nn, nn::OptimizerConfig};

use crate::{
    config::Platform,
    model::{Model, ModelKind},
    train::{
        checkpoint::save_checkpoint,
        data::BinArchive,
        distributed::{Collective, LocalCollective},
    },
};

const LABEL_SHAPES: [usize; 4] = [21, 3, 33, 33];
const LABEL_OFFSETS: [usize; 4] = [0, 21, 24, 57];

#[derive(Debug, Clone)]
pub struct TrainingOptions {
    pub bin_dir: PathBuf,
    pub checkpoint: Option<PathBuf>,
    pub output_prefix: PathBuf,
    pub platform: Platform,
    pub pileup: bool,
    pub add_indel_length: bool,
    pub enable_dwell_time: bool,
    pub reuse_bin: bool,
    pub max_epochs: usize,
    pub mini_epochs: usize,
    pub learning_rate: f64,
    pub batch_size: usize,
    pub chunk_size: usize,
    pub exclude_training_samples: BTreeSet<String>,
    pub random_validation: bool,
    pub validation_dir: Option<PathBuf>,
    pub device: Device,
    pub seed: Option<u64>,
    pub scale_learning_rate: bool,
    pub log_path: PathBuf,
    pub best_checkpoint: PathBuf,
}

impl TrainingOptions {
    pub fn new(bin_dir: PathBuf, output_prefix: PathBuf, platform: Platform) -> Self {
        Self {
            bin_dir,
            checkpoint: None,
            output_prefix,
            platform,
            pileup: false,
            add_indel_length: true,
            enable_dwell_time: false,
            reuse_bin: false,
            max_epochs: 30,
            mini_epochs: 1,
            learning_rate: 1e-3,
            batch_size: 2000,
            chunk_size: 200,
            exclude_training_samples: BTreeSet::new(),
            random_validation: false,
            validation_dir: None,
            device: Device::Cpu,
            seed: None,
            scale_learning_rate: false,
            log_path: PathBuf::from("training.log"),
            best_checkpoint: PathBuf::from("best_val_loss.pt"),
        }
    }
}

#[derive(Debug, Default)]
pub struct TrainingReport {
    pub epochs: usize,
    pub final_loss: f64,
    pub best_validation_loss: Option<f64>,
    pub checkpoints: Vec<PathBuf>,
}

#[derive(Debug)]
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut value = self.0;
        value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        value ^ (value >> 31)
    }

    fn shuffle<T>(&mut self, values: &mut [T]) {
        for index in (1..values.len()).rev() {
            values.swap(index, self.next() as usize % (index + 1));
        }
    }
}

fn archives(path: &Path, excluded: &BTreeSet<String>) -> Result<Vec<BinArchive>> {
    let mut names = Vec::new();
    for entry in fs::read_dir(path).with_context(|| {
        format!(
            "failed to list training archive directory {}",
            path.display()
        )
    })? {
        let path = entry?.path();
        if !path.is_file() {
            continue;
        }
        let name = path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("");
        if name.contains("_20_") || excluded.iter().any(|prefix| name.contains(prefix)) {
            continue;
        }
        names.push(path);
    }
    names.sort();
    names
        .iter()
        .map(|path| BinArchive::open(path))
        .collect::<Result<_>>()
}

fn chunks(archives: &[BinArchive], batch_size: usize, chunk_size: usize) -> Vec<(usize, usize)> {
    archives
        .iter()
        .enumerate()
        .flat_map(|(archive, dataset)| {
            (0..dataset.chunk_count(batch_size, chunk_size))
                .map(move |chunk| (archive, chunk))
                .collect::<Vec<_>>()
        })
        .collect()
}

type ChunkIndices = (usize, usize);

fn split_validation(
    archives: &[BinArchive],
    batch_size: usize,
    chunk_size: usize,
    enabled: bool,
) -> (Vec<ChunkIndices>, Vec<ChunkIndices>) {
    let mut training = Vec::new();
    let mut validation = Vec::new();
    for (archive, dataset) in archives.iter().enumerate() {
        let count = dataset.chunk_count(batch_size, chunk_size);
        if !enabled || count < 10 {
            training.extend((0..count).map(|chunk| (archive, chunk)));
            continue;
        }
        let train_count = ((count - 10) as f64 * 0.9) as usize;
        let validation_count = count - 10 - train_count;
        training.extend((0..train_count).map(|chunk| (archive, chunk)));
        if validation_count > 0 {
            validation.extend((count - validation_count..count).map(|chunk| (archive, chunk)));
        }
    }
    (training, validation)
}

fn batches(
    chunks: &[(usize, usize)],
    chunks_per_batch: usize,
    world_size: usize,
    rank: usize,
) -> Vec<Vec<(usize, usize)>> {
    if world_size == 1 {
        return chunks
            .chunks(chunks_per_batch)
            .filter(|batch| batch.len() == chunks_per_batch)
            .map(<[(usize, usize)]>::to_vec)
            .collect();
    }
    let samples = chunks.len().div_ceil(world_size);
    let total = samples * world_size;
    let indices: Vec<usize> = (rank..total)
        .step_by(world_size)
        .take(samples)
        .map(|index| index % chunks.len())
        .collect();
    indices
        .chunks(chunks_per_batch)
        .filter(|batch| batch.len() == chunks_per_batch)
        .map(|batch| batch.iter().map(|index| chunks[*index]).collect())
        .collect()
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
fn run_epoch(
    archives: &[BinArchive],
    batches: &[Vec<(usize, usize)>],
    model: &Model,
    mut optimizer: Option<&mut nn::Optimizer>,
    collective: &dyn Collective,
    random_offset: usize,
    chunk_size: usize,
    input_channels: usize,
    add_indel_length: bool,
) -> Result<f64> {
    let mut total_loss = 0.0;
    let mut count = 0;
    for batch in batches {
        let mut inputs = Vec::new();
        let mut labels = Vec::new();
        for (archive, chunk) in batch {
            let start = random_offset + chunk * chunk_size;
            if start >= archives[*archive].len() {
                continue;
            }
            let end = (start + chunk_size).min(archives[*archive].len());
            let (input, label) = archives[*archive].read_batch(start, end, input_channels)?;
            inputs.push(input);
            labels.push(label);
        }
        if inputs.is_empty() {
            continue;
        }
        let input = Tensor::cat(&inputs, 0).to_device(model.var_store.device());
        let label = Tensor::cat(&labels, 0).to_device(model.var_store.device());
        let task_count = if add_indel_length { 4 } else { 2 };
        let train = optimizer.is_some();
        let loss = if train {
            let outputs = model.forward_tasks(&input, true);
            let mut losses = Vec::new();
            for task in 0..task_count {
                losses.push(focal_loss(
                    &label.narrow(1, LABEL_OFFSETS[task] as i64, LABEL_SHAPES[task] as i64),
                    &outputs[task],
                ));
            }
            Tensor::stack(&losses, 0).sum(Kind::Float)
        } else {
            let outputs = model.forward_tasks(&input, false);
            let mut losses = Vec::new();
            for task in 0..task_count {
                losses.push(focal_loss(
                    &label.narrow(1, LABEL_OFFSETS[task] as i64, LABEL_SHAPES[task] as i64),
                    &outputs[task],
                ));
            }
            Tensor::stack(&losses, 0).sum(Kind::Float)
        };
        if let Some(optimizer) = optimizer.as_deref_mut() {
            optimizer.zero_grad();
            loss.backward();
            if collective.world_size() > 1 {
                let mut variables: Vec<_> = model.var_store.variables().into_iter().collect();
                variables.sort_by(|left, right| left.0.cmp(&right.0));
                for (_, variable) in variables {
                    let mut gradient = variable.grad();
                    if !gradient.defined() {
                        gradient = variable.zeros_like();
                    }
                    let averaged = collective.all_reduce_mean(&gradient)?;
                    tch::no_grad(|| {
                        gradient.copy_(&averaged);
                    });
                }
            }
            optimizer.step();
        }
        total_loss += collective
            .all_reduce_mean(&loss.detach())?
            .double_value(&[]);
        count += 1;
    }
    Ok(total_loss / count.max(1) as f64)
}

fn write_checkpoint(model: &Model, path: &Path) -> Result<()> {
    if let Some(parent) = path.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    save_checkpoint(&model.var_store, path)
}

pub fn train(options: &TrainingOptions) -> Result<TrainingReport> {
    let collective = LocalCollective;
    train_with_collective(options, &collective)
}

pub fn train_with_collective(
    options: &TrainingOptions,
    collective: &dyn Collective,
) -> Result<TrainingReport> {
    if options.mini_epochs == 0 || options.chunk_size == 0 {
        bail!("mini-epochs and chunk size must be positive");
    }
    if options.batch_size == 0 || !options.batch_size.is_multiple_of(options.chunk_size) {
        bail!("batch size must be a positive multiple of chunk size");
    }
    if options.max_epochs == 0 {
        bail!("maximum epochs must be positive");
    }
    let archives = archives(&options.bin_dir, &options.exclude_training_samples)?;
    if archives.is_empty() {
        bail!(
            "no training archives were found in {}",
            options.bin_dir.display()
        );
    }
    let validation_archives = options
        .validation_dir
        .as_ref()
        .map(|path| self::archives(path, &options.exclude_training_samples))
        .transpose()?;
    let (mut training_chunks, validation_chunks) = if let Some(validation) = &validation_archives {
        (
            chunks(&archives, options.batch_size, options.chunk_size),
            chunks(validation, options.batch_size, options.chunk_size),
        )
    } else {
        split_validation(
            &archives,
            options.batch_size,
            options.chunk_size,
            options.random_validation,
        )
    };
    let input_channels = if options.pileup {
        18
    } else if options.enable_dwell_time {
        9
    } else {
        8
    };
    for archive in &archives {
        let available = *archive.tensor_shape().last().unwrap_or(&0);
        if available != input_channels + usize::from(options.reuse_bin) {
            bail!(
                "archive has {available} channels; model requires {input_channels} and reuse_bin is {}",
                options.reuse_bin
            );
        }
    }
    if let Some(seed) = options.seed {
        tch::manual_seed(seed as i64);
    }
    let model = Model::new(
        if options.pileup {
            ModelKind::Pileup
        } else {
            ModelKind::FullAlignment
        },
        options.device,
        Some(input_channels as i64),
        options.add_indel_length,
        options.checkpoint.as_deref(),
    )?;
    let learning_rate = options.learning_rate
        * if options.scale_learning_rate {
            collective.world_size() as f64
        } else {
            1.0
        };
    let mut optimizer = nn::adamw(0.9, 0.999, 1e-7).build(&model.var_store, learning_rate)?;
    let chunks_per_batch = options.batch_size / options.chunk_size;
    let chunks_per_epoch = (training_chunks.len() / options.mini_epochs).max(1);
    let total_epochs = options.max_epochs * options.mini_epochs;
    let main = collective.rank() == 0;
    let mut log = if main {
        let mut log = std::io::BufWriter::new(fs::File::create(&options.log_path)?);
        use std::io::Write;
        writeln!(log, "epoch\tloss\tval_loss")?;
        Some(log)
    } else {
        None
    };
    let mut rng = Rng(options.seed.unwrap_or(0));
    let mut best_validation = None;
    let mut patience = 0;
    let mut report = TrainingReport::default();
    for epoch in 0..total_epochs {
        if epoch % options.mini_epochs == 0 {
            rng.shuffle(&mut training_chunks);
        }
        let mut epoch_chunks = training_chunks.clone();
        let shuffle_seed = collective.broadcast_i64(if main { rng.next() as i64 } else { 0 }, 0)?;
        Rng(shuffle_seed as u64).shuffle(&mut epoch_chunks);
        let random_offset = collective.broadcast_i64(
            if main {
                (rng.next() % options.chunk_size as u64) as i64
            } else {
                0
            },
            0,
        )? as usize;
        let start = (epoch % options.mini_epochs) * chunks_per_epoch;
        let end = (start + chunks_per_epoch).min(epoch_chunks.len());
        let epoch_chunks = if start < end {
            &epoch_chunks[start..end]
        } else {
            &[]
        };
        let training_batches = batches(
            epoch_chunks,
            chunks_per_batch,
            collective.world_size(),
            collective.rank(),
        );
        let loss = run_epoch(
            &archives,
            &training_batches,
            &model,
            Some(&mut optimizer),
            collective,
            random_offset,
            options.chunk_size,
            input_channels,
            options.add_indel_length,
        )?;
        let validation_loss = if options.random_validation || validation_archives.is_some() {
            let validation = validation_archives.as_deref().unwrap_or(&archives);
            let validation_batches = batches(
                &validation_chunks,
                chunks_per_batch,
                collective.world_size(),
                collective.rank(),
            );
            Some(tch::no_grad(|| {
                run_epoch(
                    validation,
                    &validation_batches,
                    &model,
                    None,
                    collective,
                    0,
                    options.chunk_size,
                    input_channels,
                    options.add_indel_length,
                )
            })?)
        } else {
            None
        };
        report.epochs = epoch + 1;
        report.final_loss = loss;
        if main {
            let checkpoint = PathBuf::from(format!(
                "{}.{:02}.pt",
                options.output_prefix.display(),
                epoch + 1
            ));
            write_checkpoint(&model, &checkpoint)?;
            report.checkpoints.push(checkpoint);
            if let Some(log) = &mut log {
                use std::io::Write;
                match validation_loss {
                    Some(value) => writeln!(log, "{}\t{loss:.6}\t{value:.6}", epoch + 1)?,
                    None => writeln!(log, "{}\t{loss:.6}\t", epoch + 1)?,
                }
                log.flush()?;
            }
        }
        if let Some(value) = validation_loss {
            if best_validation.is_none_or(|best| value < best) {
                best_validation = Some(value);
                patience = 0;
                if main {
                    write_checkpoint(&model, &options.best_checkpoint)?;
                }
            } else {
                patience += 1;
                if patience >= 10 * options.mini_epochs {
                    collective.barrier()?;
                    report.best_validation_loss = best_validation;
                    return Ok(report);
                }
            }
        }
        collective.barrier()?;
    }
    report.best_validation_loss = best_validation;
    Ok(report)
}
