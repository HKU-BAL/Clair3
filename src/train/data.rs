use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use hdf5_metno::{
    Dataset, Extent, File, SimpleExtents,
    types::{FixedAscii, VarLenAscii},
};
use ndarray::{Array2, ArrayD, IxDyn};
use tch::Tensor;

use crate::{
    config::{FLANKING_BASES, LABEL_COUNT, POSITION_COUNT, Platform},
    io::{compression::open_text, interval::IntervalMap},
    tensor::labels::{labels_from_reference, labels_from_vcf_columns},
};

const POSITION_WIDTH: usize = POSITION_COUNT + 50;
const ALT_INFO_WIDTH: usize = 5000;

#[derive(Debug, Clone, PartialEq)]
pub struct TensorExample {
    pub position: String,
    pub matrix: Vec<i32>,
    pub alt_info: String,
    pub labels: Vec<i8>,
}

pub struct BinArchive {
    _file: File,
    matrix: Dataset,
    labels: Dataset,
    positions: Dataset,
    alt_info: Dataset,
    shape: Vec<usize>,
    rows: usize,
}

impl BinArchive {
    pub fn open(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("failed to open HDF5 training archive {}", path.display()))?;
        let matrix = file.dataset("position_matrix")?;
        let labels = file.dataset("label")?;
        let positions = file.dataset("position")?;
        let alt_info = file.dataset("alt_info")?;
        let shape = matrix.shape();
        let rows = *shape
            .first()
            .context("position_matrix must have a sample dimension")?;
        if labels.shape().first() != Some(&rows)
            || positions.shape().first() != Some(&rows)
            || alt_info.shape().first() != Some(&rows)
            || labels.shape().get(1) != Some(&LABEL_COUNT)
        {
            bail!("inconsistent datasets in {}", path.display());
        }
        Ok(Self {
            _file: file,
            matrix,
            labels,
            positions,
            alt_info,
            shape: shape[1..].to_vec(),
            rows,
        })
    }

    pub fn len(&self) -> usize {
        self.rows
    }

    pub fn is_empty(&self) -> bool {
        self.rows == 0
    }

    pub fn tensor_shape(&self) -> &[usize] {
        &self.shape
    }

    pub fn chunk_count(&self, batch_size: usize, chunk_size: usize) -> usize {
        self.rows.saturating_sub(batch_size) / chunk_size
    }

    pub fn read_batch(
        &self,
        start: usize,
        end: usize,
        input_channels: usize,
    ) -> Result<(Tensor, Tensor)> {
        if start > end || end > self.rows {
            bail!(
                "training batch {start}..{end} is outside {} rows",
                self.rows
            );
        }
        let values = match self.shape.len() {
            2 => self
                .matrix
                .read_slice::<f32, _, IxDyn>((start..end, .., ..))?,
            3 => self
                .matrix
                .read_slice::<f32, _, IxDyn>((start..end, .., .., ..))?,
            rank => bail!("unsupported tensor rank {rank}"),
        };
        let values: Vec<f32> = values.iter().copied().collect();
        let mut shape = vec![(end - start) as i64];
        shape.extend(self.shape.iter().map(|dimension| *dimension as i64));
        let matrix = Tensor::from_slice(&values).reshape(&shape);
        let matrix = if input_channels < *self.shape.last().unwrap_or(&0) {
            matrix.narrow(-1, 0, input_channels as i64).contiguous()
        } else if input_channels == *self.shape.last().unwrap_or(&0) {
            matrix
        } else {
            bail!("model requires more channels than the training archive contains");
        };
        let labels = self.labels.read_slice::<f32, _, IxDyn>((start..end, ..))?;
        let labels: Vec<f32> = labels.iter().copied().collect();
        let labels =
            Tensor::from_slice(&labels).reshape([(end - start) as i64, LABEL_COUNT as i64]);
        Ok((matrix, labels))
    }

    pub fn read_examples(&self, start: usize, end: usize) -> Result<Vec<TensorExample>> {
        let (matrix, labels) = self.read_batch(start, end, *self.shape.last().unwrap_or(&0))?;
        let matrix = Vec::<f32>::try_from(matrix.flatten(0, -1))?;
        let labels = Vec::<f32>::try_from(labels.flatten(0, -1))?;
        let positions = read_strings::<POSITION_WIDTH>(&self.positions, start, end)?;
        let alternate = read_strings::<ALT_INFO_WIDTH>(&self.alt_info, start, end)?;
        let matrix_size: usize = self.shape.iter().product();
        Ok(positions
            .into_iter()
            .zip(alternate)
            .enumerate()
            .map(|(index, (position, alt_info))| TensorExample {
                position,
                matrix: matrix[index * matrix_size..(index + 1) * matrix_size]
                    .iter()
                    .map(|value| *value as i32)
                    .collect(),
                alt_info,
                labels: labels[index * LABEL_COUNT..(index + 1) * LABEL_COUNT]
                    .iter()
                    .map(|value| *value as i8)
                    .collect(),
            })
            .collect())
    }
}

fn read_strings<const N: usize>(
    dataset: &Dataset,
    start: usize,
    end: usize,
) -> Result<Vec<String>> {
    match dataset.read_slice::<FixedAscii<N>, _, IxDyn>((start..end, ..)) {
        Ok(values) => Ok(values.iter().map(ToString::to_string).collect()),
        Err(_) => {
            let values = dataset.read_slice::<VarLenAscii, _, IxDyn>((start..end, ..))?;
            Ok(values.iter().map(ToString::to_string).collect())
        }
    }
}

pub struct BinWriter {
    _file: File,
    matrix: Dataset,
    labels: Dataset,
    positions: Dataset,
    alt_info: Dataset,
    shape: Vec<usize>,
    pileup: bool,
    rows: usize,
}

impl BinWriter {
    pub fn create(path: &Path, shape: &[usize], pileup: bool) -> Result<Self> {
        if !matches!(shape.len(), 2 | 3) || shape.contains(&0) {
            bail!("invalid training tensor shape {shape:?}");
        }
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        let file = File::create(path)?;
        let matrix_shape = std::iter::once(Extent::resizable(0))
            .chain(shape.iter().map(|dimension| Extent::fixed(*dimension)))
            .collect::<Vec<_>>();
        let matrix = if pileup {
            file.new_dataset::<i32>()
                .shape(SimpleExtents::new(matrix_shape))
                .chunk(chunk_shape(shape))
                .blosc_lz4hc(5, true)
                .create("position_matrix")?
        } else {
            file.new_dataset::<i8>()
                .shape(SimpleExtents::new(matrix_shape))
                .chunk(chunk_shape(shape))
                .blosc_lz4hc(5, true)
                .create("position_matrix")?
        };
        let labels = if pileup {
            file.new_dataset::<i32>()
                .shape(SimpleExtents::new([
                    Extent::resizable(0),
                    Extent::fixed(LABEL_COUNT),
                ]))
                .chunk((500, LABEL_COUNT))
                .blosc_lz4hc(5, true)
                .create("label")?
        } else {
            file.new_dataset::<i8>()
                .shape(SimpleExtents::new([
                    Extent::resizable(0),
                    Extent::fixed(LABEL_COUNT),
                ]))
                .chunk((500, LABEL_COUNT))
                .blosc_lz4hc(5, true)
                .create("label")?
        };
        let positions = file
            .new_dataset::<FixedAscii<POSITION_WIDTH>>()
            .shape(SimpleExtents::new([Extent::resizable(0), Extent::fixed(1)]))
            .chunk((500, 1))
            .blosc_lz4hc(5, true)
            .create("position")?;
        let alt_info = file
            .new_dataset::<FixedAscii<ALT_INFO_WIDTH>>()
            .shape(SimpleExtents::new([Extent::resizable(0), Extent::fixed(1)]))
            .chunk((500, 1))
            .blosc_lz4hc(5, true)
            .create("alt_info")?;
        Ok(Self {
            _file: file,
            matrix,
            labels,
            positions,
            alt_info,
            shape: shape.to_vec(),
            pileup,
            rows: 0,
        })
    }

    pub fn append(&mut self, examples: &[TensorExample]) -> Result<()> {
        if examples.is_empty() {
            return Ok(());
        }
        let count = examples.len();
        let start = self.rows;
        let end = start + count;
        let matrix_size: usize = self.shape.iter().product();
        if examples.iter().any(|example| {
            example.matrix.len() != matrix_size || example.labels.len() != LABEL_COUNT
        }) {
            bail!("tensor or label dimensions do not match the training archive");
        }
        let matrix: Vec<i32> = examples
            .iter()
            .flat_map(|example| example.matrix.iter().copied())
            .collect();
        let labels: Vec<i8> = examples
            .iter()
            .flat_map(|example| example.labels.iter().copied())
            .collect();
        let mut matrix_shape = vec![count];
        matrix_shape.extend(self.shape.iter().copied());
        let resized_shape = std::iter::once(end)
            .chain(self.shape.iter().copied())
            .collect::<Vec<_>>();
        self.matrix.resize(resized_shape.as_slice())?;
        if self.pileup {
            let values = ArrayD::from_shape_vec(IxDyn(&matrix_shape), matrix)?;
            self.matrix.write_slice(&values, (start..end, .., ..))?;
        } else {
            let values: Vec<i8> = matrix.into_iter().map(|value| value as i8).collect();
            let values = ArrayD::from_shape_vec(IxDyn(&matrix_shape), values)?;
            self.matrix.write_slice(&values, (start..end, .., .., ..))?;
        }
        self.labels.resize((end, LABEL_COUNT))?;
        if self.pileup {
            let values: Vec<i32> = labels.into_iter().map(i32::from).collect();
            self.labels.write_slice(
                &Array2::from_shape_vec((count, LABEL_COUNT), values)?,
                (start..end, ..),
            )?;
        } else {
            self.labels.write_slice(
                &Array2::from_shape_vec((count, LABEL_COUNT), labels)?,
                (start..end, ..),
            )?;
        }
        append_strings::<POSITION_WIDTH>(
            &self.positions,
            start,
            examples.iter().map(|example| example.position.as_str()),
        )?;
        append_strings::<ALT_INFO_WIDTH>(
            &self.alt_info,
            start,
            examples.iter().map(|example| example.alt_info.as_str()),
        )?;
        self.rows = end;
        Ok(())
    }
}

fn chunk_shape(shape: &[usize]) -> Vec<usize> {
    std::iter::once(500).chain(shape.iter().copied()).collect()
}

fn append_strings<'a, const N: usize>(
    dataset: &Dataset,
    start: usize,
    values: impl Iterator<Item = &'a str>,
) -> Result<()> {
    let values: Vec<_> = values
        .map(|value| {
            let bytes = value.as_bytes();
            FixedAscii::<N>::from_ascii(&bytes[..bytes.len().min(N)])
        })
        .collect::<Result<_, _>>()?;
    let count = values.len();
    let end = start + count;
    dataset.resize((end, 1))?;
    dataset.write_slice(
        &Array2::from_shape_vec((count, 1), values)?,
        (start..end, ..),
    )?;
    Ok(())
}

#[derive(Debug, Clone)]
pub struct TensorBinOptions {
    pub platform: Platform,
    pub pileup: bool,
    pub enable_dwell_time: bool,
    pub allow_duplicate_positions: bool,
    pub maximum_non_variant_ratio: Option<f64>,
    pub shuffle: bool,
    pub seed: u64,
}

impl Default for TensorBinOptions {
    fn default() -> Self {
        Self {
            platform: Platform::Ont,
            pileup: false,
            enable_dwell_time: false,
            allow_duplicate_positions: false,
            maximum_non_variant_ratio: None,
            shuffle: false,
            seed: 0,
        }
    }
}

struct TruthVariant {
    labels: Vec<i8>,
    references: Vec<String>,
    alternates: Vec<String>,
}

fn decode_alternate(reference: &str, alternate: &str) -> (Vec<String>, Vec<String>) {
    let mut references = Vec::new();
    let mut alternates = Vec::new();
    let minimum_length = alternate
        .split(',')
        .map(|allele| allele.len().saturating_sub(1))
        .min()
        .unwrap_or(0)
        .min(reference.len().saturating_sub(1));
    for alternate in alternate.split(',') {
        let common = reference
            .chars()
            .rev()
            .zip(alternate.chars().rev())
            .take(reference.len().min(alternate.len()).saturating_sub(1))
            .take_while(|(left, right)| left == right)
            .count()
            .min(minimum_length);
        let ref_len = reference.len() - common;
        let alt_len = alternate.len() - common;
        references.push(reference[..ref_len].to_string());
        alternates.push(alternate[..alt_len].to_string());
    }
    (references, alternates)
}

fn truth_variants(
    path: &Path,
    regions: Option<&IntervalMap>,
) -> Result<(BTreeMap<String, TruthVariant>, BTreeSet<String>)> {
    let mut variants = BTreeMap::new();
    let mut missing = BTreeSet::new();
    for line in crate::io::compression::read_text_lines(path)? {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let columns: Vec<_> = line.split_whitespace().collect();
        if columns.len() < 6 {
            bail!("invalid truth variant row '{line}'");
        }
        let position: i64 = columns[1].parse()?;
        let key = format!("{}:{position}", columns[0]);
        if columns[4] == "-1" || columns[5] == "-1" {
            missing.insert(key);
            continue;
        }
        if regions.is_some_and(|regions| !regions.overlaps(columns[0], position - 1, position)) {
            continue;
        }
        let (references, alternates) = decode_alternate(columns[2], columns[3]);
        variants.insert(
            key,
            TruthVariant {
                labels: labels_from_vcf_columns(&columns)?,
                references,
                alternates,
            },
        );
    }
    Ok((variants, missing))
}

fn reference_base(base: u8) -> char {
    match base.to_ascii_uppercase() {
        b'A' | b'R' | b'W' | b'M' | b'D' | b'H' | b'V' | b'N' => 'A',
        b'C' | b'Y' | b'S' | b'B' => 'C',
        b'G' | b'K' => 'G',
        _ => 'T',
    }
}

fn has_read_support(truth: &TruthVariant, alt_info: &str) -> bool {
    let entries: Vec<_> = alt_info
        .split_once('-')
        .map(|(_, entries)| entries)
        .unwrap_or("")
        .split_whitespace()
        .collect();
    let mut found = 0;
    for entry in entries.chunks(2) {
        let allele = entry[0];
        if allele.contains(['*', '#']) || allele.starts_with('R') {
            continue;
        }
        let matches = match allele.as_bytes()[0] {
            b'X' => truth
                .alternates
                .iter()
                .any(|alternate| alternate == &allele[1..]),
            b'I' => truth
                .alternates
                .iter()
                .any(|alternate| alternate == &allele[1..]),
            b'D' => truth
                .references
                .iter()
                .zip(&truth.alternates)
                .any(|(reference, alternate)| {
                    alternate.len() == 1 && reference.get(1..) == Some(&allele[1..])
                }),
            _ => false,
        };
        found += usize::from(matches);
    }
    found >= truth.alternates.len()
}

#[derive(Debug)]
struct SplitMix64(u64);

impl SplitMix64 {
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

#[allow(clippy::too_many_arguments)]
fn parse_example(
    line: &str,
    options: &TensorBinOptions,
    truth: &BTreeMap<String, TruthVariant>,
    missing: &BTreeSet<String>,
    regions: Option<&IntervalMap>,
) -> Result<Option<TensorExample>> {
    let fields: Vec<_> = line.split('\t').collect();
    if fields.len() != 5 {
        bail!("expected five tensor columns, got {}", fields.len());
    }
    let position: i64 = fields[1].parse()?;
    if regions.is_some_and(|regions| !regions.overlaps(fields[0], position - 1, position)) {
        return Ok(None);
    }
    let key = format!("{}:{position}", fields[0]);
    if missing.contains(&key) {
        return Ok(None);
    }
    let sequence = fields[2].to_ascii_uppercase();
    let reference = *sequence.as_bytes().get(FLANKING_BASES).unwrap_or(&b'N');
    if !matches!(reference, b'A' | b'C' | b'G' | b'T') {
        return Ok(None);
    }
    if let Some(truth) = truth.get(&key)
        && !has_read_support(truth, fields[4])
    {
        return Ok(None);
    }
    let channels = if options.pileup {
        18
    } else if options.enable_dwell_time {
        9
    } else {
        8
    };
    let mut matrix: Vec<i32> = fields[3]
        .split_whitespace()
        .map(|value| value.parse::<i32>().map_err(anyhow::Error::from))
        .collect::<Result<_>>()?;
    if options.pileup {
        if matrix.len() != POSITION_COUNT * channels {
            bail!("invalid pileup tensor size {}", matrix.len());
        }
        let depth = fields[4].split('-').next().unwrap_or("0").parse::<f64>()?;
        if depth > 144.0 * 1.5 {
            let factor = depth / 144.0;
            matrix = matrix
                .into_iter()
                .map(|value| (value as f64 / factor) as i32)
                .collect();
        }
    } else {
        let row_size = POSITION_COUNT * channels;
        if !matrix.len().is_multiple_of(row_size)
            || matrix.len() / row_size > options.platform.matrix_depth()
        {
            bail!("invalid full-alignment tensor size {}", matrix.len());
        }
        let depth = matrix.len() / row_size;
        let padding = options.platform.matrix_depth() - depth;
        let mut padded = vec![0; (padding / 2) * row_size];
        padded.extend(matrix);
        padded.extend(std::iter::repeat_n(0, padding.div_ceil(2) * row_size));
        matrix = padded;
    }
    let labels = if let Some(truth) = truth.get(&key) {
        truth.labels.clone()
    } else {
        labels_from_reference(&reference_base(reference).to_string())?
    };
    Ok(Some(TensorExample {
        position: format!("{key}:{sequence}"),
        matrix,
        alt_info: fields[4].to_string(),
        labels,
    }))
}

fn write_examples(
    writer: &mut BinWriter,
    examples: &mut Vec<TensorExample>,
    reference_indices: &mut Vec<usize>,
    options: &TensorBinOptions,
    random: &mut SplitMix64,
) -> Result<()> {
    if let Some(ratio) = options.maximum_non_variant_ratio {
        let variants = examples.len() - reference_indices.len();
        let keep =
            (ratio * variants as f64 / reference_indices.len().max(1) as f64).clamp(0.0, 1.0);
        let remove: BTreeSet<usize> = reference_indices
            .iter()
            .copied()
            .filter(|_| random.next() as f64 / u64::MAX as f64 > keep)
            .collect();
        let mut index = 0;
        examples.retain(|_| {
            let keep = !remove.contains(&index);
            index += 1;
            keep
        });
    }
    if options.shuffle {
        random.shuffle(examples);
    }
    writer.append(examples)?;
    examples.clear();
    reference_indices.clear();
    Ok(())
}

pub fn tensor_to_bin(
    tensor_path: &Path,
    truth_path: &Path,
    bed_path: Option<&Path>,
    output: &Path,
    options: &TensorBinOptions,
    chunk_id: Option<usize>,
    chunk_num: Option<usize>,
) -> Result<usize> {
    if chunk_id.is_some() != chunk_num.is_some() {
        bail!("chunk ID and chunk count must be supplied together");
    }
    let regions = bed_path
        .map(|path| IntervalMap::from_bed(path, None))
        .transpose()?;
    let (truth, missing) = truth_variants(truth_path, regions.as_ref())?;
    let mut paths = Vec::new();
    if tensor_path == Path::new("PIPE") || tensor_path == Path::new("-") || tensor_path.is_file() {
        paths.push(tensor_path.to_path_buf());
    } else {
        let parent = tensor_path.parent().unwrap_or_else(|| Path::new("."));
        let prefix = tensor_path
            .file_name()
            .context("tensor path has no filename")?
            .to_string_lossy();
        for entry in fs::read_dir(parent)? {
            let path = entry?.path();
            let name = path
                .file_name()
                .and_then(|name| name.to_str())
                .unwrap_or("");
            if name.starts_with(&format!("{prefix}_")) || name.starts_with(&format!("{prefix}.")) {
                paths.push(path);
            }
        }
        paths.sort();
        if let (Some(id), Some(count)) = (chunk_id, chunk_num) {
            if count == 0 || id == 0 || id > count {
                bail!("chunk ID must be between 1 and the chunk count");
            }
            let size = paths.len().div_ceil(count);
            paths = paths.into_iter().skip((id - 1) * size).take(size).collect();
        }
    }
    let shape = if options.pileup {
        vec![POSITION_COUNT, 18]
    } else {
        vec![
            options.platform.matrix_depth(),
            POSITION_COUNT,
            if options.enable_dwell_time { 9 } else { 8 },
        ]
    };
    let mut writer = BinWriter::create(output, &shape, options.pileup)?;
    let mut examples = Vec::new();
    let mut reference_indices = Vec::new();
    let mut random = SplitMix64(options.seed);
    let mut used_positions = BTreeSet::new();
    let mut total = 0;
    for path in paths {
        let reader = open_text(&path)?;
        for line in std::io::BufRead::lines(reader) {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let Some(example) = parse_example(&line, options, &truth, &missing, regions.as_ref())?
            else {
                continue;
            };
            let position_key = example
                .position
                .rsplit_once(':')
                .context("invalid tensor position")?
                .0
                .to_string();
            if !used_positions.insert(position_key.clone()) && !options.allow_duplicate_positions {
                continue;
            }
            if !truth.contains_key(&position_key) {
                reference_indices.push(examples.len());
            }
            examples.push(example);
            total += 1;
            if examples.len() == 50_000 {
                write_examples(
                    &mut writer,
                    &mut examples,
                    &mut reference_indices,
                    options,
                    &mut random,
                )?;
            }
        }
    }
    write_examples(
        &mut writer,
        &mut examples,
        &mut reference_indices,
        options,
        &mut random,
    )?;
    Ok(total)
}

pub fn merge_bins(inputs: &[PathBuf], output: &Path, pileup: bool) -> Result<usize> {
    let mut archives: Vec<BinArchive> = inputs
        .iter()
        .map(|path| BinArchive::open(path))
        .collect::<Result<_>>()?;
    let first = archives
        .first()
        .context("no training archives were provided")?;
    let shape = first.tensor_shape().to_vec();
    if archives
        .iter()
        .any(|archive| archive.tensor_shape() != shape)
    {
        bail!("training archives have different tensor shapes");
    }
    let mut writer = BinWriter::create(output, &shape, pileup)?;
    let mut count = 0;
    for archive in &mut archives {
        for start in (0..archive.len()).step_by(500) {
            let end = (start + 500).min(archive.len());
            let examples = archive.read_examples(start, end)?;
            count += examples.len();
            writer.append(&examples)?;
        }
    }
    Ok(count)
}
