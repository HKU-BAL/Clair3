use std::{
    fmt::Write as _,
    fs,
    io::BufRead,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use tch::Tensor;

use crate::{
    call::{
        caller::{TensorInput, call_tensor_batch},
        decode::DecodeOptions,
    },
    config::{FLANKING_BASES, PREDICT_BATCH_SIZE, Platform},
    io::{compression::read_text_lines, reference::Reference, variant::VcfRecord},
    model::{Model, ModelKind},
    pipeline::{
        call::{CallConfig, full_alignment_batch, pileup_batch},
        gvcf::nonvariant_records,
    },
};

#[derive(Debug, Default)]
pub struct TensorData {
    pub shape: Vec<usize>,
    pub values: Vec<i32>,
    pub positions: Vec<String>,
    pub alt_info: Vec<String>,
    pub nonvariants: Vec<VcfRecord>,
}

impl TensorData {
    pub fn positions_with_reference(&self, reference: &Reference) -> Result<Vec<String>> {
        self.positions
            .iter()
            .map(|position| {
                let mut fields = position.rsplitn(3, ':');
                let reference_base = fields.next().context("missing reference base")?;
                let position_one_based: i64 = fields.next().context("missing position")?.parse()?;
                let contig = fields.next().context("missing contig")?;
                let sequence = if reference_base.len() == 1 {
                    String::from_utf8(reference.fetch_padded(
                        contig,
                        position_one_based - 1 - FLANKING_BASES as i64,
                        position_one_based + FLANKING_BASES as i64,
                        b'A',
                    )?)?
                } else {
                    reference_base.to_string()
                };
                if sequence.len() != 2 * FLANKING_BASES + 1 {
                    bail!("invalid reference sequence length for {contig}:{position_one_based}");
                }
                Ok(format!("{contig}:{position_one_based}:{sequence}"))
            })
            .collect()
    }

    pub fn to_text(&self, reference: &Reference) -> Result<String> {
        let positions = self.positions_with_reference(reference)?;
        let sample_size = self.shape.get(1..).unwrap_or(&[]).iter().product::<usize>();
        if sample_size == 0
            || self.values.len() != positions.len() * sample_size
            || positions.len() != self.alt_info.len()
        {
            bail!("tensor shape, values, positions, and alternate information are inconsistent");
        }
        let mut output = String::new();
        for (index, (position, alt_info)) in positions.iter().zip(&self.alt_info).enumerate() {
            let (contig_position, sequence) = position
                .rsplit_once(':')
                .context("invalid tensor reference")?;
            let (contig, position) = contig_position
                .rsplit_once(':')
                .context("invalid tensor position")?;
            let values = &self.values[index * sample_size..(index + 1) * sample_size];
            write!(output, "{contig}\t{position}\t{sequence}\t")?;
            for (index, value) in values.iter().enumerate() {
                if index > 0 {
                    output.push(' ');
                }
                write!(output, "{value}")?;
            }
            writeln!(output, "\t{alt_info}")?;
        }
        Ok(output)
    }

    pub fn write_npy(&self, path: &Path, _reference: &Reference) -> Result<()> {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        let shape: Vec<i64> = self.shape.iter().map(|value| *value as i64).collect();
        let values: Vec<i8> = self.values.iter().map(|value| *value as i8).collect();
        Tensor::from_slice(&values)
            .reshape(&shape)
            .write_npy(path)?;
        let mut info = String::new();
        for (position, alt_info) in self.positions.iter().zip(&self.alt_info) {
            writeln!(info, "{position}\t{alt_info}")?;
        }
        fs::write(PathBuf::from(format!("{}.info", path.display())), info)?;
        Ok(())
    }

    pub fn from_npy(path: &Path) -> Result<Self> {
        let tensor = Tensor::read_npy(path)?;
        let shape = tensor
            .size()
            .into_iter()
            .map(|value| value as usize)
            .collect::<Vec<_>>();
        let tensor = tensor.flatten(0, -1);
        let values = match tensor.kind() {
            tch::Kind::Int8 => Vec::<i8>::try_from(&tensor)?
                .into_iter()
                .map(i32::from)
                .collect(),
            tch::Kind::Int => Vec::<i32>::try_from(&tensor)?,
            kind => bail!("unsupported tensor dtype {kind:?} in {}", path.display()),
        };
        let mut positions = Vec::new();
        let mut alt_info = Vec::new();
        for line in read_text_lines(&PathBuf::from(format!("{}.info", path.display())))? {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let (position, alternate) = line.split_once('\t').context("invalid tensor info row")?;
            positions.push(position.to_string());
            alt_info.push(alternate.to_string());
        }
        Ok(Self {
            shape,
            values,
            positions,
            alt_info,
            nonvariants: Vec::new(),
        })
    }

    pub fn from_text(
        reader: impl BufRead,
        pileup: bool,
        platform: Platform,
        dwell_time: bool,
    ) -> Result<Self> {
        let mut values = Vec::new();
        let mut positions = Vec::new();
        let mut alt_info = Vec::new();
        let channels = if pileup {
            18
        } else if dwell_time {
            9
        } else {
            8
        };
        for line in reader.lines() {
            let line = line?;
            if line.is_empty() {
                continue;
            }
            let fields: Vec<_> = line.trim_end_matches('\r').split('\t').collect();
            if fields.len() != 5 {
                bail!("expected five tensor columns, got {}", fields.len());
            }
            positions.push(format!("{}:{}:{}", fields[0], fields[1], fields[2]));
            alt_info.push(fields[4].to_string());
            let row: Vec<i32> = fields[3]
                .split_whitespace()
                .map(|value| value.parse::<i32>().map_err(anyhow::Error::from))
                .collect::<Result<_>>()?;
            if pileup {
                if row.len() != 33 * channels {
                    bail!("invalid pileup tensor size {}", row.len());
                }
                values.extend(row);
            } else {
                let row_size = 33 * channels;
                if !row.len().is_multiple_of(row_size)
                    || row.len() / row_size > platform.matrix_depth()
                {
                    bail!("invalid full-alignment tensor size {}", row.len());
                }
                let padding = platform.matrix_depth() - row.len() / row_size;
                values.extend(std::iter::repeat_n(0, (padding / 2) * row_size));
                values.extend(row);
                values.extend(std::iter::repeat_n(0, padding.div_ceil(2) * row_size));
            }
        }
        let shape = if pileup {
            vec![positions.len(), 33, channels]
        } else {
            vec![positions.len(), platform.matrix_depth(), 33, channels]
        };
        Ok(Self {
            shape,
            values,
            positions,
            alt_info,
            nonvariants: Vec::new(),
        })
    }
}

pub fn create_tensor_data(config: &CallConfig, pileup: bool) -> Result<TensorData> {
    let reference = Reference::open(&config.reference)?;
    if pileup {
        let (batch, region) = pileup_batch(config, &reference)?;
        let nonvariants = if config.gvcf && region.start != region.end {
            nonvariant_records(
                &reference,
                &crate::io::interval::Region::new(
                    &region.contig,
                    (region.start - 33).max(0),
                    (region.end + 33).min(reference.contig_len(&region.contig) as i64),
                )?,
                &batch,
                &config.gvcf_options,
            )?
            .into_iter()
            .filter(|record| record.position > region.start && record.position < region.end)
            .collect()
        } else {
            Vec::new()
        };
        let values = batch.values;
        Ok(TensorData {
            shape: batch.shape,
            values,
            positions: batch.positions,
            alt_info: batch.alt_info,
            nonvariants,
        })
    } else {
        let batch = full_alignment_batch(config, &reference)?;
        Ok(TensorData {
            shape: batch.shape,
            values: batch.values.into_iter().map(i32::from).collect(),
            positions: batch.positions,
            alt_info: batch.alt_info,
            nonvariants: Vec::new(),
        })
    }
}

pub fn call_tensor_data(
    checkpoint: &Path,
    kind: ModelKind,
    data: &TensorData,
    options: &DecodeOptions,
    device: tch::Device,
) -> Result<Vec<VcfRecord>> {
    if data.positions.is_empty() {
        return Ok(Vec::new());
    }
    let channels = *data
        .shape
        .last()
        .context("tensor has no channel dimension")?;
    let model = Model::new(
        kind,
        device,
        Some(channels as i64),
        options.add_indel_length,
        Some(checkpoint),
    )?;
    tch::no_grad(|| {
        call_tensor_batch(
            &model,
            TensorInput {
                shape: &data.shape,
                values: &data.values,
                positions: &data.positions,
                alt_info: &data.alt_info,
            },
            options,
            PREDICT_BATCH_SIZE,
        )
    })
}
