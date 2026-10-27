use std::{
    fs::File,
    io::{BufWriter, Write},
    path::Path,
};

use anyhow::{Context, Result, bail};
use tch::{Device, Tensor};

use crate::{
    io::{
        reference::Reference,
        variant::{VcfRecord, vcf_header},
    },
    model::Model,
};

use super::decode::{DecodeOptions, Prediction, decode_variant};

#[derive(Debug, Clone, Copy)]
pub struct TensorInput<'a, T = i32> {
    pub shape: &'a [usize],
    pub values: &'a [T],
    pub positions: &'a [String],
    pub alt_info: &'a [String],
}

pub fn call_tensor_batch<T: tch::kind::Element>(
    model: &Model,
    input: TensorInput<'_, T>,
    options: &DecodeOptions,
    batch_size: usize,
) -> Result<Vec<VcfRecord>> {
    if batch_size == 0 {
        bail!("prediction batch size must be positive");
    }
    if input.shape.is_empty()
        || input.shape[0] != input.positions.len()
        || input.positions.len() != input.alt_info.len()
        || input.shape.iter().product::<usize>() != input.values.len()
    {
        bail!("inconsistent tensor shape, position, and alternate-information sizes");
    }
    let count = input.shape[0];
    let shape: Vec<i64> = input.shape.iter().map(|size| *size as i64).collect();
    let device = model.var_store.device();
    let mut records = Vec::new();
    for start in (0..count).step_by(batch_size) {
        let end = (start + batch_size).min(count);
        let batch_shape: Vec<i64> = std::iter::once((end - start) as i64)
            .chain(shape[1..].iter().copied())
            .collect();
        let item_size = input.shape[1..].iter().product::<usize>();
        let mut tensor = Tensor::from_slice(&input.values[start * item_size..end * item_size])
            .reshape(&batch_shape);
        if options.pileup {
            let scales = input.alt_info[start..end]
                .iter()
                .map(|alt| {
                    let depth = alt.split('-').next().unwrap_or("0").parse::<i64>()?;
                    Ok(if depth > 0 && depth as f64 > 144.0 * 1.5 {
                        depth as f64 / 144.0
                    } else {
                        1.0
                    })
                })
                .collect::<Result<Vec<f64>>>()?;
            let mut scale_shape = vec![1_i64; tensor.dim()];
            scale_shape[0] = (end - start) as i64;
            tensor = (tensor / Tensor::from_slice(&scales).reshape(&scale_shape))
                .to_kind(tch::Kind::Int);
        }
        let tensor = tensor.to_device(device);
        let outputs = model
            .forward_tasks(&tensor, false)
            .into_iter()
            .map(|tensor| {
                let tensor = tensor.to_device(Device::Cpu).flatten(0, -1);
                Vec::<f32>::try_from(&tensor)
            })
            .collect::<Result<Vec<_>, _>>()?;
        if outputs.len() < 2 || (options.add_indel_length && outputs.len() != 4) {
            bail!("model output does not match the requested indel-length mode");
        }
        for (local_index, index) in (start..end).enumerate() {
            if options.pileup
                && input.alt_info[index]
                    .split('-')
                    .next()
                    .unwrap_or("0")
                    .parse::<i64>()?
                    <= 0
            {
                continue;
            }
            let mut parts = input.positions[index].rsplitn(3, ':');
            let reference_sequence = parts.next().context("missing reference sequence")?;
            let position: i64 = parts
                .next()
                .context("missing candidate position")?
                .parse()?;
            let contig = parts
                .next()
                .context("missing candidate contig")?
                .to_string();
            let center = if reference_sequence.len() == 1 { 0 } else { 16 };
            let prediction = Prediction {
                gt21: outputs[0][local_index * 21..local_index * 21 + 21].to_vec(),
                genotype: outputs[1][local_index * 3..local_index * 3 + 3].to_vec(),
                length1: if options.add_indel_length {
                    outputs[2][local_index * 33..local_index * 33 + 33].to_vec()
                } else {
                    Vec::new()
                },
                length2: if options.add_indel_length {
                    outputs[3][local_index * 33..local_index * 33 + 33].to_vec()
                } else {
                    Vec::new()
                },
            };
            if let Some(record) = decode_variant(
                &contig,
                position,
                reference_sequence,
                center,
                &input.alt_info[index],
                &prediction,
                options,
            )? {
                records.push(record);
            }
        }
    }
    Ok(records)
}

pub fn write_vcf(
    output: &Path,
    reference: &Reference,
    records: &[VcfRecord],
    sample_name: &str,
    command: Option<&str>,
    gvcf: bool,
) -> Result<()> {
    let contigs = reference.contigs()?;
    let header = vcf_header(reference, command, sample_name, gvcf, &contigs)?;
    let mut writer: Box<dyn Write> = if output == Path::new("PIPE") || output == Path::new("-") {
        Box::new(std::io::stdout())
    } else {
        Box::new(BufWriter::new(File::create(output)?))
    };
    writer.write_all(header.as_bytes())?;
    for record in records {
        writeln!(writer, "{}", record.to_line())?;
    }
    writer.flush()?;
    if records.is_empty() && output != Path::new("PIPE") && output != Path::new("-") {
        std::fs::remove_file(output)?;
    }
    Ok(())
}
