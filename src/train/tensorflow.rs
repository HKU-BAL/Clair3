use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use tch::{Device, Tensor};

use crate::{
    config::Platform,
    model::{Model, ModelKind},
    train::checkpoint::save_checkpoint,
};

#[derive(Clone, Copy)]
enum Transform {
    Dense,
    Lstm,
    Conv,
    Identity,
    Zeros,
}

#[derive(Clone)]
struct Mapping {
    source: String,
    target: String,
    shape: Vec<i64>,
    transform: Transform,
}

impl Mapping {
    fn new(source: String, target: String, shape: Vec<i64>, transform: Transform) -> Self {
        Self {
            source,
            target,
            shape,
            transform,
        }
    }

    fn source_shape(&self) -> Vec<i64> {
        match self.transform {
            Transform::Dense | Transform::Lstm => vec![self.shape[1], self.shape[0]],
            Transform::Conv => vec![self.shape[2], self.shape[3], self.shape[1], self.shape[0]],
            Transform::Identity | Transform::Zeros => self.shape.clone(),
        }
    }

    fn apply(&self, tensor: Tensor) -> Tensor {
        match self.transform {
            Transform::Dense | Transform::Lstm => tensor.transpose(0, 1).contiguous(),
            Transform::Conv => tensor.permute([3, 2, 0, 1]).contiguous(),
            Transform::Identity => tensor,
            Transform::Zeros => Tensor::zeros(&self.shape, (tch::Kind::Float, Device::Cpu)),
        }
    }
}

fn shape(variables: &BTreeMap<String, Tensor>, name: &str) -> Result<Vec<i64>> {
    let tensor = variables
        .get(name)
        .with_context(|| format!("missing model tensor {name}"))?;
    Ok(tensor.size())
}

fn dense_specs(variables: &BTreeMap<String, Tensor>, name: &str) -> Result<Vec<Mapping>> {
    Ok(vec![
        Mapping::new(
            format!("{name}/kernel"),
            format!("{name}.weight"),
            shape(variables, &format!("{name}.weight"))?,
            Transform::Dense,
        ),
        Mapping::new(
            format!("{name}/bias"),
            format!("{name}.bias"),
            shape(variables, &format!("{name}.bias"))?,
            Transform::Identity,
        ),
    ])
}

fn lstm_specs(variables: &BTreeMap<String, Tensor>, layer: usize) -> Result<Vec<Mapping>> {
    let suffix = if layer == 2 { "_1" } else { "" };
    let torch = format!("LSTM{layer}");
    let mut output = Vec::new();
    for (direction, reverse) in [("forward", false), ("backward", true)] {
        let suffix = format!("{direction}_lstm{suffix}");
        let torch_suffix = if reverse { "_reverse" } else { "" };
        for (tf_name, torch_name, transform) in [
            ("kernel", "weight_ih_l0", Transform::Lstm),
            ("recurrent_kernel", "weight_hh_l0", Transform::Lstm),
            ("bias", "bias_ih_l0", Transform::Identity),
            ("bias", "bias_hh_l0", Transform::Zeros),
        ] {
            let target = format!("{torch}.{torch_name}{torch_suffix}");
            output.push(Mapping::new(
                format!("LSTM{layer}/{suffix}/{tf_name}"),
                target.clone(),
                shape(variables, &target)?,
                transform,
            ));
        }
    }
    Ok(output)
}

fn bn_specs(
    variables: &BTreeMap<String, Tensor>,
    source: &str,
    target: &str,
) -> Result<Vec<Mapping>> {
    let mut output = Vec::new();
    for (source_name, target_name) in [
        ("gamma", "weight"),
        ("beta", "bias"),
        ("moving_mean", "running_mean"),
        ("moving_variance", "running_var"),
    ] {
        let target_name = format!("{target}.{target_name}");
        output.push(Mapping::new(
            format!("{source}/{source_name}"),
            target_name.clone(),
            shape(variables, &target_name)?,
            Transform::Identity,
        ));
    }
    Ok(output)
}

fn conv_specs(
    variables: &BTreeMap<String, Tensor>,
    source: &str,
    target: &str,
    transform: Transform,
) -> Result<Vec<Mapping>> {
    Ok(vec![
        Mapping::new(
            format!("{source}/kernel"),
            format!("{target}.weight"),
            shape(variables, &format!("{target}.weight"))?,
            transform,
        ),
        Mapping::new(
            format!("{source}/bias"),
            format!("{target}.bias"),
            shape(variables, &format!("{target}.bias"))?,
            Transform::Identity,
        ),
    ])
}

fn full_specs(variables: &BTreeMap<String, Tensor>) -> Result<Vec<Mapping>> {
    let mut output = Vec::new();
    output.extend(conv_specs(
        variables,
        "conv1/conv",
        "conv1.conv",
        Transform::Conv,
    )?);
    output.extend(bn_specs(variables, "conv1/bn", "conv1.bn")?);
    for (residual, conv) in [
        ("res_block1", None),
        ("res_block2", Some("conv3")),
        ("res_block3", Some("conv5")),
    ] {
        if let Some(conv) = conv {
            output.extend(conv_specs(
                variables,
                &format!("{conv}/conv"),
                &format!("{conv}.conv"),
                Transform::Conv,
            )?);
            output.extend(bn_specs(
                variables,
                &format!("{conv}/bn"),
                &format!("{conv}.bn"),
            )?);
        }
        output.extend(conv_specs(
            variables,
            &format!("{residual}/basic_block/conv1"),
            &format!("{residual}.0.conv1"),
            Transform::Conv,
        )?);
        output.extend(bn_specs(
            variables,
            &format!("{residual}/basic_block/bn1"),
            &format!("{residual}.0.bn1"),
        )?);
        output.extend(conv_specs(
            variables,
            &format!("{residual}/basic_block/conv2"),
            &format!("{residual}.0.conv2"),
            Transform::Conv,
        )?);
        output.extend(bn_specs(
            variables,
            &format!("{residual}/basic_block/bn2"),
            &format!("{residual}.0.bn2"),
        )?);
    }
    Ok(output)
}

fn mappings(
    kind: ModelKind,
    add_indel_length: bool,
    variables: &BTreeMap<String, Tensor>,
) -> Result<Vec<Mapping>> {
    let mut output = Vec::new();
    if kind == ModelKind::Pileup {
        output.extend(lstm_specs(variables, 1)?);
        output.extend(lstm_specs(variables, 2)?);
    } else if kind == ModelKind::FullAlignment {
        output.extend(full_specs(variables)?);
    } else {
        bail!("TensorFlow conversion is supported for pileup and full-alignment models");
    }
    for name in ["L4", "L5_1", "L5_2", "Y_gt21_logits", "Y_genotype_logits"] {
        output.extend(dense_specs(variables, name)?);
    }
    if add_indel_length {
        for name in [
            "L5_3",
            "L5_4",
            "Y_indel_length_logits_1",
            "Y_indel_length_logits_2",
        ] {
            output.extend(dense_specs(variables, name)?);
        }
    }
    Ok(output)
}

fn resolve_checkpoint(path: &Path) -> Result<PathBuf> {
    if path.is_dir() {
        let index = path.join("checkpoint");
        let content = fs::read_to_string(&index)
            .with_context(|| format!("no TensorFlow checkpoint index under {}", path.display()))?;
        let prefix = content
            .lines()
            .find_map(|line| {
                line.trim()
                    .strip_prefix("model_checkpoint_path:")?
                    .trim_start()
                    .strip_prefix('"')?
                    .split('"')
                    .next()
            })
            .context("TensorFlow checkpoint state has no model path")?;
        return Ok(path.join(prefix));
    }
    Ok(
        if path.extension().and_then(|value| value.to_str()) == Some("index") {
            path.with_extension("")
        } else {
            path.to_path_buf()
        },
    )
}

fn checkpoint_sources(source: &str) -> Vec<String> {
    let mut paths = vec![source.to_owned()];
    for layer in 1..=2 {
        for direction in ["forward", "backward"] {
            for source_suffix in ["", "_1"] {
                let prefix = format!("LSTM{layer}/{direction}_lstm{source_suffix}/");
                if let Some(rest) = source.strip_prefix(&prefix) {
                    for target_suffix in ["", "_1"] {
                        paths.push(format!(
                            "LSTM{layer}/{direction}_layer{target_suffix}/cell/{rest}"
                        ));
                    }
                }
            }
        }
    }
    for layer in 1..=3 {
        let prefix = format!("res_block{layer}/basic_block/");
        if let Some(rest) = source.strip_prefix(&prefix) {
            paths.push(format!("res_block{layer}/_functional/_operations/1/{rest}"));
        }
    }
    let mut sources = Vec::new();
    for path in paths {
        for path in [path.clone(), path.replacen("/kernel", "/_kernel", 1)] {
            sources.push(format!("model/{path}/.ATTRIBUTES/VARIABLE_VALUE"));
            sources.push(format!("{path}/.ATTRIBUTES/VARIABLE_VALUE"));
        }
    }
    sources.sort();
    sources.dedup();
    sources
}

struct CheckpointReader(super::tensorflow_table::Checkpoint);

impl CheckpointReader {
    fn open(path: &Path) -> Result<Self> {
        Ok(Self(super::tensorflow_table::Checkpoint::open(path)?))
    }

    fn read(&self, mapping: &Mapping) -> Result<Tensor> {
        let mut errors = Vec::new();
        for source in checkpoint_sources(&mapping.source) {
            match self.0.read_tensor(&source) {
                Ok((dimensions, values)) => {
                    if dimensions != mapping.source_shape() {
                        bail!(
                            "tensor {source} has shape {dimensions:?}, expected {:?}",
                            mapping.source_shape()
                        );
                    }
                    return Ok(Tensor::from_slice(&values).reshape(&dimensions));
                }
                Err(error) => errors.push(format!("{source}: {error:#}")),
            }
        }
        bail!(
            "failed to read tensor {}: {}",
            mapping.source,
            errors.join("; ")
        )
    }
}

fn restore(path: &Path, mappings: &[Mapping]) -> Result<Vec<(String, Tensor)>> {
    let reader = CheckpointReader::open(&resolve_checkpoint(path)?)?;
    mappings
        .iter()
        .map(|mapping| Ok((mapping.target.clone(), mapping.apply(reader.read(mapping)?))))
        .collect()
}

pub fn convert_checkpoint(
    checkpoint: &Path,
    output: &Path,
    kind: ModelKind,
    platform: Platform,
    add_indel_length: bool,
    enable_dwell_time: bool,
) -> Result<()> {
    if kind == ModelKind::Pileup && enable_dwell_time {
        bail!("dwell time is supported only for full-alignment conversion");
    }
    let model = Model::new(
        kind,
        Device::Cpu,
        Some(if kind == ModelKind::Pileup {
            18
        } else if enable_dwell_time {
            9
        } else {
            8
        }),
        add_indel_length,
        None,
    )?;
    let _ = platform;
    let variables: BTreeMap<_, _> = model.var_store.variables().into_iter().collect();
    let mappings = mappings(kind, add_indel_length, &variables)?;
    let restored = restore(checkpoint, &mappings)?;
    let mut variables = model.var_store.variables();
    tch::no_grad(|| {
        for (name, tensor) in restored {
            let variable = variables
                .get_mut(&name)
                .with_context(|| format!("missing converted tensor {name}"))?;
            if variable.size() != tensor.size() {
                bail!(
                    "converted tensor {name} has shape {:?}, expected {:?}",
                    tensor.size(),
                    variable.size()
                );
            }
            variable.f_copy_(&tensor)?;
        }
        Ok(())
    })?;
    if let Some(parent) = output.parent()
        && !parent.as_os_str().is_empty()
    {
        fs::create_dir_all(parent)?;
    }
    save_checkpoint(&model.var_store, output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolves_object_graph_names() {
        let names = checkpoint_sources("LSTM1/backward_lstm/kernel");
        assert!(
            names.contains(
                &"model/LSTM1/backward_layer/cell/kernel/.ATTRIBUTES/VARIABLE_VALUE".to_owned()
            ),
            "{names:?}"
        );
        let names = checkpoint_sources("res_block3/basic_block/conv2/kernel");
        assert!(names.contains(&"model/res_block3/_functional/_operations/1/conv2/_kernel/.ATTRIBUTES/VARIABLE_VALUE".to_owned()), "{names:?}");
    }

    #[test]
    fn resolves_checkpoint_state() {
        let directory = tempfile::tempdir().unwrap();
        fs::write(
            directory.path().join("checkpoint"),
            "model_checkpoint_path: \"model.v2\"\nall_model_checkpoint_paths: \"model.v1\"\nall_model_checkpoint_paths: \"model.v2\"\n",
        )
        .unwrap();
        assert_eq!(
            resolve_checkpoint(directory.path()).unwrap(),
            directory.path().join("model.v2")
        );
    }

    #[test]
    #[ignore = "requires the TensorFlow test checkpoint"]
    fn restores_dense_checkpoint() {
        let checkpoint =
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/tf-checkpoints/pileup-1");
        let mappings = vec![Mapping::new(
            "L4/kernel".to_string(),
            "L4.weight".to_string(),
            vec![128, 10560],
            Transform::Dense,
        )];
        let restored = restore(&checkpoint, &mappings).unwrap();
        assert_eq!(restored.len(), 1);
        assert_eq!(restored[0].1.size(), vec![128, 10560]);
        assert!(restored[0].1.isfinite().all().int64_value(&[]) == 1);
    }
}
