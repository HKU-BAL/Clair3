use std::{collections::BTreeMap, path::Path};

use anyhow::{Context, Result, bail};
use tch::{Tensor, nn::VarStore};

use super::pickle::load_python_state_dict;

pub fn load_checkpoint(var_store: &mut VarStore, path: &Path) -> Result<()> {
    let loaded = if path.extension().and_then(|value| value.to_str()) == Some("safetensors") {
        Tensor::read_safetensors(path).map_err(anyhow::Error::from)
    } else {
        load_python_state_dict(path)
    }
    .with_context(|| format!("failed to read checkpoint {}", path.display()))?;

    let mut tensors = BTreeMap::new();
    for (mut name, tensor) in loaded {
        if let Some(stripped) = name.strip_prefix("state_dict.") {
            name = stripped.to_string();
        }
        while let Some(stripped) = name.strip_prefix("module.") {
            name = stripped.to_string();
        }
        tensors.insert(name, tensor);
    }

    let mut variables = var_store.variables();
    let mut missing: Vec<&str> = variables
        .keys()
        .filter(|name| !tensors.contains_key(*name))
        .map(String::as_str)
        .collect();
    let mut unexpected: Vec<&str> = tensors
        .keys()
        .filter(|name| !variables.contains_key(*name))
        .map(String::as_str)
        .collect();
    missing.sort_unstable();
    unexpected.sort_unstable();
    if !missing.is_empty() || !unexpected.is_empty() {
        bail!(
            "checkpoint {} does not match the model; missing: {}{}",
            path.display(),
            missing.join(", "),
            if unexpected.is_empty() {
                String::new()
            } else {
                format!("; unexpected: {}", unexpected.join(", "))
            }
        )
    }

    tch::no_grad(|| {
        for (name, variable) in &mut variables {
            let source = &tensors[name];
            if source.size() != variable.size() {
                bail!(
                    "tensor {name} in {} has shape {:?}, expected {:?}",
                    path.display(),
                    source.size(),
                    variable.size()
                )
            }
            let source = source
                .to_device(var_store.device())
                .to_kind(variable.kind());
            variable
                .f_copy_(&source)
                .with_context(|| format!("failed to copy tensor {name}"))?;
        }
        Ok(())
    })
}
