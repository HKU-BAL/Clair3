use std::{fs::File, io::Write, path::Path};

use anyhow::{Context, Result, bail};
use tch::{Device, Kind, Tensor, nn::VarStore};
use zip::{CompressionMethod, ZipWriter, write::FileOptions};

fn unicode(output: &mut Vec<u8>, value: &str) {
    output.push(b'X');
    output.extend((value.len() as u32).to_le_bytes());
    output.extend(value.as_bytes());
}

fn integer(output: &mut Vec<u8>, value: i64) {
    if (0..=u8::MAX as i64).contains(&value) {
        output.extend([b'K', value as u8]);
    } else if (0..=u16::MAX as i64).contains(&value) {
        output.push(b'M');
        output.extend((value as u16).to_le_bytes());
    } else if i32::try_from(value).is_ok() {
        output.push(b'J');
        output.extend((value as i32).to_le_bytes());
    } else {
        output.extend([0x8a, 8]);
        output.extend(value.to_le_bytes());
    }
}

fn global(output: &mut Vec<u8>, module: &str, name: &str) {
    output.push(b'c');
    output.extend(module.as_bytes());
    output.push(b'\n');
    output.extend(name.as_bytes());
    output.push(b'\n');
}

fn ordered_dict(output: &mut Vec<u8>) {
    global(output, "collections", "OrderedDict");
    output.extend(b")R");
}

fn dimensions(output: &mut Vec<u8>, dimensions: &[i64]) {
    output.push(b'(');
    for dimension in dimensions {
        integer(output, *dimension);
    }
    output.push(b't');
}

fn strides(size: &[i64]) -> Vec<i64> {
    let mut strides = vec![1; size.len()];
    for index in (0..size.len().saturating_sub(1)).rev() {
        strides[index] = strides[index + 1] * size[index + 1];
    }
    strides
}

fn tensor_pickle(tensor: &Tensor, storage_index: usize, output: &mut Vec<u8>) -> Result<Vec<u8>> {
    let tensor = tensor.to_device(Device::Cpu).contiguous();
    let size = tensor.size();
    let count = size.iter().product::<i64>();
    let flat = tensor.reshape([count]);
    let (storage, bytes) = match flat.kind() {
        Kind::Float => (
            "FloatStorage",
            Vec::<f32>::try_from(&flat)?
                .into_iter()
                .flat_map(f32::to_le_bytes)
                .collect(),
        ),
        Kind::Double => (
            "DoubleStorage",
            Vec::<f64>::try_from(&flat)?
                .into_iter()
                .flat_map(f64::to_le_bytes)
                .collect(),
        ),
        Kind::Int64 => (
            "LongStorage",
            Vec::<i64>::try_from(&flat)?
                .into_iter()
                .flat_map(i64::to_le_bytes)
                .collect(),
        ),
        kind => bail!("unsupported checkpoint tensor type {kind:?}"),
    };
    global(output, "torch._utils", "_rebuild_tensor_v2");
    output.push(b'(');
    output.push(b'(');
    unicode(output, "storage");
    global(output, "torch", storage);
    unicode(output, &storage_index.to_string());
    unicode(output, "cpu");
    integer(output, count);
    output.push(b't');
    output.push(b'Q');
    integer(output, tensor.storage_offset());
    dimensions(output, &size);
    dimensions(output, &strides(&size));
    output.push(0x89);
    ordered_dict(output);
    output.push(b't');
    output.push(b'R');
    Ok(bytes)
}

fn checkpoint_pickle(tensors: &[(String, Tensor)]) -> Result<(Vec<u8>, Vec<Vec<u8>>)> {
    let mut pickle = vec![0x80, 0x02];
    ordered_dict(&mut pickle);
    pickle.push(b'(');
    let mut storages = Vec::new();
    for (index, (name, tensor)) in tensors.iter().enumerate() {
        unicode(&mut pickle, name);
        storages.push(tensor_pickle(tensor, index, &mut pickle)?);
    }
    pickle.extend(b"u.");
    Ok((pickle, storages))
}

pub fn save_checkpoint(var_store: &VarStore, path: &Path) -> Result<()> {
    let mut tensors: Vec<(String, Tensor)> = var_store.variables().into_iter().collect();
    tensors.sort_by(|left, right| left.0.cmp(&right.0));
    let (pickle, storages) = checkpoint_pickle(&tensors)?;
    let prefix = path
        .file_stem()
        .and_then(|name| name.to_str())
        .context("checkpoint filename must be valid UTF-8")?;
    let file = File::create(path)
        .with_context(|| format!("failed to create checkpoint {}", path.display()))?;
    let mut archive = ZipWriter::new(file);
    let options = FileOptions::default().compression_method(CompressionMethod::Stored);
    archive.start_file(format!("{prefix}/data.pkl"), options)?;
    archive.write_all(&pickle)?;
    archive.start_file(format!("{prefix}/byteorder"), options)?;
    archive.write_all(b"little")?;
    archive.start_file(format!("{prefix}/.format_version"), options)?;
    archive.write_all(b"1")?;
    archive.start_file(format!("{prefix}/version"), options)?;
    archive.write_all(b"1")?;
    for (index, data) in storages.iter().enumerate() {
        archive.start_file(format!("{prefix}/data/{index}"), options)?;
        archive.write_all(data)?;
    }
    archive.finish()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn checkpoint_is_a_zip_with_python_pickle_and_little_endian_storage() {
        let tensor = Tensor::from_slice(&[1.0_f32, 2.0, 3.0, 4.0]).reshape([2, 2]);
        let (pickle, storage) = checkpoint_pickle(&[("weight".to_string(), tensor)]).unwrap();
        assert_eq!(&pickle[..2], &[0x80, 0x02]);
        assert!(
            pickle
                .windows(b"weight".len())
                .any(|window| window == b"weight")
        );
        assert!(
            pickle
                .windows(b"FloatStorage".len())
                .any(|window| window == b"FloatStorage")
        );
        assert_eq!(storage[0].len(), 16);
        assert_eq!(&storage[0][..4], &1.0_f32.to_le_bytes());
    }
}
