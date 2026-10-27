use std::{
    collections::BTreeMap,
    fs,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};

#[derive(Clone, Copy, Debug)]
struct BlockHandle {
    offset: u64,
    size: u64,
}

impl BlockHandle {
    fn decode(data: &[u8]) -> Result<Self> {
        let mut cursor = 0;
        let offset = varint(data, &mut cursor)?;
        let size = varint(data, &mut cursor)?;
        Ok(Self { offset, size })
    }
}

fn varint(data: &[u8], cursor: &mut usize) -> Result<u64> {
    let mut value = 0_u64;
    for shift in (0..70).step_by(7) {
        let byte = *data.get(*cursor).context("truncated varint")?;
        *cursor += 1;
        if shift >= 64 && byte > 1 {
            bail!("varint overflow");
        }
        value |= u64::from(byte & 0x7f) << shift;
        if byte & 0x80 == 0 {
            return Ok(value);
        }
    }
    bail!("varint is too long")
}

fn block(file: &[u8], handle: BlockHandle) -> Result<Vec<u8>> {
    let start = usize::try_from(handle.offset)?;
    let size = usize::try_from(handle.size)?;
    let end = start
        .checked_add(size + 5)
        .context("SSTable block overflow")?;
    if end > file.len() {
        bail!("SSTable block extends past the end of the file");
    }
    let compressed = &file[start..start + size];
    let compression = file[start + size];
    let stored_crc = u32::from_le_bytes(file[start + size + 1..end].try_into()?);
    let mut crc_input = compressed.to_vec();
    crc_input.push(compression);
    let actual_crc = (!crc32c(&crc_input))
        .rotate_right(15)
        .wrapping_add(0xa282_ead8);
    if stored_crc != actual_crc {
        bail!("SSTable block checksum does not match");
    }
    match compression {
        0 => Ok(compressed.to_vec()),
        1 => snappy(compressed),
        other => bail!("unsupported SSTable compression {other}"),
    }
}

fn crc32c(data: &[u8]) -> u32 {
    let mut crc = !0_u32;
    for byte in data {
        crc ^= u32::from(*byte);
        for _ in 0..8 {
            crc = (crc >> 1) ^ (0x82f6_3b78 & (crc & 1).wrapping_neg());
        }
    }
    crc
}

fn snappy(data: &[u8]) -> Result<Vec<u8>> {
    let mut cursor = 0;
    let expected = varint(data, &mut cursor)? as usize;
    let mut output = Vec::with_capacity(expected);
    while cursor < data.len() {
        let tag = data[cursor];
        cursor += 1;
        match tag & 3 {
            0 => {
                let mut length = usize::from(tag >> 2);
                if length < 60 {
                    length += 1;
                } else {
                    let bytes = length - 59;
                    let mut value = 0_usize;
                    for index in 0..bytes {
                        value |= usize::from(data[cursor + index]) << (index * 8);
                    }
                    cursor += bytes;
                    length = value + 1;
                }
                let end = cursor + length;
                output
                    .extend_from_slice(data.get(cursor..end).context("truncated Snappy literal")?);
                cursor = end;
            }
            1 => {
                let length = 4 + usize::from((tag >> 2) & 7);
                let offset = (usize::from(tag & 0xe0) << 3) | usize::from(data[cursor]);
                cursor += 1;
                copy(&mut output, offset, length)?;
            }
            2 => {
                let length = usize::from(tag >> 2) + 1;
                let offset = usize::from(u16::from_le_bytes(data[cursor..cursor + 2].try_into()?));
                cursor += 2;
                copy(&mut output, offset, length)?;
            }
            3 => {
                let length = usize::from(tag >> 2) + 1;
                let offset = u32::from_le_bytes(data[cursor..cursor + 4].try_into()?) as usize;
                cursor += 4;
                copy(&mut output, offset, length)?;
            }
            _ => unreachable!(),
        }
    }
    if output.len() != expected {
        bail!("Snappy output length differs from its header");
    }
    Ok(output)
}

fn copy(output: &mut Vec<u8>, offset: usize, length: usize) -> Result<()> {
    if offset == 0 || offset > output.len() {
        bail!("invalid Snappy back-reference");
    }
    let start = output.len() - offset;
    for index in 0..length {
        let byte = output[start + index];
        output.push(byte);
    }
    Ok(())
}

fn entries(data: &[u8]) -> Result<BTreeMap<Vec<u8>, Vec<u8>>> {
    let restarts = data
        .get(data.len().saturating_sub(4)..)
        .context("SSTable block is missing its restart array")?;
    let restart_count = u32::from_le_bytes(restarts.try_into()?) as usize;
    let restart_bytes = restart_count
        .checked_add(1)
        .and_then(|count| count.checked_mul(4))
        .context("SSTable restart array size overflow")?;
    let entries_size = data
        .len()
        .checked_sub(restart_bytes)
        .context("SSTable restart array extends past the beginning of its block")?;
    let data = &data[..entries_size];
    let mut entries = BTreeMap::new();
    let mut cursor = 0;
    let mut key = Vec::new();
    while cursor < data.len() {
        let shared = varint(data, &mut cursor)? as usize;
        let unshared = varint(data, &mut cursor)? as usize;
        let value_size = varint(data, &mut cursor)? as usize;
        let delta = data
            .get(cursor..cursor + unshared)
            .context("truncated SSTable key")?;
        cursor += unshared;
        let value = data
            .get(cursor..cursor + value_size)
            .context("truncated SSTable value")?;
        cursor += value_size;
        if shared > key.len() {
            bail!("SSTable shared key length exceeds the previous key length");
        }
        key.truncate(shared);
        key.extend_from_slice(delta);
        entries.insert(key.clone(), value.to_vec());
    }
    Ok(entries)
}

fn read_index(file: &[u8]) -> Result<BTreeMap<Vec<u8>, BlockHandle>> {
    if file.len() < 48 {
        bail!("TensorFlow checkpoint index is too short");
    }
    if u64::from_le_bytes(file[file.len() - 8..].try_into()?) != 0xdb47_7524_8b80_fb57 {
        bail!("not an SSTable checkpoint index");
    }
    let trailer = &file[file.len() - 48..];
    let mut cursor = 0;
    let _metadata_offset = varint(trailer, &mut cursor)?;
    let _metadata_size = varint(trailer, &mut cursor)?;
    let index_offset = varint(trailer, &mut cursor)?;
    let index_size = varint(trailer, &mut cursor)?;
    let index = block(
        file,
        BlockHandle {
            offset: index_offset,
            size: index_size,
        },
    )?;
    let mut handles = BTreeMap::new();
    for (key, value) in entries(&index).context("failed to decode checkpoint index block")? {
        handles.insert(key, BlockHandle::decode(&value)?);
    }
    Ok(handles)
}

fn proto<'a>(data: &'a [u8], cursor: &mut usize) -> Result<Option<(u32, ProtoValue<'a>)>> {
    if *cursor >= data.len() {
        return Ok(None);
    }
    let tag = varint(data, cursor)?;
    let field = (tag >> 3) as u32;
    match tag & 7 {
        0 => Ok(Some((field, ProtoValue::Integer(varint(data, cursor)?)))),
        1 => {
            let value = u64::from_le_bytes(
                data.get(*cursor..*cursor + 8)
                    .context("truncated fixed64")?
                    .try_into()?,
            );
            *cursor += 8;
            Ok(Some((field, ProtoValue::Integer(value))))
        }
        2 => {
            let length = varint(data, cursor)? as usize;
            let value = data
                .get(*cursor..*cursor + length)
                .context("truncated protobuf bytes")?;
            *cursor += length;
            Ok(Some((field, ProtoValue::Bytes(value))))
        }
        5 => {
            let value = u32::from_le_bytes(
                data.get(*cursor..*cursor + 4)
                    .context("truncated fixed32")?
                    .try_into()?,
            );
            *cursor += 4;
            Ok(Some((field, ProtoValue::Integer(u64::from(value)))))
        }
        wire => bail!("unsupported protobuf wire type {wire}"),
    }
}

enum ProtoValue<'a> {
    Integer(u64),
    Bytes(&'a [u8]),
}

fn header(data: &[u8]) -> Result<u64> {
    let mut cursor = 0;
    while let Some((field, value)) = proto(data, &mut cursor)? {
        if field == 1 {
            let ProtoValue::Integer(shards) = value else {
                bail!("invalid shard count")
            };
            return Ok(shards);
        }
    }
    bail!("TensorFlow bundle has no shard count")
}

struct BundleEntry {
    shape: Vec<i64>,
    shard: u64,
    offset: u64,
}

fn bundle_entry(data: &[u8]) -> Result<BundleEntry> {
    let mut cursor = 0;
    let mut shape = Vec::new();
    let mut shard = 0;
    let mut offset = 0;
    let mut float = false;
    while let Some((field, value)) = proto(data, &mut cursor)? {
        match (field, value) {
            (1, ProtoValue::Integer(dtype)) => float = dtype == 1,
            (2, ProtoValue::Bytes(shape_data)) => {
                let mut shape_cursor = 0;
                while let Some((field, value)) = proto(shape_data, &mut shape_cursor)? {
                    if field == 2 {
                        let ProtoValue::Bytes(dimension) = value else {
                            bail!("invalid tensor dimension")
                        };
                        let mut dimension_cursor = 0;
                        while let Some((field, value)) = proto(dimension, &mut dimension_cursor)? {
                            if field == 1 {
                                let ProtoValue::Integer(size) = value else {
                                    bail!("invalid tensor dimension size")
                                };
                                shape.push(size as i64);
                            }
                        }
                    }
                }
            }
            (3, ProtoValue::Integer(value)) => shard = value,
            (4, ProtoValue::Integer(value)) => offset = value,
            _ => {}
        }
    }
    if !float {
        bail!("only float32 TensorFlow checkpoint tensors are supported");
    }
    Ok(BundleEntry {
        shape,
        shard,
        offset,
    })
}

fn checkpoint_file(prefix: &Path, suffix: &str) -> PathBuf {
    let mut path = prefix.as_os_str().to_os_string();
    path.push(suffix);
    PathBuf::from(path)
}

pub struct Checkpoint {
    prefix: PathBuf,
    shards: u64,
    blocks: BTreeMap<Vec<u8>, BlockHandle>,
}

impl Checkpoint {
    pub fn open(prefix: &Path) -> Result<Self> {
        let index = checkpoint_file(prefix, ".index");
        let bytes = fs::read(&index).with_context(|| {
            format!(
                "failed to read TensorFlow checkpoint index {}",
                index.display()
            )
        })?;
        let handles = read_index(&bytes)?;
        let first = handles
            .values()
            .next()
            .context("TensorFlow checkpoint index has no data blocks")?;
        let first_block = block(&bytes, *first)?;
        let metadata = entries(&first_block).context("failed to decode checkpoint header block")?;
        let shards = metadata
            .get(&Vec::new())
            .map(|data| header(data).context("failed to decode TensorFlow bundle header"))
            .transpose()?
            .unwrap_or(1);
        Ok(Self {
            prefix: prefix.to_path_buf(),
            shards,
            blocks: handles,
        })
    }

    pub fn read_tensor(&self, name: &str) -> Result<(Vec<i64>, Vec<f32>)> {
        let handle = self
            .blocks
            .range::<Vec<u8>, _>((
                std::ops::Bound::Included(name.as_bytes().to_vec()),
                std::ops::Bound::Unbounded,
            ))
            .next()
            .map(|(_, handle)| *handle)
            .context("tensor is absent from checkpoint index")?;
        let data = fs::read(checkpoint_file(&self.prefix, ".index"))?;
        let block = block(&data, handle)?;
        let values = entries(&block)
            .with_context(|| format!("failed to decode checkpoint block containing {name}"))?;
        let entry = values
            .get(name.as_bytes())
            .context("tensor is absent from its checkpoint block")?;
        let entry = bundle_entry(entry)
            .with_context(|| format!("failed to decode TensorFlow bundle entry for {name}"))?;
        if entry.shard >= self.shards.max(1) {
            bail!("tensor shard is out of range");
        }
        let path = PathBuf::from(format!(
            "{}.data-{:05}-of-{:05}",
            self.prefix.display(),
            entry.shard,
            self.shards.max(1)
        ));
        let bytes = fs::read(&path)
            .with_context(|| format!("failed to read tensor data {}", path.display()))?;
        let start = usize::try_from(entry.offset)?;
        let count = entry.shape.iter().product::<i64>() as usize;
        let byte_count = count.checked_mul(4).context("tensor size overflow")?;
        let end = start
            .checked_add(byte_count)
            .context("tensor offset overflow")?;
        let bytes = bytes
            .get(start..end)
            .context("tensor data is shorter than expected")?;
        let values = bytes
            .chunks_exact(4)
            .map(|chunk| f32::from_le_bytes(chunk.try_into().unwrap()))
            .collect();
        Ok((entry.shape, values))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn put_varint(mut value: u64, output: &mut Vec<u8>) {
        while value >= 128 {
            output.push((value as u8 & 127) | 128);
            value >>= 7;
        }
        output.push(value as u8);
    }

    fn entry(key: &[u8], value: &[u8]) -> Vec<u8> {
        let mut output = Vec::new();
        put_varint(0, &mut output);
        put_varint(key.len() as u64, &mut output);
        put_varint(value.len() as u64, &mut output);
        output.extend_from_slice(key);
        output.extend_from_slice(value);
        output
    }

    fn handle(handle: BlockHandle) -> Vec<u8> {
        let mut output = Vec::new();
        put_varint(handle.offset, &mut output);
        put_varint(handle.size, &mut output);
        output
    }

    fn finish_block(entries: &mut Vec<u8>) {
        entries.extend_from_slice(&0_u32.to_le_bytes());
        entries.extend_from_slice(&1_u32.to_le_bytes());
    }

    fn put_block(file: &mut Vec<u8>, contents: &[u8]) -> BlockHandle {
        let handle = BlockHandle {
            offset: file.len() as u64,
            size: contents.len() as u64,
        };
        file.extend_from_slice(contents);
        file.push(0);
        let mut crc_input = contents.to_vec();
        crc_input.push(0);
        let crc = (!crc32c(&crc_input))
            .rotate_right(15)
            .wrapping_add(0xa282_ead8);
        file.extend_from_slice(&crc.to_le_bytes());
        handle
    }

    #[test]
    fn reads_sstable_entries_and_tensor_data() {
        let mut data = entry(b"", &[8, 1]);
        data.extend(entry(b"tensor", &[8, 1, 18, 4, 18, 2, 8, 2, 24, 0, 32, 0]));
        finish_block(&mut data);
        let mut metadata = Vec::new();
        finish_block(&mut metadata);
        let mut file = Vec::new();
        let data_handle = put_block(&mut file, &data);
        let metadata_handle = put_block(&mut file, &metadata);
        let mut index = entry(b"tensor", &handle(data_handle));
        finish_block(&mut index);
        let index_handle = put_block(&mut file, &index);
        let mut footer = Vec::new();
        for value in [
            metadata_handle.offset,
            metadata_handle.size,
            index_handle.offset,
            index_handle.size,
        ] {
            put_varint(value, &mut footer);
        }
        footer.resize(40, 0);
        footer.extend_from_slice(&0xdb47_7524_8b80_fb57_u64.to_le_bytes());
        file.extend_from_slice(&footer);
        let directory = tempfile::tempdir().unwrap();
        let prefix = directory.path().join("model.v1");
        fs::write(checkpoint_file(&prefix, ".index"), file).unwrap();
        fs::write(
            checkpoint_file(&prefix, ".data-00000-of-00001"),
            [1.0_f32.to_le_bytes(), 2.0_f32.to_le_bytes()].concat(),
        )
        .unwrap();
        let checkpoint = Checkpoint::open(&prefix).unwrap();
        let (shape, values) = checkpoint.read_tensor("tensor").unwrap();
        assert_eq!(shape, vec![2]);
        assert_eq!(values, vec![1.0, 2.0]);
    }

    #[test]
    fn ignores_restart_data_and_decodes_overlapping_snappy_copies() {
        let mut data = entry(b"key", b"value");
        finish_block(&mut data);
        assert_eq!(
            entries(&data).unwrap().get(b"key".as_slice()).unwrap(),
            b"value"
        );
        let compressed = [12, 8, b'a', b'b', b'c', 0x15, 3];
        assert_eq!(snappy(&compressed).unwrap(), b"abcabcabcabc");
    }
}
