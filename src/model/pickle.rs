use std::{
    collections::HashMap,
    fs::File,
    io::{Read, Seek},
    path::Path,
};

use anyhow::{Context, Result, bail};
use half::{bf16, f16};
use tch::{Kind, Tensor};
use zip::ZipArchive;

#[derive(Debug, Clone)]
struct Global {
    module: String,
    name: String,
}

#[derive(Debug, Clone)]
struct StorageSpec {
    dtype: String,
    key: String,
    numel: usize,
}

#[derive(Debug, Clone)]
struct TensorSpec {
    storage: StorageSpec,
    offset: i64,
    shape: Vec<i64>,
    stride: Vec<i64>,
}

#[derive(Debug, Clone)]
enum Value {
    None,
    Bool(bool),
    Int(i64),
    Float(f64),
    String(String),
    Global(Global),
    Tuple(Vec<Value>),
    List(Vec<Value>),
    Dict(Vec<(Value, Value)>),
    Tensor(TensorSpec),
    Mark,
}

impl Value {
    fn string(self) -> Result<String> {
        match self {
            Self::String(value) => Ok(value),
            other => bail!("expected string, got {other:?}"),
        }
    }

    fn int(self) -> Result<i64> {
        match self {
            Self::Int(value) => Ok(value),
            other => bail!("expected integer, got {other:?}"),
        }
    }

    fn bool(self) -> Result<bool> {
        match self {
            Self::Bool(value) => Ok(value),
            Self::Int(value) => Ok(value != 0),
            other => bail!("expected boolean, got {other:?}"),
        }
    }

    fn float(self) -> Result<f64> {
        match self {
            Self::Float(value) => Ok(value),
            Self::Int(value) => Ok(value as f64),
            Self::Bool(value) => Ok(if value { 1.0 } else { 0.0 }),
            other => bail!("expected float, got {other:?}"),
        }
    }

    fn usize(self) -> Result<usize> {
        let value = self.int()?;
        if value < 0 {
            bail!("expected non-negative integer, got {value}")
        }
        Ok(value as usize)
    }

    fn sequence(self) -> Result<Vec<Value>> {
        match self {
            Self::Tuple(values) | Self::List(values) => Ok(values),
            other => bail!("expected sequence, got {other:?}"),
        }
    }

    fn ints(self) -> Result<Vec<i64>> {
        self.sequence()?.into_iter().map(Self::int).collect()
    }
}

struct Parser<'a> {
    data: &'a [u8],
    position: usize,
    stack: Vec<Value>,
    memo: HashMap<u32, Value>,
    next_memo: u32,
}

impl<'a> Parser<'a> {
    fn new(data: &'a [u8]) -> Self {
        Self {
            data,
            position: 0,
            stack: vec![],
            memo: HashMap::new(),
            next_memo: 0,
        }
    }

    fn byte(&mut self) -> Result<u8> {
        let value = self
            .data
            .get(self.position)
            .copied()
            .context("unexpected end of pickle")?;
        self.position += 1;
        Ok(value)
    }

    fn bytes(&mut self, length: usize) -> Result<&'a [u8]> {
        let end = self
            .position
            .checked_add(length)
            .context("pickle offset overflow")?;
        let value = self
            .data
            .get(self.position..end)
            .context("truncated pickle value")?;
        self.position = end;
        Ok(value)
    }

    fn line(&mut self) -> Result<String> {
        let remainder = &self.data[self.position..];
        let length = remainder
            .iter()
            .position(|byte| *byte == b'\n')
            .context("unterminated pickle line")?;
        let line = String::from_utf8(self.bytes(length)?.to_vec())?;
        self.byte()?;
        Ok(line)
    }

    fn u16(&mut self) -> Result<u16> {
        Ok(u16::from_le_bytes(self.bytes(2)?.try_into().unwrap()))
    }

    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.bytes(4)?.try_into().unwrap()))
    }

    fn u64(&mut self) -> Result<u64> {
        Ok(u64::from_le_bytes(self.bytes(8)?.try_into().unwrap()))
    }

    fn string(&mut self, length: usize) -> Result<Value> {
        Ok(Value::String(String::from_utf8(
            self.bytes(length)?.to_vec(),
        )?))
    }

    fn long(&mut self, length: usize) -> Result<Value> {
        if length == 0 || length > 8 {
            bail!("unsupported pickle LONG size {length}")
        }
        let bytes = self.bytes(length)?;
        let mut value = if bytes[length - 1] & 0x80 == 0 {
            0_i64
        } else {
            -1_i64
        };
        for byte in bytes {
            value = (value << 8) | i64::from(*byte);
        }
        Ok(Value::Int(value))
    }

    fn push(&mut self, value: Value) {
        self.stack.push(value);
    }

    fn pop(&mut self) -> Result<Value> {
        self.stack.pop().context("pickle stack underflow")
    }

    fn memo_put(&mut self, index: u32) -> Result<()> {
        let value = self
            .stack
            .last()
            .context("memo write with empty stack")?
            .clone();
        self.memo.insert(index, value);
        self.next_memo = self.next_memo.max(index.saturating_add(1));
        Ok(())
    }

    fn memo_get(&mut self, index: u32) -> Result<()> {
        let value = self
            .memo
            .get(&index)
            .cloned()
            .with_context(|| format!("missing pickle memo entry {index}"))?;
        self.push(value);
        Ok(())
    }

    fn pop_mark(&mut self) -> Result<Vec<Value>> {
        let index = self
            .stack
            .iter()
            .rposition(|value| matches!(value, Value::Mark))
            .context("pickle MARK not found")?;
        let values = self.stack.split_off(index + 1);
        self.stack.pop();
        Ok(values)
    }

    fn reduce(&mut self) -> Result<()> {
        let args = self.pop()?.sequence()?;
        let callable = self.pop()?;
        let Value::Global(global) = callable else {
            bail!("unsupported pickle callable {callable:?}")
        };
        let value = match (global.module.as_str(), global.name.as_str()) {
            ("collections", "OrderedDict") | ("builtins", "dict") => {
                let mut entries = vec![];
                if let Some(first) = args.into_iter().next() {
                    for pair in first.sequence()? {
                        let mut pair = pair.sequence()?;
                        if pair.len() != 2 {
                            bail!("dictionary entry is not a pair")
                        }
                        let value = pair.pop().unwrap();
                        let key = pair.pop().unwrap();
                        entries.push((key, value));
                    }
                }
                Value::Dict(entries)
            }
            ("torch._utils", "_rebuild_tensor_v2") | ("torch._utils", "_rebuild_tensor") => {
                if args.len() < 4 {
                    bail!("tensor rebuild requires at least four arguments")
                }
                let storage = parse_storage(args[0].clone())?;
                let offset = args[1].clone().int()?;
                let shape = args[2].clone().ints()?;
                let stride = args[3].clone().ints()?;
                Value::Tensor(TensorSpec {
                    storage,
                    offset,
                    shape,
                    stride,
                })
            }
            ("torch._utils", "_rebuild_parameter") => args
                .into_iter()
                .next()
                .context("missing parameter tensor")?,
            ("builtins", "bool") => Value::Bool(
                args.into_iter()
                    .next()
                    .context("missing bool argument")?
                    .bool()?,
            ),
            ("builtins", "float") => Value::Float(
                args.into_iter()
                    .next()
                    .context("missing float argument")?
                    .float()?,
            ),
            ("torch", "device") => Value::String(
                args.into_iter()
                    .next()
                    .context("missing torch device")?
                    .string()?,
            ),
            _ => bail!(
                "unsupported pickle global {}.{}",
                global.module,
                global.name
            ),
        };
        self.push(value);
        Ok(())
    }

    fn parse(mut self) -> Result<Value> {
        loop {
            match self.byte()? {
                0x80 => {
                    self.byte()?;
                }
                0x95 => {
                    self.u64()?;
                }
                0x2e => {
                    if self.position != self.data.len() {
                        bail!("trailing data after pickle STOP")
                    }
                    return self.pop();
                }
                0x28 => self.push(Value::Mark),
                0x29 => self.push(Value::Tuple(vec![])),
                0x74 => {
                    let values = self.pop_mark()?;
                    self.push(Value::Tuple(values));
                }
                0x85..=0x87 => {
                    let count = match self.data[self.position - 1] {
                        0x85 => 1,
                        0x86 => 2,
                        _ => 3,
                    };
                    let mut values = vec![];
                    for _ in 0..count {
                        values.push(self.pop()?);
                    }
                    values.reverse();
                    self.push(Value::Tuple(values));
                }
                0x5d => self.push(Value::List(vec![])),
                0x6c => {
                    let values = self.pop_mark()?;
                    self.push(Value::List(values));
                }
                0x61 => {
                    let value = self.pop()?;
                    match self.stack.last_mut() {
                        Some(Value::List(values)) => values.push(value),
                        other => bail!("APPEND outside list: {other:?}"),
                    }
                }
                0x65 => {
                    let values = self.pop_mark()?;
                    match self.stack.last_mut() {
                        Some(Value::List(target)) => target.extend(values),
                        other => bail!("APPENDS outside list: {other:?}"),
                    }
                }
                0x7d => self.push(Value::Dict(vec![])),
                0x64 => {
                    let values = self.pop_mark()?;
                    self.push(Value::Dict(pairs(values)?));
                }
                0x73 => {
                    let value = self.pop()?;
                    let key = self.pop()?;
                    match self.stack.last_mut() {
                        Some(Value::Dict(entries)) => entries.push((key, value)),
                        other => bail!("SETITEM outside dictionary: {other:?}"),
                    }
                }
                0x75 => {
                    let values = self.pop_mark()?;
                    let entries = pairs(values)?;
                    match self.stack.last_mut() {
                        Some(Value::Dict(target)) => target.extend(entries),
                        other => bail!("SETITEMS outside dictionary: {other:?}"),
                    }
                }
                0x4e => self.push(Value::None),
                0x88 => self.push(Value::Bool(true)),
                0x89 => self.push(Value::Bool(false)),
                0x4a => {
                    let value = i32::from_le_bytes(self.bytes(4)?.try_into().unwrap());
                    self.push(Value::Int(value as i64));
                }
                0x4b => {
                    let value = self.byte()?;
                    self.push(Value::Int(value as i64));
                }
                0x4d => {
                    let value = self.u16()?;
                    self.push(Value::Int(value as i64));
                }
                0x8a => {
                    let length = self.byte()? as usize;
                    let value = self.long(length)?;
                    self.push(value);
                }
                0x8b => {
                    let length = self.u32()? as usize;
                    let value = self.long(length)?;
                    self.push(value);
                }
                0x49 => {
                    let value = self.line()?.parse()?;
                    self.push(Value::Int(value));
                }
                0x4c => {
                    let value = self.line()?.trim_end_matches('L').parse()?;
                    self.push(Value::Int(value));
                }
                0x46 => {
                    let value = self.line()?.parse()?;
                    self.push(Value::Float(value));
                }
                0x47 => {
                    let value = f64::from_be_bytes(self.bytes(8)?.try_into().unwrap());
                    self.push(Value::Float(value));
                }
                0x58 => {
                    let length = self.u32()? as usize;
                    let value = self.string(length)?;
                    self.push(value);
                }
                0x8c => {
                    let length = self.byte()? as usize;
                    let value = self.string(length)?;
                    self.push(value);
                }
                0x8d => {
                    let length = self.u64()? as usize;
                    let value = self.string(length)?;
                    self.push(value);
                }
                0x54 => {
                    let length = self.u32()? as usize;
                    let value = self.string(length)?;
                    self.push(value);
                }
                0x55 => {
                    let length = self.byte()? as usize;
                    let value = self.string(length)?;
                    self.push(value);
                }
                0x56 => {
                    let value = self.line()?;
                    self.push(Value::String(value));
                }
                0x63 => {
                    let module = self.line()?;
                    let name = self.line()?;
                    self.push(Value::Global(Global { module, name }));
                }
                0x93 => {
                    let name = self.pop()?.string()?;
                    let module = self.pop()?.string()?;
                    self.push(Value::Global(Global { module, name }));
                }
                0x52 => self.reduce()?,
                0x62 => {
                    self.pop()?;
                }
                0x71 => {
                    let index = self.byte()? as u32;
                    self.memo_put(index)?;
                }
                0x72 => {
                    let index = self.u32()?;
                    self.memo_put(index)?;
                }
                0x68 => {
                    let index = self.byte()? as u32;
                    self.memo_get(index)?;
                }
                0x6a => {
                    let index = self.u32()?;
                    self.memo_get(index)?;
                }
                0x94 => {
                    let index = self.next_memo;
                    self.memo_put(index)?;
                }
                0x51 => {
                    let index = self.next_memo;
                    self.memo_put(index)?;
                }
                0x30 => {
                    self.pop()?;
                }
                0x31 => {
                    self.pop_mark()?;
                }
                0x32 => {
                    let value = self.stack.last().context("DUP with empty stack")?.clone();
                    self.push(value);
                }
                opcode => bail!(
                    "unsupported pickle opcode 0x{opcode:02x} at offset {}",
                    self.position - 1
                ),
            }
        }
    }
}

fn pairs(values: Vec<Value>) -> Result<Vec<(Value, Value)>> {
    if !values.len().is_multiple_of(2) {
        bail!("dictionary contains an unpaired item")
    }
    let mut entries = vec![];
    let mut values = values.into_iter();
    while let Some(key) = values.next() {
        entries.push((key, values.next().unwrap()));
    }
    Ok(entries)
}

fn parse_storage(value: Value) -> Result<StorageSpec> {
    let values = value.sequence()?;
    if values.len() < 5 || values[0].clone().string()? != "storage" {
        bail!("invalid tensor storage descriptor")
    }
    let Value::Global(global) = &values[1] else {
        bail!("invalid tensor storage dtype")
    };
    Ok(StorageSpec {
        dtype: global.name.clone(),
        key: values[2].clone().string()?,
        numel: values[4].clone().usize()?,
    })
}

fn read_record<R: Read + Seek>(archive: &mut ZipArchive<R>, name: &str) -> Result<Vec<u8>> {
    let mut record = archive
        .by_name(name)
        .with_context(|| format!("missing archive record {name}"))?;
    let mut data = Vec::with_capacity(record.size() as usize);
    record.read_to_end(&mut data)?;
    Ok(data)
}

fn values<const N: usize>(
    bytes: &[u8],
    expected: usize,
    convert: impl Fn([u8; N]) -> f32,
) -> Result<Vec<f32>> {
    if bytes.len() != expected * N {
        bail!(
            "storage has {} bytes, expected {}",
            bytes.len(),
            expected * N
        )
    }
    Ok(bytes
        .chunks_exact(N)
        .map(|chunk| convert(chunk.try_into().unwrap()))
        .collect())
}

fn tensor_from_spec(
    archive: &mut ZipArchive<impl Read + Seek>,
    prefix: &str,
    spec: &TensorSpec,
) -> Result<Tensor> {
    let record = format!("{prefix}data/{}", spec.storage.key);
    let bytes = read_record(archive, &record)?;
    let count = spec.storage.numel;
    let tensor = match spec.storage.dtype.as_str() {
        "FloatStorage" => Tensor::from_slice(&values::<4>(&bytes, count, f32::from_le_bytes)?),
        "DoubleStorage" => {
            if bytes.len() != count * 8 {
                bail!("invalid float64 storage size")
            }
            let doubles: Vec<f64> = bytes
                .chunks_exact(8)
                .map(|chunk| f64::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            Tensor::from_slice(&doubles)
        }
        "LongStorage" => {
            if bytes.len() != count * 8 {
                bail!("invalid int64 storage size")
            }
            let values: Vec<i64> = bytes
                .chunks_exact(8)
                .map(|chunk| i64::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            Tensor::from_slice(&values)
        }
        "IntStorage" => {
            if bytes.len() != count * 4 {
                bail!("invalid int32 storage size")
            }
            let values: Vec<i32> = bytes
                .chunks_exact(4)
                .map(|chunk| i32::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            Tensor::from_slice(&values)
        }
        "ShortStorage" => {
            if bytes.len() != count * 2 {
                bail!("invalid int16 storage size")
            }
            let values: Vec<i16> = bytes
                .chunks_exact(2)
                .map(|chunk| i16::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            Tensor::from_slice(&values)
        }
        "CharStorage" => {
            if bytes.len() != count {
                bail!("invalid int8 storage size")
            }
            let values: Vec<i8> = bytes.iter().map(|value| *value as i8).collect();
            Tensor::from_slice(&values)
        }
        "ByteStorage" => {
            if bytes.len() != count {
                bail!("invalid uint8 storage size")
            }
            Tensor::from_slice(&bytes)
        }
        "BoolStorage" => {
            if bytes.len() != count {
                bail!("invalid bool storage size")
            }
            let values: Vec<bool> = bytes.iter().map(|value| *value != 0).collect();
            Tensor::from_slice(&values)
        }
        "HalfStorage" => {
            if bytes.len() != count * 2 {
                bail!("invalid float16 storage size")
            }
            let values: Vec<f16> = bytes
                .chunks_exact(2)
                .map(|chunk| f16::from_le_bytes(chunk.try_into().unwrap()))
                .collect();
            Tensor::from_slice(&values)
        }
        "BFloat16Storage" => {
            if bytes.len() != count * 2 {
                bail!("invalid bfloat16 storage size")
            }
            let values: Vec<f32> = bytes
                .chunks_exact(2)
                .map(|chunk| bf16::from_le_bytes(chunk.try_into().unwrap()).to_f32())
                .collect();
            Tensor::from_slice(&values).to_kind(Kind::BFloat16)
        }
        dtype => bail!("unsupported tensor storage type {dtype}"),
    };
    if spec.shape.len() != spec.stride.len() {
        bail!("tensor shape and stride ranks differ")
    }
    Ok(tensor.as_strided(&spec.shape, &spec.stride, spec.offset))
}

pub fn load_python_state_dict(path: &Path) -> Result<Vec<(String, Tensor)>> {
    let file = File::open(path)?;
    let mut archive = ZipArchive::new(std::io::BufReader::new(file))?;
    let pickle_name = archive
        .file_names()
        .find(|name| name.ends_with("data.pkl"))
        .map(str::to_string)
        .context("checkpoint does not contain data.pkl")?;
    let pickle = read_record(&mut archive, &pickle_name)?;
    let root = Parser::new(&pickle).parse()?;
    let mut entries = match root {
        Value::Dict(entries) => entries,
        other => bail!("checkpoint root is not a state dictionary: {other:?}"),
    };
    if let Some(position) = entries
        .iter()
        .position(|(key, _)| matches!(key, Value::String(name) if name == "state_dict"))
    {
        let (_, value) = entries.remove(position);
        entries = match value {
            Value::Dict(entries) => entries,
            other => bail!("state_dict entry is not a dictionary: {other:?}"),
        };
    }

    let prefix = pickle_name.trim_end_matches("data.pkl").to_string();
    let mut tensors = vec![];
    for (key, value) in entries {
        let key = key.string()?;
        let Value::Tensor(spec) = value else {
            bail!("checkpoint entry {key} is not a tensor")
        };
        tensors.push((key, tensor_from_spec(&mut archive, &prefix, &spec)?));
    }
    Ok(tensors)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_ordered_dict_checkpoint() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/models/fb.pt");
        let tensors = load_python_state_dict(&path).unwrap();
        assert!(
            tensors
                .iter()
                .any(|(name, _)| name == "attn_mha.in_proj_weight")
        );
        assert!(
            tensors
                .iter()
                .any(|(name, _)| name == "conv1.bn.num_batches_tracked")
        );
    }
}
