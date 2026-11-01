use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    io::{BufRead, BufReader},
    path::PathBuf,
};

use anyhow::{Context, Result, bail};
use clap::Args;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead, Record, record::Aux};

use crate::io::reference::Reference;

#[derive(Args, Debug)]
#[command(rename_all = "snake_case")]
pub struct PhaseReadsArgs {
    #[arg(long = "bam_fn")]
    bam: PathBuf,
    #[arg(long = "ref_fn")]
    reference: PathBuf,
    #[arg(long = "variants_fn")]
    variants: PathBuf,
    #[arg(long = "ctg_name")]
    contig: String,
    #[arg(long = "output_fn")]
    output: PathBuf,
    #[arg(long = "meta_fn")]
    metadata: PathBuf,
    #[arg(long = "region")]
    regions: Vec<String>,
}

#[derive(Debug)]
struct PhaseVariant {
    position: i64,
    reference: Vec<u8>,
    alternate: Vec<u8>,
}

#[derive(Debug)]
struct NativeRead {
    name: Vec<u8>,
    mapq: u8,
    source: i32,
    sample: i32,
    start: i64,
    end: i64,
    hp: i32,
    ps: i32,
    supplementary: bool,
    reverse: bool,
    bx: Vec<u8>,
    chromosome: Vec<u8>,
    subalignment: Vec<u8>,
    variants: Vec<(i32, i32, i32)>,
}

type Cigar = (u8, i64);

fn read_variants(path: &PathBuf) -> Result<Vec<PhaseVariant>> {
    let mut output = Vec::new();
    for line in BufReader::new(fs::File::open(path)?).lines() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() != 3 {
            bail!("invalid native phase variant line: {line}");
        }
        output.push(PhaseVariant {
            position: fields[0].parse()?,
            reference: fields[1].as_bytes().to_vec(),
            alternate: fields[2].as_bytes().to_vec(),
        });
    }
    Ok(output)
}

fn parse_region(value: &str) -> Result<(i64, i64, usize, usize)> {
    let fields: Vec<_> = value.split(',').collect();
    if fields.len() != 4 {
        bail!("invalid native phase region: {value}");
    }
    let start = fields[0].parse()?;
    let end = fields[1].parse()?;
    let first = fields[2].parse()?;
    let last = fields[3].parse()?;
    Ok((start, end, first, last))
}

fn split_left(cigar: &[Cigar], index: usize, consumed: i64) -> impl Iterator<Item = Cigar> + '_ {
    let (operation, _) = cigar[index];
    std::iter::once((operation, consumed))
        .filter(|(_, length)| *length > 0)
        .chain(cigar[..index].iter().rev().copied())
}

fn split_right(cigar: &[Cigar], index: usize, consumed: i64) -> impl Iterator<Item = Cigar> + '_ {
    let (operation, length) = cigar[index];
    std::iter::once((operation, length - consumed))
        .filter(|(_, length)| *length > 0)
        .chain(cigar[index + 1..].iter().copied())
}

fn prefix_length(cigar: impl Iterator<Item = Cigar>, target: i64) -> (i64, i64) {
    let mut reference = 0;
    let mut query = 0;
    for (operation, length) in cigar {
        match operation {
            0 | 7 | 8 => {
                reference += length;
                query += length;
                if reference >= target {
                    return (target, query + target - reference);
                }
            }
            2 => {
                reference += length;
                if reference >= target {
                    return (target, query);
                }
            }
            1 => query += length,
            4 | 5 => {}
            3 => return (target, query),
            _ => unreachable!(),
        }
    }
    assert!(reference < target);
    (reference, query)
}

fn edit_distance(first: &[u8], second: &[u8]) -> usize {
    if second.len() <= 21 {
        let mut previous = [0_usize; 22];
        let mut current = [0_usize; 22];
        for (index, value) in previous.iter_mut().take(second.len() + 1).enumerate() {
            *value = index;
        }
        for (i, left) in first.iter().enumerate() {
            current[0] = i + 1;
            for (j, right) in second.iter().enumerate() {
                current[j + 1] = (previous[j + 1] + 1)
                    .min(current[j] + 1)
                    .min(previous[j] + usize::from(left != right));
            }
            std::mem::swap(&mut previous, &mut current);
        }
        return previous[second.len()];
    }
    let mut previous: Vec<usize> = (0..=second.len()).collect();
    let mut current = vec![0; second.len() + 1];
    for (i, left) in first.iter().enumerate() {
        current[0] = i + 1;
        for (j, right) in second.iter().enumerate() {
            current[j + 1] = (previous[j + 1] + 1)
                .min(current[j] + 1)
                .min(previous[j] + usize::from(left != right));
        }
        std::mem::swap(&mut previous, &mut current);
    }
    previous[second.len()]
}

fn realign(
    variant: &PhaseVariant,
    cigar: &[Cigar],
    index: usize,
    consumed: i64,
    query_position: i64,
    query: &[u8],
    reference: &[u8],
) -> Result<Option<i32>> {
    let (left_reference, left_query) = prefix_length(split_left(cigar, index, consumed), 10);
    let (right_reference, right_query) = prefix_length(
        split_right(cigar, index, consumed),
        variant.reference.len() as i64 + 10,
    );
    let query_start = query_position - left_query;
    let query_end = query_position + right_query;
    if query_start < 0 || query_end > query.len() as i64 {
        bail!("read realignment exceeds query sequence");
    }
    let reference_start = variant.position - left_reference;
    let reference_end = variant.position + right_reference;
    if reference_start < 0 || reference_end > reference.len() as i64 {
        bail!("variant realignment exceeds reference sequence");
    }
    let query = &query[query_start as usize..query_end as usize];
    let reference_allele = &reference[reference_start as usize..reference_end as usize];
    let position = variant.position as usize;
    let left = reference_start as usize;
    let right = reference_end as usize;
    if query == reference_allele {
        return Ok(Some(0));
    }
    let reference_length = variant.reference.len();
    let mut alternate =
        Vec::with_capacity(left_query as usize + variant.alternate.len() + right_query as usize);
    alternate.extend_from_slice(&reference[left..position]);
    alternate.extend_from_slice(&variant.alternate);
    alternate.extend_from_slice(&reference[position + reference_length..right]);
    if query == alternate {
        return Ok(Some(1));
    }
    let reference_distance = edit_distance(query, reference_allele);
    let alternate_distance = edit_distance(query, &alternate);
    Ok(if reference_distance < alternate_distance {
        Some(0)
    } else if alternate_distance < reference_distance {
        Some(1)
    } else {
        None
    })
}

fn detect_alleles(
    record: &Record,
    variants: &[PhaseVariant],
    reference: &[u8],
) -> Result<Vec<(i32, i32, i32)>> {
    let cigar: Vec<Cigar> = record
        .cigar()
        .iter()
        .map(|operation| {
            let code = match operation {
                bam::record::Cigar::Match(_) => 0,
                bam::record::Cigar::Ins(_) => 1,
                bam::record::Cigar::Del(_) => 2,
                bam::record::Cigar::RefSkip(_) => 3,
                bam::record::Cigar::SoftClip(_) => 4,
                bam::record::Cigar::HardClip(_) => 5,
                bam::record::Cigar::Pad(_) => 6,
                bam::record::Cigar::Equal(_) => 7,
                bam::record::Cigar::Diff(_) => 8,
            };
            (code, operation.len() as i64)
        })
        .collect();
    if cigar.is_empty() {
        return Ok(Vec::new());
    }
    let query = record.seq().as_bytes();
    let mut output = Vec::new();
    let mut variant_index = 0;
    while variant_index < variants.len() && variants[variant_index].position < record.pos() {
        variant_index += 1;
    }
    let mut reference_position = record.pos();
    let mut query_position = 0_i64;
    for (index, (operation, length)) in cigar.iter().copied().enumerate() {
        match operation {
            0 | 7 | 8 => {
                while variant_index < variants.len()
                    && variants[variant_index].position < reference_position + length
                {
                    let variant = &variants[variant_index];
                    let consumed = variant.position - reference_position;
                    if let Some(allele) = realign(
                        variant,
                        &cigar,
                        index,
                        consumed,
                        query_position + consumed,
                        &query,
                        reference,
                    )? {
                        output.push((variant.position as i32, allele, 30));
                    }
                    variant_index += 1;
                }
                query_position += length;
                reference_position += length;
            }
            1 => {
                if variant_index < variants.len()
                    && variants[variant_index].position == reference_position
                {
                    let variant = &variants[variant_index];
                    if let Some(allele) =
                        realign(variant, &cigar, index, 0, query_position, &query, reference)?
                    {
                        output.push((variant.position as i32, allele, 30));
                    }
                    variant_index += 1;
                }
                query_position += length;
            }
            2 => {
                while variant_index < variants.len()
                    && variants[variant_index].position < reference_position + length
                {
                    let variant = &variants[variant_index];
                    if let Some(allele) = realign(
                        variant,
                        &cigar,
                        index,
                        variant.position - reference_position,
                        query_position,
                        &query,
                        reference,
                    )? {
                        output.push((variant.position as i32, allele, 30));
                    }
                    variant_index += 1;
                }
                reference_position += length;
            }
            3 => {
                while variant_index < variants.len()
                    && variants[variant_index].position < reference_position + length
                {
                    variant_index += 1;
                }
                reference_position += length;
            }
            4 => query_position += length,
            5 | 6 => {}
            _ => bail!("unsupported CIGAR operation {operation}"),
        }
    }
    Ok(output)
}

fn tag_integer(record: &Record, tag: &[u8; 2], default: i32) -> Result<i32> {
    match record.aux(tag) {
        Ok(Aux::I8(value)) => Ok(value as i32),
        Ok(Aux::I16(value)) => Ok(value as i32),
        Ok(Aux::U8(value)) => Ok(value as i32),
        Ok(Aux::U16(value)) => Ok(value as i32),
        Ok(Aux::U32(value)) => Ok(value as i32),
        Ok(Aux::I32(value)) => Ok(value),
        Ok(Aux::String(value)) => Ok(value.parse()?),
        Ok(_) => bail!("unsupported tag type for {}", String::from_utf8_lossy(tag)),
        Err(_) => Ok(default),
    }
}

fn native_read(
    record: &Record,
    contig: &str,
    variants: &[PhaseVariant],
    reference: &[u8],
) -> Result<Option<NativeRead>> {
    let variants = detect_alleles(record, variants, reference)?;
    if variants.is_empty() {
        return Ok(None);
    }
    let bx = match record.aux(b"BX") {
        Ok(Aux::String(value)) => value.as_bytes().to_vec(),
        Ok(_) => bail!("invalid BX tag"),
        Err(_) => Vec::new(),
    };
    Ok(Some(NativeRead {
        name: record.qname().to_vec(),
        mapq: record.mapq(),
        source: 0,
        sample: 0,
        start: record.pos(),
        end: record.cigar().end_pos(),
        hp: tag_integer(record, b"HP", -1)?,
        ps: tag_integer(record, b"PS", -1)?,
        supplementary: record.flags() & 0x800 != 0,
        reverse: record.flags() & 0x10 != 0,
        bx,
        chromosome: contig.as_bytes().to_vec(),
        subalignment: Vec::new(),
        variants,
    }))
}

fn merge_group(group: Vec<NativeRead>) -> Option<NativeRead> {
    let mut primary = None;
    let mut primary_count = 0;
    let mut start = i64::MAX;
    for read in &group {
        start = start.min(read.start);
        if !read.supplementary {
            primary_count += 1;
            primary = Some(read);
        }
    }
    if primary_count == 0 || primary_count > 2 {
        return None;
    }
    let mut primary = primary.expect("a primary read is present").clone_metadata();
    let mut variants = BTreeMap::new();
    let mut conflicting = BTreeSet::new();
    for read in group {
        for (position, allele, quality) in read.variants {
            if let Some((_, old_allele)) = variants.get(&position)
                && *old_allele != allele
            {
                conflicting.insert(position);
            } else {
                variants.insert(position, (allele, quality));
            }
        }
    }
    primary.start = start;
    primary.variants = variants
        .into_iter()
        .filter(|(position, _)| !conflicting.contains(position))
        .map(|(position, (allele, quality))| (position, allele, quality))
        .collect();
    Some(primary)
}

impl NativeRead {
    fn clone_metadata(&self) -> Self {
        Self {
            name: self.name.clone(),
            mapq: self.mapq,
            source: self.source,
            sample: self.sample,
            start: self.start,
            end: self.end,
            hp: self.hp,
            ps: self.ps,
            supplementary: self.supplementary,
            reverse: self.reverse,
            bx: self.bx.clone(),
            chromosome: self.chromosome.clone(),
            subalignment: self.subalignment.clone(),
            variants: Vec::new(),
        }
    }
}

fn write_u32(output: &mut Vec<u8>, value: u32) {
    output.extend_from_slice(&value.to_le_bytes());
}

fn write_i32(output: &mut Vec<u8>, value: i32) {
    write_u32(output, value as u32);
}

fn write_text(output: &mut Vec<u8>, value: &[u8]) {
    write_u32(output, value.len() as u32);
    output.extend_from_slice(value);
}

fn write_read(output: &mut Vec<u8>, read: &NativeRead) {
    write_text(output, &read.name);
    write_u32(output, 1);
    write_i32(output, read.mapq as i32);
    write_i32(output, read.source);
    write_i32(output, read.sample);
    write_i32(output, read.start as i32);
    write_i32(output, read.end as i32);
    write_i32(output, read.hp);
    write_i32(output, read.ps);
    write_i32(output, i32::from(read.supplementary));
    write_i32(output, i32::from(read.reverse));
    write_text(output, &read.bx);
    write_text(output, &read.chromosome);
    write_text(output, &read.subalignment);
    write_u32(output, read.variants.len() as u32);
    for (position, allele, quality) in &read.variants {
        write_i32(output, *position);
        write_i32(output, *allele);
        write_i32(output, *quality);
    }
}

pub fn run(args: PhaseReadsArgs) -> Result<()> {
    let reference = Reference::open(&args.reference)?;
    let sequence = reference.fetch(&args.contig, 0, reference.contig_len(&args.contig) as i64)?;
    let variants = read_variants(&args.variants)?;
    let mut output = b"C3R1\0\0\0\0".to_vec();
    let mut read_count = 0_u32;
    let mut total = 0_usize;
    let mut covered = BTreeSet::new();
    for value in &args.regions {
        let (start, end, first, last) = parse_region(value)?;
        let variants = variants
            .get(first..last)
            .context("native phase region exceeds variant list")?;
        let mut reader = IndexedReader::from_path(&args.bam)?;
        reader.set_reference(&args.reference)?;
        reader.set_threads(1)?;
        reader.fetch((args.contig.as_bytes(), start, end))?;
        let mut groups: BTreeMap<Vec<u8>, Vec<NativeRead>> = BTreeMap::new();
        let mut record = Record::new();
        while let Some(result) = reader.read(&mut record) {
            result?;
            if record.flags() & (0x4 | 0x100 | 0x400 | 0x800) != 0
                || record.mapq() < 20
                || record.cigar().is_empty()
            {
                continue;
            }
            if let Some(read) = native_read(&record, &args.contig, variants, &sequence)? {
                groups.entry(read.name.clone()).or_default().push(read);
            }
        }
        for (_, group) in groups {
            let Some(read) = merge_group(group) else {
                continue;
            };
            if read.variants.is_empty() {
                continue;
            }
            total += 1;
            covered.extend(read.variants.iter().map(|(position, _, _)| *position));
            if read.variants.len() >= 2 {
                read_count += 1;
                write_read(&mut output, &read);
            }
        }
    }
    output[4..8].copy_from_slice(&read_count.to_le_bytes());
    fs::write(&args.output, output)?;
    let covered: Vec<String> = covered.iter().map(i32::to_string).collect();
    fs::write(&args.metadata, format!("{total}\n{}\n", covered.join(",")))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unit_edit_distance_matches_simple_cases() {
        assert_eq!(edit_distance(b"abc", b"abc"), 0);
        assert_eq!(edit_distance(b"abc", b"adc"), 1);
        assert_eq!(edit_distance(b"abc", b"ac"), 1);
        assert_eq!(edit_distance(b"a", b""), 1);
        assert_eq!(edit_distance(b"", b""), 0);
    }

    #[test]
    fn cigar_prefixes_include_split_operations() {
        let cigar = vec![(0, 3), (2, 1), (0, 6), (1, 2), (0, 4)];
        assert_eq!(prefix_length(split_left(&cigar, 2, 5), 10), (9, 8));
        assert_eq!(prefix_length(split_right(&cigar, 2, 5), 11), (5, 7));
    }

    #[test]
    fn snp_realignment_uses_unique_edit_distance() {
        let variant = PhaseVariant {
            position: 4,
            reference: b"C".to_vec(),
            alternate: b"T".to_vec(),
        };
        let cigar = vec![(0, 10)];
        let reference = b"AAAACAAAAA";
        assert_eq!(
            realign(&variant, &cigar, 0, 4, 4, b"AAAATAAAAA", reference).unwrap(),
            Some(1)
        );
        assert_eq!(
            realign(&variant, &cigar, 0, 4, 4, b"AAAACAAAAA", reference).unwrap(),
            Some(0)
        );
        assert_eq!(
            realign(&variant, &cigar, 0, 4, 4, b"AAAAGAAAAA", reference).unwrap(),
            None
        );
    }
}
