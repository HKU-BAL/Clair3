use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};

use anyhow::Result;

use crate::io::variant::{VcfReader, VcfRecord};

#[derive(Debug, Default)]
pub struct VcfFile {
    pub header: Vec<String>,
    pub records: Vec<VcfRecord>,
}

pub fn read_vcf(path: &Path) -> Result<VcfFile> {
    let mut reader = VcfReader::open(path)?;
    let mut output = VcfFile {
        header: reader.header().to_vec(),
        ..Default::default()
    };
    while let Some(record) = reader.next_record()? {
        output.records.push(record);
    }
    Ok(output)
}

#[derive(Debug, Clone)]
pub struct MergeOptions {
    pub print_ref_calls: bool,
    pub haploid_precise: bool,
    pub haploid_sensitive: bool,
    pub qual: Option<f64>,
}

impl Default for MergeOptions {
    fn default() -> Self {
        Self {
            print_ref_calls: false,
            haploid_precise: false,
            haploid_sensitive: false,
            qual: Some(2.0),
        }
    }
}

fn is_reference(record: &VcfRecord) -> bool {
    record.alternate == "." || record.reference == record.alternate
}

fn update_haploid(record: &mut VcfRecord, options: &MergeOptions) -> bool {
    let genotype = record.samples[0]
        .split(':')
        .next()
        .unwrap_or("")
        .replace('|', "/");
    if options.haploid_precise {
        let genotype = match genotype.as_str() {
            "1/1" => "1",
            "0/0" => "0",
            _ => return false,
        };
        let fields: Vec<_> = record.samples[0].split(':').collect();
        record.samples[0] = std::iter::once(genotype)
            .chain(fields[1..].iter().copied())
            .collect::<Vec<_>>()
            .join(":");
    }
    if options.haploid_sensitive {
        if record.alternate.contains(',') {
            return false;
        }
        let genotype = if matches!(genotype.as_str(), "0/1" | "1/0" | "1/1") {
            "1"
        } else {
            "0"
        };
        let fields: Vec<_> = record.samples[0].split(':').collect();
        record.samples[0] = std::iter::once(genotype)
            .chain(fields[1..].iter().copied())
            .collect::<Vec<_>>()
            .join(":");
    }
    true
}

fn mark_low_qual(record: &mut VcfRecord, options: &MergeOptions) {
    if let (Some(threshold), Ok(quality)) = (options.qual, record.quality.parse::<f64>())
        && threshold != 0.0
        && quality <= threshold
    {
        record.filter = "LowQual".to_string();
    }
}

fn process(record: &mut VcfRecord, options: &MergeOptions) -> bool {
    if !update_haploid(record, options) {
        return false;
    }
    if !is_reference(record) {
        mark_low_qual(record, options);
    }
    options.print_ref_calls || !is_reference(record)
}

pub fn merge_records(
    pileup: &[VcfRecord],
    full_alignment: &[VcfRecord],
    options: &MergeOptions,
) -> Vec<VcfRecord> {
    let full_positions: BTreeSet<(&str, i64)> = full_alignment
        .iter()
        .map(|record| (record.contig.as_str(), record.position))
        .collect();
    let mut output = Vec::new();
    for record in full_alignment {
        let mut record = record.clone();
        if process(&mut record, options) {
            output.push(record);
        }
    }
    for record in pileup {
        if full_positions.contains(&(record.contig.as_str(), record.position)) {
            continue;
        }
        let mut record = record.clone();
        if process(&mut record, options) {
            output.push(record);
        }
    }
    output.sort_by_key(|record| (record.contig.clone(), record.position));
    output
}

pub fn merge_illumina_records(
    pileup: &[VcfRecord],
    full_alignment: &[VcfRecord],
    full_regions: &crate::io::interval::IntervalMap,
    contig: &str,
    options: &MergeOptions,
) -> Vec<VcfRecord> {
    let mut output = BTreeMap::new();
    for record in pileup {
        if full_regions.contains_point_shrunk(contig, record.position, 33) {
            continue;
        }
        let mut record = record.clone();
        if process(&mut record, options) {
            output.insert(record.position, record);
        }
    }
    for record in full_alignment {
        if !full_regions.contains_point_shrunk(contig, record.position, 33) {
            continue;
        }
        let mut record = record.clone();
        if process(&mut record, options) {
            output.insert(record.position, record);
        }
    }
    output.into_values().collect()
}

pub fn sort_records(records: &mut Vec<VcfRecord>, contigs: &[String]) {
    let major_order = |contig: &str| {
        let name = contig.strip_prefix("chr").unwrap_or(contig);
        let suffix = match name.parse::<usize>() {
            Ok(value) => value,
            Err(_) if name == "X" => 23,
            Err(_) if name == "Y" => 24,
            Err(_) => return usize::MAX,
        };
        if contig.starts_with("chr") {
            suffix - 1
        } else {
            suffix + 23
        }
    };
    let mut ordered = BTreeMap::new();
    for record in records.drain(..) {
        let contig_order = contigs
            .iter()
            .position(|candidate| *candidate == record.contig)
            .unwrap_or(usize::MAX);
        ordered.insert(
            (
                major_order(&record.contig),
                contig_order,
                record.contig.clone(),
                record.position,
            ),
            record,
        );
    }
    records.extend(ordered.into_values());
}
