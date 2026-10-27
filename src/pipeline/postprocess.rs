use std::{
    collections::{BTreeMap, BTreeSet},
    path::Path,
};

use anyhow::{Context, Result, bail};
use rayon::prelude::*;
use rust_htslib::bam::{self, IndexedReader, Read as BamRead};

use crate::{
    io::{interval::IntervalMap, variant::VcfRecord},
    pipeline::merge::read_vcf,
};

const NT16: &[u8; 16] = b"=ACMGRSVTWYHKDBN";
const FILTER_FLAGS: u16 = 0x4 | 0x8 | 0x100 | 0x800;

pub fn add_back_missing_variants(
    known: &[VcfRecord],
    called: &[VcfRecord],
    switch_genotype: bool,
) -> Result<Vec<VcfRecord>> {
    let called: BTreeMap<_, _> = called
        .iter()
        .map(|record| ((record.contig.clone(), record.position), record.clone()))
        .collect();
    let mut output = Vec::new();
    for known in known {
        if let Some(record) = called.get(&(known.contig.clone(), known.position)) {
            output.push(record.clone());
        } else if switch_genotype {
            let mut record = known.clone();
            record.reference = record.reference.chars().next().unwrap_or('.').to_string();
            record.alternate = ".".to_string();
            record.quality = ".".to_string();
            record.filter = ".".to_string();
            record.info = ".".to_string();
            record.format = Some("GT".to_string());
            record.samples = vec!["./.".to_string()];
            output.push(record);
        } else {
            output.push(known.clone());
        }
    }
    Ok(output)
}

fn parallel_map<T, F>(records: &[VcfRecord], threads: usize, function: F) -> Result<Vec<T>>
where
    T: Send,
    F: Fn(&VcfRecord) -> Result<T> + Sync + Send,
{
    if threads == 0 {
        bail!("threads must be positive");
    }
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()?
        .install(|| records.par_iter().map(function).collect())
}

fn pair_end_depth(bam: &Path, record: &VcfRecord, min_mq: u8, min_bq: u8) -> Result<usize> {
    let mut reader = IndexedReader::from_path(bam)?;
    reader.fetch((
        record.contig.as_bytes(),
        record.position - 1,
        record.position,
    ))?;
    let alternate = record.alternate.split(',').next().unwrap_or("");
    let reference = record.reference.as_bytes();
    let alternate = alternate.as_bytes();
    let mut names = BTreeSet::new();
    for pileup in reader.pileup() {
        let pileup = pileup?;
        if pileup.pos() as i64 != record.position - 1 {
            continue;
        }
        for alignment in pileup.alignments() {
            let read = alignment.record();
            if read.mapq() < min_mq || read.flags() & FILTER_FLAGS != 0 {
                continue;
            }
            let Some(qpos) = alignment.qpos() else {
                continue;
            };
            if read.qual().get(qpos).copied().unwrap_or(0) < min_bq {
                continue;
            }
            let base = NT16[read.seq().encoded_base(qpos) as usize];
            let matches = if reference.len() == 1 && alternate.len() == 1 {
                base == alternate[0] && matches!(alignment.indel(), bam::pileup::Indel::None)
            } else if reference.len() == 1 && alternate.len() > 1 {
                match alignment.indel() {
                    bam::pileup::Indel::Ins(length) => {
                        length as usize + 1 == alternate.len()
                            && base == reference[0]
                            && read
                                .seq()
                                .as_bytes()
                                .get(qpos + 1..qpos + length as usize + 1)
                                == Some(&alternate[1..])
                    }
                    _ => false,
                }
            } else {
                match alignment.indel() {
                    bam::pileup::Indel::Del(length) => {
                        length as usize + 1 == reference.len() && base == reference[0]
                    }
                    _ => false,
                }
            };
            if matches {
                names.insert(read.qname().to_vec());
            }
        }
    }
    Ok(names.len())
}

pub fn add_pair_end_depth(
    bam: &Path,
    records: &[VcfRecord],
    contig: Option<&str>,
    filter_tag: Option<&str>,
    threads: usize,
    min_mq: u8,
    min_bq: u8,
) -> Result<Vec<VcfRecord>> {
    let candidates: Vec<VcfRecord> = records
        .iter()
        .filter(|record| {
            contig.is_none_or(|name| name == record.contig)
                && filter_tag.is_none_or(|tags| tags.split(',').any(|tag| tag == record.filter))
                && !matches!(
                    record
                        .samples
                        .first()
                        .and_then(|sample| sample.split(':').next()),
                    Some("0/0" | "0|0")
                )
        })
        .cloned()
        .collect();
    let depths = parallel_map(&candidates, threads, |record| {
        pair_end_depth(bam, record, min_mq, min_bq)
    })?;
    Ok(candidates
        .into_iter()
        .zip(depths)
        .map(|(mut record, depth)| {
            if let Some(format) = &mut record.format {
                format.push_str(":PEAD");
            }
            if !record.samples.is_empty() {
                record.samples[0].push_str(&format!(":{depth}"));
            }
            record
        })
        .collect())
}

#[derive(Debug)]
struct SvInterval {
    contig: String,
    start: i64,
    end: i64,
    quality: Option<f64>,
}

impl SvInterval {
    fn contains(&self, contig: &str, position: i64) -> bool {
        self.contig == contig && self.start <= position && position < self.end
    }
}

fn sv_intervals(
    vcf: Option<&Path>,
    bed: Option<&Path>,
    filter: &str,
    alt_tag: &str,
) -> Result<Vec<SvInterval>> {
    if let Some(path) = bed {
        let mut output = Vec::new();
        for contig in IntervalMap::from_bed(path, None)?
            .contigs()
            .map(str::to_string)
            .collect::<Vec<_>>()
        {
            let map = IntervalMap::from_bed(path, Some(&contig))?;
            for region in map.regions(&contig) {
                output.push(SvInterval {
                    contig: contig.clone(),
                    start: region.start,
                    end: region.end,
                    quality: None,
                });
            }
        }
        return Ok(output);
    }
    let Some(path) = vcf else {
        return Ok(Vec::new());
    };
    let mut output = Vec::new();
    for record in read_vcf(path)?.records {
        if record.filter != filter || record.alternate != alt_tag {
            continue;
        }
        let end = record
            .info
            .split(';')
            .find_map(|field| field.strip_prefix("END="))
            .context("SV record has no END")?
            .parse::<i64>()?;
        if end < record.position {
            continue;
        }
        output.push(SvInterval {
            contig: record.contig,
            start: record.position,
            end,
            quality: record.quality.parse().ok(),
        });
    }
    Ok(output)
}

fn long_deletion_supported(
    bam: &Path,
    record: &VcfRecord,
    flanking: i64,
    min_mq: u8,
    min_bq: u8,
    min_proportion: f64,
) -> Result<bool> {
    let mut reader = IndexedReader::from_path(bam)?;
    let start = record.position - flanking;
    let end = record.position + flanking;
    reader.fetch((record.contig.as_bytes(), start - 1, end))?;
    for pileup in reader.pileup() {
        let pileup = pileup?;
        if (pileup.pos() as i64) < start - 1 || pileup.pos() as i64 >= end {
            continue;
        }
        let mut depth = 0_usize;
        let mut deleted = 0_usize;
        for alignment in pileup.alignments() {
            let record = alignment.record();
            if record.mapq() < min_mq || record.flags() & FILTER_FLAGS != 0 {
                continue;
            }
            if let Some(qpos) = alignment.qpos()
                && record.qual().get(qpos).copied().unwrap_or(0) < min_bq
            {
                continue;
            }
            depth += 1;
            deleted += usize::from(alignment.is_del());
        }
        if deleted as f64 / depth.max(1) as f64 >= min_proportion {
            return Ok(true);
        }
    }
    Ok(false)
}

#[allow(clippy::too_many_arguments)]
pub fn switch_zygosity(
    bam: &Path,
    records: &[VcfRecord],
    contig: Option<&str>,
    sv_vcf: Option<&Path>,
    sv_bed: Option<&Path>,
    use_sv_quality: bool,
    filter_tag: Option<&str>,
    sv_filter_tag: &str,
    sv_alt_tag: &str,
    max_af: f64,
    min_deletion_proportion: f64,
    flanking: i64,
    threads: usize,
    min_mq: u8,
    min_bq: u8,
) -> Result<Vec<VcfRecord>> {
    let intervals = sv_intervals(sv_vcf, sv_bed, sv_filter_tag, sv_alt_tag)?;
    let has_sv_source = sv_vcf.is_some() || sv_bed.is_some();
    let candidates: Vec<VcfRecord> = records
        .iter()
        .filter(|record| {
            if contig.is_some_and(|name| name != record.contig)
                || filter_tag.is_some_and(|tag| tag != record.filter)
                || record.samples.is_empty()
                || record.filter != "PASS"
                || record.reference.len() != 1
                || record.alternate.split(',').next().unwrap_or("").len() != 1
            {
                return false;
            }
            let alleles = record.samples[0]
                .split(':')
                .next()
                .unwrap_or("")
                .replace('/', "|")
                .replace('.', "0");
            let mut alleles = alleles
                .split('|')
                .filter_map(|allele| allele.parse::<i32>().ok())
                .collect::<Vec<_>>();
            alleles.sort_unstable();
            let af = record
                .format
                .as_deref()
                .and_then(|format| format.split(':').position(|field| field == "AF"))
                .and_then(|index| record.samples[0].split(':').nth(index))
                .and_then(|value| value.parse::<f64>().ok())
                .unwrap_or(0.0);
            alleles.iter().sum::<i32>() == 2
                && af <= max_af
                && (!has_sv_source
                    || intervals
                        .iter()
                        .any(|interval| interval.contains(&record.contig, record.position)))
        })
        .cloned()
        .collect();
    let switched = parallel_map(&candidates, threads, |record| {
        long_deletion_supported(
            bam,
            record,
            flanking,
            min_mq,
            min_bq,
            min_deletion_proportion,
        )
    })?;
    let switched: BTreeSet<(String, i64)> = candidates
        .iter()
        .zip(switched)
        .filter(|(_, supported)| *supported)
        .map(|(record, _)| (record.contig.clone(), record.position))
        .collect();
    let use_quality = use_sv_quality && sv_bed.is_none();
    let mut output = Vec::new();
    for record in records {
        if contig.is_some_and(|name| name != record.contig)
            || filter_tag.is_some_and(|tag| tag != record.filter)
        {
            continue;
        }
        let mut record = record.clone();
        if switched.contains(&(record.contig.clone(), record.position)) {
            let sample = record.samples[0].replacen("1/1", "0/1", 1);
            record.samples[0] = sample;
            if use_quality {
                let quality = intervals
                    .iter()
                    .find(|interval| interval.contains(&record.contig, record.position))
                    .and_then(|interval| interval.quality);
                let Some(quality) = quality else {
                    continue;
                };
                let original = record.quality.clone();
                record.quality = format!("{quality:.2}");
                record
                    .info
                    .push_str(&format!(";SVBASEDHET;ORG_CLAIR3_SCORE={original}"));
            }
        }
        output.push(record);
    }
    Ok(output)
}

pub fn before_chrom_header(header: &[String], additions: &[&str]) -> String {
    let mut lines = header.to_vec();
    let index = lines
        .iter()
        .position(|line| line.starts_with("#CHROM"))
        .unwrap_or(lines.len());
    lines.splice(
        index..index,
        additions.iter().map(|line| (*line).to_string()),
    );
    lines.join("\n") + "\n"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pair_end_depth_counts_supporting_reads() {
        let bam = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/data/toy.bam");
        let record = VcfRecord::from_line("chrToy\t5\t.\tA\tC\t50\tPASS\t.\tGT\t0/1").unwrap();
        assert_eq!(pair_end_depth(&bam, &record, 5, 0).unwrap(), 6);
        let record = VcfRecord::from_line("chrToy\t1100\t.\tA\tC,G\t50\tPASS\t.\tGT\t1/2").unwrap();
        assert_eq!(pair_end_depth(&bam, &record, 5, 0).unwrap(), 4);
        let record =
            VcfRecord::from_line("chrToy\t900\t.\tGCA\tG\t50\tPASS\t.\tGT:AF\t1/1:0.6").unwrap();
        assert!(long_deletion_supported(&bam, &record, 16, 5, 0, 0.3).unwrap());
    }
}
