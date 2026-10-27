#![allow(clippy::approx_constant)]

use std::path::Path;

use anyhow::{Context, Result, bail};

use crate::{
    io::{interval::Region, reference::Reference, variant::VcfRecord},
    tensor::pileup::PileupBatch,
};

#[derive(Debug, Clone)]
pub struct GvcfOptions {
    pub base_error: f64,
    pub gq_bin_size: i32,
    pub bp_resolution: bool,
}

impl Default for GvcfOptions {
    fn default() -> Self {
        Self {
            base_error: 0.001,
            gq_bin_size: 5,
            bp_resolution: false,
        }
    }
}

#[derive(Debug, Clone)]
struct ReferenceSite {
    position: i64,
    reference: u8,
    genotype: String,
    raw_gq: i32,
    binned_gq: i32,
    depth: u64,
    pl: [i32; 3],
}

fn normalize_log10(values: [f64; 3]) -> [f64; 3] {
    const LOG_10: f64 = 2.3025;
    let maximum = values.into_iter().fold(f64::NEG_INFINITY, f64::max);
    let sum: f64 = values
        .into_iter()
        .map(|value| 10_f64.powf(value - maximum))
        .sum();
    let normalizer = (maximum + sum.ln() / LOG_10).round_to(6);
    values.map(|value| (value - normalizer).min(0.0))
}

trait RoundTo {
    fn round_to(self, digits: i32) -> f64;
}

impl RoundTo for f64 {
    fn round_to(self, digits: i32) -> f64 {
        let multiplier = 10_f64.powi(digits);
        (self * multiplier).round() / multiplier
    }
}

fn reference_site(
    contig: &str,
    position: i64,
    reference: u8,
    n_ref: u64,
    n_total: u64,
    options: &GvcfOptions,
) -> ReferenceSite {
    let n_alts = n_total.saturating_sub(n_ref);
    let log_error = options.base_error.ln() / 2.3025;
    let log_no_error = (-options.base_error).ln_1p() / 2.3025;
    let probabilities = if n_total == 0 {
        normalize_log10([-1.0, -1.0, -1.0])
    } else {
        normalize_log10([
            n_ref as f64 * log_no_error + n_alts as f64 * log_error,
            -(n_total as f64) * 0.3010,
            n_ref as f64 * log_error + n_alts as f64 * log_no_error,
        ])
    };
    let valid = probabilities[0] == probabilities.into_iter().fold(f64::NEG_INFINITY, f64::max);
    let reference_probability = 10_f64.powf(probabilities[0]);
    let raw_gq = if reference_probability == 1.0 {
        50
    } else {
        ((-10.0 * (-reference_probability).ln_1p() / 2.3025).round_to(6) as i32).min(50)
    };
    let binned_gq = if raw_gq >= 1 {
        (raw_gq - 1).div_euclid(options.gq_bin_size) * options.gq_bin_size + 1
    } else {
        0
    };
    let phred = probabilities.map(|value| -10.0 * value);
    let minimum = phred.into_iter().fold(f64::INFINITY, f64::min);
    let mut pl = phred.map(|value| (value - minimum) as i32);
    let mut raw_gq = raw_gq;
    let mut binned_gq = binned_gq;
    let mut reference = reference.to_ascii_uppercase();
    if !matches!(reference, b'A' | b'C' | b'G' | b'T') {
        reference = b'N';
        raw_gq = 1;
        binned_gq = 1;
        pl = [0, 0, 0];
    }
    let _ = contig;
    ReferenceSite {
        position,
        reference,
        genotype: if valid { "0/0" } else { "./." }.to_string(),
        raw_gq,
        binned_gq,
        depth: n_total,
        pl,
    }
}

fn site_record(contig: &str, site: &ReferenceSite, end: i64, contig_length: i64) -> VcfRecord {
    let end = if end == contig_length - 1 {
        contig_length
    } else {
        end
    };
    let genotype = if site.reference == b'N' {
        "./."
    } else {
        &site.genotype
    };
    VcfRecord {
        contig: contig.to_string(),
        position: site.position,
        id: ".".to_string(),
        reference: (site.reference as char).to_string(),
        alternate: "<NON_REF>".to_string(),
        quality: "0".to_string(),
        filter: ".".to_string(),
        info: format!("END={end}"),
        format: Some("GT:GQ:MIN_DP:PL".to_string()),
        samples: vec![format!(
            "{}:{}:{}:{},{},{}",
            genotype, site.binned_gq, site.depth, site.pl[0], site.pl[1], site.pl[2]
        )],
    }
}

fn write_block(
    output: &mut Vec<VcfRecord>,
    contig: &str,
    block: &[ReferenceSite],
    min_raw_gq: i32,
    min_depth: u64,
    bp_resolution: bool,
    contig_length: i64,
) {
    if block.is_empty() {
        return;
    }
    let first = &block[0];
    if (bp_resolution || first.genotype == "./.") && first.reference != b'N' {
        output.extend(
            block
                .iter()
                .map(|site| site_record(contig, site, site.position, contig_length)),
        );
        return;
    }
    let mut site = first.clone();
    site.binned_gq = if first.reference == b'N' {
        1
    } else {
        min_raw_gq
    };
    site.depth = min_depth;
    output.push(site_record(
        contig,
        &site,
        block[block.len() - 1].position,
        contig_length,
    ));
}

pub fn nonvariant_records(
    reference: &Reference,
    region: &Region,
    pileup: &PileupBatch,
    options: &GvcfOptions,
) -> Result<Vec<VcfRecord>> {
    if options.gq_bin_size <= 0 {
        bail!("GVCF GQ bin size must be positive");
    }
    if !(0.0..1.0).contains(&options.base_error) {
        bail!("GVCF base error must be between 0 and 1");
    }
    if pileup.gvcf_ref_counts.len() != (region.end - region.start) as usize
        || pileup.gvcf_total_counts.len() != pileup.gvcf_ref_counts.len()
    {
        bail!("pileup batch has no GVCF counts for the requested region");
    }
    let contig_length = reference.contig_len(&region.contig) as i64;
    if !pileup.has_pileup_rows {
        let site = ReferenceSite {
            position: (region.start + 1).max(1),
            reference: b'N',
            genotype: "./.".to_string(),
            raw_gq: 1,
            binned_gq: 1,
            depth: 0,
            pl: [0, 0, 0],
        };
        return Ok(vec![site_record(
            &region.contig,
            &site,
            region.end,
            contig_length,
        )]);
    }
    let sequence = reference.fetch_region(region)?;
    let mut output = Vec::new();
    let mut block: Vec<ReferenceSite> = Vec::new();
    let mut min_raw_gq = 0;
    let mut min_depth = 0_u64;
    let mut max_depth = 0_u64;
    for (index, reference_base) in sequence.iter().enumerate() {
        let site = reference_site(
            &region.contig,
            region.start + index as i64 + 1,
            *reference_base,
            pileup.gvcf_ref_counts[index],
            pileup.gvcf_total_counts[index],
            options,
        );
        let split = block.last().is_some_and(|previous: &ReferenceSite| {
            previous.binned_gq != site.binned_gq
                || previous.genotype != site.genotype
                || ((previous.reference == b'N') != (site.reference == b'N'))
                || (site.depth < min_depth && max_depth > (site.depth as f64 * 1.3).ceil() as u64)
                || (site.depth > max_depth && site.depth > (min_depth as f64 * 1.3).ceil() as u64)
        });
        if split {
            write_block(
                &mut output,
                &region.contig,
                &block,
                min_raw_gq,
                min_depth,
                options.bp_resolution,
                contig_length,
            );
            block.clear();
        }
        if block.is_empty() {
            min_depth = site.depth;
            max_depth = site.depth;
            min_raw_gq = site.raw_gq;
        } else {
            min_depth = min_depth.min(site.depth);
            max_depth = max_depth.max(site.depth);
            min_raw_gq = min_raw_gq.min(site.raw_gq);
        }
        block.push(site);
    }
    write_block(
        &mut output,
        &region.contig,
        &block,
        min_raw_gq,
        min_depth,
        options.bp_resolution,
        contig_length,
    );
    Ok(output)
}

fn gvcf_variant(record: &VcfRecord) -> Result<VcfRecord> {
    if record.alternate == "." || record.reference == record.alternate {
        bail!("reference records must not be converted to GVCF variants");
    }
    let mut record = record.clone();
    let format = record.format.as_deref().unwrap_or("");
    let fields: Vec<_> = format.split(':').collect();
    let index = |name: &str| {
        fields
            .iter()
            .position(|field| *field == name)
            .with_context(|| format!("missing {name} in variant FORMAT"))
    };
    let ad = index("AD")?;
    let af = index("AF")?;
    let pl = index("PL")?;
    let alt_count = record.alternate.split(',').count();
    let mut sample: Vec<_> = record.samples[0].split(':').map(str::to_string).collect();
    for field in [ad, af, pl] {
        if field >= sample.len() {
            bail!("variant sample is shorter than its FORMAT");
        }
    }
    sample[ad].push_str(",0");
    sample[af].push_str(",0.0000");
    sample[pl].push_str(&",990".repeat(alt_count + 2));
    record.alternate.push_str(",<NON_REF>");
    record.samples[0] = sample.join(":");
    Ok(record)
}

fn set_field(record: &mut VcfRecord, name: &str, value: &str) {
    if let Some(format) = &record.format {
        let fields: Vec<_> = format.split(':').collect();
        if let Some(index) = fields.iter().position(|field| *field == name) {
            let mut sample: Vec<_> = record.samples[0].split(':').map(str::to_string).collect();
            if index < sample.len() {
                sample[index] = value.to_string();
                record.samples[0] = sample.join(":");
            }
        }
    }
}

pub fn merge_nonvariants(
    variants: &[VcfRecord],
    nonvariants: &[VcfRecord],
    reference: &Reference,
    haploid: bool,
) -> Result<Vec<VcfRecord>> {
    let mut variants: Vec<VcfRecord> = variants
        .iter()
        .filter(|record| record.alternate != "." && record.reference != record.alternate)
        .map(gvcf_variant)
        .collect::<Result<_>>()?;
    variants.sort_by_key(|record| (record.contig.clone(), record.position));
    let mut variant_index = 0;
    let mut output = Vec::new();
    let mut pending_start: Option<(String, i64)> = None;
    for block in nonvariants {
        let end = block
            .info
            .strip_prefix("END=")
            .context("non-variant record has no END")?
            .parse::<i64>()?;
        let contig_length = reference.contig_len(&block.contig) as i64;
        let mut start = pending_start
            .take()
            .filter(|(contig, _)| *contig == block.contig)
            .map(|(_, start)| start)
            .unwrap_or(block.position);
        if start > end {
            pending_start = Some((block.contig.clone(), start));
            continue;
        }
        while variant_index < variants.len()
            && (variants[variant_index].contig < block.contig
                || (variants[variant_index].contig == block.contig
                    && variants[variant_index].position
                        + variants[variant_index].reference.len() as i64
                        - 1
                        < start))
        {
            output.push(variants[variant_index].clone());
            variant_index += 1;
        }
        while variant_index < variants.len()
            && variants[variant_index].contig == block.contig
            && variants[variant_index].position <= end
        {
            let variant = &variants[variant_index];
            let variant_end = variant.position + variant.reference.len() as i64 - 1;
            if variant.position > start {
                let mut left = block.clone();
                left.info = format!("END={}", variant.position - 1);
                output.push(left);
            }
            output.push(variant.clone());
            start = start.max(variant_end + 1);
            variant_index += 1;
        }
        if start <= end {
            let mut right = block.clone();
            if start != block.position {
                right.position = start;
                right.reference =
                    String::from_utf8(reference.fetch(&block.contig, start - 1, start)?)?;
            }
            if end == contig_length - 1 {
                right.info = format!("END={contig_length}");
            }
            output.push(right);
        } else {
            pending_start = Some((block.contig.clone(), start));
        }
    }
    output.extend(variants.into_iter().skip(variant_index));
    if haploid {
        for record in &mut output {
            if record.alternate != "<NON_REF>" {
                continue;
            }
            let genotype = record.samples[0]
                .split(':')
                .next()
                .unwrap_or("./.")
                .to_string();
            let genotype = match genotype.as_str() {
                "0/0" => "0",
                "./." => ".",
                other => other,
            };
            set_field(record, "GT", genotype);
        }
    }
    Ok(output)
}

pub fn write_gvcf(path: &Path, header: &str, records: &[VcfRecord]) -> Result<()> {
    use std::io::Write;
    let mut writer = std::io::BufWriter::new(std::fs::File::create(path)?);
    writer.write_all(header.as_bytes())?;
    for record in records {
        writeln!(writer, "{}", record.to_line())?;
    }
    writer.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn likelihoods_use_legacy_phred_and_gq_bins() {
        let options = GvcfOptions::default();
        let reference = reference_site("chr1", 1, b'A', 20, 20, &options);
        assert_eq!(reference.genotype, "0/0");
        assert_eq!(reference.pl, [0, 60, 599]);
        assert_eq!(reference.binned_gq, 46);
        let missing = reference_site("chr1", 1, b'N', 0, 0, &options);
        assert_eq!(missing.binned_gq, 1);
        assert_eq!(missing.pl, [0, 0, 0]);
    }

    #[test]
    fn gvcf_blocks_use_minimum_depth_and_bp_resolution_preserves_sites() {
        let reference = Reference::open(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures/data/toy.fa")
                .as_path(),
        )
        .unwrap();
        let region = Region::new("chrToy", 1000, 1003).unwrap();
        let batch = PileupBatch {
            gvcf_ref_counts: vec![30, 25, 28],
            gvcf_total_counts: vec![30, 25, 28],
            has_pileup_rows: true,
            ..Default::default()
        };
        let records =
            nonvariant_records(&reference, &region, &batch, &GvcfOptions::default()).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].samples[0].split(':').nth(2), Some("25"));

        let records = nonvariant_records(
            &reference,
            &region,
            &batch,
            &GvcfOptions {
                bp_resolution: true,
                ..Default::default()
            },
        )
        .unwrap();
        let depths: Vec<_> = records
            .iter()
            .map(|record| record.samples[0].split(':').nth(2).unwrap())
            .collect();
        assert_eq!(depths, ["30", "25", "28"]);
    }

    #[test]
    fn merging_skips_blocks_inside_deletions() {
        let reference = Reference::open(
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("tests/fixtures/data/toy.fa")
                .as_path(),
        )
        .unwrap();
        let variant = VcfRecord::from_line(
            "chrToy\t900\t.\tGCA\tG\t50\tPASS\t.\tGT:GQ:DP:AD:AF:PL\t1/1:50:12:0,12:1.0000:54,14,0",
        )
        .unwrap();
        let block = |position, end| {
            VcfRecord::from_line(&format!("chrToy\t{position}\t.\tA\t<NON_REF>\t0\t.\tEND={end}\tGT:GQ:MIN_DP:PL\t0/0:36:12:0,36,359")).unwrap()
        };
        let merged = merge_nonvariants(
            &[variant],
            &[block(900, 900), block(901, 902), block(903, 1000)],
            &reference,
            true,
        )
        .unwrap();
        assert_eq!(merged.len(), 2);
        assert_eq!(merged[0].position, 900);
        assert!(merged[0].samples[0].ends_with(",990,990,990"));
        assert_eq!(merged[1].position, 903);
        assert_eq!(merged[1].samples[0].split(':').next(), Some("0"));
    }
}
