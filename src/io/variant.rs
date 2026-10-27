use std::{collections::BTreeSet, io::BufRead, path::Path};

use anyhow::{Result, bail};

use super::{compression::open_text, reference::Reference};

#[derive(Debug, Clone, PartialEq)]
pub struct VcfRecord {
    pub contig: String,
    pub position: i64,
    pub id: String,
    pub reference: String,
    pub alternate: String,
    pub quality: String,
    pub filter: String,
    pub info: String,
    pub format: Option<String>,
    pub samples: Vec<String>,
}

impl VcfRecord {
    pub fn from_line(line: &str) -> Result<Self> {
        let columns: Vec<&str> = line.trim_end().split('\t').collect();
        if columns.len() < 8 {
            bail!("invalid VCF record with {} columns", columns.len())
        }
        Ok(Self {
            contig: columns[0].to_string(),
            position: columns[1].parse()?,
            id: columns[2].to_string(),
            reference: columns[3].to_string(),
            alternate: columns[4].to_string(),
            quality: columns[5].to_string(),
            filter: columns[6].to_string(),
            info: columns[7].to_string(),
            format: columns.get(8).map(|value| value.to_string()),
            samples: columns
                .get(9..)
                .unwrap_or(&[])
                .iter()
                .map(|value| value.to_string())
                .collect(),
        })
    }

    pub fn position_zero_based(&self) -> i64 {
        self.position - 1
    }

    pub fn to_line(&self) -> String {
        let mut columns = vec![
            self.contig.clone(),
            self.position.to_string(),
            self.id.clone(),
            self.reference.clone(),
            self.alternate.clone(),
            self.quality.clone(),
            self.filter.clone(),
            self.info.clone(),
        ];
        if let Some(format) = &self.format {
            columns.push(format.clone());
            columns.extend(self.samples.clone());
        }
        columns.join("\t")
    }

    pub fn sample_field(&self, sample: usize, field: usize) -> Option<&str> {
        self.samples.get(sample)?.split(':').nth(field)
    }
}

pub struct VcfReader {
    reader: Box<dyn BufRead>,
    header: Vec<String>,
    pending: Option<String>,
}

impl VcfReader {
    pub fn open(path: &Path) -> Result<Self> {
        let mut reader = open_text(path)?;
        let mut header = vec![];
        let mut pending = None;
        loop {
            let mut line = String::new();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            let line = line.trim_end_matches(['\r', '\n']).to_string();
            if line.starts_with('#') {
                header.push(line);
            } else if !line.is_empty() {
                pending = Some(line);
                break;
            }
        }
        Ok(Self {
            reader,
            header,
            pending,
        })
    }

    pub fn header(&self) -> &[String] {
        &self.header
    }

    pub fn next_record(&mut self) -> Result<Option<VcfRecord>> {
        loop {
            let line = match self.pending.take() {
                Some(line) => line,
                None => {
                    let mut line = String::new();
                    if self.reader.read_line(&mut line)? == 0 {
                        return Ok(None);
                    }
                    line.trim_end_matches(['\r', '\n']).to_string()
                }
            };
            if line.is_empty() {
                continue;
            }
            if line.starts_with('#') {
                bail!("VCF header line found after records")
            }
            return Ok(Some(VcfRecord::from_line(&line)?));
        }
    }
}

pub fn candidate_positions(path: &Path, contig: Option<&str>) -> Result<Vec<i64>> {
    let mut reader = VcfReader::open(path)?;
    let mut positions = BTreeSet::new();
    while let Some(record) = reader.next_record()? {
        if contig.is_some_and(|name| name != record.contig) {
            continue;
        }
        positions.insert(record.position);
    }
    Ok(positions.into_iter().collect())
}

pub fn vcf_header(
    reference: &Reference,
    command: Option<&str>,
    sample_name: &str,
    gvcf: bool,
    contigs: &[String],
) -> Result<String> {
    let mut lines = vec![
        "##fileformat=VCFv4.2".to_string(),
        "##source=Clair3".to_string(),
        "##clair3_version=1.2.0".to_string(),
    ];
    if let Some(command) = command {
        lines.push(format!("##cmdline={command}"));
    }
    lines.push(format!("##reference={}", reference.path().display()));
    lines.extend([
        "##FILTER=<ID=PASS,Description=\"All filters passed\">".to_string(),
        "##FILTER=<ID=LowQual,Description=\"Low quality variant\">".to_string(),
        "##FILTER=<ID=RefCall,Description=\"Reference call\">".to_string(),
        "##INFO=<ID=P,Number=0,Type=Flag,Description=\"Result from pileup calling\">".to_string(),
        "##INFO=<ID=F,Number=0,Type=Flag,Description=\"Result from full-alignment calling\">"
            .to_string(),
    ]);
    if gvcf {
        lines.extend([
            "##ALT=<ID=NON_REF,Description=\"Represents any possible alternative allele at this location\">".to_string(),
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position (for use with symbolic alleles)\">".to_string(),
        ]);
    }
    lines.extend([
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">".to_string(),
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">".to_string(),
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads 1. with MQ below 5 or an user-specified threshold, or 2. selected by 'samtools view -F 2316', are filtered)\">".to_string(),
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">".to_string(),
    ]);
    if gvcf {
        lines.push("##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP observed within the GVCF block\">".to_string());
    }
    lines.extend([
        "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">".to_string(),
        "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Observed allele frequency in reads, for each ALT allele, in the same order as listed, or the REF allele for a RefCall\">".to_string(),
    ]);
    for contig in contigs {
        lines.push(format!(
            "##contig=<ID={contig},length={}>",
            reference.contig_len(contig)
        ));
    }
    lines.push(format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}"
    ));
    Ok(lines.join("\n") + "\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::interval::Region;

    #[test]
    fn streams_vcf_and_preserves_records() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/data/known.vcf.gz");
        let mut reader = VcfReader::open(&path).unwrap();
        assert!(
            reader
                .header()
                .iter()
                .any(|line| line == "##fileformat=VCFv4.2")
        );
        let first = reader.next_record().unwrap().unwrap();
        assert_eq!(first.contig, "chrToy");
        assert_eq!(first.position, 5);
        assert_eq!(first.to_line(), "chrToy\t5\t.\tA\tC\t50\tPASS\t.\tGT\t0/1");
        assert_eq!(first.sample_field(0, 0), Some("0/1"));
        let positions = candidate_positions(&path, Some("chr:colon")).unwrap();
        assert_eq!(positions, vec![100]);
    }

    #[test]
    fn builds_legacy_header_in_fasta_order() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/data/toy.fa");
        let reference = Reference::open(&path).unwrap();
        let contigs = reference.contigs().unwrap();
        let header = vcf_header(&reference, Some("clair3 test"), "TOY", true, &contigs).unwrap();
        assert!(header.contains("##cmdline=clair3 test\n"));
        assert!(header.contains("##ALT=<ID=NON_REF"));
        assert!(header.contains("##contig=<ID=chrToy,length=4096>\n##contig=<ID=chr:colon,length=1024>\n##contig=<ID=chrEmpty,length=512>"));
        assert!(header.ends_with("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTOY\n"));
        let region = Region::new("chrToy", 0, 1).unwrap();
        assert_eq!(region.to_one_based_string(), "chrToy:1-1");
    }
}
