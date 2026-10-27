use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result};
use rust_htslib::faidx;

use super::interval::Region;

pub struct Reference {
    path: PathBuf,
    reader: faidx::Reader,
}

impl Reference {
    pub fn open(path: &Path) -> Result<Self> {
        let reader = faidx::Reader::from_path(path)
            .with_context(|| format!("failed to open indexed FASTA {}", path.display()))?;
        Ok(Self {
            path: path.to_path_buf(),
            reader,
        })
    }

    pub fn path(&self) -> &Path {
        &self.path
    }

    pub fn contigs(&self) -> Result<Vec<String>> {
        self.reader.seq_names().map_err(anyhow::Error::from)
    }

    pub fn contig_len(&self, contig: &str) -> u64 {
        self.reader.fetch_seq_len(contig)
    }

    pub fn contig_lengths(&self) -> Result<BTreeMap<String, u64>> {
        let mut lengths = BTreeMap::new();
        for contig in self.contigs()? {
            lengths.insert(contig.clone(), self.contig_len(&contig));
        }
        Ok(lengths)
    }

    pub fn fetch(&self, contig: &str, start: i64, end: i64) -> Result<Vec<u8>> {
        if start < 0 || end < start {
            anyhow::bail!("invalid FASTA interval {contig}:{start}-{end}")
        }
        if end == start {
            return Ok(Vec::new());
        }
        let sequence = self
            .reader
            .fetch_seq(contig, start as usize, end as usize - 1)?;
        Ok(sequence
            .into_iter()
            .map(|base| base.to_ascii_uppercase())
            .collect())
    }

    pub fn fetch_region(&self, region: &Region) -> Result<Vec<u8>> {
        self.fetch(&region.contig, region.start, region.end)
    }

    pub fn fetch_padded(
        &self,
        contig: &str,
        start: i64,
        end: i64,
        padding_base: u8,
    ) -> Result<Vec<u8>> {
        if end < start {
            anyhow::bail!("invalid padded FASTA interval {contig}:{start}-{end}")
        }
        let length = self.contig_len(contig) as i64;
        let clipped_start = start.clamp(0, length);
        let clipped_end = end.clamp(0, length);
        let mut sequence = vec![padding_base; (clipped_start - start).max(0) as usize];
        if clipped_end > clipped_start {
            sequence.extend(self.fetch(contig, clipped_start, clipped_end)?);
        }
        sequence.extend(std::iter::repeat_n(
            padding_base,
            (end - clipped_end).max(0) as usize,
        ));
        Ok(sequence)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_fixture_reference_and_boundaries() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/data/toy.fa");
        let reference = Reference::open(&path).unwrap();
        assert_eq!(reference.contig_len("chrToy"), 4096);
        assert_eq!(reference.contig_len("chr:colon"), 1024);
        assert_eq!(reference.fetch("chrToy", 0, 1).unwrap().len(), 1);
        assert_eq!(reference.fetch("chrToy", 4095, 4096).unwrap().len(), 1);
        let padded = reference.fetch_padded("chrToy", -2, 2, b'A').unwrap();
        assert_eq!(padded.len(), 4);
        assert_eq!(&padded[..2], b"AA");
        let region = crate::io::interval::Region::new("chrToy", 0, 10).unwrap();
        assert_eq!(reference.fetch_region(&region).unwrap().len(), 10);
    }
}
