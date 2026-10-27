use std::{collections::BTreeMap, path::Path};

use anyhow::{Result, bail};

use super::compression::read_text_lines;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub contig: String,
    pub start: i64,
    pub end: i64,
}

impl Region {
    pub fn new(contig: impl Into<String>, start: i64, end: i64) -> Result<Self> {
        let contig = contig.into();
        if start < 0 || end < start {
            bail!("invalid region {contig}:{start}-{end}")
        }
        Ok(Self { contig, start, end })
    }

    pub fn from_bed(line: &str) -> Result<(Self, Option<String>)> {
        let columns: Vec<&str> = line.split_whitespace().collect();
        if columns.len() < 3 {
            bail!("invalid BED row '{line}'")
        }
        let start = columns[1].parse()?;
        let end = columns[2].parse()?;
        let region = Self::new(columns[0], start, end)?;
        let metadata = (columns.len() > 3).then(|| columns[3..].join("\t"));
        Ok((region, metadata))
    }

    pub fn from_one_based(contig: &str, start: i64, end: i64) -> Result<Self> {
        if start < 1 || end < start {
            bail!("invalid one-based region {contig}:{start}-{end}")
        }
        Self::new(contig, start - 1, end)
    }

    pub fn from_region_string(value: &str) -> Result<Self> {
        let Some((contig, coordinates)) = value.rsplit_once(':') else {
            bail!("region '{value}' has no contig coordinates")
        };
        let Some((start, end)) = coordinates.split_once('-') else {
            bail!("region '{value}' has no start/end coordinates")
        };
        Self::from_one_based(contig, start.parse()?, end.parse()?)
    }

    pub fn to_one_based_string(&self) -> String {
        format!("{}:{}-{}", self.contig, self.start + 1, self.end)
    }

    pub fn overlaps(&self, start: i64, end: i64) -> bool {
        self.start < end && start < self.end
    }
}

#[derive(Debug, Default)]
pub struct IntervalMap {
    contigs: BTreeMap<String, Vec<Region>>,
    unmerged: BTreeMap<String, Vec<Region>>,
    prefix_end: BTreeMap<String, Vec<i64>>,
}

impl IntervalMap {
    pub fn insert(&mut self, region: Region) {
        self.contigs
            .entry(region.contig.clone())
            .or_default()
            .push(region);
    }

    pub fn merge_overlaps(&mut self) {
        for intervals in self.contigs.values_mut() {
            intervals.sort_by_key(|region| (region.start, region.end));
            let mut merged: Vec<Region> = Vec::with_capacity(intervals.len());
            for region in intervals.drain(..) {
                if let Some(last) = merged.last_mut()
                    && region.start <= last.end
                {
                    last.end = last.end.max(region.end);
                    continue;
                }
                merged.push(region);
            }
            *intervals = merged;
        }
    }

    pub fn contains_point(&self, contig: &str, position: i64) -> bool {
        let Some(intervals) = self.contigs.get(contig) else {
            return false;
        };
        let index = intervals.partition_point(|region| region.start <= position);
        index > 0 && intervals[index - 1].end > position
    }

    pub fn contains_point_shrunk(&self, contig: &str, position: i64, padding: i64) -> bool {
        let Some(intervals) = self.unmerged.get(contig) else {
            return self.contigs.get(contig).is_some_and(|intervals| {
                let index = intervals.partition_point(|region| region.start + padding <= position);
                index > 0 && intervals[index - 1].end - padding > position
            });
        };
        let index = intervals.partition_point(|region| region.start + padding <= position);
        index > 0 && self.prefix_end[contig][index - 1] > position + padding
    }

    pub fn overlaps(&self, contig: &str, start: i64, end: i64) -> bool {
        let Some(intervals) = self.contigs.get(contig) else {
            return false;
        };
        let index = intervals.partition_point(|region| region.start < end);
        index > 0
            && intervals[..index]
                .iter()
                .rev()
                .take_while(|region| region.end > start)
                .any(|region| region.overlaps(start, end))
    }

    pub fn from_bed(path: &Path, contig: Option<&str>) -> Result<Self> {
        let mut map = Self::default();
        for line in read_text_lines(path)? {
            let line = line?;
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let (mut region, _) = Region::from_bed(&line)?;
            if contig.is_some_and(|name| name != region.contig) {
                continue;
            }
            if region.start == region.end {
                region.end += 1;
            }
            map.unmerged
                .entry(region.contig.clone())
                .or_default()
                .push(region.clone());
            map.insert(region);
        }
        map.merge_overlaps();
        for (contig, intervals) in &mut map.unmerged {
            intervals.sort_by_key(|region| (region.start, region.end));
            let mut maximum_end = i64::MIN;
            map.prefix_end.insert(
                contig.clone(),
                intervals
                    .iter()
                    .map(|region| {
                        maximum_end = maximum_end.max(region.end);
                        maximum_end
                    })
                    .collect(),
            );
        }
        Ok(map)
    }

    pub fn contigs(&self) -> impl Iterator<Item = &str> {
        self.contigs.keys().map(String::as_str)
    }

    pub fn regions(&self, contig: &str) -> &[Region] {
        if let Some(regions) = self.unmerged.get(contig) {
            return regions;
        }
        self.contigs.get(contig).map_or(&[], Vec::as_slice)
    }

    pub fn bounds(&self, contig: &str) -> Option<(i64, i64)> {
        let intervals = self.contigs.get(contig)?;
        Some((
            intervals.first()?.start,
            intervals.iter().map(|region| region.end).max()?,
        ))
    }

    pub fn is_empty(&self) -> bool {
        self.contigs.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_regions_and_merges_intervals() {
        let (region, metadata) = Region::from_bed("chr:colon\t10\t20\tpayload").unwrap();
        assert_eq!(region.contig, "chr:colon");
        assert_eq!(metadata.as_deref(), Some("payload"));
        assert_eq!(
            Region::from_region_string("chr:colon:11-20").unwrap(),
            region
        );
        assert_eq!(region.to_one_based_string(), "chr:colon:11-20");

        let mut map = IntervalMap::default();
        map.insert(Region::new("chr1", 0, 10).unwrap());
        map.insert(Region::new("chr1", 10, 20).unwrap());
        map.insert(Region::new("chr1", 30, 40).unwrap());
        map.merge_overlaps();
        assert!(map.contains_point("chr1", 15));
        assert!(!map.contains_point("chr1", 25));
        assert!(map.overlaps("chr1", 19, 31));
        assert!(!map.overlaps("chr1", 20, 30));
    }

    #[test]
    fn reads_fixture_bed() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures/data/small.bed");
        let map = IntervalMap::from_bed(&path, Some("chrToy")).unwrap();
        assert!(map.contains_point("chrToy", 499));
        assert!(!map.contains_point("chrToy", 1000));
    }
}
