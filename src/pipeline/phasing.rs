use std::{
    fs,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    process::Command,
};

use anyhow::{Result, bail};
use rust_htslib::bam::{IndexedReader, Read as BamRead, Record};

const OPTIMIZED_PHASE_CORE: &[u8] = include_bytes!(env!("CLAIR3_PHASE_CORE"));

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PhaseRegion {
    pub start: i64,
    pub end: i64,
    pub first_variant: usize,
    pub end_variant: usize,
}

fn merge_span(positions: &[i64], groups: &mut Vec<PhaseRegion>, start: i64, end: i64) {
    let first = positions.partition_point(|position| *position < start);
    let last = positions.partition_point(|position| *position < end);
    if first == last {
        return;
    }
    let first_group =
        groups.partition_point(|group| positions[group.end_variant - 1] < positions[first]);
    let last_group = groups.partition_point(|group| group.start <= positions[last - 1]);
    let mut first_variant = first;
    let mut end_variant = last;
    let mut region_start = positions[first];
    let mut region_end = positions[last - 1] + 1;
    for group in &groups[first_group..last_group] {
        first_variant = first_variant.min(group.first_variant);
        end_variant = end_variant.max(group.end_variant);
        region_start = region_start.min(group.start);
        region_end = region_end.max(group.end);
    }
    groups.splice(
        first_group..last_group,
        [PhaseRegion {
            start: region_start,
            end: region_end,
            first_variant,
            end_variant,
        }],
    );
}

pub fn scan_phase_regions(bam: &Path, contig: &str, positions: &[i64]) -> Result<Vec<PhaseRegion>> {
    if positions.is_empty() || !positions.is_sorted() {
        bail!("phasing positions must be nonempty and sorted");
    }
    let mut reader = IndexedReader::from_path(bam)?;
    reader.fetch(contig.as_bytes())?;
    let mut groups = Vec::new();
    let mut record = Record::new();
    while let Some(result) = reader.read(&mut record) {
        result?;
        let flags = record.flags();
        if flags & 0x904 != 0 || record.mapq() < 20 || record.cigar().is_empty() {
            continue;
        }
        merge_span(
            positions,
            &mut groups,
            record.pos(),
            record.cigar().end_pos(),
        );
    }
    Ok(groups)
}

pub fn capture_alignment_spans(
    bam: &Path,
    reference: &Path,
    contig: &str,
    output: &Path,
) -> Result<()> {
    let mut reader = IndexedReader::from_path(bam)?;
    reader.set_reference(reference)?;
    reader.fetch(contig.as_bytes())?;
    let mut writer = BufWriter::new(fs::File::create(output)?);
    let mut record = Record::new();
    while let Some(result) = reader.read(&mut record) {
        result?;
        if record.flags() & 0x904 != 0 || record.mapq() < 20 || record.cigar().is_empty() {
            continue;
        }
        writeln!(writer, "{}\t{}", record.pos(), record.cigar().end_pos())?;
    }
    writer.flush()?;
    Ok(())
}

pub fn phase_regions_from_files(positions: &[i64], paths: &[PathBuf]) -> Result<Vec<PhaseRegion>> {
    if positions.is_empty() || !positions.is_sorted() {
        bail!("phasing positions must be nonempty and sorted");
    }
    let mut groups = Vec::new();
    for path in paths {
        for line in BufReader::new(fs::File::open(path)?).lines() {
            let line = line?;
            let Some((start, end)) = line.split_once('\t') else {
                bail!("invalid alignment span in {}", path.display());
            };
            merge_span(positions, &mut groups, start.parse()?, end.parse()?);
        }
    }
    Ok(groups)
}

fn whatshap_python(executable: &str) -> String {
    let path = Path::new(executable);
    let resolved = if path.components().count() == 1 {
        std::env::var_os("PATH")
            .into_iter()
            .flat_map(|paths| std::env::split_paths(&paths).collect::<Vec<_>>())
            .map(|directory| directory.join(path))
            .find(|candidate| candidate.is_file())
    } else {
        Some(path.to_path_buf())
    };
    let Some(header) = resolved.and_then(|path| fs::read(path).ok()) else {
        return "python3".to_string();
    };
    let first_line = String::from_utf8_lossy(&header);
    let Some(shebang) = first_line
        .lines()
        .next()
        .filter(|line| line.starts_with("#!"))
    else {
        return "python3".to_string();
    };
    shebang[2..]
        .split_whitespace()
        .find(|part| part.contains("python"))
        .unwrap_or("python3")
        .to_string()
}

pub fn parallel_phase_command(whatshap: &str, threads: usize, regions: Option<&Path>) -> Command {
    let mut command = Command::new(whatshap_python(whatshap));
    command
        .arg("-c")
        .arg(include_str!("parallel_phase.py"))
        .env("CLAIR3_PHASE_THREADS", threads.max(1).to_string());
    if let Some(regions) = regions {
        command.env("CLAIR3_PHASE_REGIONS", regions);
        if !OPTIMIZED_PHASE_CORE.is_empty()
            && let Some(directory) = regions.parent()
        {
            let core = directory.join("clair3_whatshap_core.so");
            if fs::write(&core, OPTIMIZED_PHASE_CORE).is_ok() {
                command.env("CLAIR3_FAST_CORE", core);
            }
        }
    }
    command
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn span_components_keep_all_connected_variants_together() {
        let positions = vec![100, 200, 300, 400];
        let mut groups = Vec::new();
        merge_span(&positions, &mut groups, 100, 201);
        merge_span(&positions, &mut groups, 300, 401);
        merge_span(&positions, &mut groups, 150, 360);
        assert_eq!(
            groups,
            vec![PhaseRegion {
                start: 100,
                end: 401,
                first_variant: 0,
                end_variant: 4,
            }]
        );
    }

    #[test]
    fn captured_spans_match_direct_alignment_scan() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR"));
        let bam = root.join("tests/fixtures/data/toy.bam");
        let reference = root.join("tests/fixtures/data/toy.fa");
        let directory = tempfile::tempdir().unwrap();
        let path = directory.path().join("toy.spans");
        capture_alignment_spans(&bam, &reference, "chrToy", &path).unwrap();
        let positions: Vec<_> = (0..4096).step_by(37).collect();
        let expected = scan_phase_regions(&bam, "chrToy", &positions).unwrap();
        let actual = phase_regions_from_files(&positions, &[path]).unwrap();
        assert_eq!(actual, expected);
    }

    #[test]
    fn disconnected_spans_stay_independent() {
        let positions = vec![100, 200, 300, 400];
        let mut groups = Vec::new();
        merge_span(&positions, &mut groups, 200, 201);
        merge_span(&positions, &mut groups, 400, 401);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].first_variant, 1);
        assert_eq!(groups[1].first_variant, 3);
    }
}
