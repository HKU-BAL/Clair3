use std::{fs, path::PathBuf};

use clair3::{
    io::reference::Reference,
    tensor::full_alignment::{FullAlignmentOptions, create_full_alignment},
};
use tch::Tensor;

fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

#[test]
fn python_illumina_tensors_match_rust() {
    let data = root().join("tests/fixtures/data");
    let reference = Reference::open(&data.join("toy.fa")).unwrap();
    let region = clair3::io::interval::Region::new("chrToy", 0, 4096).unwrap();
    let candidates: Vec<i64> = serde_json::from_str(
        &fs::read_to_string(root().join("tests/fixtures/tensors/pure_ilmn_candidates.json"))
            .unwrap(),
    )
    .unwrap();
    let options = FullAlignmentOptions {
        matrix_depth: 55,
        illumina: true,
        ..Default::default()
    };
    let batch = create_full_alignment(
        &data.join("toy.bam"),
        &reference,
        &region,
        &candidates,
        &[],
        &options,
    )
    .unwrap();
    let expected = Tensor::read_npy(root().join("tests/fixtures/tensors/pure_ilmn.npy")).unwrap();
    let shape = expected
        .size()
        .into_iter()
        .map(|value| value as usize)
        .collect::<Vec<_>>();
    assert_eq!(batch.shape, shape);
    let expected_values = Vec::<i8>::try_from(&expected.flatten(0, -1)).unwrap();
    let differences: Vec<(usize, i8, i8)> = batch
        .values
        .iter()
        .zip(&expected_values)
        .enumerate()
        .filter(|(_, (a, e))| a != e)
        .take(20)
        .map(|(index, (a, e))| (index, *a, *e))
        .collect();
    assert!(
        differences.is_empty(),
        "first differences (index, actual, expected): {differences:?}"
    );
    let expected_info =
        fs::read_to_string(root().join("tests/fixtures/tensors/pure_ilmn.info")).unwrap();
    let mut positions = vec![];
    let mut alt_info = vec![];
    for line in expected_info.lines() {
        let (position, alt) = line.split_once('\t').unwrap();
        positions.push(position.to_string());
        alt_info.push(alt.to_string());
    }
    assert_eq!(batch.positions, positions);
    assert_eq!(batch.alt_info, alt_info);
}

#[test]
fn illumina_n_bases_use_legacy_strand_encoding() {
    use rust_htslib::bam::{
        self, Header, Writer,
        record::{Cigar, CigarString, Record},
    };

    let directory = tempfile::tempdir().unwrap();
    let bam_path = directory.path().join("ambiguous.bam");
    let mut header = Header::new();
    let mut contig = bam::header::HeaderRecord::new(b"SQ");
    contig.push_tag(b"SN", "chrToy").push_tag(b"LN", 4096);
    header.push_record(&contig);
    let mut writer = Writer::from_path(&bam_path, &header, bam::Format::Bam).unwrap();
    for (name, ambiguous) in [("ambiguous", true), ("forward", false)] {
        let mut sequence = vec![b'A'; 40];
        if ambiguous {
            sequence[19] = b'N';
        }
        let mut record = Record::new();
        record.set(
            name.as_bytes(),
            Some(&CigarString(vec![Cigar::Match(40)])),
            &sequence,
            &[40; 40],
        );
        record.set_tid(0);
        record.set_pos(480);
        record.set_mapq(60);
        record.set_flags(0);
        record.set_mtid(-1);
        record.set_mpos(-1);
        record.set_insert_size(0);
        writer.write(&record).unwrap();
    }
    drop(writer);
    bam::index::build(&bam_path, None, bam::index::Type::Bai, 0).unwrap();
    let data = root().join("tests/fixtures/data");
    let reference = Reference::open(&data.join("toy.fa")).unwrap();
    let batch = create_full_alignment(
        &bam_path,
        &reference,
        &clair3::io::interval::Region::new("chrToy", 0, 4096).unwrap(),
        &[499],
        &[],
        &FullAlignmentOptions {
            matrix_depth: 55,
            illumina: true,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(batch.values[((26 * 33 + 15) * 8) + 2], 100);
    assert_eq!(batch.values[((26 * 33 + 16) * 8) + 2], 50);
    assert_eq!(batch.values[((27 * 33 + 16) * 8) + 2], 100);
}
