use std::{fs, path::PathBuf};

use clair3::{
    config::{LABEL_COUNT, Platform},
    tensor::labels::{labels_from_reference, labels_from_vcf_columns},
    train::data::{
        BinArchive, BinWriter, TensorBinOptions, TensorExample, merge_bins, tensor_to_bin,
    },
};

fn directory() -> PathBuf {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/verify/training");
    fs::create_dir_all(&path).unwrap();
    path
}

fn examples() -> Vec<TensorExample> {
    (0..2)
        .map(|index| TensorExample {
            position: format!("chrToy:{}:{}", index + 5, "A".repeat(33)),
            matrix: (0..33 * 18)
                .map(|value| value + index * 1000 - 20)
                .collect(),
            alt_info: format!("{}-XC 1 RA 1", index + 2),
            labels: labels_from_reference("A").unwrap(),
        })
        .collect()
}

#[test]
fn hdf5_roundtrip_and_merge_preserve_datasets() {
    let directory = directory();
    let expected = examples();
    for (index, example) in expected.iter().enumerate() {
        let path = directory.join(format!("part{index}.bin"));
        let mut writer = BinWriter::create(&path, &[33, 18], true).unwrap();
        writer.append(std::slice::from_ref(example)).unwrap();
        drop(writer);
    }
    let first = BinArchive::open(&directory.join("part0.bin")).unwrap();
    assert_eq!(first.len(), 1);
    assert_eq!(first.tensor_shape(), &[33, 18]);
    let (matrix, labels) = first.read_batch(0, 1, 18).unwrap();
    assert_eq!(matrix.size(), vec![1, 33, 18]);
    assert_eq!(labels.size(), vec![1, LABEL_COUNT as i64]);
    assert_eq!(first.read_examples(0, 1).unwrap(), expected[..1]);

    let merged = directory.join("merged.bin");
    assert_eq!(
        merge_bins(
            &[directory.join("part0.bin"), directory.join("part1.bin")],
            &merged,
            true,
        )
        .unwrap(),
        2
    );
    let merged = BinArchive::open(&merged).unwrap();
    assert_eq!(merged.len(), 2);
    assert_eq!(merged.read_examples(0, 2).unwrap(), expected);
}

#[test]
fn tensor_conversion_labels_filters_and_deduplicates() {
    let directory = directory();
    let tensor = directory.join("tensor.txt");
    let truth = directory.join("truth.txt");
    fs::write(&truth, "chrToy 5 A C 0 1\n").unwrap();
    let values = (0..33 * 18)
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    fs::write(
        &tensor,
        format!(
            "chrToy\t5\t{}\t{}\t2-XC 1 RA 1\nchrToy\t6\t{}\t{}\t2-RA 2\n",
            "A".repeat(33),
            values,
            "A".repeat(33),
            (0..33 * 18)
                .map(|value| (value + 594).to_string())
                .collect::<Vec<_>>()
                .join(" ")
        ),
    )
    .unwrap();
    let options = TensorBinOptions {
        platform: Platform::Ont,
        pileup: true,
        ..Default::default()
    };
    let output = directory.join("tensors.bin");
    assert_eq!(
        tensor_to_bin(&tensor, &truth, None, &output, &options, None, None).unwrap(),
        2
    );
    let archive = BinArchive::open(&output).unwrap();
    let examples = archive.read_examples(0, 2).unwrap();
    assert_eq!(
        examples[0].labels,
        labels_from_vcf_columns(&["chrToy", "5", "A", "C", "0", "1"]).unwrap()
    );
    assert_eq!(examples[1].labels, labels_from_reference("A").unwrap());
    assert_eq!(examples[0].matrix, (0..33 * 18).collect::<Vec<i32>>());

    let filtered = directory.join("filtered.bin");
    let options = TensorBinOptions {
        maximum_non_variant_ratio: Some(0.0),
        ..options.clone()
    };
    tensor_to_bin(&tensor, &truth, None, &filtered, &options, None, None).unwrap();
    assert_eq!(BinArchive::open(&filtered).unwrap().len(), 1);

    let duplicate = directory.join("duplicate.txt");
    let first = fs::read_to_string(&tensor)
        .unwrap()
        .lines()
        .next()
        .unwrap()
        .to_string();
    fs::write(&duplicate, format!("{first}\n{first}\n")).unwrap();
    let unique = directory.join("unique.bin");
    tensor_to_bin(&duplicate, &truth, None, &unique, &options, None, None).unwrap();
    assert_eq!(BinArchive::open(&unique).unwrap().len(), 1);
    let duplicated = directory.join("duplicated.bin");
    let options = TensorBinOptions {
        allow_duplicate_positions: true,
        maximum_non_variant_ratio: None,
        ..options
    };
    tensor_to_bin(&duplicate, &truth, None, &duplicated, &options, None, None).unwrap();
    assert_eq!(BinArchive::open(&duplicated).unwrap().len(), 2);
}

#[test]
fn high_depth_pileup_counts_survive_text_and_training_conversion() {
    use clair3::{
        io::reference::Reference,
        pipeline::{
            call::CallConfig,
            tensor::{TensorData, create_tensor_data},
        },
    };

    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let data = root.join("tests/fixtures/data");
    let mut config = CallConfig::new(
        data.join("toy.bam"),
        data.join("toy.fa"),
        PathBuf::new(),
        "chrToy".to_string(),
        Platform::Ont,
    );
    config.start = Some(2484);
    config.end = Some(2516);
    let reference = Reference::open(&data.join("toy.fa")).unwrap();
    let tensor = create_tensor_data(&config, true).unwrap();
    assert!(tensor.values.contains(&-150));
    let text = tensor.to_text(&reference).unwrap();
    let restored = TensorData::from_text(text.as_bytes(), true, Platform::Ont, false).unwrap();
    assert!(restored.values.contains(&-150));

    let directory = directory();
    let tensor_path = directory.join("high-depth.txt");
    fs::write(&tensor_path, text).unwrap();
    let truth = directory.join("high-depth.var");
    fs::write(&truth, "chrToy 2500 T A 0 1\n").unwrap();
    let output = directory.join("high-depth.bin");
    let options = TensorBinOptions {
        platform: Platform::Ont,
        pileup: true,
        maximum_non_variant_ratio: Some(0.0),
        ..Default::default()
    };
    tensor_to_bin(&tensor_path, &truth, None, &output, &options, None, None).unwrap();
    let archive = BinArchive::open(&output).unwrap();
    assert_eq!(archive.len(), 1);
    let examples = archive.read_examples(0, 1).unwrap();
    let index = restored
        .positions
        .iter()
        .position(|position| position.starts_with("chrToy:2500:"))
        .unwrap();
    let depth = restored.alt_info[index]
        .split('-')
        .next()
        .unwrap()
        .parse::<f64>()
        .unwrap();
    let expected: Vec<i32> = restored.values[index * 33 * 18..(index + 1) * 33 * 18]
        .iter()
        .map(|value| (*value as f64 / (depth / 144.0)) as i32)
        .collect();
    assert_eq!(examples[0].matrix, expected);
}
