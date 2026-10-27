use anyhow::{Result, bail};

use crate::config::LABEL_COUNT;

pub const GT21_LABELS: [&str; 21] = [
    "AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT", "DelDel", "ADel", "CDel", "GDel",
    "TDel", "InsIns", "AIns", "CIns", "GIns", "TIns", "InsDel",
];
pub const VARIANT_LENGTH_COUNT: usize = 33;
pub const VARIANT_LENGTH_OFFSET: isize = 16;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(usize)]
pub enum Genotype {
    HomoReference = 0,
    HomoVariant = 1,
    HeteroVariant = 2,
    HeteroVariantMulti = 3,
}

impl Genotype {
    pub fn from_alleles(first: i32, second: i32) -> Self {
        match (first, second) {
            (0, 0) => Self::HomoReference,
            (a, b) if a == b => Self::HomoVariant,
            (a, b) if a != 0 && b != 0 => Self::HeteroVariantMulti,
            _ => Self::HeteroVariant,
        }
    }

    pub const fn for_task(self) -> Self {
        match self {
            Self::HeteroVariantMulti => Self::HeteroVariant,
            other => other,
        }
    }

    pub const fn as_vcf(self) -> &'static str {
        match self {
            Self::HomoReference => "0/0",
            Self::HomoVariant => "1/1",
            Self::HeteroVariant => "0/1",
            Self::HeteroVariantMulti => "1/2",
        }
    }
}

pub fn gt21_index_from_label(label: &str) -> Option<usize> {
    GT21_LABELS.iter().position(|candidate| *candidate == label)
}

pub fn partial_label(reference: &str, alternate: &str) -> Result<String> {
    if reference.len() > alternate.len() {
        return Ok("Del".to_string());
    }
    if reference.len() < alternate.len() {
        return Ok("Ins".to_string());
    }
    let Some(base) = alternate.chars().next() else {
        bail!("alternate allele cannot be empty")
    };
    Ok(base.to_string())
}

pub fn mix_partial_labels(first: &str, second: &str) -> String {
    if first.len() == 1 && second.len() == 1 {
        return if first <= second {
            format!("{first}{second}")
        } else {
            format!("{second}{first}")
        };
    }
    if first.len() > 1 && second.len() == 1 {
        return format!("{second}{first}");
    }
    if second.len() > 1 && first.len() == 1 {
        return format!("{first}{second}");
    }
    if first == second {
        return format!("{first}{second}");
    }
    "InsDel".to_string()
}

pub fn gt21_from_alleles(reference: &str, alternates: &[String]) -> Result<usize> {
    if alternates.len() < 2 {
        bail!("at least two allele labels are required")
    }
    let first = partial_label(reference, &alternates[0])?;
    let second = partial_label(reference, &alternates[1])?;
    let mixed = mix_partial_labels(&first, &second);
    let Some(index) = gt21_index_from_label(&mixed) else {
        bail!("unsupported GT21 label '{mixed}'")
    };
    Ok(index)
}

pub fn gt21_from_variant(
    reference: &str,
    alternate: &str,
    genotype_1: i32,
    genotype_2: i32,
) -> Result<usize> {
    let mut alternates: Vec<String> = alternate.split(',').map(str::to_string).collect();
    if alternates.len() == 1 {
        let first = if genotype_1 == 0 || genotype_2 == 0 {
            reference
        } else {
            &alternates[0]
        };
        alternates.insert(0, first.to_string());
    }
    gt21_from_alleles(reference, &alternates)
}

fn clamp_variant_length(length: isize) -> usize {
    length
        .clamp(-VARIANT_LENGTH_OFFSET, VARIANT_LENGTH_OFFSET)
        .saturating_add(VARIANT_LENGTH_OFFSET) as usize
}

pub fn labels_from_reference(reference_base: &str) -> Result<Vec<i8>> {
    let mut labels = vec![0; LABEL_COUNT];
    let gt21 = gt21_index_from_label(&format!("{reference_base}{reference_base}"));
    let Some(gt21) = gt21 else {
        bail!("unsupported reference base '{reference_base}'")
    };
    labels[gt21] = 1;
    labels[21] = 1;
    labels[24 + VARIANT_LENGTH_OFFSET as usize] = 1;
    labels[24 + VARIANT_LENGTH_COUNT + VARIANT_LENGTH_OFFSET as usize] = 1;
    Ok(labels)
}

pub fn labels_from_vcf_columns(columns: &[&str]) -> Result<Vec<i8>> {
    if columns.len() < 6 {
        bail!("expected contig, position, REF, ALT, and two genotype alleles")
    }
    let reference = columns[2];
    let alternate = columns[3];
    let genotype_1: i32 = columns[4].parse()?;
    let genotype_2: i32 = columns[5].parse()?;
    let mut alternates: Vec<String> = alternate.split(',').map(str::to_string).collect();
    if alternates.len() == 1 {
        let first = if genotype_1 == 0 || genotype_2 == 0 {
            reference
        } else {
            &alternates[0]
        };
        alternates.insert(0, first.to_string());
    }

    let mut labels = vec![0; LABEL_COUNT];
    labels[gt21_from_alleles(reference, &alternates)?] = 1;
    let genotype = Genotype::from_alleles(genotype_1, genotype_2).for_task();
    labels[21 + genotype as usize] = 1;

    let reference_length = reference.len() as isize;
    let mut lengths: Vec<usize> = alternates
        .iter()
        .map(|alt| clamp_variant_length(alt.len() as isize - reference_length))
        .collect();
    lengths.sort_unstable();
    labels[24 + lengths[0]] = 1;
    labels[24 + VARIANT_LENGTH_COUNT + lengths[1]] = 1;
    Ok(labels)
}

#[cfg(test)]
mod tests {
    use serde_json::Value;

    use super::*;

    fn strings(value: &Value) -> Vec<&str> {
        value
            .as_array()
            .unwrap()
            .iter()
            .map(|item| item.as_str().unwrap())
            .collect()
    }

    #[test]
    fn matches_legacy_golden_labels() {
        let golden: Value =
            serde_json::from_str(include_str!("../../tests/fixtures/labels.json")).unwrap();
        assert_eq!(GT21_LABELS.to_vec(), strings(&golden["gt21_labels"]));

        let reference = golden["reference_labels"].as_object().unwrap();
        for (base, expected) in reference {
            let expected: Vec<i8> = expected
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap() as i8)
                .collect();
            assert_eq!(
                labels_from_reference(base).unwrap(),
                expected,
                "reference {base}"
            );
        }

        for case in golden["variant_cases"].as_array().unwrap() {
            let columns = strings(&case["columns"]);
            let expected: Vec<i8> = case["label"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_i64().unwrap() as i8)
                .collect();
            assert_eq!(
                labels_from_vcf_columns(&columns).unwrap(),
                expected,
                "columns {columns:?}"
            );
            let gt21 = gt21_from_variant(
                columns[2],
                columns[3],
                columns[4].parse().unwrap(),
                columns[5].parse().unwrap(),
            )
            .unwrap();
            assert_eq!(gt21, case["gt21"].as_u64().unwrap() as usize);
        }

        for case in golden["partial_cases"].as_array().unwrap() {
            let reference = case["ref"].as_str().unwrap();
            let alt1 = case["alt1"].as_str().unwrap();
            let alt2 = case["alt2"].as_str().unwrap();
            let first = partial_label(reference, alt1).unwrap();
            let second = partial_label(reference, alt2).unwrap();
            assert_eq!(first, case["partial1"].as_str().unwrap());
            assert_eq!(second, case["partial2"].as_str().unwrap());
            assert_eq!(
                mix_partial_labels(&first, &second),
                case["mixed"].as_str().unwrap()
            );
        }
    }

    #[test]
    fn genotype_matches_legacy_encoding() {
        assert_eq!(Genotype::from_alleles(0, 0), Genotype::HomoReference);
        assert_eq!(Genotype::from_alleles(1, 1), Genotype::HomoVariant);
        assert_eq!(Genotype::from_alleles(0, 1), Genotype::HeteroVariant);
        assert_eq!(Genotype::from_alleles(1, 2), Genotype::HeteroVariantMulti);
        assert_eq!(
            Genotype::HeteroVariantMulti.for_task(),
            Genotype::HeteroVariant
        );
    }
}
