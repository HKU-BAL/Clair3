use std::cmp::Ordering;

use anyhow::{Result, bail};

use crate::io::variant::VcfRecord;

const GT21_LABELS: [&str; 21] = [
    "AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "GT", "TT", "DelDel", "ADel", "CDel", "GDel",
    "TDel", "InsIns", "AIns", "CIns", "GIns", "TIns", "InsDel",
];
const HOMO_SNP: [usize; 4] = [0, 4, 7, 9];
const HETERO_SNP: [usize; 6] = [1, 2, 3, 5, 6, 8];
const LENGTH_OFFSET: usize = 16;

#[derive(Debug, Clone)]
pub struct Prediction {
    pub gt21: Vec<f32>,
    pub genotype: Vec<f32>,
    pub length1: Vec<f32>,
    pub length2: Vec<f32>,
}

#[derive(Debug, Clone)]
pub struct DecodeOptions {
    pub show_reference: bool,
    pub haploid_precise: bool,
    pub haploid_sensitive: bool,
    pub quality_for_pass: Option<f64>,
    pub add_indel_length: bool,
    pub gvcf: bool,
    pub pileup: bool,
    pub enable_long_indel: bool,
    pub maximum_variant_length: usize,
    pub keep_iupac: bool,
}

impl Default for DecodeOptions {
    fn default() -> Self {
        Self {
            show_reference: false,
            haploid_precise: false,
            haploid_sensitive: false,
            quality_for_pass: Some(2.0),
            add_indel_length: false,
            gvcf: false,
            pileup: false,
            enable_long_indel: false,
            maximum_variant_length: 50,
            keep_iupac: false,
        }
    }
}

#[derive(Debug)]
struct AltInfo {
    depth: i64,
    entries: Vec<(String, i64)>,
}

impl AltInfo {
    fn parse(text: &str) -> Result<Self> {
        let mut parts = text.trim_end().split('-');
        let depth = parts.next().unwrap_or("0").parse()?;
        let entries_text = parts.next().unwrap_or("");
        let fields: Vec<&str> = entries_text.split(' ').collect();
        let mut entries: Vec<(String, i64)> = Vec::new();
        for pair in fields.chunks(2) {
            if pair.len() != 2 {
                continue;
            }
            let count = pair[1].parse()?;
            if let Some(entry) = entries.iter_mut().find(|entry| entry.0 == pair[0]) {
                entry.1 = count;
            } else {
                entries.push((pair[0].to_string(), count));
            }
        }
        Ok(Self { depth, entries })
    }

    fn of_type(&self, kind: char) -> Vec<(String, i64)> {
        self.entries
            .iter()
            .filter(|(key, _)| key.starts_with(kind))
            .map(|(key, count)| (key[kind.len_utf8()..].to_string(), *count))
            .collect()
    }

    fn reference_count(&self) -> i64 {
        self.entries
            .iter()
            .filter(|(key, _)| key.starts_with('R'))
            .fold(0, |_, (_, count)| (*count).max(0))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CallKind {
    Reference,
    HomoSnp,
    HeteroSnp,
    HomoInsertion,
    HeteroAcgtInsertion,
    HeteroInsertionInsertion,
    HomoDeletion,
    HeteroAcgtDeletion,
    HeteroDeletionDeletion,
    InsertionDeletion,
}

struct Outcome {
    kind: CallKind,
    reference: String,
    alternate: String,
    probability: f32,
    hetero: bool,
}

fn base_from_iupac(base: char) -> char {
    match base.to_ascii_uppercase() {
        'A' | 'R' | 'W' | 'M' | 'D' | 'H' | 'V' | 'N' => 'A',
        'C' | 'Y' | 'S' | 'B' => 'C',
        'G' | 'K' => 'G',
        'T' | 'U' => 'T',
        _ => 'N',
    }
}

fn convert_iupac_to_n(text: &str) -> String {
    if text == "." {
        return text.to_string();
    }
    text.chars()
        .map(|base| {
            if "ACGTN,.".contains(base.to_ascii_uppercase()) {
                base
            } else {
                'N'
            }
        })
        .collect()
}

fn sorted_alt_bases(alt: &AltInfo, alternate: Option<&str>) -> (Vec<String>, Option<String>) {
    let mut bases: Vec<(String, i64)> = alt
        .entries
        .iter()
        .filter(|(key, _)| key.starts_with('X') && key.len() > 1)
        .map(|(key, count)| (key[1..].to_string(), *count))
        .collect();
    bases.sort_by_key(|entry| std::cmp::Reverse(entry.1));
    let Some(first) = bases.first() else {
        return (bases.into_iter().map(|(base, _)| base).collect(), None);
    };
    let requested_count = alternate.and_then(|base| {
        bases
            .iter()
            .find(|(candidate, _)| candidate == base)
            .map(|(_, count)| *count)
    });
    let selected = if requested_count.is_none_or(|count| first.1 - count >= 9) {
        first.0.clone()
    } else {
        alternate.unwrap().to_string()
    };
    (
        bases.into_iter().map(|(base, _)| base).collect(),
        Some(selected),
    )
}

fn select_indel(
    entries: &[(String, i64)],
    proposed_length: Option<usize>,
    insertion: bool,
    ignore: &str,
    maximum_length: usize,
    multiple: bool,
) -> Vec<String> {
    let mut proposed = Vec::new();
    let mut normal = Vec::new();
    let target = proposed_length.map(|length| length + usize::from(insertion));
    for (bases, count) in entries {
        if bases == ignore {
            continue;
        }
        if target.is_some_and(|target| bases.len() == target) {
            proposed.push((bases.clone(), *count));
        } else if !bases.is_empty() && bases.len() <= maximum_length {
            normal.push((bases.clone(), *count));
        }
    }
    proposed.sort_by_key(|entry| std::cmp::Reverse(entry.1));
    if target.is_some() && !proposed.is_empty() {
        return vec![proposed[0].0.clone()];
    }
    if multiple {
        normal.sort_by_key(|entry| entry.1);
        normal.reverse();
    } else {
        normal.sort_by_key(|entry| std::cmp::Reverse(entry.1));
        return normal.into_iter().take(1).map(|(base, _)| base).collect();
    }
    if !insertion && normal.len() <= 1 {
        return Vec::new();
    }
    if !insertion && normal[0].0.len() <= normal[1].0.len() {
        normal.swap(0, 1);
    }
    normal.into_iter().take(2).map(|(base, _)| base).collect()
}

fn one_indel(
    entries: &[(String, i64)],
    proposed_length: Option<usize>,
    insertion: bool,
    ignore: &str,
    maximum_length: usize,
) -> Option<String> {
    select_indel(
        entries,
        proposed_length,
        insertion,
        ignore,
        maximum_length,
        false,
    )
    .into_iter()
    .next()
}

fn reference_outcome(reference: String, probability: f32) -> Outcome {
    Outcome {
        kind: CallKind::Reference,
        reference: reference.clone(),
        alternate: reference,
        probability,
        hetero: false,
    }
}

#[allow(unused_parens)]
fn decode_outcome(
    reference_sequence: &str,
    center: usize,
    prediction: &Prediction,
    options: &DecodeOptions,
    alt: &AltInfo,
) -> Result<Outcome> {
    let reference_raw = reference_sequence[center..].chars().next().unwrap();
    let reference_acgt = base_from_iupac(reference_raw);
    let reference_gt21 = match reference_acgt {
        'A' => 0,
        'C' => 4,
        'G' => 7,
        'T' => 9,
        _ => bail!("unsupported reference base {reference_acgt}"),
    };
    let homo_reference = prediction.genotype[0];
    let homo_variant = prediction.genotype[1];
    let hetero_variant = prediction.genotype[2];
    let length0 = if options.add_indel_length {
        prediction.length1[LENGTH_OFFSET] * prediction.length2[LENGTH_OFFSET]
    } else {
        1.0
    };
    let reference_probability = length0 * homo_reference * prediction.gt21[reference_gt21];
    if (!options.add_indel_length
        && homo_reference >= 0.5
        && prediction.gt21[reference_gt21] >= 0.5)
        || (options.add_indel_length
            && prediction.length1[LENGTH_OFFSET] >= 0.5
            && prediction.length2[LENGTH_OFFSET] >= 0.5
            && homo_reference >= 0.5
            && prediction.gt21[reference_gt21] >= 0.5)
    {
        return Ok(reference_outcome(
            reference_acgt.to_string(),
            reference_probability,
        ));
    }
    let mut homo_snp: Vec<f32> = HOMO_SNP
        .iter()
        .map(|index| length0 * homo_variant * prediction.gt21[*index])
        .collect();
    let mut hetero_snp: Vec<f32> = HETERO_SNP
        .iter()
        .map(|index| length0 * hetero_variant * prediction.gt21[*index])
        .collect();
    let mut homo_insertion = if options.add_indel_length {
        (1..=LENGTH_OFFSET)
            .map(|length| {
                prediction.length1[LENGTH_OFFSET + length]
                    * prediction.length2[LENGTH_OFFSET + length]
                    * homo_variant
                    * prediction.gt21[15]
            })
            .collect()
    } else {
        vec![homo_variant * prediction.gt21[15]]
    };
    let homo_insertion_lengths: Vec<usize> = if options.add_indel_length {
        (1..=LENGTH_OFFSET).collect()
    } else {
        Vec::new()
    };
    let mut hetero_ins_ins = if options.add_indel_length {
        let mut result = Vec::new();
        for first in 1..=LENGTH_OFFSET {
            for second in first..=LENGTH_OFFSET {
                result.push((
                    (first, second),
                    prediction.length1[LENGTH_OFFSET + first]
                        * prediction.length2[LENGTH_OFFSET + second]
                        * hetero_variant
                        * prediction.gt21[15],
                ));
            }
        }
        result
    } else {
        vec![((0, 0), hetero_variant * prediction.gt21[15])]
    };
    let mut hetero_acgt_insertion: Vec<(usize, usize, f32)> = if options.add_indel_length {
        (1..=LENGTH_OFFSET)
            .flat_map(|length| {
                (0..4).map(move |base| {
                    (
                        base,
                        length,
                        prediction.length1[LENGTH_OFFSET]
                            * prediction.length2[LENGTH_OFFSET + length]
                            * prediction.gt21[16 + base]
                            * hetero_variant,
                    )
                })
            })
            .collect()
    } else {
        (0..4)
            .map(|base| (base, 0, prediction.gt21[16 + base] * hetero_variant))
            .collect()
    };
    let mut homo_deletion = if options.add_indel_length {
        (1..=LENGTH_OFFSET)
            .map(|length| {
                prediction.length1[LENGTH_OFFSET - length]
                    * prediction.length2[LENGTH_OFFSET - length]
                    * homo_variant
                    * prediction.gt21[10]
            })
            .collect()
    } else {
        vec![homo_variant * prediction.gt21[10]]
    };
    let homo_deletion_lengths: Vec<usize> = if options.add_indel_length {
        (1..=LENGTH_OFFSET).collect()
    } else {
        Vec::new()
    };
    let mut hetero_del_del = if options.add_indel_length {
        let mut result = Vec::new();
        for first in 1..=LENGTH_OFFSET {
            for second in 1..=LENGTH_OFFSET {
                if first == second && first != LENGTH_OFFSET && second != LENGTH_OFFSET {
                    continue;
                }
                result.push((
                    (first.min(second), first.max(second)),
                    prediction.length1[LENGTH_OFFSET - first]
                        * prediction.length2[LENGTH_OFFSET - second]
                        * hetero_variant
                        * prediction.gt21[10],
                ));
            }
        }
        result
    } else {
        vec![((0, 0), hetero_variant * prediction.gt21[10])]
    };
    let mut hetero_acgt_deletion: Vec<(usize, usize, f32)> = if options.add_indel_length {
        (1..=LENGTH_OFFSET)
            .flat_map(|length| {
                (0..4).map(move |base| {
                    (
                        base,
                        length,
                        prediction.length1[LENGTH_OFFSET - length]
                            * prediction.length2[LENGTH_OFFSET]
                            * prediction.gt21[11 + base]
                            * hetero_variant,
                    )
                })
            })
            .collect()
    } else {
        (0..4)
            .map(|base| (base, 0, prediction.gt21[11 + base] * hetero_variant))
            .collect()
    };
    let mut insertion_deletion = if options.add_indel_length {
        let mut result = Vec::new();
        for deletion in 1..=LENGTH_OFFSET {
            for insertion in 1..=LENGTH_OFFSET {
                result.push((
                    (deletion, insertion),
                    prediction.length1[LENGTH_OFFSET - deletion]
                        * prediction.length2[LENGTH_OFFSET + insertion]
                        * hetero_variant
                        * prediction.gt21[20],
                ));
            }
        }
        result
    } else {
        vec![((0, 0), hetero_variant * prediction.gt21[20])]
    };
    let insertions = alt.of_type('I');
    let deletions = alt.of_type('D');
    let mut outcome = loop {
        let maximum = homo_snp
            .iter()
            .chain(&hetero_snp)
            .chain(&homo_insertion)
            .chain(&homo_deletion)
            .chain(hetero_ins_ins.iter().map(|(_, probability)| probability))
            .chain(
                hetero_acgt_insertion
                    .iter()
                    .map(|(_, _, probability)| probability),
            )
            .chain(hetero_del_del.iter().map(|(_, probability)| probability))
            .chain(
                hetero_acgt_deletion
                    .iter()
                    .map(|(_, _, probability)| probability),
            )
            .chain(
                insertion_deletion
                    .iter()
                    .map(|(_, probability)| probability),
            )
            .fold(reference_probability, |maximum, probability| {
                maximum.max(*probability)
            });
        if maximum == reference_probability {
            return Ok(reference_outcome(reference_acgt.to_string(), maximum));
        }
        let index =
            |probabilities: &[f32]| probabilities.iter().position(|p| *p == maximum).unwrap();
        if homo_snp.contains(&maximum) {
            let index = index(&homo_snp);
            let labels = ["AA", "CC", "GG", "TT"];
            let pair = labels[index];
            let proposed = if pair.chars().next().unwrap() == reference_raw {
                pair.chars().nth(1).unwrap().to_string()
            } else {
                pair.chars().next().unwrap().to_string()
            };
            let (_, Some(alternate)) = sorted_alt_bases(alt, Some(&proposed)) else {
                homo_snp[index] = 0.0;
                continue;
            };
            if alternate.chars().next().unwrap() == reference_raw {
                break (Outcome {
                    kind: CallKind::HomoSnp,
                    reference: reference_raw.to_string(),
                    alternate,
                    probability: maximum,
                    hetero: false,
                });
            }
            break (Outcome {
                kind: CallKind::HomoSnp,
                reference: reference_raw.to_string(),
                alternate,
                probability: maximum,
                hetero: false,
            });
        }
        if hetero_snp.contains(&maximum) {
            let index = index(&hetero_snp);
            let pair = ["AC", "AG", "AT", "CG", "CT", "GT"][index];
            let first = pair.chars().next().unwrap();
            let second = pair.chars().nth(1).unwrap();
            let selected = if first != reference_raw && second != reference_raw {
                let (bases, _) = sorted_alt_bases(alt, None);
                if bases.len() < 2 {
                    hetero_snp[index] = 0.0;
                    continue;
                }
                Some(bases[..2].join(","))
            } else {
                let proposed = if first == reference_raw {
                    second
                } else {
                    first
                };
                let (_, selected) = sorted_alt_bases(alt, Some(&proposed.to_string()));
                match selected {
                    None => {
                        hetero_snp[index] = 0.0;
                        continue;
                    }
                    Some(base) if base.chars().next().unwrap() == reference_raw => {
                        break (Outcome {
                            kind: CallKind::HeteroSnp,
                            reference: reference_raw.to_string(),
                            alternate: base,
                            probability: maximum,
                            hetero: false,
                        });
                    }
                    Some(base) => Some(base),
                }
            };
            break (Outcome {
                kind: CallKind::HeteroSnp,
                reference: reference_raw.to_string(),
                alternate: selected.unwrap(),
                probability: maximum,
                hetero: false,
            });
        }
        if homo_insertion.contains(&maximum) {
            let index = index(&homo_insertion);
            let proposed = options
                .add_indel_length
                .then(|| homo_insertion_lengths[index])
                .filter(|length| *length < LENGTH_OFFSET);
            let Some(alternate) = one_indel(
                &insertions,
                proposed,
                true,
                "",
                options.maximum_variant_length,
            ) else {
                homo_insertion[index] = 0.0;
                continue;
            };
            break (Outcome {
                kind: CallKind::HomoInsertion,
                reference: reference_raw.to_string(),
                alternate,
                probability: maximum,
                hetero: false,
            });
        }
        if hetero_acgt_insertion
            .iter()
            .any(|(_, _, probability)| *probability == maximum)
        {
            let index = hetero_acgt_insertion
                .iter()
                .position(|(_, _, probability)| *probability == maximum)
                .unwrap();
            let (base, length, _) = hetero_acgt_insertion[index];
            let proposed = options
                .add_indel_length
                .then_some(length)
                .filter(|length| *length < LENGTH_OFFSET);
            let Some(mut alternate) = one_indel(
                &insertions,
                proposed,
                true,
                "",
                options.maximum_variant_length,
            ) else {
                hetero_acgt_insertion[index].2 = 0.0;
                continue;
            };
            let hetero_base = if options.add_indel_length {
                b"ACGT"[base] as char
            } else {
                b"ACGT"[index] as char
            };
            if hetero_base != reference_raw {
                let (bases, _) = sorted_alt_bases(alt, None);
                if bases.is_empty() {
                    break (Outcome {
                        kind: CallKind::HeteroAcgtInsertion,
                        reference: reference_raw.to_string(),
                        alternate,
                        probability: maximum,
                        hetero: false,
                    });
                }
                alternate = format!("{},{}", bases[0], alternate);
            }
            break (Outcome {
                kind: CallKind::HeteroAcgtInsertion,
                reference: reference_raw.to_string(),
                alternate,
                probability: maximum,
                hetero: false,
            });
        }
        if hetero_ins_ins
            .iter()
            .any(|(_, probability)| *probability == maximum)
        {
            let index = hetero_ins_ins
                .iter()
                .position(|(_, probability)| *probability == maximum)
                .unwrap();
            let (first_length, second_length) = hetero_ins_ins[index].0;
            let first = options
                .add_indel_length
                .then_some(first_length)
                .filter(|length| *length < LENGTH_OFFSET);
            let second = options
                .add_indel_length
                .then_some(second_length)
                .filter(|length| *length < LENGTH_OFFSET);
            let mut bases = Vec::new();
            if options.add_indel_length
                && let Some(first) =
                    one_indel(&insertions, first, true, "", options.maximum_variant_length)
                && let Some(second) = one_indel(
                    &insertions,
                    second,
                    true,
                    &first,
                    options.maximum_variant_length,
                )
            {
                bases = vec![first, second];
            }
            if bases.len() < 2 {
                bases = select_indel(
                    &insertions,
                    None,
                    true,
                    "",
                    options.maximum_variant_length,
                    true,
                );
            }
            if bases.len() < 2 {
                hetero_ins_ins[index].1 = 0.0;
                continue;
            }
            if bases[0] == bases[1] {
                break (Outcome {
                    kind: CallKind::HeteroInsertionInsertion,
                    reference: reference_raw.to_string(),
                    alternate: bases[0].clone(),
                    probability: maximum,
                    hetero: false,
                });
            }
            break (Outcome {
                kind: CallKind::HeteroInsertionInsertion,
                reference: reference_raw.to_string(),
                alternate: format!("{},{}", bases[1], bases[0]),
                probability: maximum,
                hetero: false,
            });
        }
        if homo_deletion.contains(&maximum) {
            let index = index(&homo_deletion);
            let proposed = options
                .add_indel_length
                .then(|| homo_deletion_lengths[index])
                .filter(|length| *length < LENGTH_OFFSET);
            let Some(deleted) = one_indel(
                &deletions,
                proposed,
                false,
                "",
                options.maximum_variant_length,
            ) else {
                homo_deletion[index] = 0.0;
                continue;
            };
            let reference = format!("{reference_raw}{deleted}");
            let alternate = reference.chars().next().unwrap().to_string();
            break (Outcome {
                kind: CallKind::HomoDeletion,
                reference,
                alternate,
                probability: maximum,
                hetero: false,
            });
        }
        if hetero_acgt_deletion
            .iter()
            .any(|(_, _, probability)| *probability == maximum)
        {
            let index = hetero_acgt_deletion
                .iter()
                .position(|(_, _, probability)| *probability == maximum)
                .unwrap();
            let (base, length, _) = hetero_acgt_deletion[index];
            let proposed = options
                .add_indel_length
                .then_some(length)
                .filter(|length| *length < LENGTH_OFFSET);
            let Some(deleted) = one_indel(
                &deletions,
                proposed,
                false,
                "",
                options.maximum_variant_length,
            ) else {
                hetero_acgt_deletion[index].2 = 0.0;
                continue;
            };
            let reference = format!("{reference_raw}{deleted}");
            let first = reference.chars().next().unwrap();
            let mut alternate = first.to_string();
            let hetero_base = if options.add_indel_length {
                b"ACGT"[base] as char
            } else {
                b"ACGT"[index] as char
            };
            if hetero_base != first {
                let second = format!("{}{}", hetero_base, &reference[1..]);
                alternate = format!("{},{}", first, second);
            }
            break (Outcome {
                kind: CallKind::HeteroAcgtDeletion,
                reference,
                alternate,
                probability: maximum,
                hetero: false,
            });
        }
        if hetero_del_del
            .iter()
            .any(|(_, probability)| *probability == maximum)
        {
            let index = hetero_del_del
                .iter()
                .position(|(_, probability)| *probability == maximum)
                .unwrap();
            let (mut first_length, mut second_length) = hetero_del_del[index].0;
            if options.add_indel_length && first_length < second_length {
                std::mem::swap(&mut first_length, &mut second_length);
            }
            let mut bases = Vec::new();
            if options.add_indel_length
                && let Some(first) = one_indel(
                    &deletions,
                    Some(first_length).filter(|length| *length < LENGTH_OFFSET),
                    false,
                    "",
                    options.maximum_variant_length,
                )
                && let Some(second) = one_indel(
                    &deletions,
                    Some(second_length).filter(|length| *length < LENGTH_OFFSET),
                    false,
                    &first,
                    options.maximum_variant_length,
                )
            {
                bases = if first.len() > second.len() {
                    vec![first, second]
                } else {
                    vec![second, first]
                };
            }
            if bases.len() < 2 {
                bases = select_indel(
                    &deletions,
                    None,
                    false,
                    "",
                    options.maximum_variant_length,
                    true,
                );
            }
            if bases.len() < 2 {
                hetero_del_del[index].1 = 0.0;
                continue;
            }
            let reference = format!("{reference_raw}{}", bases[0]);
            let first = reference.chars().next().unwrap().to_string();
            let second = format!("{}{}", first, &reference[bases[1].len() + 1..]);
            if first == second || reference == first || reference == second {
                break (Outcome {
                    kind: CallKind::HeteroDeletionDeletion,
                    reference,
                    alternate: first,
                    probability: maximum,
                    hetero: false,
                });
            }
            break (Outcome {
                kind: CallKind::HeteroDeletionDeletion,
                reference,
                alternate: format!("{first},{second}"),
                probability: maximum,
                hetero: false,
            });
        }
        let index = insertion_deletion
            .iter()
            .position(|(_, probability)| *probability == maximum)
            .unwrap();
        let (deletion_length, insertion_length) = insertion_deletion[index].0;
        let insertion = options
            .add_indel_length
            .then_some(insertion_length)
            .filter(|length| *length < LENGTH_OFFSET);
        let deletion = options
            .add_indel_length
            .then_some(deletion_length)
            .filter(|length| *length < LENGTH_OFFSET);
        let insertion = one_indel(
            &insertions,
            insertion,
            true,
            "",
            options.maximum_variant_length,
        );
        let deletion = one_indel(
            &deletions,
            deletion,
            false,
            "",
            options.maximum_variant_length,
        );
        let (Some(insertion), Some(deletion)) = (insertion, deletion) else {
            insertion_deletion[index].1 = 0.0;
            continue;
        };
        let reference = format!("{reference_raw}{deletion}");
        let alternate = format!(
            "{},{}",
            reference.chars().next().unwrap(),
            insertion.to_string() + &reference[1..]
        );
        break (Outcome {
            kind: CallKind::InsertionDeletion,
            reference,
            alternate,
            probability: maximum,
            hetero: false,
        });
    };
    outcome.hetero = matches!(
        outcome.kind,
        CallKind::HeteroSnp
            | CallKind::HeteroAcgtInsertion
            | CallKind::HeteroInsertionInsertion
            | CallKind::HeteroAcgtDeletion
            | CallKind::HeteroDeletionDeletion
            | CallKind::InsertionDeletion
    ) || hetero_snp.contains(&outcome.probability)
        || hetero_acgt_insertion
            .iter()
            .any(|(_, _, probability)| *probability == outcome.probability)
        || hetero_ins_ins
            .iter()
            .any(|(_, probability)| *probability == outcome.probability)
        || hetero_acgt_deletion
            .iter()
            .any(|(_, _, probability)| *probability == outcome.probability)
        || hetero_del_del
            .iter()
            .any(|(_, probability)| *probability == outcome.probability)
        || insertion_deletion
            .iter()
            .any(|(_, probability)| *probability == outcome.probability);
    Ok(outcome)
}

fn long_indel_count(
    entries: &[(String, i64)],
    proposed_insertion: &str,
    proposed_deletion_length: usize,
    deletion: bool,
) -> i64 {
    if proposed_insertion.len() <= 50 && proposed_deletion_length <= 50 {
        return 0;
    }
    let proposed_length = if deletion {
        proposed_deletion_length
    } else {
        proposed_insertion.len().saturating_sub(1)
    };
    let minimum = (proposed_length as f64 * 0.9).max(50.0);
    let maximum = proposed_length as f64 * 1.1;
    entries
        .iter()
        .filter(|(bases, _)| {
            !(deletion && bases.len() == proposed_deletion_length)
                && !(!deletion && bases == proposed_insertion)
                && bases.len() as f64 >= minimum
                && bases.len() as f64 <= maximum
        })
        .map(|(_, count)| *count)
        .sum()
}

fn partial_label(reference: &str, alternate: &str) -> String {
    match reference.len().cmp(&alternate.len()) {
        Ordering::Greater => "Del".to_string(),
        Ordering::Less => "Ins".to_string(),
        Ordering::Equal => alternate.chars().next().unwrap_or('.').to_string(),
    }
}

fn mixed_label(first: &str, second: &str) -> String {
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
    if !first.is_empty() && first == second {
        return format!("{first}{second}");
    }
    "InsDel".to_string()
}

fn genotype_likelihoods(
    genotype: &str,
    prediction: &Prediction,
    reference: &str,
    alternate: &str,
) -> Vec<i32> {
    let _ = genotype;
    if alternate == "." {
        return vec![990];
    }
    let alternates: Vec<&str> = alternate.split(',').collect();
    let pairs: &[(usize, usize)] = if alternates.len() == 1 {
        &[(0, 0), (0, 1), (1, 1)]
    } else {
        &[(0, 0), (0, 1), (1, 1), (0, 2), (1, 2), (2, 2)]
    };
    let reference = if reference.len() == 1 {
        base_from_iupac(reference.chars().next().unwrap()).to_string()
    } else {
        reference.to_string()
    };
    let mut likelihoods = Vec::new();
    for (first, second) in pairs {
        let first_base = if *first == 0 {
            reference.as_str()
        } else {
            alternates[*first - 1]
        };
        let second_base = if *second == 0 {
            reference.as_str()
        } else {
            alternates[*second - 1]
        };
        let label = mixed_label(
            &partial_label(&reference, first_base),
            &partial_label(&reference, second_base),
        );
        let Some(gt21_index) = GT21_LABELS.iter().position(|candidate| *candidate == label) else {
            return vec![990; pairs.len()];
        };
        let zygosity = if *first == 0 && *second == 0 {
            0
        } else if first == second {
            1
        } else {
            2
        };
        likelihoods.push(prediction.gt21[gt21_index] * prediction.genotype[zygosity]);
    }
    let sum: f32 = likelihoods.iter().sum();
    let mut phred: Vec<f64> = likelihoods
        .iter()
        .map(|probability| {
            let probability = (*probability / sum) as f64 + 1e-8;
            -10.0 * probability.log10()
        })
        .collect();
    let minimum = phred.iter().copied().fold(f64::INFINITY, f64::min);
    phred
        .drain(..)
        .map(|value| (value - minimum).ceil() as i32)
        .collect()
}

pub fn decode_variant(
    contig: &str,
    position: i64,
    reference_sequence: &str,
    center: usize,
    alternate_info: &str,
    prediction: &Prediction,
    options: &DecodeOptions,
) -> Result<Option<VcfRecord>> {
    if prediction.gt21.len() != 21 || prediction.genotype.len() != 3 {
        bail!("expected 21 GT21 and 3 genotype probabilities");
    }
    if options.add_indel_length
        && (prediction.length1.len() != 33 || prediction.length2.len() != 33)
    {
        bail!("expected two vectors of 33 indel-length probabilities");
    }
    if center >= reference_sequence.len() || !reference_sequence.is_ascii() {
        bail!("invalid reference sequence or center");
    }
    let alt = AltInfo::parse(alternate_info)?;
    let outcome = decode_outcome(reference_sequence, center, prediction, options, &alt)?;
    if (!options.show_reference && outcome.kind == CallKind::Reference)
        || (outcome.kind != CallKind::Reference && outcome.reference == outcome.alternate)
    {
        return Ok(None);
    }
    let is_multi = outcome.alternate.contains(',');
    if options.haploid_precise && outcome.hetero {
        return Ok(None);
    }
    if options.haploid_sensitive && is_multi {
        return Ok(None);
    }
    let mut genotype = match outcome.kind {
        CallKind::Reference => "0/0",
        CallKind::HomoSnp | CallKind::HomoInsertion | CallKind::HomoDeletion => "1/1",
        _ => "0/1",
    };
    if is_multi {
        genotype = "1/2";
    }
    if options.haploid_precise || options.haploid_sensitive {
        genotype = if genotype.contains('1') { "1" } else { "0" };
    }
    let reference_count = alt.reference_count();
    let snp = alt.of_type('X');
    let insertions = alt.of_type('I');
    let deletions = alt.of_type('D');
    let count = |entries: &[(String, i64)], base: &str| {
        entries
            .iter()
            .find(|(candidate, _)| candidate == base)
            .map(|(_, count)| *count)
            .unwrap_or(0)
    };
    let long_count = |entries: &[(String, i64)], insertion: &str, deletion: usize| {
        if options.enable_long_indel {
            long_indel_count(entries, insertion, deletion, false)
        } else {
            0
        }
    };
    let mut supported = 0;
    let mut alternate_counts = Vec::new();
    match outcome.kind {
        CallKind::Reference => supported = reference_count,
        CallKind::HomoSnp | CallKind::HeteroSnp => {
            for base in outcome.alternate.chars().filter(|base| *base != ',') {
                let value = count(&snp, &base.to_string());
                supported += value;
                alternate_counts.push(value);
            }
        }
        CallKind::HomoInsertion | CallKind::HeteroInsertionInsertion => {
            for base in outcome.alternate.split(',') {
                let value = count(&insertions, base) + long_count(&insertions, base, 0);
                supported += value;
                alternate_counts.push(value);
            }
        }
        CallKind::HeteroAcgtInsertion => {
            let alleles: Vec<&str> = outcome.alternate.split(',').collect();
            let (snp_base, insertion) = if is_multi {
                (alleles[0].chars().next(), alleles[1])
            } else {
                (None, alleles[0])
            };
            if let Some(base) = snp_base {
                let value = count(&snp, &base.to_string());
                supported += value;
                alternate_counts.push(value);
            }
            let value = count(&insertions, insertion) + long_count(&insertions, insertion, 0);
            supported += value;
            alternate_counts.push(value);
        }
        CallKind::HomoDeletion => {
            if !deletions.is_empty() {
                let deleted = &outcome.reference[1..];
                let value = count(&deletions, deleted) + long_count(&deletions, "", deleted.len());
                supported = value;
                alternate_counts.push(value);
            }
        }
        CallKind::HeteroDeletionDeletion => {
            if deletions.len() > 1 {
                for base in outcome.alternate.split(',') {
                    let length = outcome.reference.len() - base.len();
                    let value = deletions
                        .iter()
                        .find(|(candidate, _)| candidate.len() == length)
                        .map(|(_, count)| *count)
                        .unwrap_or(0)
                        + long_count(&deletions, "", length);
                    supported += value;
                    alternate_counts.push(value);
                }
            }
        }
        CallKind::HeteroAcgtDeletion => {
            let alleles: Vec<&str> = outcome.alternate.split(',').collect();
            if is_multi {
                let base = alleles[1].chars().next().unwrap();
                let value = count(&snp, &base.to_string());
                supported += value;
                alternate_counts.push(value);
            }
            let deleted = &outcome.reference[1..];
            let value = count(&deletions, deleted) + long_count(&deletions, "", deleted.len());
            supported += value;
            alternate_counts.push(value);
        }
        CallKind::InsertionDeletion => {
            for base in outcome.alternate.split(',') {
                let length = outcome.reference.len() as i64 - base.len() as i64;
                let value = if length < 0 {
                    let inserted = if outcome.reference.len() > 1 {
                        &base[..base.len() - (outcome.reference.len() - 1)]
                    } else {
                        base
                    };
                    count(&insertions, inserted) + long_count(&insertions, inserted, 0)
                } else {
                    deletions
                        .iter()
                        .find(|(candidate, _)| candidate.len() == length as usize)
                        .map(|(_, count)| *count)
                        .unwrap_or(0)
                        + long_count(&deletions, "", length as usize)
                };
                supported += value;
                alternate_counts.push(value);
            }
        }
    }
    let alternate = if outcome.kind == CallKind::Reference {
        ".".to_string()
    } else {
        outcome.alternate.clone()
    };
    let allele_frequency = |count: i64| {
        if alt.depth == 0 {
            0.0
        } else {
            (count as f64 / alt.depth as f64).min(1.0)
        }
    };
    let allele_frequency_string = if alternate_counts.len() <= 1 {
        format!("{:.4}", allele_frequency(supported))
    } else {
        alternate_counts
            .iter()
            .map(|count| format!("{:.4}", allele_frequency(*count)))
            .collect::<Vec<_>>()
            .join(",")
    };
    let probability = outcome.probability as f64;
    let mut quality = (-10.0_f64 / 10.0_f64.ln()
        * (((1.0 - probability) + 1e-10) / (probability + 1e-10)).ln()
        + 10.0)
        .max(0.0);
    quality = (quality * 100.0).round_ties_even() / 100.0;
    let filter = if outcome.kind == CallKind::Reference {
        "RefCall"
    } else if options
        .quality_for_pass
        .is_none_or(|threshold| quality >= threshold)
    {
        "PASS"
    } else {
        "LowQual"
    };
    let (reference, alternate) = if options.keep_iupac {
        (outcome.reference, alternate)
    } else {
        (
            convert_iupac_to_n(&outcome.reference),
            convert_iupac_to_n(&alternate),
        )
    };
    let likelihoods = options.gvcf.then(|| {
        genotype_likelihoods(genotype, prediction, &reference, &alternate)
            .iter()
            .map(i32::to_string)
            .collect::<Vec<_>>()
            .join(",")
    });
    let allele_depth = if alternate_counts.is_empty() {
        reference_count.to_string()
    } else {
        format!(
            "{},{}",
            reference_count,
            alternate_counts
                .iter()
                .map(|count| count.to_string())
                .collect::<Vec<_>>()
                .join(",")
        )
    };
    let mut sample = format!(
        "{}:{}:{}:{}:{}",
        genotype, quality as i32, alt.depth, allele_depth, allele_frequency_string
    );
    let mut format = "GT:GQ:DP:AD:AF".to_string();
    if let Some(likelihoods) = likelihoods {
        format.push_str(":PL");
        sample.push(':');
        sample.push_str(&likelihoods);
    }
    Ok(Some(VcfRecord {
        contig: contig.to_string(),
        position,
        id: ".".to_string(),
        reference,
        alternate,
        quality: format!("{quality:.2}"),
        filter: filter.to_string(),
        info: if options.pileup { "P" } else { "F" }.to_string(),
        format: Some(format),
        samples: vec![sample],
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multiallelic_indels_reverse_equal_count_groups() {
        let entries = ["TAA", "TA", "TAAAAA"]
            .into_iter()
            .map(|bases| (bases.to_string(), 11))
            .collect::<Vec<_>>();
        assert_eq!(
            select_indel(&entries, None, true, "", 50, true),
            ["TAAAAA", "TA"]
        );
        assert_eq!(
            one_indel(&entries, None, true, "", 50).as_deref(),
            Some("TAA")
        );
    }
}
