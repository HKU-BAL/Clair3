use std::path::Path;

use anyhow::{Context, Result};

use crate::{io::variant::VcfRecord, pipeline::merge::read_vcf};

#[derive(Debug, Clone, Copy, Default)]
struct Counts {
    false_positive: usize,
    false_negative: usize,
    truth_positive: usize,
    query_positive: usize,
}

#[derive(Debug, Clone)]
pub struct MetricRow {
    pub variant_type: &'static str,
    pub query_fp: usize,
    pub truth_fn: usize,
    pub truth_tp: usize,
    pub query_tp: usize,
    pub precision: f64,
    pub recall: f64,
    pub f1: f64,
}

fn field<'a>(record: &'a VcfRecord, format: &[&str], sample: usize, name: &str) -> Result<&'a str> {
    let index = format
        .iter()
        .position(|field| *field == name)
        .with_context(|| format!("missing hap.py field {name}"))?;
    record
        .samples
        .get(sample)
        .and_then(|sample| sample.split(':').nth(index))
        .with_context(|| format!("missing hap.py sample field {name}"))
}

fn metrics(counts: Counts) -> (f64, f64, f64) {
    let precision = if counts.query_positive + counts.false_positive == 0 {
        0.0
    } else {
        counts.query_positive as f64 / (counts.query_positive + counts.false_positive) as f64
    };
    let recall = if counts.truth_positive + counts.false_negative == 0 {
        0.0
    } else {
        counts.truth_positive as f64 / (counts.truth_positive + counts.false_negative) as f64
    };
    let f1 = if precision + recall == 0.0 {
        0.0
    } else {
        2.0 * precision * recall / (precision + recall)
    };
    let round = |value: f64| (value * 1_000_000.0).round() / 1_000_000.0;
    (round(precision), round(recall), round(f1))
}

fn row(name: &'static str, counts: Counts) -> MetricRow {
    let (precision, recall, f1) = metrics(counts);
    MetricRow {
        variant_type: name,
        query_fp: counts.false_positive,
        truth_fn: counts.false_negative,
        truth_tp: counts.truth_positive,
        query_tp: counts.query_positive,
        precision,
        recall,
        f1,
    }
}

fn count_decision(counts: &mut Counts, query: bool, decision: &str) {
    match (query, decision) {
        (true, "TP") => counts.query_positive += 1,
        (false, "TP") => counts.truth_positive += 1,
        (true, "FP") => counts.false_positive += 1,
        (false, "FN") => counts.false_negative += 1,
        _ => {}
    }
}

pub fn overall_metrics(path: &Path, contig: Option<&str>) -> Result<Vec<MetricRow>> {
    let file = read_vcf(path)?;
    let mut overall = Counts::default();
    let mut snp = Counts::default();
    let mut indel = Counts::default();
    let mut insertion = Counts::default();
    let mut deletion = Counts::default();
    for record in file.records {
        if contig.is_some_and(|name| name != record.contig) {
            continue;
        }
        let format: Vec<_> = record.format.as_deref().unwrap_or("").split(':').collect();
        let truth_decision = field(&record, &format, 0, "BD")?;
        let truth_location = field(&record, &format, 0, "BI")?;
        let truth_type = field(&record, &format, 0, "BVT")?;
        let query_decision = field(&record, &format, 1, "BD")?;
        let query_location = field(&record, &format, 1, "BI")?;
        let query_type = field(&record, &format, 1, "BVT")?;
        for (kind, query) in [(query_type, true), (truth_type, false)] {
            let decision = if query {
                query_decision
            } else {
                truth_decision
            };
            let location = if query {
                query_location
            } else {
                truth_location
            };
            let counts = match kind {
                "SNP" => &mut snp,
                "INDEL" => &mut indel,
                _ => continue,
            };
            count_decision(counts, query, decision);
            if kind == "INDEL" && location.starts_with('i') {
                count_decision(&mut insertion, query, decision);
            } else if kind == "INDEL" && location.starts_with('d') {
                count_decision(&mut deletion, query, decision);
            }
        }
    }
    overall.truth_positive = snp.truth_positive + indel.truth_positive;
    overall.query_positive = snp.query_positive + indel.query_positive;
    overall.false_positive = snp.false_positive + indel.false_positive;
    overall.false_negative = snp.false_negative + indel.false_negative;
    Ok(vec![
        row("Overall", overall),
        row("SNP", snp),
        row("INDEL", indel),
        row("INS", insertion),
        row("DEL", deletion),
    ])
}

fn format_float(value: f64) -> String {
    if value.fract() == 0.0 {
        format!("{value:.1}")
    } else {
        value.to_string()
    }
}

pub fn format_metrics(rows: &[MetricRow]) -> String {
    let mut lines = vec![
        [
            "VariantType",
            "TRUTH.FP",
            "TRUTH.FN",
            "TRUTH.TP",
            "QUERY.TP",
            "METRIC.Precision",
            "METRIC.Recall",
            "METRIC.F1_Score",
        ]
        .iter()
        .map(|value| format!("{value:<20}"))
        .collect::<String>(),
    ];
    for row in rows {
        let values = [
            row.variant_type.to_string(),
            row.query_fp.to_string(),
            row.truth_fn.to_string(),
            row.truth_tp.to_string(),
            row.query_tp.to_string(),
            format_float(row.precision),
            format_float(row.recall),
            format_float(row.f1),
        ];
        lines.push(
            values
                .iter()
                .map(|value| format!("{value:<20}"))
                .collect::<String>(),
        );
    }
    lines.push(String::new());
    lines.join("\n") + "\n\n"
}

pub fn overall_report(
    path: &Path,
    contig: Option<&str>,
    happy_log: Option<&Path>,
) -> Result<String> {
    let rows = overall_metrics(path, contig)?;
    let mut report = format_metrics(&rows);
    if let Some(path) = happy_log.filter(|path| path.is_file()) {
        let mut snp = Vec::new();
        let mut indel = Vec::new();
        for line in std::fs::read_to_string(path)?
            .lines()
            .filter(|line| line.contains("PASS"))
        {
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len() < 14 {
                continue;
            }
            let values = vec![
                fields[6].to_string(),
                fields[4].to_string(),
                fields[3].to_string(),
                fields[11].to_string(),
                fields[10].to_string(),
                fields[13].to_string(),
            ];
            if line.contains("INDEL") {
                indel = values;
            } else if line.contains("SNP") {
                snp = values;
            }
        }
        let values = [&rows[1], &rows[2]]
            .iter()
            .flat_map(|row| [row.precision, row.recall, row.f1])
            .map(|value| format!("{value:.6}"))
            .chain([
                rows[1].query_fp.to_string(),
                rows[1].truth_fn.to_string(),
                rows[1].truth_tp.to_string(),
                rows[2].query_fp.to_string(),
                rows[2].truth_fn.to_string(),
                rows[2].truth_tp.to_string(),
            ])
            .collect::<Vec<_>>();
        report.push_str("Double check with happy log:\n");
        report.push_str(&values.join(" "));
        report.push('\n');
        report.push_str(
            &[
                snp[3..].to_vec(),
                indel[3..].to_vec(),
                snp[..3].to_vec(),
                indel[..3].to_vec(),
            ]
            .concat()
            .join(" "),
        );
        report.push_str("\n\n\n");
    }
    let probabilities = rows
        .iter()
        .flat_map(|row| [row.precision, row.recall, row.f1])
        .map(format_float)
        .collect::<Vec<_>>();
    report.push_str(&probabilities.join(" "));
    report.push('\n');
    let counts = [
        rows[0].query_fp,
        rows[0].truth_fn,
        rows[0].truth_tp,
        rows[1].query_fp,
        rows[1].truth_fn,
        rows[1].truth_tp,
        rows[2].query_fp,
        rows[2].truth_fn,
        rows[2].truth_tp,
        rows[3].query_tp,
        rows[3].truth_tp,
        rows[3].truth_tp,
        rows[4].query_fp,
        rows[4].truth_fn,
        rows[4].truth_tp,
    ];
    report.push_str(
        &counts
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(" "),
    );
    report.push('\n');
    Ok(report)
}
