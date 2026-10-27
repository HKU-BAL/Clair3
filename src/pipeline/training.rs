use std::{
    collections::BTreeSet,
    fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use clap::Args;
use flate2::{Compression, write::GzEncoder};

use crate::{
    io::{compression::read_text_lines, variant::candidate_positions},
    pipeline::{commands::ensure_output_parent, merge::read_vcf},
};

#[derive(Args)]
#[command(rename_all = "snake_case")]
pub struct GetTruthArgs {
    #[arg(long = "vcf_fn")]
    vcf: PathBuf,
    #[arg(long = "var_fn", default_value = "PIPE")]
    output: PathBuf,
    #[arg(long = "ctgName")]
    contig: Option<String>,
    #[arg(long = "ctgStart")]
    start: Option<i64>,
    #[arg(long = "ctgEnd")]
    end: Option<i64>,
    #[arg(long = "truth_vcf_fn")]
    truth_vcf: Option<PathBuf>,
}

fn write_truth<W: Write>(args: &GetTruthArgs, writer: &mut W) -> Result<()> {
    let input = read_vcf(&args.vcf)?;
    let truth_positions = args
        .truth_vcf
        .as_ref()
        .map(|path| candidate_positions(path, args.contig.as_deref()))
        .transpose()?
        .unwrap_or_default();
    let mut positions = BTreeSet::new();
    let mut last_contig = args.contig.clone();
    for record in input.records {
        if args
            .contig
            .as_deref()
            .is_some_and(|name| name != record.contig)
        {
            continue;
        }
        if args
            .start
            .zip(args.end)
            .is_some_and(|(start, end)| record.position < start || record.position > end)
        {
            continue;
        }
        last_contig = Some(record.contig.clone());
        let sample = record
            .samples
            .last()
            .context("VCF record has no sample column")?;
        let genotype = sample.split(':').next().unwrap_or("0");
        let mut alleles = genotype
            .replace('/', "|")
            .replace('.', "0")
            .split('|')
            .map(str::parse::<i32>)
            .collect::<std::result::Result<Vec<_>, _>>()?;
        if alleles.len() != 2 {
            bail!(
                "expected a diploid genotype at {}:{}",
                record.contig,
                record.position
            );
        }
        alleles.sort_unstable();
        let mut alternate = record.alternate;
        if alternate.contains('*') {
            let alternate_alleles: Vec<_> = alternate.split(',').collect();
            if alleles[0] + alleles[1] != 3 || alternate_alleles.len() != 2 {
                eprintln!("error with variant representation");
                continue;
            }
            alternate = alternate_alleles
                .into_iter()
                .filter(|base| *base != "*")
                .collect::<Vec<_>>()
                .join("");
            alleles = vec![0, 1];
        }
        positions.insert(record.position);
        writeln!(
            writer,
            "{} {} {} {} {} {}",
            record.contig, record.position, record.reference, alternate, alleles[0], alleles[1]
        )?;
    }
    for position in truth_positions {
        if positions.contains(&position) {
            continue;
        }
        let contig = last_contig
            .as_deref()
            .context("cannot determine contig for missing truth variants")?;
        writeln!(writer, "{contig} {position} None None -1 -1")?;
    }
    Ok(())
}

pub fn get_truth(args: GetTruthArgs) -> Result<()> {
    if args.output == Path::new("PIPE") || args.output == Path::new("-") {
        let mut writer = BufWriter::new(std::io::stdout().lock());
        write_truth(&args, &mut writer)?;
        writer.flush()?;
    } else {
        ensure_output_parent(&args.output)?;
        let mut writer = GzEncoder::new(fs::File::create(&args.output)?, Compression::default());
        write_truth(&args, &mut writer)?;
        writer.finish()?;
    }
    Ok(())
}

#[derive(Args)]
#[command(rename_all = "snake_case")]
pub struct SplitExtendBedArgs {
    #[arg(long = "bed_fn")]
    bed: Option<PathBuf>,
    #[arg(long = "output_fn")]
    output: PathBuf,
    #[arg(long = "expand_region_size", default_value_t = 33)]
    expand: i64,
    #[arg(long = "ctgName")]
    contig: Option<String>,
    #[arg(long = "ctgStart")]
    start: Option<i64>,
    #[arg(long = "ctgEnd")]
    end: Option<i64>,
}

pub fn split_extend_bed(args: SplitExtendBedArgs) -> Result<()> {
    let Some(bed) = args.bed else {
        return Ok(());
    };
    if args.expand < 0 {
        bail!("expand_region_size must be nonnegative");
    }
    if args
        .start
        .zip(args.end)
        .is_some_and(|(start, end)| end < start)
    {
        bail!("ctgEnd must be greater than or equal to ctgStart");
    }
    let mut regions = Vec::new();
    for line in read_text_lines(&bed)? {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() < 3 {
            bail!("invalid BED row: {line}");
        }
        let contig = fields[0];
        if args.contig.as_deref().is_some_and(|name| name != contig) {
            continue;
        }
        let start: i64 = fields[1].parse()?;
        let end: i64 = fields[2].parse()?;
        if start < 0 || end < start {
            bail!("invalid BED interval: {line}");
        }
        if args.start.is_some_and(|minimum| end < minimum) {
            continue;
        }
        if args.end.is_some_and(|maximum| start > maximum) {
            continue;
        }
        regions.push((contig.to_string(), start, end));
    }
    regions.sort_by(|left, right| {
        left.0
            .cmp(&right.0)
            .then(left.1.cmp(&right.1))
            .then(left.2.cmp(&right.2))
    });
    let mut output = Vec::new();
    let mut current: Option<(String, i64, i64)> = None;
    for (contig, start, end) in regions {
        let start = (start - args.expand).max(0);
        let end = end + args.expand;
        if let Some((previous_contig, previous_start, previous_end)) = &mut current {
            if *previous_contig == contig && start <= *previous_end {
                *previous_end = (*previous_end).max(end);
            } else {
                output.push(format!("{previous_contig} {previous_start} {previous_end}"));
                current = Some((contig, start, end));
            }
        } else {
            current = Some((contig, start, end));
        }
    }
    if let Some((contig, start, end)) = current {
        output.push(format!("{contig} {start} {end}"));
    }
    ensure_output_parent(&args.output)?;
    fs::write(&args.output, output.join("\n"))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn truth_and_bed_commands_work() {
        let directory = tempfile::tempdir().unwrap();
        let vcf = directory.path().join("input.vcf");
        let header =
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
        fs::write(
            &vcf,
            format!(
                "{header}chrToy\t5\t.\tA\tC\t.\t.\t.\tGT\t0/1\nchrToy\t100\t.\tG\tC,*\t.\t.\t.\tGT\t1/2\n"
            ),
        )
        .unwrap();
        let truth_vcf = directory.path().join("truth.vcf");
        fs::write(
            &truth_vcf,
            format!(
                "{header}chrToy\t5\t.\tA\tC\t.\t.\t.\tGT\t0/1\nchrToy\t1500\t.\tT\tG\t.\t.\t.\tGT\t0/1\n"
            ),
        )
        .unwrap();
        let output = directory.path().join("truth.var.gz");
        get_truth(GetTruthArgs {
            vcf,
            output: output.clone(),
            contig: Some("chrToy".to_string()),
            start: Some(1),
            end: Some(2000),
            truth_vcf: Some(truth_vcf),
        })
        .unwrap();
        let rows: Vec<_> = read_text_lines(&output)
            .unwrap()
            .map(|line| line.unwrap())
            .collect();
        assert_eq!(
            rows,
            [
                "chrToy 5 A C 0 1",
                "chrToy 100 G C 0 1",
                "chrToy 1500 None None -1 -1"
            ]
        );
        let bed = directory.path().join("input.bed");
        fs::write(
            &bed,
            "chrToy\t0\t10\nchrToy\t20\t30\nchrToy\t50\t60\nchr:colon\t0\t5\n",
        )
        .unwrap();
        let output = directory.path().join("extended.bed");
        split_extend_bed(SplitExtendBedArgs {
            bed: Some(bed),
            output: output.clone(),
            expand: 5,
            contig: None,
            start: Some(0),
            end: Some(100),
        })
        .unwrap();
        assert_eq!(
            fs::read_to_string(output).unwrap(),
            "chr:colon 0 10\nchrToy 0 35\nchrToy 45 65"
        );
    }
}
