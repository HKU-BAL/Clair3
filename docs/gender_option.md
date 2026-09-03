# New Feature: `--gender` / `--par_regions_bed`


Clair3 now supports a `--gender` option so that variant calling handles the sex chromosomes (chrX/chrY)
with the correct ploidy. Previously, chrX/chrY of male samples were treated as diploid by default,
producing false heterozygous (`0/1`) calls on these haploid chromosomes (issue #66 of ClairS, https://github.com/HKU-BAL/ClairS/issues/66).
This update lets the user specify the sample gender so that chrX/chrY are called at the correct ploidy.

---

## Parameters

### `--gender {unknown, male, female}` (default: `unknown`)

| Value | Effect |
|---|---|
| `unknown` | Existing behavior: chrX and chrY are both called as diploid (backward compatible, no change). |
| `male` | Automatically applies haploid precise mode to chrX/chrY: `0/1` is dropped and `1/1` becomes haploid `1`. 
| `female` | Removes chrY from the calling contig list; chrX remains diploid. |

### `--par_regions_bed=FILE` (optional)

A BED file of pseudo-autosomal regions (PAR), 0-based half-open intervals. PAR are the segments at
both ends of chrX/chrY that are homologous to each other and recombine during meiosis; these regions
are in fact diploid. When this BED is provided, haploid mode will NOT be applied to sites on chrX/chrY inside these intervals. 
This is only meaningful when `--gender male` is set.

**Note**: PAR coordinates differ between reference builds. Please use the BED file that matches
your reference:

| Reference | Official PAR BED|
|---|---|
| GRCh38 (hg38) | https://www.ncbi.nlm.nih.gov/grc/human|
| CHM13 (T2T) | https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0_PAR.bed |

---

## Interaction with existing `--haploid_precise` / `--haploid_sensitive` parameters

- These two haploid flags are unchanged: they are **global** and apply to all contigs.
- If a global haploid flag is set, it **takes precedence** over `--gender`.
- `--gender male` is a convenience: it automatically applies haploid (precise) mode to the
  sex chromosomes only, with no extra flags required.

---

## Examples

```bash
# Male sample: chrX/chrY haploid, autosomes diploid
--gender male

# Male sample with PAR exemption: chrX/chrY haploid, PAR regions kept diploid
--gender male --par_regions_bed=GRCh38_PAR.bed

# Female sample: chrY removed, chrX remains diploid
--gender female
```

---

## Accuracy benefit (GRCh38 HG002 chrX/chrY, 20x sample (ONT platform), against the GIAB chrXY truth set)

| Condition | SNP Recall | SNP Precision | SNP F1 |
|---|---|---|
| `--gender unknown`  | 0.967 | 0.810 | 0.882 |
| `--gender male` | 0.935 | 0.979 | 0.973 |
| `--gender male --par_regions_bed` | 0.967 | 0.978 | 0.973 |

`--gender male` eliminates the false heterozygous calls on chrX/chrY (genotype false positives
`FP.gt` drop from 1516 to 7); `--par_regions_bed` keeps the PAR
regions diploid and recovers recall.
