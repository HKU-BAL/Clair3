# Clair3's performance on indels > 50bp

(Document revision 0, benchmarks based on Clair3 v0.1-r8)

Inspired by @**[HenrivdGeest](https://github.com/HenrivdGeest)** in [#64](https://github.com/HKU-BAL/Clair3/issues/64), we questioned: __Is it still sensible to limit to calling indels ≤50bp in short variant callers using long-reads? How good is Clair3 on finding indels longer than 50bp?__

It's a continuous debate on whether to categorize indels > 50bp as SV, or those > 1000bp. But back in the days when NGS short-read is the only choice for SV calling, "> 50bp" is a rule of thumb. However, with long-reads, the average length of a reliable gap opened in alignment is much longer than 50bp.

Clair3‘s design does not have a length limit on the indels called, although, in the current implementation (v0.1-r8), we ignored those indel signals > 50bp. We made a few changes in Clair3 to make it output all I(nsertion) and D(eletion) decisions made by the 21-genotype task, and benchmarked those > 50bp called in [50x of HG002](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=NHGRI_UCSC_panel/HG002/nanopore/Guppy_4.2.2/) against the [GIAB HG002 SV v0.6 truth set](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NIST_SV_v0.6/). The results are as follows.

Clair3 achieved 0.963 precision and 0.529 recall without any filtering. Compare to SV callers including cuteSV (Precision: 0.928, Recall: 0.966) and Sniffle (Precision: 0.913, Recall: 0.926), the recall is low, but the precision looks good. It makes sense because while the SV callers are using multiple techniques to capture SV signals, Clair3 (or basically all small variant callers) is only effectively relying on the "split-read" technique to find SV. But compared to the recall achievable by short reads, which I don't have a handy number but I doubt it could exceed 0.1 (let me know if you have a better estimate), 0.529 is pretty good as a starting point. We looked into a few SVs we missed and we think the recall has a potential to go beyond 0.7 if we relieve the requirement for a consistent starting position to call an indel, and pay more attention to the indels between 50 to 200bp, at which the recall is particularly low.

Back to the question __Is it still sensible to limit to calling indels ≤50bp in short variant callers using long-reads?__ I think it yes and no. Yes because performance-wise, short variant caller cannot take over SV callers' job to call insertions and deletions between 50bp and 1000bp (yet). No because the implicit indel length cap 50bp is no longer suitable for short variant callers using long-reads, and should be continuously revised in the future. In the next Clair3 release (v0.1-r9), we will provide an option for outputting indel calls longer than 50bp.

## Compare to SV callers

| Callers  | Precision | Recall | F1-score |
| -------- | --------- | ------ | -------- |
| Clair3   | 0.963     | 0.529  | 0.684    |
| cuteSV   | 0.928     | 0.966  | 0.947    |
| Sniffles | 0.913     | 0.926  | 0.920    |

_Using truvari --includebed HG002_SVs_Tier1_v0.6.bed --passonly -r 1000 -p 0_

## Deletions - ONT HG002 50x Guppy4.2.2 GRCh37 

| Length (bp) | Precision | Recall | FP | TP  | FN  |
| ----------- | --------- | ------ | -- | --- | --- |
| 50-100      | 0.959     | 0.449  | 28 | 654 | 803 |
| 100-200     | 0.964     | 0.513  | 13 | 353 | 335 |
| 200-300     | 0.979     | 0.648  | 4  | 186 | 101 |
| 300-400     | 0.992     | 0.860  | 7  | 826 | 135 |
| 400-500     | 1.000     | 0.763  | 0  | 58  | 18  |
| 500-600     | 0.963     | 0.578  | 1  | 26  | 19  |
| 600-700     | 0.962     | 0.625  | 1  | 25  | 15  |
| 700-800     | 1.000     | 0.766  | 0  | 36  | 11  |
| 800-900     | 1.000     | 0.788  | 0  | 26  | 7   |
| 900-1000    | 1.000     | 0.730  | 0  | 27  | 10  |
| \>1000      | 0.983     | 0.748  | 7  | 395 | 133 |


## Insertions - ONT HG002 50x Guppy4.2.2 GRCh37
| Length (bp) | Precision | Recall | FP | TP  | FN  |
| ----------- | --------- | ------ | -- | --- | --- |
| 50-100      | 0.919     | 0.399  | 53 | 603 | 909 |
| 100-200     | 0.906     | 0.340  | 33 | 317 | 614 |
| 200-300     | 0.926     | 0.487  | 18 | 226 | 238 |
| 300-400     | 0.986     | 0.645  | 10 | 725 | 399 |
| 400-500     | 0.958     | 0.391  | 3  | 68  | 106 |
| 500-600     | 0.953     | 0.477  | 3  | 61  | 67  |
| 600-700     | 0.942     | 0.395  | 3  | 49  | 75  |
| 700-800     | 0.981     | 0.464  | 1  | 51  | 59  |
| 800-900     | 1.000     | 0.403  | 0  | 27  | 40  |
| 900-1000    | 0.963     | 0.433  | 1  | 26  | 34  |
| \>1000      | 0.983     | 0.459  | 6  | 343 | 405 |

