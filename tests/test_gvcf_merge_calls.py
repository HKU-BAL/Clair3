"""Regression tests for gvcfGenerator.mergeCalls interval merging."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from preprocess import utils

VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
)


def _var_line(pos, ref="A", alt="G"):
    return (
        f"chr1\t{pos}\t.\t{ref}\t{alt}\t30\tPASS\tP\t"
        "GT:GQ:DP:AD:AF:PL\t0/1:30:20:10,10:0.5:30,0,30\n"
    )


def _nonvar_line(start, end, ref="C"):
    return (
        f"chr1\t{start}\t.\t{ref}\t<NON_REF>\t0\t.\tEND={end}\t"
        "GT:GQ:MIN_DP:PL\t0/0:40:20:0,60,900\n"
    )


class GvcfMergeCallsTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.ref_seq = "ACGT" * 100
        compression = patch.object(utils, "COMPRESS_GVCF", False)
        compression.start()
        self.addCleanup(compression.stop)

    def _ref_base(self, pos_spec):
        pos = int(pos_spec.split(":")[1].split("-")[0])
        return self.ref_seq[pos - 1]

    def _merge(self, var_lines, nonvar_lines):
        var_path = self.root / "var.vcf"
        nonvar_path = self.root / "nonvar.gvcf"
        out_path = self.root / "out.gvcf"
        var_path.write_text(VCF_HEADER + "".join(var_lines))
        nonvar_path.write_text(VCF_HEADER + "".join(nonvar_lines))
        generator = utils.gvcfGenerator(ref_path="dummy.fa")
        with patch.object(generator, "readReferenceBaseAtPos", side_effect=self._ref_base):
            generator.mergeCalls(
                str(var_path), str(nonvar_path), str(out_path), "SAMPLE", "chr1"
            )
        records = []
        for line in out_path.read_text().splitlines():
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            pos = int(cols[1])
            ref = cols[3]
            alt = cols[4]
            info = cols[7]
            if info.startswith("END="):
                end = int(info.split("=")[1])
                kind = "nonvar"
            else:
                end = pos + len(ref) - 1
                kind = "var"
            records.append((kind, pos, end, ref, alt))
        return records

    def test_variant_before_nonvar_block_is_preserved(self):
        records = self._merge(
            [_var_line(50)],
            [_nonvar_line(100, 200)],
        )
        self.assertEqual(
            records,
            [
                ("var", 50, 50, "A", "G,<NON_REF>"),
                ("nonvar", 100, 200, "C", "<NON_REF>"),
            ],
        )

    def test_last_variant_inside_nonvar_block_splits_right_tail(self):
        records = self._merge(
            [_var_line(150)],
            [_nonvar_line(100, 200)],
        )
        self.assertEqual(
            records,
            [
                ("nonvar", 100, 149, "C", "<NON_REF>"),
                ("var", 150, 150, "A", "G,<NON_REF>"),
                ("nonvar", 151, 200, self.ref_seq[150], "<NON_REF>"),
            ],
        )

    def test_last_variant_at_nonvar_start_trims_right_tail(self):
        records = self._merge(
            [_var_line(100, ref="C", alt="G")],
            [_nonvar_line(100, 200, ref="C")],
        )
        self.assertEqual(
            records,
            [
                ("var", 100, 100, "C", "G,<NON_REF>"),
                ("nonvar", 101, 200, self.ref_seq[100], "<NON_REF>"),
            ],
        )

    def test_empty_variant_file_keeps_first_nonvar_block(self):
        records = self._merge(
            [],
            [_nonvar_line(100, 200), _nonvar_line(201, 300, ref="G")],
        )
        self.assertEqual(
            records,
            [
                ("nonvar", 100, 200, "C", "<NON_REF>"),
                ("nonvar", 201, 300, "G", "<NON_REF>"),
            ],
        )

    def test_multiple_variants_within_single_nonvar_block(self):
        records = self._merge(
            [_var_line(120), _var_line(160)],
            [_nonvar_line(100, 200)],
        )
        self.assertEqual(
            records,
            [
                ("nonvar", 100, 119, "C", "<NON_REF>"),
                ("var", 120, 120, "A", "G,<NON_REF>"),
                ("nonvar", 121, 159, self.ref_seq[120], "<NON_REF>"),
                ("var", 160, 160, "A", "G,<NON_REF>"),
                ("nonvar", 161, 200, self.ref_seq[160], "<NON_REF>"),
            ],
        )

    def test_deletion_spanning_block_and_gap_before_next_block(self):
        records = self._merge(
            [_var_line(110, ref="ACGTACGTACGT", alt="A"), _var_line(210)],  # del spans 110..121
            [_nonvar_line(100, 115), _nonvar_line(150, 200)],
        )
        self.assertEqual(
            records,
            [
                ("nonvar", 100, 109, "C", "<NON_REF>"),
                ("var", 110, 121, "ACGTACGTACGT", "A,<NON_REF>"),
                ("nonvar", 150, 200, "C", "<NON_REF>"),
                ("var", 210, 210, "A", "G,<NON_REF>"),
            ],
        )

    def test_deletion_covering_subsequent_overlapping_variant(self):
        records = self._merge(
            [_var_line(100, ref="ACGTACGTACG", alt="A"), _var_line(105)],
            [_nonvar_line(100, 200)],
        )
        self.assertEqual(
            records,
            [
                ("var", 100, 110, "ACGTACGTACG", "A,<NON_REF>"),
                ("var", 105, 105, "A", "G,<NON_REF>"),
                ("nonvar", 111, 200, self.ref_seq[110], "<NON_REF>"),
            ],
        )


if __name__ == "__main__":
    unittest.main()
