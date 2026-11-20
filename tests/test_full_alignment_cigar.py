"""Native tensor regressions. Build libclair3, then run unittest discovery.

Requires the usual numpy dependency and pysam to construct indexed test BAMs.
The fixtures are generated locally; no reference genome or model is downloaded.
"""
from array import array
from pathlib import Path
import tempfile
import unittest

import libclair3
import numpy as np
import pysam


class FullAlignmentCigarTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.root = Path(self.directory.name)
        self.fasta = self.root / "reference.fa"
        self.fasta.write_text(">test\n" + "A" * 4000 + "\n")
        pysam.faidx(str(self.fasta))

    def tearDown(self):
        self.directory.cleanup()

    def write_read(self, cigar, sequence, reverse=False, signals=None, start=100):
        bam = self.root / "reads.bam"
        read = pysam.AlignedSegment()
        read.query_name = "read"
        read.flag = 16 if reverse else 0
        read.reference_id = 0
        read.reference_start = start
        read.mapping_quality = 60
        read.cigartuples = cigar
        read.query_sequence = sequence
        read.query_qualities = array("B", [40] * len(sequence))
        if signals is not None:
            moves = [5]
            for length in signals:
                moves.extend([1] + [0] * (length - 1))
            read.set_tag("mv", array("b", moves))
        with pysam.AlignmentFile(bam, "wb", header={
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "test", "LN": 4000}],
        }) as out:
            out.write(read)
        pysam.index(str(bam))
        return bam

    def write_reads(self, reads):
        """Write forward-strand all-A reads given as (start, cigar) pairs."""
        bam = self.root / "reads.bam"
        with pysam.AlignmentFile(bam, "wb", header={
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "test", "LN": 4000}],
        }) as out:
            for index, (start, cigar) in enumerate(
                    sorted(reads, key=lambda read: read[0])):
                length = sum(n for op, n in cigar if op in (0, 1, 4, 7, 8))
                read = pysam.AlignedSegment()
                read.query_name = "read%d" % index
                read.flag = 0
                read.reference_id = 0
                read.reference_start = start
                read.mapping_quality = 60
                read.cigartuples = cigar
                read.query_sequence = "A" * length
                read.query_qualities = array("B", [40] * length)
                out.write(read)
        pysam.index(str(bam))
        return bam

    @staticmethod
    def depth_and_counts(alt_info):
        """Return (depth, {allele: count}) from a native alt_info string."""
        _, depth, _, alleles = alt_info.decode().split("-", 3)
        fields = alleles.split()
        return int(depth), dict(zip(fields[::2], map(int, fields[1::2])))

    def tensors(self, bam, candidates, dwell=False):
        ffi, lib = libclair3.ffi, libclair3.lib
        result = lib.calculate_clair3_full_alignment(
            b"test:1-4000", str(bam).encode(), str(self.fasta).encode(),
            ffi.new("struct Variant *[]", 0), 0,
            ffi.new("size_t[]", candidates), len(candidates),
            False, 5, 0, 1, 50, dwell, False,
        )
        try:
            channels = 9 if dwell else 8
            tensor = np.frombuffer(ffi.buffer(
                result.matrix, len(candidates) * 33 * channels,
            ), dtype=np.int8).reshape(len(candidates), 33, channels).copy()
            alleles = [ffi.string(result.all_alt_info[i])
                       for i in range(len(candidates))]
            return tensor, alleles
        finally:
            lib.destroy_fa_data(result)

    def test_window_extending_beyond_both_read_ends(self):
        bam = self.write_read([(0, 5)], "AAAAA")
        tensor, _ = self.tensors(bam, [102])
        expected = np.zeros((33, 8), dtype=np.int8)
        # Window 86..118 contains this read only at 100..104.
        expected[14:19] = [100, 0, 100, 100, 100, 0, 0, 60]
        np.testing.assert_array_equal(tensor[0], expected)

    def test_mixed_cigar_coordinates_and_dwell_on_both_strands(self):
        # 2S 5M 2I 4= 3D 2X 5N 5M 2S
        cigar = [(4, 2), (0, 5), (1, 2), (7, 4), (2, 3),
                 (8, 2), (3, 5), (0, 5), (4, 2)]
        sequence = "TT" + "ACGTA" + "CC" + "TGCA" + "AT" + "CGTAC" + "GG"
        signals = [1 + i % 7 for i in range(len(sequence))]
        # Reference positions and their observed query offsets, independently
        # enumerated for the CIGAR above. Inserted and clipped bases consume
        # query coordinates; deletions and reference skips do not.
        observed = dict(zip(
            [100, 101, 102, 103, 104, 105, 106, 107, 108, 112, 113,
             119, 120, 121, 122, 123],
            [2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        ))
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                bam = self.write_read(cigar, sequence, reverse, signals)
                tensor, _ = self.tensors(bam, [109], dwell=True)
                aligned_signals = signals[::-1] if reverse else signals
                for position in range(93, 126):
                    cell = tensor[0, position - 93]
                    if position not in observed:
                        np.testing.assert_array_equal(cell, 0)
                        continue
                    query_index = observed[position]
                    self.assertEqual(cell[0], 100)
                    self.assertEqual(cell[2], 50 if reverse else 100)
                    self.assertEqual(cell[3], 100)
                    self.assertEqual(cell[4], 100)
                    self.assertEqual(cell[7], 60)
                    expected_signal = aligned_signals[query_index]
                    if position == 104:
                        expected_signal += sum(aligned_signals[7:9])
                    self.assertEqual(cell[8], expected_signal)

    def test_disjoint_windows_keep_query_coordinates(self):
        sequence = "A" * 30 + "C" * 50 + "G" * 700 + "T" * 20
        bam = self.write_read([(0, len(sequence))], sequence)
        candidates = [125, 155, 700, 890]
        combined, combined_alleles = self.tensors(bam, candidates)
        for i, candidate in enumerate(candidates):
            isolated, alleles = self.tensors(bam, [candidate])
            np.testing.assert_array_equal(combined[i], isolated[0])
            self.assertEqual(combined_alleles[i], alleles[0])
        # Central mismatch encodings A/reference=0, C=25, G=75, T=50.
        np.testing.assert_array_equal(combined[:, 16, 1], [0, 25, 75, 50])

    def test_missing_move_tag_and_empty_candidates(self):
        bam = self.write_read([(0, 40)], "A" * 40)
        ordinary, alleles = self.tensors(bam, [120])
        with_dwell, dwell_alleles = self.tensors(bam, [120], dwell=True)
        np.testing.assert_array_equal(with_dwell[:, :, :8], ordinary)
        np.testing.assert_array_equal(with_dwell[:, :, 8], 0)
        self.assertEqual(alleles, dwell_alleles)
        empty, empty_alleles = self.tensors(bam, [])
        self.assertEqual(empty.shape, (0, 33, 8))
        self.assertEqual(empty_alleles, [])

    def test_reads_clipped_at_deletion_breakpoint_are_not_reference(self):
        # Candidate 200 is the anchor of a 20 bp deletion. Reads that stop
        # within a few bases after the anchor and are then clipped cannot show
        # whether the deletion is present and must not count as reference.
        reads = (
            [(150, [(0, 51), (2, 20), (0, 50)])] * 2      # deletion carriers
            + [(150, [(0, 52), (4, 50)])] * 3             # clipped at breakpoint
            + [(150, [(0, 120)])]                         # spans as reference
            + [(150, [(0, 60), (4, 30)])]                 # clipped 9 bp later
            + [(150, [(0, 51), (2, 2), (0, 2), (4, 50)])]  # 2 bp deletion carrier
        )
        _, alleles = self.tensors(self.write_reads(reads), [200])
        depth, counts = self.depth_and_counts(alleles[0])
        # Previously depth 8 and reference 5 (8 reference bases - 3 deletions).
        self.assertEqual(depth, 5)
        self.assertEqual(counts, {"D" + "A" * 20: 2, "DAA": 1, "RA": 2})

    def test_clipped_read_at_snv_only_candidate_is_reference(self):
        reads = [(550, [(0, 51), (4, 40)]), (550, [(0, 100)])]
        _, alleles = self.tensors(self.write_reads(reads), [600])
        self.assertEqual(alleles[0], b"601-2-A-RA 2 ")


if __name__ == "__main__":
    unittest.main()
