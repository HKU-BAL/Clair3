"""Count width, buffer ownership, indel reuse, and ambiguous-base regressions."""
from array import array
import ctypes
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import libclair3
import numpy as np
import pysam
from preprocess import CreateTensorPileupFromCffi as pileup
from preprocess.medaka_utils import Region


class PileupNativeTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.reference = self.root / 'reference.fa'
        self.reference.write_text('>test\n' + 'A' * 5000 + '\n')
        pysam.faidx(str(self.reference))

    def bam(self, reads):
        path = self.root / 'reads.bam'
        with pysam.AlignmentFile(path, 'wb', header={
            'HD': {'SO': 'coordinate'}, 'SQ': [{'SN': 'test', 'LN': 5000}],
        }) as output:
            for index, (start, cigar, sequence, reverse) in enumerate(reads):
                read = pysam.AlignedSegment()
                read.query_name = 'read' + str(index)
                read.reference_id = 0
                read.reference_start = start
                read.flag = 16 if reverse else 0
                read.mapping_quality = 60
                read.cigartuples = cigar
                read.query_sequence = sequence
                read.query_qualities = array('B', [40] * len(sequence))
                output.write(read)
        pysam.index(str(path))
        return path

    def counts(self, bam):
        ffi, lib = libclair3.ffi, libclair3.lib
        handle = lib.create_bam_fset(str(bam).encode(), str(self.reference).encode())
        result = lib.calculate_clair3_pileup(
            b'test:1-5000', handle, str(self.reference).encode(),
            1, 0.08, 0.15, 5, 2000, False, 144, True, True,
        )
        try:
            return pileup._plp_data_to_numpy(result, 18, gvcf=True)
        finally:
            lib.destroy_plp_data(result, True)
            lib.destroy_bam_fset(handle)

    def test_array_copy_uses_the_loaded_native_count_width(self):
        ffi, lib = libclair3.ffi, libclair3.lib
        data = lib.create_plp_data(2, 2, 18, 1, 1, 0)
        width = ffi.sizeof(ffi.typeof(data.matrix).item)
        data.matrix[0] = -17 if width == 4 else (1 << 64) - 17
        data.matrix[19] = 70000
        data.major[0], data.major[1] = 100, 101
        counts, positions, _, _ = pileup._plp_data_to_numpy(data, 18)
        lib.destroy_plp_data(data, False)
        self.assertEqual(counts.dtype.itemsize, width)
        self.assertEqual(counts[0, 0], -17)
        self.assertEqual(counts[1, 1], 70000)
        self.assertTrue(counts.flags.owndata)
        np.testing.assert_array_equal(positions['major'], [100, 101])

    def test_single_chunk_does_not_copy_and_adjacent_chunks_still_join(self):
        counts = np.arange(72, dtype=np.int32).reshape(4, 18)
        positions = np.zeros(4, dtype=[('major', int), ('minor', int)])
        positions['major'] = [100, 101, 102, 103]
        stitch = getattr(pileup, '__enforce_pileup_chunk_contiguity')
        chunks, _, _ = stitch([(counts, positions, [], [])])
        self.assertTrue(np.shares_memory(chunks[0][0], counts))
        joined, _, _ = stitch([(counts[:2], positions[:2], [], []),
                              (counts[2:], positions[2:], [], [])])
        np.testing.assert_array_equal(joined[0][0], counts)
        positions['major'][2:] += 10
        split, _, _ = stitch([(counts, positions, [], [])])
        self.assertEqual([len(c) for c, p in split], [2, 2])

    def test_fixed_allocation_and_growth_use_compact_elements(self):
        ffi, lib = libclair3.ffi, libclair3.lib
        data = lib.create_plp_data(2, 2, 18, 1, 1, 18)
        try:
            self.assertEqual(ffi.sizeof(ffi.typeof(data.matrix).item), 4)
            data.matrix[0], data.matrix[35] = -17, 70000
            # This existing C helper is not part of the generated CFFI API.
            native = ctypes.CDLL(libclair3.__file__)
            grow = native.enlarge_plp_data
            grow.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_size_t]
            grow.restype = None
            grow(int(ffi.cast('uintptr_t', data)), 5, 18)
            self.assertEqual(data.buffer_cols, 5)
            values = np.frombuffer(ffi.buffer(data.matrix, 5 * 18 * 4), dtype=np.int32)
            self.assertEqual((values[0], values[35]), (-17, 70000))
            np.testing.assert_array_equal(values[36:], 0)
        finally:
            lib.destroy_plp_data(data, False)

    def test_indel_counts_reset_after_long_deletions_and_coverage_gaps(self):
        reads = [
            (100, [(0, 5), (2, 70), (0, 5), (1, 2), (0, 5), (2, 1), (0, 10)], 'A' * 27, False),
            (100, [(0, 5), (2, 70), (0, 5)], 'A' * 10, False),
            (100, [(0, 5), (2, 34), (0, 5)], 'A' * 10, False),
            (220, [(0, 5), (2, 33), (0, 20)], 'A' * 25, True),
            (300, [(0, 5), (2, 1), (0, 20)], 'A' * 25, False),
        ]
        counts, positions, _, _ = self.counts(self.bam(reads))
        rows = dict(zip(positions['major'], counts))
        self.assertEqual((rows[104][6], rows[104][7]), (3, 2))
        self.assertEqual((rows[224][15], rows[224][16]), (1, 1))
        self.assertEqual((rows[179][4], rows[179][5]), (1, 1))
        for position, row in rows.items():
            if position not in [104, 184, 304]:
                self.assertEqual(row[6], 0)
            if position != 224:
                self.assertEqual(row[15], 0)
            if position != 179:
                self.assertEqual(row[4], 0)

    def test_ambiguous_query_base_does_not_write_before_matrix(self):
        counts, _, _, _ = self.counts(self.bam([
            (100, [(0, 1), (1, 2), (0, 20)], 'N' + 'A' * 22, False),
        ]))
        expected = np.zeros(18, dtype=np.int32)
        expected[4:6] = 1  # An insertion anchored on N still contributes.
        np.testing.assert_array_equal(counts[0], expected)
        self.assertEqual(counts[1, 0], -1)

    def test_counts_exceeding_int16_remain_exact(self):
        bam = self.bam((100, [(0, 20)], 'A' * 20, bool(i % 2)) for i in range(70000))
        counts, _, _, gvcf = self.counts(bam)
        self.assertEqual(counts.dtype, np.dtype(np.int32))
        np.testing.assert_array_equal(counts[:, 0], -35000)
        np.testing.assert_array_equal(counts[:, 9], -35000)
        self.assertEqual(gvcf[0][100], 70000)

    def test_file_handle_pool_matches_worker_count(self):
        bam = self.bam([(100, [(0, 40)], 'A' * 40, False)])
        original_init = pileup.BAMHandler.__init__
        sizes = []
        def record_init(instance, bam, fasta, size=16):
            sizes.append(size)
            original_init(instance, bam, fasta, size)
        with patch.object(pileup.BAMHandler, '__init__', record_init):
            pileup.pileup_counts_clair3(Region('test', 0, 5000), str(bam), str(self.reference),
                                        1, 0.08, 0.15, 5, False, 50, 144, workers=2)
        self.assertEqual(sizes, [2])


if __name__ == '__main__':
    unittest.main()
