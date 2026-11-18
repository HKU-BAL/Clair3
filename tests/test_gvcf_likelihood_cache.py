"""Exact-output and ownership checks for the bounded gVCF likelihood cache."""
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from preprocess import utils


class GvcfLikelihoodCacheTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.reference = self.root / "reference.fa"
        self.reference.write_text(">test\n" + "A" * 10000 + "\n")
        Path(str(self.reference) + ".fai").write_text("test\t10000\t6\t10000\t10001\n")
        self.calculators = []
        compression = patch.object(utils, "COMPRESS_GVCF", False)
        compression.start()
        self.addCleanup(compression.stop)

    def calculator(self, error=0.001, bin_size=5, bp_resolution=False):
        output = self.root / str(len(self.calculators))
        calculator = utils.variantInfoCalculator(
            str(output), str(self.reference), error, bin_size, "test",
            bp_resolution=bp_resolution, sample_name="sample",
        )
        self.calculators.append(calculator)
        self.addCleanup(calculator.close_vcf_writer)
        return calculator, output / "sample.tmp.gvcf"

    @staticmethod
    def summary(n_ref=29, n_total=30, **changes):
        summary = dict(chr="test", pos=1, ref="A", n_ref=n_ref, n_total=n_total)
        summary.update(changes)
        return summary

    def test_reuse_does_not_share_site_metadata_or_mutable_pl(self):
        calculator, _ = self.calculator()
        with patch.object(calculator, "_cal_reference_likelihood",
                          wraps=calculator._cal_reference_likelihood) as calculate:
            first = calculator.reference_likelihood(self.summary())
            expected_pl = first["pl"][:]
            first["pl"][0] = -999
            second = calculator.reference_likelihood(self.summary(pos=20, ref="C", chr="other"))
            self.assertEqual(calculate.call_count, 1)
        self.assertEqual(second["pl"], expected_pl)
        self.assertEqual((second["pos"], second["END"], second["ref"], second["chr"]),
                         (20, 20, "C", "other"))

    def test_parameters_and_ambiguous_reference_do_not_cross_caches(self):
        first, _ = self.calculator(error=0.001, bin_size=5)
        second, _ = self.calculator(error=0.1, bin_size=1)
        summary = self.summary(10, 12)
        expected = second._reference_likelihood_values(10, 12)
        first.reference_likelihood(summary)
        actual = second.reference_likelihood(summary)
        self.assertEqual(actual["pl"], list(expected[3]))
        self.assertEqual(actual["binned_gq"], expected[2])
        self.assertNotEqual(first.reference_likelihood(summary)["pl"], actual["pl"])
        ambiguous = second.reference_likelihood(dict(summary, ref="N"))
        self.assertEqual((ambiguous["gq"], ambiguous["pl"]), (1, [0, 0, 0]))
        self.assertEqual(second.reference_likelihood(summary), actual)

    def test_cache_is_bounded_and_eviction_preserves_values(self):
        calculator, _ = self.calculator()
        expected = calculator.reference_likelihood(self.summary(0, 0))
        for depth in range(1, 4200):
            calculator.reference_likelihood(self.summary(depth, depth))
        self.assertEqual(calculator._cached_reference_likelihood.cache_info().currsize, 4096)
        self.assertEqual(calculator.reference_likelihood(self.summary(0, 0)), expected)

    def test_gvcf_blocks_equal_uncached_output(self):
        for bp_resolution in (False, True):
            with self.subTest(bp_resolution=bp_resolution):
                before, before_path = self.calculator(bp_resolution=bp_resolution)
                after, after_path = self.calculator(bp_resolution=bp_resolution)
                before._cached_reference_likelihood = before._reference_likelihood_values
                for position in range(1, 1200):
                    depth = (position // 31) % 60
                    summary = self.summary(
                        max(0, depth - position % 3), depth, pos=position,
                        ref="N" if 250 <= position < 300 else "ACGT"[position % 4],
                    )
                    before.make_gvcf_online(summary)
                    after.make_gvcf_online(summary)
                for calculator in (before, after):
                    calculator.make_gvcf_online({}, push_current=True)
                    calculator.vcf_writer.flush()
                self.assertEqual(before_path.read_bytes(), after_path.read_bytes())


if __name__ == "__main__":
    unittest.main()
