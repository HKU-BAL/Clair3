"""Regression checks for the native gVCF likelihood helper."""
import math
import unittest

from preprocess import utils


class LikelihoodMathTest(unittest.TestCase):
    def setUp(self):
        self.calculator = utils.mathcalculator()
        if not self.calculator.speedUp:
            self.skipTest("The native probability helper could not be compiled")

    def test_maximum_stops_at_the_requested_length(self):
        # The fourth value belongs to the allocation, but not the requested
        # three-element range. It makes a one-past-end read deterministic.
        values = self.calculator.ffi.new(
            "double[]", [-2.0, -4.0, -5.0, 1000000.0])
        self.assertEqual(self.calculator.lib.getMyMaxItem(values, 3), -2.0)

    def test_logsumexp_handles_small_probabilities(self):
        self.assertAlmostEqual(
            self.calculator.log10sumexp([-1000.0, -1001.0, -1002.0]),
            -1000.0 + math.log(1.11) / utils.LOG_10,
            places=10,
        )


if __name__ == "__main__":
    unittest.main()
