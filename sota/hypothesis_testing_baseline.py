"""Compatibility import for the renamed Goldschmidt Welch HE baseline.

Use :mod:`sota.goldschmidt_welch_he_baseline` in new experiments.
"""

from sota.goldschmidt_welch_he_baseline import GoldschmidtWelchHEBaseline


HEHypothesisTestingBaseline = GoldschmidtWelchHEBaseline

__all__ = ["GoldschmidtWelchHEBaseline", "HEHypothesisTestingBaseline"]
