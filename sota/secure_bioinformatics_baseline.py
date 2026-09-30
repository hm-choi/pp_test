"""Compatibility import for a retired duplicate baseline.

The previous file was HE-only Goldschmidt code but was named as though it were
a complete Secure Bioinformatics protocol. It remains only for old imports.
"""

from sota.goldschmidt_welch_he_baseline import GoldschmidtWelchHEBaseline


SecureBioinformaticsBaseline = GoldschmidtWelchHEBaseline

__all__ = ["GoldschmidtWelchHEBaseline", "SecureBioinformaticsBaseline"]
