"""invSqrt input domain for the HE-DAP search.

V in [_INV_SQRT_MIN, _INV_SQRT_MAX] = [0.01, 1000] (Welch's t-test), the same
domain for every dataset comparison, as in HEHypothesisTesting.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

from coeffs.invSqrt_coeffs import _INV_SQRT_MIN, _INV_SQRT_MAX


def all_domains():

    return [{"name": "welch", "lo": _INV_SQRT_MIN, "hi": _INV_SQRT_MAX}]


def get_domains(spec):

    domains = all_domains()

    if spec == "all":
        return domains

    names = spec.split(",")
    selected = [d for d in domains if d["name"] in names]

    if len(selected) != len(names):
        raise ValueError(f"Unknown domain in {spec}")

    return selected
