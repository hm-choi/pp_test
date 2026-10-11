"""invSqrt input domain for the HE-DAP search.

x = V / V_max in [V_BAR_MIN, 1] (Welch's t-test, V normalized by its public
bound), the same domain for every dataset comparison, as in
HEHypothesisTesting.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[2]))

from operators.HEhypo_test import V_BAR_MIN


def all_domains():

    return [{"name": "welch", "lo": V_BAR_MIN, "hi": 1.0}]


def get_domains(spec):

    domains = all_domains()

    if spec == "all":
        return domains

    names = spec.split(",")
    selected = [d for d in domains if d["name"] in names]

    if len(selected) != len(names):
        raise ValueError(f"Unknown domain in {spec}")

    return selected
