"""invSqrt input domains for the HE-DAP search and measurements.

x = V in [INV_SQRT_V_MIN, V_max] (coefficients generated from the public
bounds), one domain per dataset comparison (Welch V_max), as in
HEHypothesisTesting.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))
sys.path.append(str(Path(__file__).resolve().parents[2]))

from pptest_cases import load_cases
from operators.HEhypo_test import INV_SQRT_V_MIN


def all_domains():

    domains = []

    for name, case in load_cases().items():
        domains.append({
            "name": f"{name}_welch",
            "lo": INV_SQRT_V_MIN, "hi": case["welch_v_max"],
        })

    return domains


def get_domains(spec):

    domains = all_domains()

    if spec == "all":
        return domains

    names = spec.split(",")
    selected = [d for d in domains if d["name"] in names]

    if len(selected) != len(names):
        raise ValueError(f"Unknown domain in {spec}")

    return selected
