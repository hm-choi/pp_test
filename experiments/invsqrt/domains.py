"""invSqrt input domains for the HE-DAP search and measurements.

- normalized : x = V / V_max in [v_min, 1] (fixed coefficients,
               coeffs/invSqrt_coeffs_normalized.py);
- raw        : x = V in [1e-3, V_max] (coefficients generated from the public
               bounds), one domain per dataset comparison and test
               (Welch V_max, F Var2_max), as in HEHypothesisTesting.
"""

import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1]))

from pptest_cases import load_cases


NORMALIZED_V_MIN = [1e-3, 1e-4, 1e-5]
RAW_V_MIN = 1e-3


def all_domains():

    domains = [
        {"name": f"normalized_vmin{v:g}", "method": "normalized", "lo": v, "hi": 1.0}
        for v in NORMALIZED_V_MIN
    ]

    for name, case in load_cases().items():
        domains.append({
            "name": f"raw_{name}_welch", "method": "raw",
            "lo": RAW_V_MIN, "hi": case["welch_v_max"],
        })
        domains.append({
            "name": f"raw_{name}_f", "method": "raw",
            "lo": RAW_V_MIN, "hi": case["f_var2_max"],
        })

    return domains


def get_domains(spec):

    domains = all_domains()

    if spec == "all":
        return domains

    if spec in ("normalized", "raw"):
        return [d for d in domains if d["method"] == spec]

    names = spec.split(",")
    selected = [d for d in domains if d["name"] in names]

    if len(selected) != len(names):
        raise ValueError(f"Unknown domain in {spec}")

    return selected
