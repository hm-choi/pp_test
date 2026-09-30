"""Goldschmidt oracle-bound diagnostic for the insurance Welch instance.

The bounds below were selected after inspecting this plaintext instance.  This
is useful for diagnosing the best-case Goldschmidt arithmetic path, but is NOT
a privacy-valid or paper-faithful performance comparison.  A valid comparison
must instead receive positive denominator lower bounds as public protocol
assumptions before encryption, or use a protocol that establishes them.
"""

from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from pp_test.engine.HEengine import HEEngine
from operators.hypothesis_testing import HEHypothesisTesting
from operators.operator import HEOperator
from sota.goldschmidt_welch_he_baseline import (
    GoldschmidtWelchHEBaseline,
    PublicWelchBounds,
)


ALPHA = 0.05
CHEBYSHEV_DEGREE = 15
GOLDSCHMIDT_ITERATIONS = 6

# ORACLE DIAGNOSTIC ONLY: these values were selected from the insurance
# plaintext reference V and D.  Do not cite runtimes from this configuration as
# a private-protocol or paper-to-paper comparison.
ORACLE_INSURANCE_BOUNDS = PublicWelchBounds(
    v_lower=0.50,
    v_upper=0.55,
    d_lower=0.00080,
    d_upper=0.00095,
    t2_upper=1200.0,
    df_upper=350.0,
)


def plaintext_welch(x1: np.ndarray, x2: np.ndarray) -> dict[str, float]:
    var1 = float(np.var(x1, ddof=1))
    var2 = float(np.var(x2, ddof=1))
    v = var1 / len(x1) + var2 / len(x2)
    d = ((var1 / len(x1)) ** 2) / (len(x1) - 1) + (
        (var2 / len(x2)) ** 2
    ) / (len(x2) - 1)
    t = float((np.mean(x1) - np.mean(x2)) / np.sqrt(v))
    return {"T": t, "T2": t * t, "V": v, "D": d, "DF": v * v / d}


def assert_instance_is_in_public_range(
    reference: dict[str, float], bounds: PublicWelchBounds
) -> None:
    """Evaluation-only check; it does not supply values to HE operations."""

    checks = (
        ("V", reference["V"], bounds.v_lower, bounds.v_upper),
        ("D", reference["D"], bounds.d_lower, bounds.d_upper),
        ("T2", reference["T2"], 0.0, bounds.t2_upper),
        ("DF", reference["DF"], 0.0, bounds.df_upper),
    )
    for name, value, lower, upper in checks:
        if not lower <= value <= upper:
            raise ValueError(
                f"Public Goldschmidt {name} bound is invalid for this "
                f"evaluation instance: {value} not in [{lower}, {upper}]."
            )


def main() -> None:
    insurance = pd.read_csv(PROJECT_ROOT / "datasets" / "insurance.csv")
    x1 = insurance.loc[insurance["smoker"] == "yes", "charges"].to_numpy(
        dtype=np.float64
    ) / 1_000.0
    x2 = insurance.loc[insurance["smoker"] == "no", "charges"].to_numpy(
        dtype=np.float64
    ) / 1_000.0
    reference = plaintext_welch(x1, x2)
    ORACLE_INSURANCE_BOUNDS.validate()
    assert_instance_is_in_public_range(reference, ORACLE_INSURANCE_BOUNDS)

    input_upper_bound = float(max(np.max(x1), np.max(x2)))
    plaintext_critical = float(stats.t.ppf(1.0 - ALPHA / 2.0, reference["DF"]))
    plaintext_score = reference["T2"] - plaintext_critical**2

    engine = HEEngine(warmup_bootstrap=True)
    operator = HEOperator(engine)
    proposed = HEHypothesisTesting(engine)
    goldschmidt = GoldschmidtWelchHEBaseline(engine, verbose=True)
    encrypted_x1, encrypted_x2 = operator.encrypt(x1), operator.encrypt(x2)

    proposed.reset_bootstrap_count()
    started_at = perf_counter()
    proposed_result = proposed.HE_Welch_T_Test(
        encrypted_x1,
        encrypted_x2,
        R=input_upper_bound,
        alpha=ALPHA,
        critical_degree=CHEBYSHEV_DEGREE,
        score_bound=1_200.0,
        debug=True,
    )
    proposed_welch_seconds = perf_counter() - started_at
    proposed_welch_bootstraps = proposed.bootstrap_count()
    proposed_t = float(operator.decrypt(proposed_result["T"])[0])
    proposed_inv_df = float(operator.decrypt(proposed_result["INV_DF"])[0])
    # This reciprocal is reporting-only. No plaintext reciprocal is fed back
    # into the encrypted computation.
    proposed_df_for_report = 1.0 / proposed_inv_df

    goldschmidt.reset_bootstrap_count()
    started_at = perf_counter()
    baseline_result = goldschmidt.welch_t2_df(
        encrypted_x1,
        encrypted_x2,
        ORACLE_INSURANCE_BOUNDS,
        num_iter=GOLDSCHMIDT_ITERATIONS,
        return_debug=True,
    )
    baseline_seconds = perf_counter() - started_at
    baseline_bootstraps = goldschmidt.bootstrap_count()
    baseline_t2 = float(operator.decrypt(baseline_result["T2"])[0])
    baseline_df = float(operator.decrypt(baseline_result["DF"])[0])

    encrypted_score = float(operator.decrypt(proposed_result["score"])[0])
    encrypted_step = float(operator.decrypt(proposed_result["step"])[0])

    print("========== ORACLE-BOUND Goldschmidt diagnostic ==========")
    print("WARNING: Goldschmidt bounds were selected from plaintext reference")
    print("         This is not a privacy-valid or paper-faithful speed comparison.")
    print("engine                 : HEaaN FGb, log_slots=15, CPU")
    print(f"dataset                : insurance charges by smoker ({len(x1)}, {len(x2)})")
    print(f"Goldschmidt iterations : {GOLDSCHMIDT_ITERATIONS}")
    print("Goldschmidt depth/iter : 1 (q update and residual square are parallel)")
    print(f"Goldschmidt bounds     : {ORACLE_INSURANCE_BOUNDS}")
    print("\nplaintext T2, df      : " f"{reference['T2']}, {reference['DF']}")
    print("\n[Proposed: InvSqrt Welch arithmetic]")
    print(f"encrypted T            : {proposed_t}")
    print(f"reported T2            : {proposed_t * proposed_t}")
    print(f"reported df (=1/InvDF) : {proposed_df_for_report}")
    print(f"T absolute error       : {abs(proposed_t - reference['T'])}")
    print(f"df absolute error      : {abs(proposed_df_for_report - reference['DF'])}")
    print(f"runtime (seconds)      : {proposed_welch_seconds}")
    print(f"bootstraps             : {proposed_welch_bootstraps}")
    print("\n[Goldschmidt HE arithmetic baseline]")
    print(f"encrypted T2           : {baseline_t2}")
    print(f"encrypted df           : {baseline_df}")
    print(f"T2 absolute error      : {abs(baseline_t2 - reference['T2'])}")
    print(f"df absolute error      : {abs(baseline_df - reference['DF'])}")
    print(f"runtime (seconds)      : {baseline_seconds}")
    print(f"bootstraps             : {baseline_bootstraps}")
    print("\n[Proposed encrypted decision extension]")
    print(f"alpha                  : {ALPHA}")
    print(f"plaintext score        : {plaintext_score}")
    print(f"encrypted score        : {encrypted_score}")
    print(f"score absolute error   : {abs(encrypted_score - plaintext_score)}")
    print(f"encrypted step          : {encrypted_step}")
    print("runtime, bootstraps are included in the proposed path above.")
    print(f"decision match         : {(plaintext_score > 0.0) == (encrypted_score > 0.0)}")
    print("Goldschmidt decision   : N/A (no paper-faithful HE-only lookup/decision)")


if __name__ == "__main__":
    main()
