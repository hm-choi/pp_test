"""Compare df and 1/df parameterizations of the Welch critical value.

The experiment has two parts.

1. A plaintext approximation sweep compares equal-degree Chebyshev
   polynomials under three parameterizations.
2. An optional HEaaN run compares the legacy encrypted DF and InvDF
   arithmetic on the same insurance data, including the degree-15 critical
   value evaluation.

Run the complete experiment on the HEaaN server with::

    python3 experiments/experiment3/test.py

The approximation-only part can run without HEaaN with::

    python3 experiments/experiment3/test.py --plaintext-only
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd
from numpy.polynomial import chebyshev
from scipy import stats


PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

ALPHA = 0.05
DEGREES = (7, 15, 31, 63, 127)
DF_GLOBAL_MIN = 1.0
DF_GLOBAL_MAX = 2048.0
GRID_SIZE = 200_000
INPUT_SCALE = 1_000.0
DATASET_PATH = PROJECT_ROOT / "datasets" / "insurance.csv"
APPROXIMATION_RESULT_PATH = Path(__file__).with_name(
    "approximation_results.csv"
)
HE_RESULT_PATH = Path(__file__).with_name("he_results.csv")


def critical_value(df: np.ndarray | float) -> np.ndarray:
    return stats.t.ppf(1.0 - ALPHA / 2.0, df)


def inv_df_target(z: np.ndarray) -> np.ndarray:
    """Student-t critical value on z = 2/df - 1."""
    inv_df = (np.asarray(z, dtype=np.float64) + 1.0) / 2.0
    result = np.empty_like(inv_df)
    finite_df = inv_df > np.finfo(np.float64).eps
    result[finite_df] = critical_value(1.0 / inv_df[finite_df])
    result[~finite_df] = stats.norm.ppf(1.0 - ALPHA / 2.0)
    return result


def df_target(z: np.ndarray, lower: float, upper: float) -> np.ndarray:
    """Student-t critical value after mapping df in [lower, upper] to z."""
    df = lower + (np.asarray(z, dtype=np.float64) + 1.0) * (
        upper - lower
    ) / 2.0
    return critical_value(df)


def make_coefficients(
    representation: str,
    degree: int,
    lower: float,
    upper: float,
) -> np.ndarray:
    if representation == "inv_df":
        return chebyshev.chebinterpolate(inv_df_target, degree)
    return chebyshev.chebinterpolate(
        lambda z: df_target(z, lower, upper), degree
    )


def approximation_grid(lower: float, upper: float) -> np.ndarray:
    """Cover both the small-df endpoint and the full interval densely."""
    linear = np.linspace(lower, upper, GRID_SIZE)
    if lower > 0.0:
        logarithmic = np.geomspace(lower, upper, GRID_SIZE)
        return np.unique(np.concatenate((linear, logarithmic)))
    return linear


def run_approximation_sweep() -> pd.DataFrame:
    insurance = pd.read_csv(DATASET_PATH)
    n1 = int((insurance["smoker"] == "yes").sum())
    n2 = int((insurance["smoker"] == "no").sum())
    # For positive group variances, the Welch-Satterthwaite df is bounded by
    # min(n1 - 1, n2 - 1) and n1 + n2 - 2.
    group_df_min = float(min(n1 - 1, n2 - 1))
    group_df_max = float(n1 + n2 - 2)

    configurations = (
        ("inv_df", "universal InvDF [0,1]", DF_GLOBAL_MIN, DF_GLOBAL_MAX),
        ("df", "global DF [1,2048]", DF_GLOBAL_MIN, DF_GLOBAL_MAX),
        (
            "df",
            f"group-aware DF [{group_df_min:g},{group_df_max:g}]",
            group_df_min,
            group_df_max,
        ),
    )
    rows: list[dict[str, float | int | str]] = []
    for representation, label, lower, upper in configurations:
        df_grid = approximation_grid(lower, upper)
        expected = critical_value(df_grid)
        if representation == "inv_df":
            z_grid = 2.0 / df_grid - 1.0
        else:
            z_grid = 2.0 * (df_grid - lower) / (upper - lower) - 1.0

        for degree in DEGREES:
            coefficients = make_coefficients(
                representation, degree, lower, upper
            )
            predicted = chebyshev.chebval(z_grid, coefficients)
            error = np.abs(predicted - expected)
            rows.append(
                {
                    "parameterization": label,
                    "degree": degree,
                    "df_min": lower,
                    "df_max": upper,
                    "max_absolute_error": float(np.max(error)),
                    "rmse": float(np.sqrt(np.mean(error**2))),
                    "p99_absolute_error": float(np.quantile(error, 0.99)),
                    "worst_df": float(df_grid[np.argmax(error)]),
                    "multiplicative_depth": int(math.ceil(math.log2(degree))),
                }
            )

    result = pd.DataFrame(rows)
    result.to_csv(APPROXIMATION_RESULT_PATH, index=False)
    print("\n========== Plaintext Chebyshev parameterization sweep ==========")
    print(result.to_string(index=False))
    print(f"\nCSV written to: {APPROXIMATION_RESULT_PATH}")
    return result


def evaluate_df_critical_he(
    engine,
    operator,
    encrypted_df,
    coefficients: np.ndarray,
    lower: float,
    upper: float,
):
    """Evaluate a direct-DF Chebyshev polynomial on an HEData scalar."""
    import heaan as hn

    from hedata.data import HEData

    required_levels = math.ceil(math.log2(len(coefficients) - 1)) + 3

    # HEaaN requires at least level 3 *before* bootstrapping.  DF commonly
    # arrives here at level 3, while the constant multiplication in the
    # affine map consumes one level.  Refresh DF first; otherwise z would be
    # level 2 and could no longer be bootstrapped.
    df_for_evaluation = operator.copy_new(encrypted_df)
    levels_before_affine_map = required_levels + 1
    if df_for_evaluation.level() < levels_before_affine_map:
        if df_for_evaluation.level() < 3:
            raise RuntimeError(
                "Direct-DF critical evaluation received a ciphertext below "
                f"the HEaaN bootstrap minimum: level={df_for_evaluation.level()}."
            )
        print(
            "DF critical pre-bootstrap level:",
            df_for_evaluation.level(),
        )
        df_for_evaluation = operator.do_bootstrapping(df_for_evaluation, 11)

    # z = 2(df - lower)/(upper - lower) - 1.
    z = operator.mult_const(
        df_for_evaluation, 2.0 / (upper - lower)
    )
    z = operator.sub_const(z, (upper + lower) / (upper - lower))
    if z.level() < required_levels:
        # This is a defensive check for parameter presets whose refreshed
        # level is lower than expected.  At this point bootstrap is legal only
        # while z still has at least three levels.
        if z.level() < 3:
            raise RuntimeError(
                "Insufficient level for the direct-DF polynomial after "
                f"normalization: level={z.level()}, required={required_levels}."
            )
        z = operator.do_bootstrapping(z, 11)

    he_coefficients = hn.math.approx.ChebyshevCoefficients(
        coefficients, len(coefficients)
    )
    ciphertexts = [
        hn.math.approx.evaluate_chebyshev_expansion(
            engine.evaluator(),
            engine.bootstrapping(),
            ciphertext,
            he_coefficients,
            1.0,
        )
        for ciphertext in z.ciphertexts()
    ]
    return HEData(ciphertexts, z.size(), ciphertexts[0].level, z.scale())


def run_he_comparison() -> pd.DataFrame:
    """Compare encrypted DF and InvDF paths on identical observations."""
    from pp_test.engine.HEengine import HEEngine
    from operators.hypothesis_testing_backup import (
        HEHypothesisTesting as LegacyHEHypothesisTesting,
    )
    from operators.operator import HEOperator

    data = pd.read_csv(DATASET_PATH)
    x1 = data.loc[data["smoker"] == "yes", "charges"].to_numpy(
        dtype=np.float64
    ) / INPUT_SCALE
    x2 = data.loc[data["smoker"] == "no", "charges"].to_numpy(
        dtype=np.float64
    ) / INPUT_SCALE
    n1, n2 = len(x1), len(x2)
    r = float(max(np.max(x1), np.max(x2)))

    # Compute the reference explicitly for compatibility with older SciPy
    # versions whose Ttest_indResult does not expose a ``df`` attribute.
    variance1 = float(np.var(x1, ddof=1))
    variance2 = float(np.var(x2, ddof=1))
    a = variance1 / n1
    b = variance2 / n2
    v = a + b
    plaintext_t = float((np.mean(x1) - np.mean(x2)) / np.sqrt(v))
    plaintext_df = float(
        v**2 / (a**2 / (n1 - 1) + b**2 / (n2 - 1))
    )
    plaintext_inv_df = 1.0 / plaintext_df
    plaintext_critical = float(critical_value(plaintext_df))

    engine = HEEngine(warmup_bootstrap=True)
    io_operator = HEOperator(engine)
    encrypted_x1 = io_operator.encrypt(x1)
    encrypted_x2 = io_operator.encrypt(x2)
    rows: list[dict[str, float | int | str]] = []

    # InvDF path: universal critical-value coefficients can be reused for all
    # public group sizes.
    inv_test = LegacyHEHypothesisTesting(engine)
    inv_test.reset_bootstrap_count()
    started = perf_counter()
    inv_result = inv_test.HE_Welch_T_Test_With_InvDF(
        encrypted_x1, encrypted_x2, R=r, return_debug=True
    )
    arithmetic_seconds = perf_counter() - started
    arithmetic_bootstraps = inv_test.bootstrap_count()
    started = perf_counter()
    inv_critical_ct = inv_test.HE_T_Critical_Value_From_InvDF(
        inv_result["INV_DF"], ALPHA, degree=15
    )
    critical_seconds = perf_counter() - started
    total_bootstraps = inv_test.bootstrap_count()
    encrypted_t = float(io_operator.decrypt(inv_result["T"])[0])
    encrypted_inv_df = float(io_operator.decrypt(inv_result["INV_DF"])[0])
    encrypted_critical = float(io_operator.decrypt(inv_critical_ct)[0])
    rows.append(
        {
            "path": "InvDF (universal [0,1])",
            "plaintext_parameter": plaintext_inv_df,
            "encrypted_parameter": encrypted_inv_df,
            "parameter_absolute_error": abs(
                encrypted_inv_df - plaintext_inv_df
            ),
            "parameter_relative_error": abs(
                encrypted_inv_df - plaintext_inv_df
            ) / plaintext_inv_df,
            "t_absolute_error": abs(encrypted_t - plaintext_t),
            "critical_absolute_error": abs(
                encrypted_critical - plaintext_critical
            ),
            "arithmetic_seconds": arithmetic_seconds,
            "critical_seconds": critical_seconds,
            "total_seconds": arithmetic_seconds + critical_seconds,
            "arithmetic_bootstraps": arithmetic_bootstraps,
            "critical_bootstraps": total_bootstraps - arithmetic_bootstraps,
            "total_bootstraps": total_bootstraps,
        }
    )

    # Direct DF path: use the tight, public group-size-dependent interval.
    df_lower = float(min(n1 - 1, n2 - 1))
    df_upper = float(n1 + n2 - 2)
    df_coefficients = make_coefficients("df", 15, df_lower, df_upper)
    df_test = LegacyHEHypothesisTesting(engine)
    df_test.reset_bootstrap_count()
    critical_operator = HEOperator(engine)
    critical_operator.reset_bootstrap_count()
    started = perf_counter()
    df_result = df_test.HE_Welch_T_Test_With_DF(
        encrypted_x1, encrypted_x2, R=r, return_debug=True
    )
    arithmetic_seconds = perf_counter() - started
    arithmetic_bootstraps = df_test.bootstrap_count()
    started = perf_counter()
    df_critical_ct = evaluate_df_critical_he(
        engine,
        critical_operator,
        df_result["DF"],
        df_coefficients,
        df_lower,
        df_upper,
    )
    critical_seconds = perf_counter() - started
    critical_bootstraps = critical_operator.bootstrap_count()
    encrypted_t = float(io_operator.decrypt(df_result["T"])[0])
    encrypted_df = float(io_operator.decrypt(df_result["DF"])[0])
    encrypted_critical = float(io_operator.decrypt(df_critical_ct)[0])
    rows.append(
        {
            "path": f"DF (group-aware [{df_lower:g},{df_upper:g}])",
            "plaintext_parameter": plaintext_df,
            "encrypted_parameter": encrypted_df,
            "parameter_absolute_error": abs(encrypted_df - plaintext_df),
            "parameter_relative_error": abs(encrypted_df - plaintext_df)
            / plaintext_df,
            "t_absolute_error": abs(encrypted_t - plaintext_t),
            "critical_absolute_error": abs(
                encrypted_critical - plaintext_critical
            ),
            "arithmetic_seconds": arithmetic_seconds,
            "critical_seconds": critical_seconds,
            "total_seconds": arithmetic_seconds + critical_seconds,
            "arithmetic_bootstraps": arithmetic_bootstraps,
            "critical_bootstraps": critical_bootstraps,
            "total_bootstraps": arithmetic_bootstraps + critical_bootstraps,
        }
    )

    result = pd.DataFrame(rows)
    result.to_csv(HE_RESULT_PATH, index=False)
    print("\n========== Encrypted DF versus InvDF ==========")
    print(f"dataset                : insurance charges by smoker ({n1}, {n2})")
    print(f"alpha / degree         : {ALPHA} / 15")
    print(f"plaintext T / df       : {plaintext_t} / {plaintext_df}")
    print(f"plaintext critical     : {plaintext_critical}")
    print(result.to_string(index=False))
    print(f"\nCSV written to: {HE_RESULT_PATH}")
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--plaintext-only",
        action="store_true",
        help="skip the HEaaN comparison and run only the coefficient sweep",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_approximation_sweep()
    if not args.plaintext_only:
        run_he_comparison()


if __name__ == "__main__":
    main()
