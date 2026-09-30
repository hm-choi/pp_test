"""Boundary stress benchmark for encrypted Welch critical-value polynomials."""

from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from pp_test.engine.HEengine import HEEngine
from pp_test.operators.HEhypo_test import HEHypothesisTesting
from operators.operator import HEOperator


ALPHA = 0.05
DEGREES = (7, 15, 31, 63, 127)
MARGINS = (-0.01, -0.005, -0.002, -0.001, 0.001, 0.002, 0.005, 0.01)
GROUP_SIZE = 128
# Public normalization bound for this synthetic stress family.  With the
# configured margins, the plaintext score lies in approximately [-0.04, 0.04].
# The sign approximation requires score / SCORE_BOUND to stay in [-1, 1].
SCORE_BOUND = 0.05
RESULT_PATH = Path(__file__).with_name("test3_results.csv")


def main() -> None:
    noise = 5.0 * np.linspace(-1.0, 1.0, GROUP_SIZE)
    variance = float(np.var(noise, ddof=1))
    df_target = float(2 * GROUP_SIZE - 2)
    critical_target = float(stats.t.ppf(1.0 - ALPHA / 2.0, df_target))
    standard_error = np.sqrt(2.0 * variance / GROUP_SIZE)

    engine = HEEngine(warmup_bootstrap=False)
    operator = HEOperator(engine)
    hypothesis_test = HEHypothesisTesting(engine)
    rows = []

    for margin in MARGINS:
        mean_gap = (critical_target + margin) * standard_error
        x1 = 10.0 + noise + mean_gap / 2.0
        x2 = 10.0 + noise - mean_gap / 2.0

        var1 = float(np.var(x1, ddof=1))
        var2 = float(np.var(x2, ddof=1))
        v = var1 / GROUP_SIZE + var2 / GROUP_SIZE
        plaintext_t = float((np.mean(x1) - np.mean(x2)) / np.sqrt(v))
        plaintext_df = float(
            v**2
            / (
                ((var1 / GROUP_SIZE) ** 2) / (GROUP_SIZE - 1)
                + ((var2 / GROUP_SIZE) ** 2) / (GROUP_SIZE - 1)
            )
        )
        plaintext_critical = float(stats.t.ppf(1.0 - ALPHA / 2.0, plaintext_df))
        plaintext_score = plaintext_t**2 - plaintext_critical**2
        if abs(plaintext_score) >= SCORE_BOUND:
            raise ValueError(
                "Stress score exceeds SCORE_BOUND; enlarge the public "
                "normalization bound before evaluating the sign polynomial."
            )
        r = float(max(np.max(x1), np.max(x2)))

        hypothesis_test.reset_bootstrap_count()
        started_at = perf_counter()
        encrypted = hypothesis_test.HE_Welch_T_Test(
            operator.encrypt(x1),
            operator.encrypt(x2),
            R=r,
            alpha=ALPHA,
            critical_degree=15,
            score_bound=SCORE_BOUND,
            debug=True,
        )
        welch_seconds = perf_counter() - started_at
        welch_bootstraps = hypothesis_test.bootstrap_count()
        encrypted_t = float(operator.decrypt(encrypted["T"])[0])
        encrypted_inv_df = float(operator.decrypt(encrypted["INV_DF"])[0])
        encrypted_t_squared = encrypted["T_SQUARED"]
        encrypted_step_decision = float(operator.decrypt(encrypted["step"])[0]) > 0.5
        plaintext_decision = plaintext_score > 0.0

        for degree in DEGREES:
            hypothesis_test.reset_bootstrap_count()
            started_at = perf_counter()
            critical_ct = hypothesis_test._critical_value_from_inv_df(
                encrypted["INV_DF"], ALPHA, degree
            )
            score_ct = operator.sub(
                encrypted_t_squared,
                operator.mult(critical_ct, critical_ct),
            )
            critical_seconds = perf_counter() - started_at
            critical_bootstraps = hypothesis_test.bootstrap_count()
            encrypted_critical = float(operator.decrypt(critical_ct)[0])
            encrypted_score = float(operator.decrypt(score_ct)[0])

            rows.append(
                {
                    "degree": degree,
                    "margin": margin,
                    "score_bound": SCORE_BOUND,
                    "plaintext_t": plaintext_t,
                    "encrypted_t_error": abs(encrypted_t - plaintext_t),
                    "plaintext_df": plaintext_df,
                    "encrypted_inv_df_error": abs(
                        encrypted_inv_df - 1.0 / plaintext_df
                    ),
                    "critical_error": abs(encrypted_critical - plaintext_critical),
                    "plaintext_score": plaintext_score,
                    "encrypted_score": encrypted_score,
                    "score_error": abs(encrypted_score - plaintext_score),
                    # This sweep measures critical-value/score arithmetic.
                    # The integrated HE step is evaluated once per margin
                    # (degree 15) above; evaluating it for every degree would
                    # add the same costly sign-polynomial work five times.
                    "score_sign_match": plaintext_decision
                    == (encrypted_score > 0.0),
                    "step_decision_match": (
                        plaintext_decision == encrypted_step_decision
                        if degree == 15
                        else None
                    ),
                    "welch_seconds": welch_seconds,
                    "critical_seconds": critical_seconds,
                    "welch_bootstraps": welch_bootstraps,
                    "critical_bootstraps": critical_bootstraps,
                    "score_level": score_ct.level(),
                }
            )

    result = pd.DataFrame(rows)
    result.to_csv(RESULT_PATH, index=False)
    print("\n========== Degree stress-test summary ==========")
    print(result.to_string(index=False))
    print(f"\nCSV written to: {RESULT_PATH}")
    print("\n========== Degree summary ==========")
    print(
        result.groupby("degree", as_index=False).agg(
            score_sign_match_rate=("score_sign_match", "mean"),
            step_match_rate=("step_decision_match", "mean"),
            max_critical_error=("critical_error", "max"),
            max_score_error=("score_error", "max"),
            mean_critical_seconds=("critical_seconds", "mean"),
            total_critical_bootstraps=("critical_bootstraps", "sum"),
            minimum_score_level=("score_level", "min"),
        ).to_string(index=False)
    )


if __name__ == "__main__":
    main()
