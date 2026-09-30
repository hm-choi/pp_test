"""Simple real-data runs of the unified encrypted Welch decision API."""

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
CRITICAL_VALUE_DEGREE = 15
ADULT_SAMPLES_PER_GROUP = 1_024


def select_group(dataframe, group_column, group_value, value_column, scale=1.0):
    return dataframe.loc[
        dataframe[group_column] == group_value, value_column
    ].to_numpy(dtype=np.float64) / scale


def plaintext_welch(x1: np.ndarray, x2: np.ndarray) -> tuple[float, float, float]:
    var1, var2 = np.var(x1, ddof=1), np.var(x2, ddof=1)
    v = var1 / len(x1) + var2 / len(x2)
    t = float((np.mean(x1) - np.mean(x2)) / np.sqrt(v))
    df = float(v * v / (
        ((var1 / len(x1)) ** 2) / (len(x1) - 1)
        + ((var2 / len(x2)) ** 2) / (len(x2) - 1)
    ))
    critical = float(stats.t.ppf(1.0 - ALPHA / 2.0, df))
    return t, 1.0 / df, t * t - critical * critical


def main() -> None:
    insurance = pd.read_csv(PROJECT_ROOT / "datasets" / "insurance.csv")
    adult = pd.read_csv(PROJECT_ROOT / "datasets" / "adult_dataset.csv")
    cases = [
        (
            "Insurance: charges by smoker status",
            select_group(insurance, "smoker", "yes", "charges", 1_000.0),
            select_group(insurance, "smoker", "no", "charges", 1_000.0),
            1_200.0,
        ),
        (
            "Adult: educational-num by income",
            select_group(adult, "income", ">50K", "educational-num")[:ADULT_SAMPLES_PER_GROUP],
            select_group(adult, "income", "<=50K", "educational-num")[:ADULT_SAMPLES_PER_GROUP],
            400.0,
        ),
        (
            "Adult: age by gender",
            select_group(adult, "gender", "Female", "age")[:ADULT_SAMPLES_PER_GROUP],
            select_group(adult, "gender", "Male", "age")[:ADULT_SAMPLES_PER_GROUP],
            100.0,
        ),
    ]
    engine = HEEngine(warmup_bootstrap=True)
    operator = HEOperator(engine)
    hypothesis_test = HEHypothesisTesting(engine)

    for title, x1, x2, score_bound in cases:
        plaintext_t, plaintext_inv_df, plaintext_score = plaintext_welch(x1, x2)
        hypothesis_test.reset_bootstrap_count()
        started_at = perf_counter()
        encrypted = hypothesis_test.HE_Welch_T_Test(
            operator.encrypt(x1), operator.encrypt(x2),
            R=float(max(np.max(x1), np.max(x2))),
            alpha=ALPHA, critical_degree=CRITICAL_VALUE_DEGREE,
            score_bound=score_bound, debug=True,
        )
        runtime = perf_counter() - started_at
        t = float(operator.decrypt(encrypted["T"])[0])
        inv_df = float(operator.decrypt(encrypted["INV_DF"])[0])
        score = float(operator.decrypt(encrypted["score"])[0])
        step = float(operator.decrypt(encrypted["step"])[0])
        print(f"\n========== {title} ==========")
        print(f"plaintext T, InvDF     : {plaintext_t}, {plaintext_inv_df}")
        print(f"encrypted T, InvDF     : {t}, {inv_df}")
        print(f"plaintext score         : {plaintext_score}")
        print(f"encrypted score, step   : {score}, {step}")
        print(f"T absolute error        : {abs(t - plaintext_t)}")
        print(f"InvDF absolute error    : {abs(inv_df - plaintext_inv_df)}")
        print(f"score absolute error    : {abs(score - plaintext_score)}")
        print(f"step decision match     : {(plaintext_score > 0.0) == (step >= 0.5)}")
        print(f"runtime, bootstraps     : {runtime}, {hypothesis_test.bootstrap_count()}")


if __name__ == "__main__":
    main()
