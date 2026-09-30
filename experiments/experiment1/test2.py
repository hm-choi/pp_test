"""Debug the unified encrypted Welch decision path on real and boundary data."""

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
BOUNDARY_GROUP_SIZE = 128


def plaintext_welch(x1: np.ndarray, x2: np.ndarray) -> dict[str, float]:
    var1, var2 = float(np.var(x1, ddof=1)), float(np.var(x2, ddof=1))
    v = var1 / len(x1) + var2 / len(x2)
    df = v * v / (
        ((var1 / len(x1)) ** 2) / (len(x1) - 1)
        + ((var2 / len(x2)) ** 2) / (len(x2) - 1)
    )
    t = float((np.mean(x1) - np.mean(x2)) / np.sqrt(v))
    critical = float(stats.t.ppf(1.0 - ALPHA / 2.0, df))
    return {
        "T": t,
        "INV_DF": 1.0 / df,
        "critical_value": critical,
        "score": t * t - critical * critical,
    }


def run_case(
    title: str,
    x1: np.ndarray,
    x2: np.ndarray,
    score_bound: float,
    operator: HEOperator,
    hypothesis_test: HEHypothesisTesting,
) -> None:
    reference = plaintext_welch(x1, x2)
    r = float(max(np.max(x1), np.max(x2)))
    hypothesis_test.reset_bootstrap_count()
    started_at = perf_counter()
    encrypted = hypothesis_test.HE_Welch_T_Test(
        operator.encrypt(x1), operator.encrypt(x2), R=r,
        alpha=ALPHA, critical_degree=CRITICAL_VALUE_DEGREE,
        score_bound=score_bound, debug=True,
    )
    runtime_seconds = perf_counter() - started_at
    t = float(operator.decrypt(encrypted["T"])[0])
    inv_df = float(operator.decrypt(encrypted["INV_DF"])[0])
    critical = float(operator.decrypt(encrypted["critical_value"])[0])
    score = float(operator.decrypt(encrypted["score"])[0])
    step = float(operator.decrypt(encrypted["step"])[0])

    print(f"\n========== {title} ==========")
    print(f"n1, n2                 : {len(x1)}, {len(x2)}")
    print(f"alpha, degree          : {ALPHA}, {CRITICAL_VALUE_DEGREE}")
    print(f"score bound            : {score_bound}")
    print(f"plaintext T             : {reference['T']}")
    print(f"encrypted T             : {t}")
    print(f"T absolute error        : {abs(t - reference['T'])}")
    print(f"plaintext InvDF         : {reference['INV_DF']}")
    print(f"encrypted InvDF         : {inv_df}")
    print(f"InvDF absolute error    : {abs(inv_df - reference['INV_DF'])}")
    print(f"plaintext critical t    : {reference['critical_value']}")
    print(f"encrypted critical t    : {critical}")
    print(f"plaintext score         : {reference['score']}")
    print(f"encrypted score         : {score}")
    print(f"encrypted step(score)   : {step}")
    print(f"score absolute error    : {abs(score - reference['score'])}")
    print(f"runtime (seconds)       : {runtime_seconds}")
    print(f"bootstraps              : {hypothesis_test.bootstrap_count()}")
    print(f"score-sign match        : {(reference['score'] > 0.0) == (score > 0.0)}")
    print(f"step decision match     : {(reference['score'] > 0.0) == (step >= 0.5)}")


def main() -> None:
    insurance = pd.read_csv(PROJECT_ROOT / "datasets" / "insurance.csv")
    x1 = insurance.loc[insurance["smoker"] == "yes", "charges"].to_numpy(
        dtype=np.float64
    ) / 1_000.0
    x2 = insurance.loc[insurance["smoker"] == "no", "charges"].to_numpy(
        dtype=np.float64
    ) / 1_000.0
    noise = 5.0 * np.linspace(-1.0, 1.0, BOUNDARY_GROUP_SIZE)
    variance = float(np.var(noise, ddof=1))
    df = float(2 * BOUNDARY_GROUP_SIZE - 2)
    critical = float(stats.t.ppf(1.0 - ALPHA / 2.0, df))
    standard_error = np.sqrt(2.0 * variance / BOUNDARY_GROUP_SIZE)
    mean_gap = (critical + 0.001) * standard_error
    boundary_x1 = 10.0 + noise + mean_gap / 2.0
    boundary_x2 = 10.0 + noise - mean_gap / 2.0

    engine = HEEngine(warmup_bootstrap=False)
    operator = HEOperator(engine)
    hypothesis_test = HEHypothesisTesting(engine)
    run_case("Insurance Welch debug", x1, x2, 1_200.0, operator, hypothesis_test)
    run_case("Boundary decision debug", boundary_x1, boundary_x2, 0.01, operator, hypothesis_test)


if __name__ == "__main__":
    main()
