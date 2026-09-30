from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from utils.welch import (
    welch_t_test_from_samples,
    welch_t_test_from_sufficient_statistics,
)
from pp_test.engine.HEengine import HEEngine
from operators.hypothesis_testing import HEHypothesisTesting
from operators.operator import HEOperator


DATASET_PATH = PROJECT_ROOT / "datasets" / "insurance.csv"
ALPHA = 0.05
INPUT_SCALE = 1_000.0


def print_result(title: str, result: dict[str, float | int]) -> None:
    print(f"\n========== {title} ==========")
    for key, value in result.items():
        print(f"{key:10s}: {value}")


def decrypt_scalar(operator: HEOperator, ciphertext) -> float:
    """Decrypt the replicated scalar stored in the first slot."""
    return float(operator.decrypt(ciphertext)[0])


def run_he_welch_test(
    x1: np.ndarray, x2: np.ndarray, input_upper_bound: float
) -> tuple[float, float, float]:
    """Run encrypted Welch T and df computation and return decrypted scalars."""
    engine = HEEngine(warmup_bootstrap=False)
    operator = HEOperator(engine)
    hypothesis_test = HEHypothesisTesting(engine)

    encrypted_x1 = operator.encrypt(x1)
    encrypted_x2 = operator.encrypt(x2)

    started_at = perf_counter()
    encrypted_result = hypothesis_test.HE_Welch_T_Test_With_DF(
        encrypted_x1,
        encrypted_x2,
        R=input_upper_bound,
        return_debug=True,
    )
    runtime_seconds = perf_counter() - started_at

    return (
        decrypt_scalar(operator, encrypted_result["T"]),
        decrypt_scalar(operator, encrypted_result["DF"]),
        runtime_seconds,
    )


def main() -> None:
    df = pd.read_csv(DATASET_PATH)
    # Keep all HE inputs in the public interval [0, R] with modest magnitude.
    x1 = df.loc[df["smoker"] == "yes", "charges"].to_numpy(dtype=np.float64) / INPUT_SCALE
    x2 = df.loc[df["smoker"] == "no", "charges"].to_numpy(dtype=np.float64) / INPUT_SCALE
    input_upper_bound = float(max(np.max(x1), np.max(x2)))

    print("Group sizes")
    print(f"smoker = yes: n1 = {len(x1)}")
    print(f"smoker = no : n2 = {len(x2)}")

    reference = welch_t_test_from_samples(x1, x2)
    he_style = welch_t_test_from_sufficient_statistics(x1, x2)
    scipy_result = stats.ttest_ind(x1, x2, equal_var=False)

    print_result("Plaintext Welch t-test", reference)
    print_result("HE sufficient-statistics reference", he_style)

    print("\n========== Reference checks ==========")
    print("SciPy t-statistic:", scipy_result.statistic)
    print("SciPy p-value    :", scipy_result.pvalue)
    print("t diff, SciPy    :", abs(reference["t_stat"] - scipy_result.statistic))
    print("p diff, SciPy    :", abs(reference["p_value"] - scipy_result.pvalue))
    print("t diff, HE-style :", abs(reference["t_stat"] - he_style["t_stat"]))
    print("df diff, HE-style:", abs(reference["df"] - he_style["df"]))
    print("p diff, HE-style :", abs(reference["p_value"] - he_style["p_value"]))

    he_t_stat, he_df, he_runtime = run_he_welch_test(x1, x2, input_upper_bound)
    he_p_value = 2.0 * stats.t.sf(abs(he_t_stat), he_df)
    reference_decision = reference["p_value"] < ALPHA
    he_decision = he_p_value < ALPHA

    print("\n========== HEaaN Welch test ==========")
    print("input scale       :", INPUT_SCALE)
    print("R                 :", input_upper_bound)
    print("runtime (seconds) :", he_runtime)
    print("encrypted T       :", he_t_stat)
    print("encrypted df      :", he_df)
    print("T absolute error  :", abs(he_t_stat - reference["t_stat"]))
    print("df absolute error :", abs(he_df - reference["df"]))
    print("decision match    :", he_decision == reference_decision)

    print("\n========== Conclusion ==========")
    print(f"alpha = {ALPHA}")
    if reference_decision:
        print("Reject H0: smoker groups have significantly different mean charges.")
    else:
        print("Fail to reject H0: the mean charges do not differ significantly.")


if __name__ == "__main__":
    main()
