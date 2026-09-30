"""Simple real-data validation for the existing HE Z-test and F-test APIs.

Z-test uses public/known group variances. F-test returns an encrypted
variance ratio; its two-sided p-value and decision are computed only after
decryption, matching the scope of the current HE_F_Test implementation.
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
from pp_test.operators.HEhypo_test import HEHypothesisTesting
from operators.operator import HEOperator


ALPHA = 0.05
SAMPLES_PER_GROUP = 1_024
# Public bound for this Adult-age Z-test experiment.  The measured plaintext
# score is checked below before it is used to normalize the encrypted score.
Z_SCORE_BOUND = 100.0
F_SCORE_BOUND = 0.05
RESULT_PATH = Path(__file__).with_name("test_results.csv")


def two_sided_f_p_value(f_statistic: float, df1: int, df2: int) -> float:
    """Two-sided F-test p-value for the fixed orientation F=s1^2/s2^2."""
    cdf = float(stats.f.cdf(f_statistic, df1, df2))
    return min(1.0, 2.0 * min(cdf, 1.0 - cdf))


def main() -> None:
    adult = pd.read_csv(PROJECT_ROOT / "datasets" / "adult_dataset.csv")
    female_age = adult.loc[adult["gender"] == "Female", "age"].to_numpy(
        dtype=np.float64
    )[:SAMPLES_PER_GROUP]
    male_age = adult.loc[adult["gender"] == "Male", "age"].to_numpy(
        dtype=np.float64
    )[:SAMPLES_PER_GROUP]
    n1, n2 = len(female_age), len(male_age)
    if n1 != SAMPLES_PER_GROUP or n2 != SAMPLES_PER_GROUP:
        raise ValueError("Adult dataset does not contain enough samples per group.")

    # Experimental Z-test setting: these values are supplied as public,
    # known variances. They are shown explicitly to avoid treating an
    # encrypted sample variance as a known population variance.
    known_sigma1_sq = float(np.var(female_age, ddof=1))
    known_sigma2_sq = float(np.var(male_age, ddof=1))
    z_standard_error = np.sqrt(known_sigma1_sq / n1 + known_sigma2_sq / n2)
    plaintext_z = float((np.mean(female_age) - np.mean(male_age)) / z_standard_error)
    plaintext_z_p = float(2.0 * stats.norm.sf(abs(plaintext_z)))
    plaintext_z_critical = float(stats.norm.ppf(1.0 - ALPHA / 2.0))
    plaintext_z_score = plaintext_z**2 - plaintext_z_critical**2
    if abs(plaintext_z_score) >= Z_SCORE_BOUND:
        raise ValueError("Increase the public Z_SCORE_BOUND for this experiment.")

    plaintext_var1 = float(np.var(female_age, ddof=1))
    plaintext_var2 = float(np.var(male_age, ddof=1))
    plaintext_f = plaintext_var1 / plaintext_var2
    df1, df2 = n1 - 1, n2 - 1
    plaintext_f_p = two_sided_f_p_value(plaintext_f, df1, df2)
    lower_f_critical = float(stats.f.ppf(ALPHA / 2.0, df1, df2))
    upper_f_critical = float(stats.f.ppf(1.0 - ALPHA / 2.0, df1, df2))
    plaintext_f_score = (plaintext_f - lower_f_critical) * (
        plaintext_f - upper_f_critical
    )
    if abs(plaintext_f_score) >= F_SCORE_BOUND:
        raise ValueError("Increase the public F_SCORE_BOUND for this experiment.")

    engine = HEEngine(warmup_bootstrap=True)
    operator = HEOperator(engine)
    hypothesis_test = HEHypothesisTesting(engine)
    encrypted_female_age = operator.encrypt(female_age)
    encrypted_male_age = operator.encrypt(male_age)
    rows = []

    started_at = perf_counter()
    encrypted_z = hypothesis_test.HE_Z_Test(
        encrypted_female_age,
        encrypted_male_age,
        known_sigma1_sq,
        known_sigma2_sq,
        return_debug=True,
        alpha=ALPHA,
        score_bound=Z_SCORE_BOUND,
    )
    z_runtime = perf_counter() - started_at
    z_value = float(operator.decrypt(encrypted_z["Z"])[0])
    z_score = float(operator.decrypt(encrypted_z["score"])[0])
    z_step = float(operator.decrypt(encrypted_z["step"])[0])
    z_p_value = float(2.0 * stats.norm.sf(abs(z_value)))
    z_match = (plaintext_z_p < ALPHA) == (z_step >= 0.5)
    rows.append(
        {
            "test": "Z",
            "comparison": "Adult age: Female vs Male",
            "plaintext_statistic": plaintext_z,
            "encrypted_statistic": z_value,
            "absolute_error": abs(z_value - plaintext_z),
            "plaintext_score": plaintext_z_score,
            "encrypted_score": z_score,
            "score_absolute_error": abs(z_score - plaintext_z_score),
            "encrypted_step": z_step,
            "plaintext_p_value": plaintext_z_p,
            "decrypted_statistic_p_value": z_p_value,
            "decision_match": z_match,
            "runtime_seconds": z_runtime,
            "bootstraps": hypothesis_test.bootstrap_count(),
        }
    )

    hypothesis_test.reset_bootstrap_count()
    started_at = perf_counter()
    encrypted_f = hypothesis_test.HE_F_Test(
        encrypted_female_age,
        encrypted_male_age,
        R=float(max(np.max(female_age), np.max(male_age))),
        return_debug=True,
        lower_critical=lower_f_critical,
        upper_critical=upper_f_critical,
        score_bound=F_SCORE_BOUND,
    )
    f_runtime = perf_counter() - started_at
    f_value = float(operator.decrypt(encrypted_f["F"])[0])
    f_score = float(operator.decrypt(encrypted_f["score"])[0])
    f_step = float(operator.decrypt(encrypted_f["step"])[0])
    f_p_value = two_sided_f_p_value(f_value, df1, df2)
    f_match = (plaintext_f_p < ALPHA) == (f_step >= 0.5)
    rows.append(
        {
            "test": "F",
            "comparison": "Adult age variance: Female vs Male",
            "plaintext_statistic": plaintext_f,
            "encrypted_statistic": f_value,
            "absolute_error": abs(f_value - plaintext_f),
            "plaintext_score": plaintext_f_score,
            "encrypted_score": f_score,
            "score_absolute_error": abs(f_score - plaintext_f_score),
            "encrypted_step": f_step,
            "plaintext_p_value": plaintext_f_p,
            "decrypted_statistic_p_value": f_p_value,
            "decision_match": f_match,
            "runtime_seconds": f_runtime,
            "bootstraps": hypothesis_test.bootstrap_count(),
        }
    )

    result = pd.DataFrame(rows)
    print("\n========== Experiment 2: Z-test and F-test ==========")
    print(f"alpha                 : {ALPHA}")
    print(f"group sizes           : {n1}, {n2}")
    print(f"public Z variances    : {known_sigma1_sq}, {known_sigma2_sq}")
    print(result.to_string(index=False))
    result.to_csv(RESULT_PATH, index=False)
    print(f"\nCSV written to: {RESULT_PATH}")


if __name__ == "__main__":
    main()
