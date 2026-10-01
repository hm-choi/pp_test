"""Dataset comparisons shared by the PP-TEST experiments.

Each case returns two groups on [0, R] with the public values used by
HEHypothesisTesting (group sizes n1, n2 and the public range bound R).
Adult groups use the first 1,024 records of each group, as in the earlier
experiments (experiments/experiment1, experiment2).
"""

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATASET_DIR = PROJECT_ROOT / "datasets"

ADULT_SAMPLES_PER_GROUP = 1_024


def _select(df, group_column, group_value, value_column, scale=1.0, limit=None):

    x = df.loc[df[group_column] == group_value, value_column].to_numpy(dtype=np.float64) / scale
    return x[:limit] if limit else x


def load_cases():

    insurance = pd.read_csv(DATASET_DIR / "insurance.csv")
    adult = pd.read_csv(DATASET_DIR / "adult_dataset.csv")
    n = ADULT_SAMPLES_PER_GROUP

    raw = {
        "insurance_charges_smoker": (
            "Insurance charges/1000: smoker yes vs no",
            _select(insurance, "smoker", "yes", "charges", 1_000.0),
            _select(insurance, "smoker", "no", "charges", 1_000.0),
        ),
        "adult_edu_income": (
            "Adult educational-num: >50K vs <=50K",
            _select(adult, "income", ">50K", "educational-num", limit=n),
            _select(adult, "income", "<=50K", "educational-num", limit=n),
        ),
        "adult_age_gender": (
            "Adult age: Female vs Male",
            _select(adult, "gender", "Female", "age", limit=n),
            _select(adult, "gender", "Male", "age", limit=n),
        ),
    }

    cases = {}

    for name, (title, x1, x2) in raw.items():
        n1, n2 = len(x1), len(x2)
        R = float(max(np.max(x1), np.max(x2)))

        cases[name] = {
            "name": name,
            "title": title,
            "x1": x1,
            "x2": x2,
            "n1": n1,
            "n2": n2,
            "R": R,
            # Public bounds used by HEHypothesisTesting.
            "welch_v_max": (R ** 2 / 4.0) * (1.0 / (n1 - 1) + 1.0 / (n2 - 1)),
            "f_var2_max": n2 / (n2 - 1) * R ** 2 / 4.0,
        }

    return cases
