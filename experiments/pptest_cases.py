"""Dataset comparisons shared by the PP-TEST experiments.

Each case returns two groups on [0, R] with the public values used by
HEHypothesisTesting (group sizes n1, n2 and the public range bound R).
Adult groups use the first 1,024 records of each group; the other datasets
use every record of the two groups.

Sources (UCI Machine Learning Repository, CC BY 4.0), files in datasets/:
- Heart Disease, processed Cleveland data (heart_disease_processed_cleveland.data);
- Diabetes 130-US Hospitals for Years 1999-2008 (diabetes_130_us_hospitals.csv);
- Default of Credit Card Clients (default_of_credit_card_clients.csv, converted
  from the original .xls without changes);
- Bank Marketing, bank-full.csv;
- Wine Quality, winequality-red.csv / winequality-white.csv.
"""

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DATASET_DIR = PROJECT_ROOT / "datasets"

ADULT_SAMPLES_PER_GROUP = 1_024

HEART_COLUMNS = [
    "age", "sex", "cp", "trestbps", "chol", "fbs", "restecg",
    "thalach", "exang", "oldpeak", "slope", "ca", "thal", "num",
]


def _select(df, group_column, group_value, value_column, scale=1.0, limit=None):

    x = df.loc[df[group_column] == group_value, value_column].to_numpy(dtype=np.float64) / scale
    return x[:limit] if limit else x


def load_cases():

    insurance = pd.read_csv(DATASET_DIR / "insurance.csv")
    adult = pd.read_csv(DATASET_DIR / "adult_dataset.csv")
    heart = pd.read_csv(
        DATASET_DIR / "heart_disease_processed_cleveland.data",
        header=None, names=HEART_COLUMNS, na_values="?",
    )
    heart["disease"] = np.where(heart["num"] > 0, "yes", "no")
    diabetes = pd.read_csv(
        DATASET_DIR / "diabetes_130_us_hospitals.csv",
        usecols=["readmitted", "num_lab_procedures"],
    )
    # Shift to start at 0 (values 1..132).
    diabetes["num_lab_procedures"] = diabetes["num_lab_procedures"] - 1
    credit = pd.read_csv(DATASET_DIR / "default_of_credit_card_clients.csv")
    bank = pd.read_csv(DATASET_DIR / "bank-full.csv", sep=";")
    wine = pd.concat([
        pd.read_csv(DATASET_DIR / "winequality-red.csv", sep=";").assign(color="red"),
        pd.read_csv(DATASET_DIR / "winequality-white.csv", sep=";").assign(color="white"),
    ])
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
        "heart_chol_disease": (
            "Heart Disease (Cleveland) cholesterol: disease vs no disease",
            _select(heart, "disease", "yes", "chol"),
            _select(heart, "disease", "no", "chol"),
        ),
        "diabetes_labs_readmitted": (
            "Diabetes 130 num_lab_procedures - 1: readmitted <30 days vs not readmitted",
            _select(diabetes, "readmitted", "<30", "num_lab_procedures"),
            _select(diabetes, "readmitted", "NO", "num_lab_procedures"),
        ),
        "credit_limit_default": (
            "Credit card LIMIT_BAL/10^4: default vs no default",
            _select(credit, "default payment next month", 1, "LIMIT_BAL", 10_000.0),
            _select(credit, "default payment next month", 0, "LIMIT_BAL", 10_000.0),
        ),
        "bank_age_subscribed": (
            "Bank Marketing age: subscribed vs not subscribed",
            _select(bank, "y", "yes", "age"),
            _select(bank, "y", "no", "age"),
        ),
        "wine_alcohol_color": (
            "Wine Quality alcohol x 10: red vs white",
            _select(wine, "color", "red", "alcohol", 0.1),
            _select(wine, "color", "white", "alcohol", 0.1),
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
        }

    return cases
