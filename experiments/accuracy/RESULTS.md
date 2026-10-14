# Accuracy and decision agreement

Script `run_accuracy.py`; raw results `results/accuracy_raw.csv` (one row per case, test, repetition and alpha),
summary `results/accuracy_summary.csv`.

HEaaN FGb (GPU, NVIDIA Tesla V100), log_slots = 15. invSqrt domain V in [0.01, 1000] (degree-63 Chebyshev + 7 Newton iterations); ApproxTCrit degree 15. Two-sided alpha in {0.01, 0.025, 0.05, 0.1}, 10 repetitions. For every repetition the encrypted
statistics are computed once and the decision is evaluated separately for every alpha. Relative errors are
against the plaintext (SciPy) values; the critical value is compared as c^2.

| case | n1 / n2 | R | V |
|---|---|---|---|
| insurance_charges_smoker | 274 / 1064 | 62.65 | 0.52 |
| adult_edu_income | 1024 / 1024 | 15 | 0.0122 |
| adult_age_gender | 1024 / 1024 | 73 | 0.377 |
| diabetes_labs_readmitted | 11357 / 54864 | 131 | 0.0399 |
| credit_limit_default | 6636 / 23364 | 99 | 0.0275 |
| bank_age_subscribed | 5289 / 39922 | 77 | 0.037 |

## Relative errors, mean (standard deviation), 6 comparisons

| test | statistic | 1/df | critical c^2 | score |
|---|---|---|---|---|
| Welch | 1.36e-07 (2.93e-07) | 2.13e-06 (2.98e-06) | 1.19e-07 (4.35e-08) | 3.91e-06 (6.20e-06) |
| F | -- | -- | public | 2.65e-06 (4.97e-06) |
| Z | 6.50e-09 (1.19e-08) | -- | public | 1.71e-08 (3.22e-08) |

## Decision agreement (step output > 0.5 vs plaintext)

| test | runs | match | reject / fail configurations |
|---|---|---|---|
| welch | 240 | 240 | 24 / 0 |
| f | 240 | 240 | 21 / 3 |
| z | 240 | 240 | 24 / 0 |
