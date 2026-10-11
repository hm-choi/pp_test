# Accuracy and decision agreement

Script `run_accuracy.py`; raw results `results/accuracy_raw.csv` (one row per case, test, repetition and alpha),
summary `results/accuracy_summary.csv`.

HEaaN FGb (CPU), log_slots = 15. Two-sided alpha in {0.001, 0.01, 0.025, 0.05, 0.1}, 2 repetitions. For every repetition
the encrypted statistics are computed once and the decision is evaluated separately for every alpha; copies of the
intermediate ciphertexts are decrypted for the measurement. Each comparison has the public input interval [L, U]
(minimum and maximum of the two groups) and R = U - L.

| case | n1 / n2 | L | U | R | Welch df | V / V_max |
|---|---|---|---|---|---|---|
| insurance_charges_smoker | 274 / 1064 | 1.122 | 63.77 | 62.65 | 311.9 | 0.115 |
| adult_edu_income | 1024 / 1024 | 1 | 16 | 15 | 2037.8 | 0.111 |
| adult_age_gender | 1024 / 1024 | 17 | 90 | 73 | 2035.0 | 0.145 |
| heart_chol_disease | 139 / 164 | 126 | 564 | 438 | 298.6 | 0.055 |
| diabetes_labs_readmitted | 11357 / 54864 | 1 | 132 | 131 | 16691.2 | 0.087 |
| credit_limit_default | 6636 / 23364 | 1 | 100 | 99 | 11982.1 | 0.058 |
| bank_age_subscribed | 5289 / 39922 | 18 | 95 | 77 | 6109.2 | 0.117 |

Decision scores and public bounds B (the step function is applied to s / B):
- Welch: s = T^2 - c^2, B = max(R^2 / (V_BAR_MIN V_max), c_max^2) with V_BAR_MIN = 0.01; c^2 from the
  degree-15 polynomial of u = 1/df on [0, 1] (ApproxTCrit), c_max^2 at u = 1/(min(n1, n2) - 1); invSqrt: HE-DAP
  selection for V / V_max in [V_BAR_MIN, 1] (`../invsqrt/RESULTS.md`).
- F: s = (s1^2 - F_L s2^2)(s1^2 - F_U s2^2), each factor divided by m_j = max(s1max^2, F_j s2max^2),
  s_imax^2 = n_i R^2 / (4 (n_i - 1)).
- Z: s = Z^2 - z^2, B = max(R^2 / V_public, z^2); the sample variances are used as the public variances.

## Decision agreement

A combination (case, test, alpha) matches when, in every repetition, the step output (> 0.5), the sign of the decrypted
score and (Welch, Z) the decision from the HE-derived p-value all equal the SciPy decision.

| test | combinations | match | reject | fail to reject |
|---|---|---|---|---|
| Welch | 35 | 35 | 30 | 5 |
| F | 35 | 35 | 25 | 10 |
| Z | 35 | 35 | 30 | 5 |

Fail to reject: Welch heart_chol_disease (alpha 0.001, 0.01, 0.025, 0.05, 0.1); F adult_age_gender (alpha 0.001, 0.01); F adult_edu_income (alpha 0.001, 0.01, 0.025); F heart_chol_disease (alpha 0.001, 0.01, 0.025, 0.05, 0.1); Z heart_chol_disease (alpha 0.001, 0.01, 0.025, 0.05, 0.1).

## Intermediate-value accuracy

Mean (standard deviation) and maximum over cases, alphas and repetitions. Statistic, 1/df and critical value c^2:
relative error; score: absolute error divided by the public bound, |s - s_HE| / B; p-value: absolute error.

| test | statistic | 1/df | critical value c^2 | score / B | p-value | plaintext \|s\| / B (min - max) |
|---|---|---|---|---|---|---|
| Welch | 1.96e-07 (2.47e-07), max 7.80e-07 | 1.75e-06 (1.98e-06), max 6.08e-06 | 1.31e-04 (2.92e-04), max 9.17e-04 | 6.16e-09 (1.81e-08), max 9.42e-08 | 4.35e-10 (1.35e-09), max 5.63e-09 | 2.82e-06 - 1.23e-02 |
| F | -- | -- | public | 1.41e-10 (9.97e-11), max 3.60e-10 | -- | 3.79e-06 - 8.31e-03 |
| Z | 5.96e-09 (1.08e-08), max 3.18e-08 | -- | public | 4.73e-11 (8.49e-11), max 2.62e-10 | 1.07e-10 (2.95e-10), max 1.11e-09 | 4.31e-05 - 1.42e-01 |

The F test computes no F statistic or p-value; the F and Z critical values are public.
