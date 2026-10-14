# Critical-value approximation: method comparison

Script `run_critical_boundary.py`; raw results `results/plaintext.csv`, `results/he_raw.csv`.

Target c^2 = t_{1-alpha/2, df}^2 (two-sided), reference SciPy `t.ppf`; relative error. Global condition:
65,536 df points log-spaced over [3, 10^4] (two ciphertexts). Inverse-DF uses the stored coefficients fitted on
u = 1/df in [0, 1]; Direct-DF fits the df polynomial on [3, 10^4] and computes df = 1/u with invSqrt.

## Encrypted, alpha = 0.05, global [3, 10^4], 10 repetitions (GPU)

| degree | Direct-DF MaxRE | Direct-DF MRE | Inverse-DF MaxRE | Inverse-DF MRE |
|---|---|---|---|---|
| 7 | 6.09e-01 (1.45e-09) | 7.68e-02 (1.51e-09) | 1.00e-03 (3.35e-08) | 6.61e-04 (1.39e-10) |
| 15 | 5.76e-01 (2.25e-08) | 6.16e-02 (2.60e-09) | 5.33e-07 (1.18e-07) | 4.34e-08 (1.23e-10) |
| 31 | 4.67e-01 (3.05e-07) | 4.09e-02 (5.03e-09) | 5.30e-07 (7.39e-08) | 4.33e-08 (1.43e-10) |
| 63 | 2.37e-01 (2.79e-06) | 2.11e-02 (4.77e-09) | 5.42e-07 (6.27e-08) | 4.35e-08 (1.41e-10) |
| 127 | 4.58e-02 (1.21e-05) | 4.46e-03 (1.27e-08) | 5.50e-07 (5.87e-08) | 4.33e-08 (1.66e-10) |

## Plaintext, global [3, 10^4]

| alpha | degree | Direct-DF MaxRE | Inverse-DF MaxRE |
|---|---|---|---|
| 0.01 | 7 | 7.96e-01 | 1.50e-01 |
| 0.01 | 15 | 7.68e-01 | 1.34e-07 |
| 0.01 | 31 | 6.64e-01 | 1.64e-08 |
| 0.01 | 63 | 3.92e-01 | 1.19e-08 |
| 0.01 | 127 | 9.43e-02 | 8.48e-09 |
| 0.025 | 7 | 7.01e-01 | 1.02e-02 |
| 0.025 | 15 | 6.69e-01 | 7.17e-09 |
| 0.025 | 31 | 5.59e-01 | 9.78e-09 |
| 0.025 | 63 | 3.04e-01 | 1.03e-08 |
| 0.025 | 127 | 6.47e-02 | 9.26e-09 |
| 0.05 | 7 | 6.09e-01 | 1.00e-03 |
| 0.05 | 15 | 5.76e-01 | 9.20e-09 |
| 0.05 | 31 | 4.67e-01 | 8.38e-09 |
| 0.05 | 63 | 2.37e-01 | 9.08e-09 |
| 0.05 | 127 | 4.57e-02 | 8.80e-09 |
| 0.1 | 7 | 5.00e-01 | 6.64e-05 |
| 0.1 | 15 | 4.68e-01 | 1.08e-08 |
| 0.1 | 31 | 3.67e-01 | 1.42e-08 |
| 0.1 | 63 | 1.74e-01 | 1.60e-08 |
| 0.1 | 127 | 3.00e-02 | 1.44e-08 |
