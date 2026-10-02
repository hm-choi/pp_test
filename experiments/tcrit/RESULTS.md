# Section 4.4: critical-value approximation (ApproxTCrit)

Exact reference: SciPy `t.ppf` (two-sided). df in [1, 2000] (dense log grid + all integers). Coefficients:
`coeffs/t_critical_coeffs.py` (generated in `coeffs/approx.ipynb`, section 4-1) on 1/df in [1/2000, 1].

## Plaintext MaxRE (alpha = 0.05)

| target | degree | levels | 1/df (ours) | 1/df extrapolation df in (2000,4096] | direct df on [1,2000] |
|---|---|---|---|---|---|
| t | 7 | 4 | 2.60e-06 | 2.68e-06 | 8.34e-01 |
| t | 15 | 5 | 4.65e-09 | 5.95e-10 | 7.97e-01 |
| t | 31 | 6 | 4.16e-09 | 2.34e-10 | 6.53e-01 |
| t | 63 | 7 | 4.97e-09 | 8.78e-10 | 3.17e-01 |
| t | 127 | 8 | 4.19e-09 | 1.35e-09 | 5.22e-02 |
| t2 | 7 | 4 | 1.01e-03 | 1.05e-03 | 9.72e-01 |
| t2 | 15 | 5 | 9.29e-09 | 9.45e-10 | 9.59e-01 |
| t2 | 31 | 6 | 8.34e-09 | 5.27e-10 | 8.86e-01 |
| t2 | 63 | 7 | 9.94e-09 | 1.60e-09 | 8.68e-01 |
| t2 | 127 | 8 | 8.38e-09 | 2.00e-09 | 2.83e-01 |

Lookup tables (MaxRE against the exact quantile). Plaintext software computes the quantile for real-valued df directly;
lookup tables are used when the quantile cannot be evaluated (printed tables, table-lookup based secure testing) and need the
real-valued Welch df mapped to a tabulated df. Textbook grid: df of `distribution_table/t_distribution_table.csv`
(1..30, 40, 50, 60, 80, 100, 120, 1000).

| target | table / mapping | df >= 1 | df >= 2 | df >= 10 | df >= 30 |
|---|---|---|---|---|---|
| t | textbook df grid, round down | 1.95e+00 | 3.52e-01 | 1.23e-02 | 1.05e-02 |
| t | all integer df, floor | 1.95e+00 | 3.52e-01 | 1.23e-02 | 1.34e-03 |
| t | all integer df, round | 1.11e+00 | 2.04e-01 | 6.45e-03 | 6.73e-04 |
| t | 1/df Chebyshev degree 15 (ours) | 4.65e-09 |  |  |  |
| t2 | textbook df grid, round down | 7.72e+00 | 8.27e-01 | 2.48e-02 | 2.11e-02 |
| t2 | all integer df, floor | 7.72e+00 | 8.27e-01 | 2.48e-02 | 2.68e-03 |
| t2 | all integer df, round | 3.46e+00 | 4.49e-01 | 1.29e-02 | 1.35e-03 |
| t2 | 1/df Chebyshev degree 15 (ours) | 9.29e-09 |  |  |  |

Other alphas (0.01, 0.025, 0.1) and the group-size-aware direct-df comparison: `results/plaintext.csv`.

## HE evaluation (10 reps, input level 12)

| target | degree | levels | time (s) | MaxRE vs exact | MaxAE vs plaintext poly |
|---|---|---|---|---|---|
| t | 7 | 12->8 | 0.204 +- 0.005 | 2.82e-06 | 3.89e-06 |
| t | 15 | 12->7 | 0.321 +- 0.026 | 4.30e-07 | 4.42e-06 |
| t | 31 | 12->6 | 0.511 +- 0.013 | 3.95e-07 | 4.59e-06 |
| t | 63 | 12->5 | 0.639 +- 0.012 | 3.64e-07 | 3.62e-06 |
| t | 127 | 12->4 | 0.884 +- 0.020 | 4.23e-07 | 4.85e-06 |
| t2 | 7 | 12->8 | 0.212 +- 0.011 | 1.01e-03 | 8.53e-05 |
| t2 | 15 | 12->7 | 0.347 +- 0.031 | 6.88e-07 | 1.04e-04 |
| t2 | 31 | 12->6 | 0.516 +- 0.036 | 7.04e-07 | 1.01e-04 |
| t2 | 63 | 12->5 | 0.756 +- 0.012 | 6.85e-07 | 1.08e-04 |
| t2 | 127 | 12->4 | 0.891 +- 0.022 | 7.67e-07 | 9.08e-05 |

Degree 15 (2^4-1) reaches the CKKS noise floor (HE MaxRE ~4e-7); higher degrees add time and levels without accuracy gain.
