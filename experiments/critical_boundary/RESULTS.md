# Critical-value approximation: method comparison

Target c^2 = t_{1-alpha/2, df}^2 (two-sided), reference SciPy `t.ppf`; relative error. Script `run_critical_boundary.py`;
raw results `results/plaintext.csv`, `results/he_raw.csv`.

Conditions (every method is evaluated on the same 65,536 df points per condition, log-spaced over the range plus
the Welch df of the cases):
- global: df in [1, 10^4];
- case range: df in [min(n1, n2) - 1, n1 + n2 - 2] for each dataset comparison (public; it always contains the Welch df):
  Insurance [273, 1336], Adult [1023, 2046], Heart [138, 301], Diabetes [11356, 66219], Credit [6635, 29998], Bank [5288, 45209] (Adult: both comparisons).

Methods:
- 1/df polynomial: Chebyshev interpolation in u = 1/df. Global condition: the stored coefficients
  (`coeffs/t_critical_coeffs.py`, fitted on u in [0, 1], i.e. every df >= 1, as used by PP-TEST); case ranges:
  refitted on [1/df_hi, 1/df_lo].
- df polynomial: Chebyshev interpolation in df, fitted on the condition range.

## Plaintext, alpha = 0.001

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-10000 |
|---|---|---|---|
| 1/df polynomial | 7 | 4.74e+01 (2.53e+01) | 8.48e-01 / 2.04e+01 / 3.95e+01 / 4.74e+01 |
| 1/df polynomial | 15 | 8.96e-04 (4.46e-04) | 1.40e-05 / 3.74e-04 / 7.63e-04 / 8.96e-04 |
| 1/df polynomial | 31 | 3.41e-08 (1.39e-08) | 9.81e-09 / 1.95e-08 / 2.87e-08 / 3.41e-08 |
| 1/df polynomial | 63 | 3.81e-07 (1.70e-07) | 8.69e-09 / 2.06e-07 / 3.14e-07 / 3.81e-07 |
| 1/df polynomial | 127 | 1.35e-07 (6.60e-08) | 9.70e-09 / 7.34e-08 / 1.13e-07 / 1.35e-07 |
| df polynomial | 7 | 1.00e+00 (2.51e-01) | 1.00e+00 / 9.88e-01 / 4.46e-01 / 1.25e-01 |
| df polynomial | 15 | 1.00e+00 (2.29e-01) | 1.00e+00 / 9.85e-01 / 3.20e-01 / 8.72e-02 |
| df polynomial | 31 | 1.00e+00 (2.34e-01) | 1.00e+00 / 9.67e-01 / 4.72e-01 / 3.27e-01 |
| df polynomial | 63 | 4.22e+00 (7.79e-01) | 9.99e-01 / 4.22e+00 / 2.75e+00 / 1.31e+00 |
| df polynomial | 127 | 4.03e+01 (7.02e+00) | 8.05e+00 / 4.03e+01 / 3.02e+01 / 1.42e+01 |

Case ranges (MaxRE):

| method | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank |
|---|---|---|---|---|---|---|---|
| 1/df polynomial | 7 | 3.10e-15 | 2.78e-15 | 2.83e-15 | 2.62e-15 | 2.95e-15 | 3.44e-15 |
| 1/df polynomial | 15 | 1.86e-14 | 1.88e-14 | 1.88e-14 | 1.89e-14 | 1.87e-14 | 1.85e-14 |
| 1/df polynomial | 31 | 5.91e-14 | 5.89e-14 | 5.92e-14 | 5.94e-14 | 5.94e-14 | 5.89e-14 |
| 1/df polynomial | 63 | 2.95e-13 | 2.94e-13 | 2.94e-13 | 2.94e-13 | 2.94e-13 | 2.95e-13 |
| 1/df polynomial | 127 | 5.31e-13 | 5.30e-13 | 5.31e-13 | 5.31e-13 | 5.30e-13 | 5.30e-13 |
| df polynomial | 7 | 1.87e-05 | 8.87e-09 | 1.88e-07 | 9.05e-07 | 5.07e-07 | 7.48e-06 |
| df polynomial | 15 | 8.12e-09 | 2.58e-14 | 4.33e-13 | 7.87e-10 | 1.44e-10 | 2.51e-08 |
| df polynomial | 31 | 5.81e-14 | 5.87e-14 | 5.92e-14 | 5.89e-14 | 5.90e-14 | 2.54e-13 |
| df polynomial | 63 | 2.94e-13 | 2.94e-13 | 2.94e-13 | 2.94e-13 | 2.95e-13 | 2.94e-13 |
| df polynomial | 127 | 5.30e-13 | 5.30e-13 | 5.30e-13 | 5.30e-13 | 5.31e-13 | 5.31e-13 |

## Plaintext, alpha = 0.01

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-10000 |
|---|---|---|---|
| 1/df polynomial | 7 | 1.50e-01 (8.45e-02) | 1.48e-02 / 9.04e-02 / 1.34e-01 / 1.50e-01 |
| 1/df polynomial | 15 | 1.34e-07 (6.98e-08) | 1.27e-08 / 7.91e-08 / 1.21e-07 / 1.34e-07 |
| 1/df polynomial | 31 | 1.64e-08 (7.05e-09) | 9.76e-09 / 9.86e-09 / 1.20e-08 / 1.64e-08 |
| 1/df polynomial | 63 | 1.19e-08 (3.39e-09) | 5.26e-09 / 8.50e-09 / 6.33e-09 / 1.19e-08 |
| 1/df polynomial | 127 | 8.48e-09 (2.66e-09) | 6.12e-09 / 7.03e-09 / 5.26e-09 / 8.48e-09 |
| df polynomial | 7 | 9.98e-01 (2.10e-01) | 9.98e-01 / 9.29e-01 / 3.07e-01 / 8.16e-02 |
| df polynomial | 15 | 9.98e-01 (1.90e-01) | 9.98e-01 / 9.19e-01 / 2.11e-01 / 5.30e-02 |
| df polynomial | 31 | 9.97e-01 (1.73e-01) | 9.97e-01 / 8.65e-01 / 2.42e-01 / 1.63e-01 |
| df polynomial | 63 | 1.28e+00 (2.68e-01) | 9.85e-01 / 1.28e+00 / 7.16e-01 / 2.95e-01 |
| df polynomial | 127 | 3.42e+00 (6.56e-01) | 2.30e+00 / 3.42e+00 / 1.92e+00 / 8.10e-01 |

Case ranges (MaxRE):

| method | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank |
|---|---|---|---|---|---|---|---|
| 1/df polynomial | 7 | 3.56e-15 | 2.80e-15 | 2.91e-15 | 3.48e-15 | 2.68e-15 | 3.08e-15 |
| 1/df polynomial | 15 | 1.85e-14 | 1.87e-14 | 1.86e-14 | 1.90e-14 | 1.86e-14 | 1.91e-14 |
| 1/df polynomial | 31 | 5.93e-14 | 5.92e-14 | 5.92e-14 | 5.90e-14 | 5.91e-14 | 5.87e-14 |
| 1/df polynomial | 63 | 2.94e-13 | 2.95e-13 | 2.95e-13 | 2.94e-13 | 2.94e-13 | 2.94e-13 |
| 1/df polynomial | 127 | 5.30e-13 | 5.30e-13 | 5.32e-13 | 5.30e-13 | 5.30e-13 | 5.31e-13 |
| df polynomial | 7 | 1.19e-05 | 5.68e-09 | 1.15e-07 | 5.84e-07 | 3.27e-07 | 4.83e-06 |
| df polynomial | 15 | 5.07e-09 | 2.32e-14 | 2.60e-13 | 5.08e-10 | 9.30e-11 | 1.61e-08 |
| df polynomial | 31 | 5.97e-14 | 5.93e-14 | 5.86e-14 | 5.88e-14 | 5.90e-14 | 1.61e-13 |
| df polynomial | 63 | 2.95e-13 | 2.94e-13 | 2.95e-13 | 2.94e-13 | 2.95e-13 | 2.94e-13 |
| df polynomial | 127 | 5.31e-13 | 5.31e-13 | 5.31e-13 | 5.31e-13 | 5.30e-13 | 5.31e-13 |

## Plaintext, alpha = 0.025

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-10000 |
|---|---|---|---|
| 1/df polynomial | 7 | 1.02e-02 (5.91e-03) | 1.83e-03 / 6.94e-03 / 9.37e-03 / 1.02e-02 |
| 1/df polynomial | 15 | 9.79e-09 (8.35e-10) | 9.79e-09 / 5.01e-09 / 1.26e-09 / 7.17e-09 |
| 1/df polynomial | 31 | 9.78e-09 (2.92e-09) | 8.47e-09 / 7.28e-09 / 6.00e-09 / 9.78e-09 |
| 1/df polynomial | 63 | 1.03e-08 (1.85e-09) | 9.86e-09 / 1.03e-08 / 3.80e-09 / 9.26e-09 |
| 1/df polynomial | 127 | 9.26e-09 (7.99e-10) | 7.24e-09 / 9.26e-09 / 7.65e-10 / 8.65e-09 |
| df polynomial | 7 | 9.92e-01 (1.87e-01) | 9.92e-01 / 8.64e-01 / 2.48e-01 / 6.45e-02 |
| df polynomial | 15 | 9.91e-01 (1.70e-01) | 9.91e-01 / 8.48e-01 / 1.67e-01 / 4.08e-02 |
| df polynomial | 31 | 9.86e-01 (1.48e-01) | 9.86e-01 / 7.75e-01 / 1.75e-01 / 1.16e-01 |
| df polynomial | 63 | 9.56e-01 (1.77e-01) | 9.56e-01 / 7.64e-01 / 4.05e-01 / 1.60e-01 |
| df polynomial | 127 | 1.64e+00 (2.75e-01) | 1.28e+00 / 1.64e+00 / 6.51e-01 / 2.65e-01 |

Case ranges (MaxRE):

| method | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank |
|---|---|---|---|---|---|---|---|
| 1/df polynomial | 7 | 3.35e-15 | 3.00e-15 | 3.15e-15 | 3.00e-15 | 3.01e-15 | 3.00e-15 |
| 1/df polynomial | 15 | 1.94e-14 | 1.92e-14 | 1.85e-14 | 1.87e-14 | 1.89e-14 | 1.89e-14 |
| 1/df polynomial | 31 | 5.96e-14 | 5.88e-14 | 6.00e-14 | 5.92e-14 | 5.94e-14 | 5.88e-14 |
| 1/df polynomial | 63 | 2.94e-13 | 2.94e-13 | 2.95e-13 | 2.94e-13 | 2.95e-13 | 2.94e-13 |
| 1/df polynomial | 127 | 5.31e-13 | 5.31e-13 | 5.31e-13 | 5.31e-13 | 5.31e-13 | 5.31e-13 |
| df polynomial | 7 | 9.32e-06 | 4.47e-09 | 8.93e-08 | 4.61e-07 | 2.58e-07 | 3.81e-06 |
| df polynomial | 15 | 3.94e-09 | 2.31e-14 | 2.01e-13 | 4.00e-10 | 7.33e-11 | 1.27e-08 |
| df polynomial | 31 | 5.93e-14 | 5.88e-14 | 5.97e-14 | 5.94e-14 | 5.90e-14 | 1.27e-13 |
| df polynomial | 63 | 2.95e-13 | 2.94e-13 | 2.94e-13 | 2.94e-13 | 2.95e-13 | 2.94e-13 |
| df polynomial | 127 | 5.31e-13 | 5.31e-13 | 5.31e-13 | 5.30e-13 | 5.31e-13 | 5.31e-13 |

## Plaintext, alpha = 0.05

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-10000 |
|---|---|---|---|
| 1/df polynomial | 7 | 1.00e-03 (5.98e-04) | 2.74e-04 / 7.45e-04 / 9.42e-04 / 1.00e-03 |
| 1/df polynomial | 15 | 9.20e-09 (9.93e-10) | 9.77e-11 / 9.20e-09 / 3.23e-09 / 5.15e-09 |
| 1/df polynomial | 31 | 8.38e-09 (3.63e-10) | 1.80e-11 / 5.78e-09 / 2.36e-10 / 8.38e-09 |
| 1/df polynomial | 63 | 9.07e-09 (5.03e-10) | 2.83e-11 / 6.99e-09 / 1.63e-09 / 9.07e-09 |
| 1/df polynomial | 127 | 8.79e-09 (1.80e-10) | 1.38e-11 / 8.79e-09 / 6.04e-10 / 6.26e-09 |
| df polynomial | 7 | 9.75e-01 (1.66e-01) | 9.75e-01 / 7.86e-01 / 2.02e-01 / 5.19e-02 |
| df polynomial | 15 | 9.73e-01 (1.50e-01) | 9.73e-01 / 7.66e-01 / 1.34e-01 / 3.23e-02 |
| df polynomial | 31 | 9.63e-01 (1.27e-01) | 9.63e-01 / 6.80e-01 / 1.32e-01 / 8.68e-02 |
| df polynomial | 63 | 9.08e-01 (1.28e-01) | 9.08e-01 / 5.01e-01 / 2.56e-01 / 9.81e-02 |
| df polynomial | 127 | 9.14e-01 (1.43e-01) | 7.76e-01 / 9.14e-01 / 2.85e-01 / 1.13e-01 |

Case ranges (MaxRE):

| method | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank |
|---|---|---|---|---|---|---|---|
| 1/df polynomial | 7 | 2.86e-15 | 2.88e-15 | 2.84e-15 | 3.12e-15 | 3.00e-15 | 3.24e-15 |
| 1/df polynomial | 15 | 1.89e-14 | 1.88e-14 | 1.86e-14 | 1.88e-14 | 1.91e-14 | 1.91e-14 |
| 1/df polynomial | 31 | 5.91e-14 | 5.93e-14 | 5.92e-14 | 5.92e-14 | 5.94e-14 | 5.87e-14 |
| 1/df polynomial | 63 | 2.94e-13 | 2.95e-13 | 2.94e-13 | 2.94e-13 | 2.94e-13 | 2.95e-13 |
| 1/df polynomial | 127 | 5.31e-13 | 5.31e-13 | 5.31e-13 | 5.31e-13 | 5.31e-13 | 5.31e-13 |
| df polynomial | 7 | 7.45e-06 | 3.59e-09 | 7.07e-08 | 3.70e-07 | 2.07e-07 | 3.06e-06 |
| df polynomial | 15 | 3.14e-09 | 2.10e-14 | 1.62e-13 | 3.22e-10 | 5.89e-11 | 1.02e-08 |
| df polynomial | 31 | 5.87e-14 | 5.94e-14 | 5.91e-14 | 5.92e-14 | 5.91e-14 | 1.03e-13 |
| df polynomial | 63 | 2.94e-13 | 2.94e-13 | 2.94e-13 | 2.94e-13 | 2.94e-13 | 2.94e-13 |
| df polynomial | 127 | 5.31e-13 | 5.30e-13 | 5.31e-13 | 5.30e-13 | 5.30e-13 | 5.30e-13 |

## Plaintext, alpha = 0.1

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-10000 |
|---|---|---|---|
| 1/df polynomial | 7 | 6.64e-05 (4.08e-05) | 2.62e-05 / 5.33e-05 / 6.35e-05 / 6.64e-05 |
| 1/df polynomial | 15 | 1.08e-08 (5.03e-10) | 1.71e-10 / 7.43e-09 / 1.08e-08 / 1.03e-08 |
| 1/df polynomial | 31 | 1.42e-08 (6.99e-10) | 6.84e-11 / 8.16e-09 / 1.42e-08 / 1.37e-08 |
| 1/df polynomial | 63 | 1.60e-08 (1.05e-09) | 1.15e-10 / 9.50e-09 / 9.71e-09 / 1.60e-08 |
| 1/df polynomial | 127 | 1.44e-08 (3.41e-10) | 4.35e-11 / 6.65e-09 / 1.33e-08 / 1.44e-08 |
| df polynomial | 7 | 9.31e-01 (1.40e-01) | 9.31e-01 / 6.75e-01 / 1.57e-01 / 3.98e-02 |
| df polynomial | 15 | 9.26e-01 (1.27e-01) | 9.26e-01 / 6.52e-01 / 1.03e-01 / 2.43e-02 |
| df polynomial | 31 | 9.05e-01 (1.05e-01) | 9.05e-01 / 5.59e-01 / 9.45e-02 / 6.19e-02 |
| df polynomial | 63 | 8.13e-01 (8.97e-02) | 8.13e-01 / 3.15e-01 / 1.55e-01 / 5.82e-02 |
| df polynomial | 127 | 5.08e-01 (7.34e-02) | 5.08e-01 / 4.85e-01 / 1.22e-01 / 4.76e-02 |

Case ranges (MaxRE):

| method | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank |
|---|---|---|---|---|---|---|---|
| 1/df polynomial | 7 | 3.26e-15 | 2.79e-15 | 2.92e-15 | 2.95e-15 | 3.12e-15 | 3.45e-15 |
| 1/df polynomial | 15 | 1.87e-14 | 1.87e-14 | 1.89e-14 | 1.85e-14 | 1.84e-14 | 1.85e-14 |
| 1/df polynomial | 31 | 5.90e-14 | 5.91e-14 | 5.94e-14 | 5.93e-14 | 5.91e-14 | 5.94e-14 |
| 1/df polynomial | 63 | 2.94e-13 | 2.95e-13 | 2.94e-13 | 2.95e-13 | 2.94e-13 | 2.94e-13 |
| 1/df polynomial | 127 | 5.31e-13 | 5.30e-13 | 5.31e-13 | 5.30e-13 | 5.31e-13 | 5.31e-13 |
| df polynomial | 7 | 5.68e-06 | 2.74e-09 | 5.34e-08 | 2.83e-07 | 1.59e-07 | 2.34e-06 |
| df polynomial | 15 | 2.38e-09 | 2.08e-14 | 1.25e-13 | 2.46e-10 | 4.51e-11 | 7.83e-09 |
| df polynomial | 31 | 5.93e-14 | 5.95e-14 | 5.93e-14 | 5.91e-14 | 5.92e-14 | 8.04e-14 |
| df polynomial | 63 | 2.94e-13 | 2.95e-13 | 2.95e-13 | 2.94e-13 | 2.94e-13 | 2.95e-13 |
| df polynomial | 127 | 5.30e-13 | 5.30e-13 | 5.30e-13 | 5.30e-13 | 5.30e-13 | 5.30e-13 |

## HE evaluation

`run_critical_boundary.py --part he`, 10 repetitions (polynomials). Every method starts from the same encrypted
u = 1/df at level 12 (65,536 points in two ciphertexts) under the conditions above. MaxRE / MRE against the exact c^2,
mean (standard deviation) over the repetitions; time excludes encryption and decryption.

- 1/df (stored): `HEHypothesisTesting._critical_value_from_inv_df` with the stored coefficients.
- 1/df (range fit): 1/df polynomial refitted on the case range.
- df (HE division): df = 1/u computed in HE (invSqrt of u * df_lo on [df_lo / df_hi, 1], squared; degree 63 with 4
  Newton iterations on the global range, degree 31 with 1 iteration on the case ranges), then the df polynomial fitted
  on the condition range. The time includes the division.

### Polynomials, global range [1, 10^4]

| input | degree | time (s) | BTS | MaxRE alpha = 0.001 | MaxRE alpha = 0.01 | MaxRE alpha = 0.025 | MaxRE alpha = 0.05 | MaxRE alpha = 0.1 |
|---|---|---|---|---|---|---|---|---|
| 1/df (stored) | 7 | 0.37 (0.03) | 0 | 4.74e+01 (1.78e-04) | 1.50e-01 (6.36e-07) | 1.02e-02 (1.31e-07) | 1.00e-03 (4.48e-08) | 6.65e-05 (2.57e-08) |
| 1/df (stored) | 15 | 0.65 (0.06) | 0 | 9.08e-04 (1.34e-06) | 1.16e-06 (8.22e-08) | 8.46e-07 (1.47e-07) | 7.07e-07 (1.01e-07) | 4.52e-07 (3.71e-08) |
| 1/df (stored) | 31 | 1.13 (0.18) | 0 | 2.81e-05 (3.21e-06) | 1.17e-06 (1.11e-07) | 9.25e-07 (2.13e-07) | 7.45e-07 (1.84e-07) | 5.30e-07 (1.06e-07) |
| 1/df (stored) | 63 | 1.75 (0.44) | 0 | 2.81e-05 (3.86e-06) | 1.14e-06 (1.41e-07) | 9.20e-07 (1.56e-07) | 6.39e-07 (7.80e-08) | 4.65e-07 (4.68e-08) |
| 1/df (stored) | 127 | 2.99 (1.14) | 0 | 2.88e-05 (2.26e-06) | 1.11e-06 (1.33e-07) | 8.83e-07 (5.94e-08) | 6.90e-07 (1.25e-07) | 4.61e-07 (5.86e-08) |
| df (HE division) | 7 | 15.80 (1.00) | 2 | 1.00e+00 (2.17e-15) | 9.98e-01 (1.54e-13) | 9.92e-01 (5.59e-13) | 9.75e-01 (4.56e-12) | 9.31e-01 (1.35e-11) |
| df (HE division) | 15 | 16.04 (1.00) | 2 | 1.00e+00 (3.51e-14) | 9.98e-01 (1.26e-12) | 9.91e-01 (4.94e-12) | 9.73e-01 (1.77e-11) | 9.26e-01 (3.46e-11) |
| df (HE division) | 31 | 20.45 (1.22) | 3 | 1.00e+00 (6.60e-11) | 9.97e-01 (1.96e-09) | 9.86e-01 (6.59e-09) | 9.63e-01 (1.21e-08) | 9.05e-01 (3.05e-08) |
| df (HE division) | 63 | 21.24 (1.07) | 3 | 4.22e+00 (1.98e-04) | 1.28e+00 (4.06e-05) | 9.56e-01 (1.68e-07) | 9.08e-01 (2.54e-07) | 8.13e-01 (5.35e-07) |
| df (HE division) | 127 | 23.13 (1.48) | 3 | 4.08e+01 (1.63e+00) | 3.42e+00 (8.92e-04) | 1.64e+00 (3.48e-04) | 9.13e-01 (8.59e-05) | 5.08e-01 (7.48e-06) |

Degree-15 1/df polynomial at alpha = 0.05: MaxRE 7.07e-07 (1.01e-07), MRE 4.62e-08 (1.79e-10).

Division stage of the df input: time, bootstraps, MaxRE of the computed df.

| condition | time (s) | BTS | MaxRE of df |
|---|---|---|---|
| global | 14.40 (0.42) | 2 | 1.70e-03 |
| Insurance | 10.02 (0.57) | 1 | 2.16e-04 |
| Adult | 1.06 (0.03) | 0 | 4.05e-04 |
| Heart | 9.30 (0.25) | 1 | 5.97e-05 |
| Diabetes | 9.97 (0.43) | 1 | 1.45e-02 |
| Credit | 10.08 (0.63) | 1 | 6.57e-03 |
| Bank | 10.32 (0.94) | 1 | 7.93e-03 |

### Polynomials, case ranges, alpha = 0.05 (MaxRE over the range, mean over repetitions)

| input | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank |
|---|---|---|---|---|---|---|---|
| 1/df (stored) | 7 | 9.18e-04 | 9.52e-04 | 6.09e-04 | 1.01e-03 | 1.01e-03 | 1.01e-03 |
| 1/df (stored) | 15 | 4.79e-07 | 4.68e-07 | 4.51e-07 | 4.73e-07 | 4.76e-07 | 4.57e-07 |
| 1/df (stored) | 31 | 4.81e-07 | 4.75e-07 | 4.67e-07 | 4.82e-07 | 4.84e-07 | 4.92e-07 |
| 1/df (range fit) | 7 | 4.82e-07 | 4.99e-07 | 4.71e-07 | 4.64e-07 | 4.82e-07 | 4.69e-07 |
| 1/df (range fit) | 15 | 4.85e-07 | 4.73e-07 | 5.10e-07 | 4.89e-07 | 4.76e-07 | 4.60e-07 |
| 1/df (range fit) | 31 | 4.60e-07 | 4.81e-07 | 4.82e-07 | 4.59e-07 | 4.61e-07 | 4.64e-07 |
| df (HE division) | 7 | 7.54e-06 | 4.68e-07 | 4.94e-07 | 6.58e-07 | 5.41e-07 | 3.20e-06 |
| df (HE division) | 15 | 4.67e-07 | 4.68e-07 | 4.72e-07 | 5.02e-07 | 4.98e-07 | 4.75e-07 |
| df (HE division) | 31 | 4.67e-07 | 4.68e-07 | 4.72e-07 | 5.02e-07 | 4.98e-07 | 4.77e-07 |

### Polynomials, case ranges, alpha = 0.001 (MaxRE over the range, mean over repetitions)

| input | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank |
|---|---|---|---|---|---|---|---|
| 1/df (stored) | 7 | 4.33e+01 | 4.49e+01 | 2.85e+01 | 4.79e+01 | 4.78e+01 | 4.79e+01 |
| 1/df (stored) | 15 | 6.15e-04 | 7.28e-04 | 8.40e-04 | 9.53e-04 | 9.43e-04 | 9.48e-04 |
| 1/df (stored) | 31 | 2.92e-05 | 3.22e-05 | 2.81e-05 | 3.04e-05 | 3.19e-05 | 2.99e-05 |
| 1/df (range fit) | 7 | 1.09e-06 | 1.18e-06 | 1.20e-06 | 1.13e-06 | 1.08e-06 | 1.15e-06 |
| 1/df (range fit) | 15 | 1.15e-06 | 1.19e-06 | 1.22e-06 | 1.16e-06 | 1.15e-06 | 1.11e-06 |
| 1/df (range fit) | 31 | 1.21e-06 | 1.12e-06 | 1.15e-06 | 1.16e-06 | 1.20e-06 | 1.14e-06 |
| df (HE division) | 7 | 1.89e-05 | 1.15e-06 | 1.22e-06 | 1.61e-06 | 1.32e-06 | 7.83e-06 |
| df (HE division) | 15 | 1.14e-06 | 1.14e-06 | 1.16e-06 | 1.23e-06 | 1.22e-06 | 1.16e-06 |
| df (HE division) | 31 | 1.14e-06 | 1.14e-06 | 1.16e-06 | 1.23e-06 | 1.22e-06 | 1.16e-06 |