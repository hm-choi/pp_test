# Critical-value approximation: method comparison

Target c^2 = t_{1-alpha/2, df}^2 (two-sided), reference SciPy `t.ppf`; relative error. Script `run_critical_boundary.py`;
raw results `results/plaintext.csv`, `results/he_raw.csv`, `results/he_lookup_raw.csv`.

Conditions (every method is evaluated on the same 65,536 df points per condition, log-spaced over the range plus
the Welch df of the cases):
- global: df in [1, 2000];
- case range: df in [min(n1, n2) - 1, n1 + n2 - 2] for each dataset comparison (public; it always contains the Welch df):
  Insurance [273, 1336], Adult [1023, 2046], Heart [138, 301], Diabetes [11356, 66219], Credit [6635, 29998], Bank [5288, 45209] (Adult: both comparisons).

Methods:
- 1/df polynomial: Chebyshev interpolation in u = 1/df. Global condition: the stored coefficients
  (`coeffs/t_critical_coeffs.py`, fitted on u in [0, 1], i.e. every df >= 1, as used by PP-TEST); case ranges:
  refitted on [1/df_hi, 1/df_lo].
- df polynomial: Chebyshev interpolation in df, fitted on the condition range.
- Lookup tables: textbook df grid (1..30, 40, 50, 60, 80, 100, 120, 1000; interpolation adds df = infinity) or every
  integer df up to the largest range, mapped by floor, round or linear interpolation of c^2 in 1/df.

## Plaintext, alpha = 0.001

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 |
|---|---|---|---|
| 1/df polynomial | 7 | 4.49e+01 (2.08e+01) | 8.48e-01 / 2.04e+01 / 3.95e+01 / 4.49e+01 |
| 1/df polynomial | 15 | 8.96e-04 (3.65e-04) | 1.40e-05 / 3.74e-04 / 7.63e-04 / 8.96e-04 |
| 1/df polynomial | 31 | 3.41e-08 (1.34e-08) | 9.82e-09 / 1.95e-08 / 2.87e-08 / 3.41e-08 |
| 1/df polynomial | 63 | 3.81e-07 (1.73e-07) | 8.71e-09 / 2.06e-07 / 3.14e-07 / 3.81e-07 |
| 1/df polynomial | 127 | 1.34e-07 (6.04e-08) | 9.61e-09 / 7.34e-08 / 1.13e-07 / 1.34e-07 |
| df polynomial | 7 | 1.00e+00 (2.73e-01) | 1.00e+00 / 9.84e-01 / 2.75e-01 / 1.13e-01 |
| df polynomial | 15 | 1.00e+00 (3.03e-01) | 1.00e+00 / 9.58e-01 / 6.48e-01 / 2.99e-01 |
| df polynomial | 31 | 6.15e+00 (1.34e+00) | 9.98e-01 / 6.15e+00 / 3.35e+00 / 1.70e+00 |
| df polynomial | 63 | 5.27e+01 (1.11e+01) | 1.14e+01 / 5.27e+01 / 4.49e+01 / 2.27e+01 |
| df polynomial | 127 | 8.72e+01 (1.84e+01) | 1.36e+01 / 8.72e+01 / 7.62e+01 / 3.91e+01 |
| textbook table, floor |  | 4.05e+02 (8.10e+00) | 4.05e+02 / 4.98e+00 / 6.87e-02 / 5.42e-02 |
| textbook table + inf, interp. in 1/df |  | 2.55e+01 (1.14e+00) | 2.55e+01 / 4.98e-01 / 7.03e-04 / 4.76e-04 |
| integer table, floor |  | 4.05e+02 (8.09e+00) | 4.05e+02 / 4.98e+00 / 6.87e-02 / 6.87e-03 |
| integer table, round |  | 5.80e+01 (8.75e-01) | 5.80e+01 / 1.99e+00 / 3.55e-02 / 3.48e-03 |
| integer table, interp. in 1/df |  | 2.55e+01 (1.14e+00) | 2.55e+01 / 4.98e-01 / 7.03e-04 / 8.02e-06 |

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
| textbook table, floor |  | 4.48e-02 | 3.03e-03 | 3.05e-02 | 5.85e-03 | 5.74e-03 | 5.81e-03 |
| textbook table + inf, interp. in 1/df |  | 3.13e-04 | 6.17e-06 | 3.37e-04 | 1.98e-06 | 3.16e-06 | 3.79e-06 |
| integer table, floor |  | 7.93e-05 | 5.64e-06 | 3.13e-04 | 4.51e-08 | 1.33e-07 | 2.09e-07 |
| integer table, round |  | 3.97e-05 | 2.82e-06 | 1.57e-04 | 2.29e-08 | 6.69e-08 | 1.05e-07 |
| integer table, interp. in 1/df |  | 1.12e-09 | 5.63e-12 | 1.72e-08 | 9.84e-16 | 3.44e-15 | 8.36e-15 |

## Plaintext, alpha = 0.01

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 |
|---|---|---|---|
| 1/df polynomial | 7 | 1.42e-01 (7.14e-02) | 1.48e-02 / 9.04e-02 / 1.34e-01 / 1.42e-01 |
| 1/df polynomial | 15 | 1.34e-07 (5.90e-08) | 1.27e-08 / 7.91e-08 / 1.21e-07 / 1.34e-07 |
| 1/df polynomial | 31 | 1.64e-08 (6.83e-09) | 9.75e-09 / 9.86e-09 / 1.20e-08 / 1.64e-08 |
| 1/df polynomial | 63 | 1.19e-08 (3.51e-09) | 5.30e-09 / 8.50e-09 / 6.33e-09 / 1.19e-08 |
| 1/df polynomial | 127 | 8.48e-09 (2.51e-09) | 6.20e-09 / 7.06e-09 / 5.26e-09 / 8.48e-09 |
| df polynomial | 7 | 9.98e-01 (2.26e-01) | 9.98e-01 / 9.15e-01 / 1.79e-01 / 6.75e-02 |
| df polynomial | 15 | 9.96e-01 (2.11e-01) | 9.96e-01 / 8.42e-01 / 3.15e-01 / 1.38e-01 |
| df polynomial | 31 | 1.64e+00 (3.79e-01) | 9.78e-01 / 1.64e+00 / 6.42e-01 / 3.07e-01 |
| df polynomial | 63 | 3.74e+00 (8.50e-01) | 2.75e+00 / 3.74e+00 / 2.46e+00 / 1.09e+00 |
| df polynomial | 127 | 2.78e+00 (5.79e-01) | 1.47e+00 / 2.78e+00 / 1.77e+00 / 7.99e-01 |
| textbook table, floor |  | 4.01e+01 (1.20e+00) | 4.01e+01 / 1.89e+00 / 4.13e-02 / 3.40e-02 |
| textbook table + inf, interp. in 1/df |  | 3.62e+00 (2.02e-01) | 3.62e+00 / 1.67e-01 / 2.74e-04 / 1.97e-04 |
| integer table, floor |  | 4.01e+01 (1.20e+00) | 4.01e+01 / 1.89e+00 / 4.13e-02 / 4.33e-03 |
| integer table, round |  | 1.18e+01 (2.59e-01) | 1.18e+01 / 9.19e-01 / 2.15e-02 / 2.20e-03 |
| integer table, interp. in 1/df |  | 3.62e+00 (2.02e-01) | 3.62e+00 / 1.67e-01 / 2.74e-04 / 3.30e-06 |

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
| textbook table, floor |  | 2.86e-02 | 1.96e-03 | 1.95e-02 | 3.77e-03 | 3.70e-03 | 3.74e-03 |
| textbook table + inf, interp. in 1/df |  | 1.32e-04 | 2.63e-06 | 1.43e-04 | 8.45e-07 | 1.35e-06 | 1.61e-06 |
| integer table, floor |  | 5.10e-05 | 3.64e-06 | 2.01e-04 | 2.91e-08 | 8.61e-08 | 1.35e-07 |
| integer table, round |  | 2.55e-05 | 1.82e-06 | 1.01e-04 | 1.48e-08 | 4.32e-08 | 6.76e-08 |
| integer table, interp. in 1/df |  | 4.74e-10 | 2.40e-12 | 7.27e-09 | 9.37e-16 | 1.87e-15 | 3.75e-15 |

## Plaintext, alpha = 0.025

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 |
|---|---|---|---|
| 1/df polynomial | 7 | 9.64e-03 (5.06e-03) | 1.83e-03 / 6.94e-03 / 9.37e-03 / 9.64e-03 |
| 1/df polynomial | 15 | 9.79e-09 (8.42e-10) | 9.79e-09 / 5.01e-09 / 1.26e-09 / 7.17e-09 |
| 1/df polynomial | 31 | 9.78e-09 (2.89e-09) | 8.48e-09 / 7.29e-09 / 6.00e-09 / 9.78e-09 |
| 1/df polynomial | 63 | 1.03e-08 (1.95e-09) | 9.85e-09 / 1.03e-08 / 3.80e-09 / 9.26e-09 |
| 1/df polynomial | 127 | 9.27e-09 (7.70e-10) | 7.22e-09 / 9.27e-09 / 7.65e-10 / 8.64e-09 |
| df polynomial | 7 | 9.91e-01 (2.01e-01) | 9.91e-01 / 8.43e-01 / 1.41e-01 / 5.17e-02 |
| df polynomial | 15 | 9.85e-01 (1.77e-01) | 9.85e-01 / 7.47e-01 / 2.23e-01 / 9.62e-02 |
| df polynomial | 31 | 9.43e-01 (2.34e-01) | 9.43e-01 / 9.40e-01 / 3.27e-01 / 1.54e-01 |
| df polynomial | 63 | 1.57e+00 (3.32e-01) | 1.44e+00 / 1.57e+00 / 7.84e-01 / 3.31e-01 |
| df polynomial | 127 | 7.41e-01 (1.53e-01) | 5.42e-01 / 7.41e-01 / 3.86e-01 / 1.66e-01 |
| textbook table, floor |  | 1.58e+01 (5.74e-01) | 1.58e+01 / 1.21e+00 / 3.16e-02 / 2.65e-02 |
| textbook table + inf, interp. in 1/df |  | 1.61e+00 (9.56e-02) | 1.61e+00 / 9.42e-02 / 1.67e-04 / 1.23e-04 |
| integer table, floor |  | 1.58e+01 (5.68e-01) | 1.58e+01 / 1.21e+00 / 3.16e-02 / 3.38e-03 |
| integer table, round |  | 5.97e+00 (1.55e-01) | 5.97e+00 / 6.28e-01 / 1.65e-02 / 1.72e-03 |
| integer table, interp. in 1/df |  | 1.61e+00 (9.56e-02) | 1.61e+00 / 9.42e-02 / 1.67e-04 / 2.06e-06 |

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
| textbook table, floor |  | 2.25e-02 | 1.54e-03 | 1.53e-02 | 2.97e-03 | 2.92e-03 | 2.95e-03 |
| textbook table + inf, interp. in 1/df |  | 8.35e-05 | 1.66e-06 | 9.01e-05 | 5.34e-07 | 8.50e-07 | 1.02e-06 |
| integer table, floor |  | 4.02e-05 | 2.87e-06 | 1.58e-04 | 2.30e-08 | 6.80e-08 | 1.07e-07 |
| integer table, round |  | 2.01e-05 | 1.43e-06 | 7.94e-05 | 1.16e-08 | 3.41e-08 | 5.34e-08 |
| integer table, interp. in 1/df |  | 2.99e-10 | 1.51e-12 | 4.58e-09 | 1.06e-15 | 1.59e-15 | 2.65e-15 |

## Plaintext, alpha = 0.05

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 |
|---|---|---|---|
| 1/df polynomial | 7 | 9.50e-04 (5.17e-04) | 2.74e-04 / 7.45e-04 / 9.42e-04 / 9.50e-04 |
| 1/df polynomial | 15 | 9.20e-09 (1.07e-09) | 9.77e-11 / 9.20e-09 / 3.23e-09 / 5.15e-09 |
| 1/df polynomial | 31 | 8.38e-09 (4.04e-10) | 1.80e-11 / 5.78e-09 / 2.36e-10 / 8.38e-09 |
| 1/df polynomial | 63 | 9.08e-09 (5.78e-10) | 2.83e-11 / 7.00e-09 / 1.63e-09 / 9.08e-09 |
| 1/df polynomial | 127 | 8.80e-09 (2.10e-10) | 1.38e-11 / 8.80e-09 / 6.04e-10 / 6.25e-09 |
| df polynomial | 7 | 9.72e-01 (1.78e-01) | 9.72e-01 / 7.59e-01 / 1.13e-01 / 4.07e-02 |
| df polynomial | 15 | 9.59e-01 (1.51e-01) | 9.59e-01 / 6.49e-01 / 1.66e-01 / 7.05e-02 |
| df polynomial | 31 | 8.86e-01 (1.61e-01) | 8.86e-01 / 5.99e-01 / 1.93e-01 / 8.92e-02 |
| df polynomial | 63 | 8.68e-01 (1.64e-01) | 8.32e-01 / 8.68e-01 / 3.28e-01 / 1.34e-01 |
| df polynomial | 127 | 2.83e-01 (5.57e-02) | 2.36e-01 / 2.83e-01 / 1.18e-01 / 4.92e-02 |
| textbook table, floor |  | 7.72e+00 (3.24e-01) | 7.72e+00 / 8.28e-01 / 2.48e-02 / 2.11e-02 |
| textbook table + inf, interp. in 1/df |  | 8.28e-01 (5.09e-02) | 8.28e-01 / 5.58e-02 / 1.06e-04 / 7.94e-05 |
| integer table, floor |  | 7.72e+00 (3.20e-01) | 7.72e+00 / 8.28e-01 / 2.48e-02 / 2.70e-03 |
| integer table, round |  | 3.46e+00 (1.02e-01) | 3.46e+00 / 4.49e-01 / 1.29e-02 / 1.37e-03 |
| integer table, interp. in 1/df |  | 8.28e-01 (5.09e-02) | 8.28e-01 / 5.58e-02 / 1.06e-04 / 1.33e-06 |

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
| textbook table, floor |  | 1.80e-02 | 1.24e-03 | 1.23e-02 | 2.39e-03 | 2.34e-03 | 2.37e-03 |
| textbook table + inf, interp. in 1/df |  | 5.45e-05 | 1.09e-06 | 5.89e-05 | 3.49e-07 | 5.57e-07 | 6.67e-07 |
| integer table, floor |  | 3.23e-05 | 2.30e-06 | 1.27e-04 | 1.85e-08 | 5.46e-08 | 8.57e-08 |
| integer table, round |  | 1.62e-05 | 1.15e-06 | 6.37e-05 | 9.36e-09 | 2.74e-08 | 4.29e-08 |
| integer table, interp. in 1/df |  | 1.95e-10 | 9.92e-13 | 2.99e-09 | 5.78e-16 | 1.04e-15 | 1.73e-15 |

## Plaintext, alpha = 0.1

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 |
|---|---|---|---|
| 1/df polynomial | 7 | 6.35e-05 (3.56e-05) | 2.62e-05 / 5.33e-05 / 6.35e-05 / 6.29e-05 |
| 1/df polynomial | 15 | 1.08e-08 (5.96e-10) | 1.71e-10 / 7.43e-09 / 1.08e-08 / 1.03e-08 |
| 1/df polynomial | 31 | 1.42e-08 (7.73e-10) | 6.84e-11 / 8.16e-09 / 1.42e-08 / 1.37e-08 |
| 1/df polynomial | 63 | 1.60e-08 (1.16e-09) | 1.15e-10 / 9.51e-09 / 9.72e-09 / 1.60e-08 |
| 1/df polynomial | 127 | 1.44e-08 (3.65e-10) | 4.35e-11 / 6.65e-09 / 1.33e-08 / 1.44e-08 |
| df polynomial | 7 | 9.24e-01 (1.50e-01) | 9.24e-01 / 6.45e-01 / 8.58e-02 / 3.05e-02 |
| df polynomial | 15 | 8.98e-01 (1.22e-01) | 8.98e-01 / 5.28e-01 / 1.17e-01 / 4.94e-02 |
| df polynomial | 31 | 7.82e-01 (1.08e-01) | 7.82e-01 / 3.65e-01 / 1.10e-01 / 5.02e-02 |
| df polynomial | 63 | 4.53e-01 (8.01e-02) | 4.47e-01 / 4.53e-01 / 1.34e-01 / 5.35e-02 |
| df polynomial | 127 | 1.07e-01 (1.97e-02) | 1.04e-01 / 1.07e-01 / 3.43e-02 / 1.39e-02 |
| textbook table, floor |  | 3.68e+00 (1.78e-01) | 3.68e+00 / 5.39e-01 / 1.85e-02 / 1.60e-02 |
| textbook table + inf, interp. in 1/df |  | 3.88e-01 (2.44e-02) | 3.88e-01 / 2.95e-02 / 6.07e-05 / 4.67e-05 |
| integer table, floor |  | 3.68e+00 (1.75e-01) | 3.68e+00 / 5.39e-01 / 1.85e-02 / 2.05e-03 |
| integer table, round |  | 1.90e+00 (6.41e-02) | 1.90e+00 / 3.03e-01 / 9.68e-03 / 1.04e-03 |
| integer table, interp. in 1/df |  | 3.88e-01 (2.44e-02) | 3.88e-01 / 2.95e-02 / 6.07e-05 / 7.82e-07 |

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
| textbook table, floor |  | 1.37e-02 | 9.48e-04 | 9.38e-03 | 1.83e-03 | 1.79e-03 | 1.81e-03 |
| textbook table + inf, interp. in 1/df |  | 3.24e-05 | 6.47e-07 | 3.50e-05 | 2.08e-07 | 3.31e-07 | 3.97e-07 |
| integer table, floor |  | 2.47e-05 | 1.76e-06 | 9.71e-05 | 1.41e-08 | 4.18e-08 | 6.56e-08 |
| integer table, round |  | 1.24e-05 | 8.82e-07 | 4.87e-05 | 7.16e-09 | 2.09e-08 | 3.28e-08 |
| integer table, interp. in 1/df |  | 1.16e-10 | 5.90e-13 | 1.77e-09 | 6.57e-16 | 8.21e-16 | 1.15e-15 |

## HE evaluation

`run_critical_boundary.py --part he`, 2 repetitions (polynomials). Every method starts from the same encrypted
u = 1/df at level 12 (65,536 points in two ciphertexts) under the conditions above. MaxRE / MRE against the exact c^2,
mean (standard deviation) over the repetitions; time excludes encryption and decryption.

- 1/df (stored): `HEHypothesisTesting._critical_value_from_inv_df` with the stored coefficients.
- 1/df (range fit): 1/df polynomial refitted on the case range.
- df (HE division): df = 1/u computed in HE (invSqrt of u * df_lo on [df_lo / df_hi, 1], squared; degree 63 with 4
  Newton iterations on the global range, degree 31 with 1 iteration on the case ranges), then the df polynomial fitted
  on the condition range. The time includes the division.

### Polynomials, global range [1, 2000]

| input | degree | time (s) | BTS | MaxRE alpha = 0.001 | MaxRE alpha = 0.01 | MaxRE alpha = 0.025 | MaxRE alpha = 0.05 | MaxRE alpha = 0.1 |
|---|---|---|---|---|---|---|---|---|
| 1/df (stored) | 7 | 0.38 (0.04) | 0 | 4.49e+01 (1.83e-04) | 1.42e-01 (4.30e-07) | 9.64e-03 (1.76e-09) | 9.50e-04 (7.53e-08) | 6.37e-05 (1.07e-08) |
| 1/df (stored) | 15 | 0.64 (0.03) | 0 | 9.07e-04 (1.79e-06) | 1.12e-06 (2.13e-08) | 8.29e-07 (1.57e-07) | 6.84e-07 (6.28e-09) | 4.45e-07 (6.44e-08) |
| 1/df (stored) | 31 | 1.09 (0.13) | 0 | 2.52e-05 (1.15e-06) | 1.17e-06 (6.22e-08) | 8.61e-07 (6.41e-08) | 7.41e-07 (2.00e-07) | 5.19e-07 (9.52e-09) |
| 1/df (stored) | 63 | 1.81 (0.48) | 0 | 2.88e-05 (2.39e-06) | 1.16e-06 (2.03e-07) | 8.94e-07 (1.26e-07) | 7.49e-07 (1.20e-07) | 5.60e-07 (8.12e-08) |
| 1/df (stored) | 127 | 3.04 (1.20) | 0 | 2.86e-05 (5.31e-06) | 1.10e-06 (1.33e-07) | 8.30e-07 (7.24e-10) | 7.20e-07 (1.87e-07) | 4.39e-07 (7.91e-08) |
| df (HE division) | 7 | 14.98 (0.71) | 2 | 1.00e+00 (1.39e-14) | 9.98e-01 (6.21e-13) | 9.91e-01 (1.85e-13) | 9.72e-01 (7.65e-12) | 9.24e-01 (1.02e-11) |
| df (HE division) | 15 | 15.22 (0.71) | 2 | 1.00e+00 (4.58e-13) | 9.96e-01 (1.28e-11) | 9.85e-01 (5.40e-12) | 9.59e-01 (9.27e-11) | 8.98e-01 (1.99e-10) |
| df (HE division) | 31 | 19.44 (0.87) | 3 | 6.15e+00 (8.28e-05) | 1.64e+00 (1.02e-05) | 9.43e-01 (5.47e-08) | 8.86e-01 (9.64e-08) | 7.82e-01 (1.56e-07) |
| df (HE division) | 63 | 20.24 (0.84) | 3 | 5.27e+01 (1.92e-04) | 3.74e+00 (1.85e-07) | 1.57e+00 (4.11e-05) | 8.68e-01 (1.31e-05) | 4.53e-01 (2.01e-05) |
| df (HE division) | 127 | 22.02 (1.10) | 3 | 8.72e+01 (5.21e-04) | 2.78e+00 (8.80e-06) | 7.41e-01 (1.35e-05) | 2.83e-01 (2.91e-06) | 1.08e-01 (6.24e-06) |

Degree-15 1/df polynomial at alpha = 0.05: MaxRE 6.84e-07 (6.28e-09), MRE 4.74e-08 (9.10e-11).

Division stage of the df input: time, bootstraps, MaxRE of the computed df.

| condition | time (s) | BTS | MaxRE of df |
|---|---|---|---|
| global | 14.49 (0.68) | 2 | 3.16e-04 |
| Insurance | 9.71 (0.10) | 1 | 2.09e-04 |
| Adult | 1.08 (0.01) | 0 | 3.38e-04 |
| Heart | 9.88 (0.12) | 1 | 4.93e-05 |
| Diabetes | 9.87 (0.10) | 1 | 9.47e-03 |
| Credit | 9.77 (0.04) | 1 | 5.06e-03 |
| Bank | 9.99 (0.21) | 1 | 6.17e-03 |

### Polynomials, case ranges, alpha = 0.05 (MaxRE over the range, mean over repetitions)

| input | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank |
|---|---|---|---|---|---|---|---|
| 1/df (stored) | 7 | 9.18e-04 | 9.52e-04 | 6.09e-04 | 1.01e-03 | 1.01e-03 | 1.01e-03 |
| 1/df (stored) | 15 | 4.59e-07 | 5.61e-07 | 4.84e-07 | 4.78e-07 | 4.88e-07 | 4.93e-07 |
| 1/df (stored) | 31 | 5.23e-07 | 4.43e-07 | 4.92e-07 | 4.63e-07 | 4.80e-07 | 4.46e-07 |
| 1/df (range fit) | 7 | 5.25e-07 | 6.17e-07 | 4.43e-07 | 4.60e-07 | 4.75e-07 | 5.35e-07 |
| 1/df (range fit) | 15 | 4.01e-07 | 4.43e-07 | 4.71e-07 | 4.24e-07 | 5.02e-07 | 5.21e-07 |
| 1/df (range fit) | 31 | 4.38e-07 | 5.20e-07 | 4.46e-07 | 4.69e-07 | 4.90e-07 | 4.23e-07 |
| df (HE division) | 7 | 7.52e-06 | 4.56e-07 | 5.56e-07 | 5.94e-07 | 5.35e-07 | 3.17e-06 |
| df (HE division) | 15 | 4.88e-07 | 4.54e-07 | 5.88e-07 | 4.48e-07 | 4.85e-07 | 4.61e-07 |
| df (HE division) | 31 | 4.87e-07 | 4.55e-07 | 5.88e-07 | 4.49e-07 | 4.85e-07 | 4.61e-07 |

### Polynomials, case ranges, alpha = 0.001 (MaxRE over the range, mean over repetitions)

| input | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank |
|---|---|---|---|---|---|---|---|
| 1/df (stored) | 7 | 4.33e+01 | 4.49e+01 | 2.85e+01 | 4.79e+01 | 4.78e+01 | 4.79e+01 |
| 1/df (stored) | 15 | 6.12e-04 | 7.30e-04 | 8.39e-04 | 9.54e-04 | 9.42e-04 | 9.50e-04 |
| 1/df (stored) | 31 | 2.66e-05 | 2.97e-05 | 3.16e-05 | 3.13e-05 | 3.09e-05 | 2.99e-05 |
| 1/df (range fit) | 7 | 1.20e-06 | 1.20e-06 | 1.05e-06 | 1.19e-06 | 1.06e-06 | 1.15e-06 |
| 1/df (range fit) | 15 | 1.20e-06 | 1.09e-06 | 1.25e-06 | 1.21e-06 | 1.11e-06 | 1.11e-06 |
| 1/df (range fit) | 31 | 1.27e-06 | 1.25e-06 | 1.20e-06 | 1.15e-06 | 1.13e-06 | 1.17e-06 |
| df (HE division) | 7 | 1.89e-05 | 1.12e-06 | 1.36e-06 | 1.45e-06 | 1.31e-06 | 7.75e-06 |
| df (HE division) | 15 | 1.20e-06 | 1.11e-06 | 1.45e-06 | 1.09e-06 | 1.18e-06 | 1.13e-06 |
| df (HE division) | 31 | 1.19e-06 | 1.11e-06 | 1.45e-06 | 1.09e-06 | 1.18e-06 | 1.13e-06 |

### Lookup tables in HE

Comparisons df >= d_j evaluated as u <= 1/d_j without division: one composite sign of (u - 1/d_j) / width for all
entries at once (packed in the slots, 32768 entries per ciphertext); floor and linear interpolation in 1/df are public
weighted slot sums of the step vector. One df per ciphertext: 20 interior points of the global range, 3 interior points
of each case range and the Welch df; 10 repetitions. MaxRE at alpha = 0.05 / 0.001.

| condition | table | entries | ciphertexts | compare (s) | BTS | floor MaxRE | interp MaxRE | max rel. diff HE vs plain |
|---|---|---|---|---|---|---|---|---|
| global | textbook | 38 | 1 | 10.26 (0.63) | 4 | 2.93e+00 / 4.05e+01 | 8.19e-01 / 1.54e+01 | 3.25e-05 |
| global | integer | 2000 | 1 | 10.20 (0.61) | 4 | 2.93e+00 / 4.05e+01 | 8.19e-01 / 1.54e+01 | 9.30e-05 |
| Insurance | textbook | 38 | 1 | 9.83 (0.77) | 4 | 1.77e-02 / 4.41e-02 | 4.98e-05 / 2.86e-04 | 3.09e-05 |
| Insurance | integer | 1064 | 1 | 10.16 (0.74) | 4 | 2.13e-05 / 5.23e-05 | 2.20e-07 / 5.46e-07 | 5.46e-07 |
| Adult | textbook | 38 | 1 | 10.13 (0.61) | 4 | 1.23e-03 / 3.04e-03 | 1.20e-06 / 3.84e-05 | 3.31e-05 |
| Adult | integer | 1024 | 1 | 10.44 (0.46) | 4 | 2.02e-06 / 4.88e-06 | 2.16e-07 / 4.73e-07 | 3.93e-06 |
| Heart | textbook | 38 | 1 | 10.25 (0.63) | 4 | 1.22e-02 / 3.04e-02 | 5.89e-05 / 3.54e-04 | 2.51e-05 |
| Heart | integer | 164 | 1 | 10.29 (0.70) | 4 | 6.16e-05 / 1.51e-04 | 1.56e-07 / 3.46e-07 | 3.44e-07 |
| Diabetes | textbook | 38 | 1 | 10.30 (0.78) | 4 | 2.37e-03 / 5.81e-03 | 4.32e-07 / 2.69e-05 | 2.78e-05 |
| Diabetes | integer | 54864 | 2 | 19.46 (1.16) | 4 | 8.85e-07 / 2.15e-06 | 1.99e-07 / 4.88e-07 | 2.14e-06 |
| Credit | textbook | 38 | 1 | 10.35 (0.65) | 4 | 2.31e-03 / 5.67e-03 | 4.96e-07 / 2.72e-05 | 2.49e-05 |
| Credit | integer | 23364 | 1 | 10.56 (0.79) | 4 | 8.98e-07 / 2.20e-06 | 1.51e-07 / 3.82e-07 | 2.20e-06 |
| Bank | textbook | 38 | 1 | 10.07 (1.39) | 4 | 2.33e-03 / 5.71e-03 | 6.65e-07 / 2.99e-05 | 2.95e-05 |
| Bank | integer | 39922 | 2 | 19.93 (1.08) | 4 | 1.10e-06 / 2.76e-06 | 2.10e-07 / 5.04e-07 | 2.75e-06 |

Weighted slot sum per output: 0.16 s on average.
