# Critical-value approximation: method comparison

Target c^2 = t_{1-alpha/2, df}^2 (two-sided), reference SciPy `t.ppf`; relative error. Script:
`experiments/tcrit/tcrit_experiment.py --part plain`, raw results `results/plaintext.csv` (all degrees, alphas, ranges, df bands).

Conditions (every method is fitted/tabulated and evaluated on the same df range, evaluation grid = 20001 log-spaced df plus every integer df):
- global: df in [1, 2000];
- case range: df in [min(n1, n2) - 1, n1 + n2 - 2], which always contains the Welch-Satterthwaite df: Insurance [273, 1336],
  Adult educational-num and Adult age [1023, 2046] (identical range, identical results).

Methods: Chebyshev interpolation with input 1/df or df (degrees 7..127; the global 1/df fit is the stored table
`coeffs/t_critical_coeffs.py`, case-range fits are refitted on the case range); lookup tables with the textbook df grid
(1..30, 40, 50, 60, 80, 100, 120, 1000; interpolation adds df = infinity) or every integer df, mapped by floor, round or
linear interpolation of c^2 in 1/df.

## alpha = 0.001

| method | degree | global MaxRE (MRE) | global MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 | Insurance MaxRE (MRE) | Adult MaxRE (MRE) |
|---|---|---|---|---|---|
| 1/df polynomial | 7 | 4.8e+01 (2.3e+01) | 8.5e-01 / 2.0e+01 / 3.9e+01 / 4.8e+01 | 3.1e-15 (9.0e-16) | 2.9e-15 (8.4e-16) |
| 1/df polynomial | 15 | 9.4e-04 (4.2e-04) | 1.4e-05 / 3.7e-04 / 7.6e-04 / 9.4e-04 | 1.9e-14 (5.4e-15) | 1.9e-14 (5.3e-15) |
| 1/df polynomial | 31 | 3.0e-08 (1.2e-08) | 9.6e-09 / 1.7e-08 / 2.6e-08 / 3.0e-08 | 5.9e-14 (6.8e-15) | 5.9e-14 (6.2e-15) |
| 1/df polynomial | 63 | 3.8e-07 (1.7e-07) | 8.6e-09 / 2.1e-07 / 3.1e-07 / 3.8e-07 | 2.9e-13 (1.7e-14) | 2.9e-13 (1.6e-14) |
| 1/df polynomial | 127 | 1.0e-07 (4.7e-08) | 9.2e-09 / 5.6e-08 / 8.7e-08 / 1.0e-07 | 5.3e-13 (5.8e-15) | 5.3e-13 (5.3e-15) |
| df polynomial | 7 | 1.0e+00 (2.5e-01) | 1.0e+00 / 9.8e-01 / 2.8e-01 / 1.1e-01 | 1.9e-05 (5.9e-06) | 8.9e-09 (4.1e-09) |
| df polynomial | 15 | 1.0e+00 (2.8e-01) | 1.0e+00 / 9.6e-01 / 6.5e-01 / 3.0e-01 | 8.1e-09 (2.6e-09) | 2.6e-14 (4.2e-15) |
| df polynomial | 31 | 6.2e+00 (1.2e+00) | 1.0e+00 / 6.2e+00 / 3.3e+00 / 1.7e+00 | 5.8e-14 (7.0e-15) | 5.9e-14 (6.2e-15) |
| df polynomial | 63 | 5.3e+01 (1.0e+01) | 1.1e+01 / 5.3e+01 / 4.5e+01 / 2.3e+01 | 2.9e-13 (1.7e-14) | 2.9e-13 (1.6e-14) |
| df polynomial | 127 | 8.7e+01 (1.7e+01) | 1.4e+01 / 8.7e+01 / 7.6e+01 / 3.9e+01 | 5.3e-13 (5.8e-15) | 5.3e-13 (5.3e-15) |
| textbook table, floor |  | 4.0e+02 (7.4e+00) | 4.0e+02 / 5.0e+00 / 6.9e-02 / 5.4e-02 | 4.5e-02 (3.1e-02) | 3.0e-03 (1.8e-03) |
| textbook table + inf, interp. in 1/df |  | 2.5e+01 (1.0e+00) | 2.5e+01 / 5.0e-01 / 7.0e-04 / 4.8e-04 | 3.1e-04 (1.2e-04) | 6.2e-06 (4.7e-06) |
| integer table, floor |  | 4.0e+02 (7.3e+00) | 4.0e+02 / 5.0e+00 / 6.9e-02 / 6.8e-03 | 8.0e-05 (1.1e-05) | 5.6e-06 (1.5e-06) |
| integer table, round |  | 5.8e+01 (7.9e-01) | 5.8e+01 / 2.0e+00 / 3.6e-02 / 3.4e-03 | 4.0e-05 (5.7e-06) | 2.8e-06 (7.3e-07) |
| integer table, interp. in 1/df |  | 2.5e+01 (1.0e+00) | 2.5e+01 / 5.0e-01 / 7.0e-04 / 8.0e-06 | 1.1e-09 (1.1e-10) | 5.6e-12 (1.2e-12) |

## alpha = 0.01

| method | degree | global MaxRE (MRE) | global MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 | Insurance MaxRE (MRE) | Adult MaxRE (MRE) |
|---|---|---|---|---|---|
| 1/df polynomial | 7 | 1.5e-01 (7.9e-02) | 1.5e-02 / 9.0e-02 / 1.3e-01 / 1.5e-01 | 3.4e-15 (9.5e-16) | 2.9e-15 (9.0e-16) |
| 1/df polynomial | 15 | 1.4e-07 (6.6e-08) | 1.3e-08 / 7.8e-08 / 1.2e-07 / 1.4e-07 | 1.8e-14 (5.4e-15) | 1.8e-14 (5.2e-15) |
| 1/df polynomial | 31 | 1.6e-08 (6.4e-09) | 9.7e-09 / 9.8e-09 / 1.3e-08 / 1.6e-08 | 5.9e-14 (6.8e-15) | 5.9e-14 (6.3e-15) |
| 1/df polynomial | 63 | 1.2e-08 (3.5e-09) | 5.2e-09 / 8.5e-09 / 6.4e-09 / 1.2e-08 | 2.9e-13 (1.7e-14) | 2.9e-13 (1.6e-14) |
| 1/df polynomial | 127 | 1.1e-08 (2.6e-09) | 5.9e-09 / 6.7e-09 / 5.3e-09 / 1.1e-08 | 5.3e-13 (5.9e-15) | 5.3e-13 (5.3e-15) |
| df polynomial | 7 | 1.0e+00 (2.1e-01) | 1.0e+00 / 9.1e-01 / 1.8e-01 / 6.7e-02 | 1.2e-05 (3.7e-06) | 5.7e-09 (2.6e-09) |
| df polynomial | 15 | 1.0e+00 (1.9e-01) | 1.0e+00 / 8.4e-01 / 3.2e-01 / 1.4e-01 | 5.1e-09 (1.6e-09) | 2.3e-14 (4.4e-15) |
| df polynomial | 31 | 1.6e+00 (3.5e-01) | 9.8e-01 / 1.6e+00 / 6.4e-01 / 3.1e-01 | 6.0e-14 (6.9e-15) | 5.9e-14 (6.2e-15) |
| df polynomial | 63 | 3.7e+00 (7.8e-01) | 2.8e+00 / 3.7e+00 / 2.5e+00 / 1.1e+00 | 2.9e-13 (1.7e-14) | 2.9e-13 (1.6e-14) |
| df polynomial | 127 | 2.8e+00 (5.3e-01) | 1.5e+00 / 2.8e+00 / 1.8e+00 / 8.0e-01 | 5.3e-13 (5.8e-15) | 5.3e-13 (5.3e-15) |
| textbook table, floor |  | 4.0e+01 (1.1e+00) | 4.0e+01 / 1.9e+00 / 4.1e-02 / 3.4e-02 | 2.9e-02 (2.0e-02) | 2.0e-03 (1.1e-03) |
| textbook table + inf, interp. in 1/df |  | 3.6e+00 (1.8e-01) | 3.6e+00 / 1.7e-01 / 2.7e-04 / 2.0e-04 | 1.3e-04 (5.2e-05) | 2.6e-06 (2.0e-06) |
| integer table, floor |  | 4.0e+01 (1.1e+00) | 4.0e+01 / 1.9e+00 / 4.1e-02 / 4.3e-03 | 5.1e-05 (7.4e-06) | 3.6e-06 (9.4e-07) |
| integer table, round |  | 1.2e+01 (2.3e-01) | 1.2e+01 / 9.2e-01 / 2.1e-02 / 2.2e-03 | 2.6e-05 (3.7e-06) | 1.8e-06 (4.7e-07) |
| integer table, interp. in 1/df |  | 3.6e+00 (1.8e-01) | 3.6e+00 / 1.7e-01 / 2.7e-04 / 3.3e-06 | 4.7e-10 (4.7e-11) | 2.4e-12 (5.2e-13) |

## alpha = 0.025

| method | degree | global MaxRE (MRE) | global MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 | Insurance MaxRE (MRE) | Adult MaxRE (MRE) |
|---|---|---|---|---|---|
| 1/df polynomial | 7 | 1.0e-02 (5.6e-03) | 1.8e-03 / 6.9e-03 / 9.3e-03 / 1.0e-02 | 3.4e-15 (8.7e-16) | 3.0e-15 (8.5e-16) |
| 1/df polynomial | 15 | 9.8e-09 (8.8e-10) | 9.8e-09 / 5.0e-09 / 1.2e-09 / 6.7e-09 | 1.9e-14 (5.5e-15) | 1.9e-14 (5.3e-15) |
| 1/df polynomial | 31 | 9.8e-09 (3.0e-09) | 8.5e-09 / 7.1e-09 / 6.1e-09 / 9.8e-09 | 6.0e-14 (6.8e-15) | 5.9e-14 (6.2e-15) |
| 1/df polynomial | 63 | 1.0e-08 (1.9e-09) | 9.8e-09 / 1.0e-08 / 3.8e-09 / 9.5e-09 | 2.9e-13 (1.7e-14) | 2.9e-13 (1.6e-14) |
| 1/df polynomial | 127 | 9.6e-09 (6.2e-10) | 7.3e-09 / 8.7e-09 / 1.3e-09 / 9.6e-09 | 5.3e-13 (5.9e-15) | 5.3e-13 (5.3e-15) |
| df polynomial | 7 | 9.9e-01 (1.8e-01) | 9.9e-01 / 8.4e-01 / 1.4e-01 / 5.2e-02 | 9.3e-06 (2.9e-06) | 4.5e-09 (2.0e-09) |
| df polynomial | 15 | 9.8e-01 (1.6e-01) | 9.8e-01 / 7.5e-01 / 2.2e-01 / 9.6e-02 | 3.9e-09 (1.2e-09) | 2.3e-14 (4.6e-15) |
| df polynomial | 31 | 9.4e-01 (2.1e-01) | 9.4e-01 / 9.4e-01 / 3.3e-01 / 1.5e-01 | 5.9e-14 (6.9e-15) | 5.9e-14 (6.2e-15) |
| df polynomial | 63 | 1.6e+00 (3.0e-01) | 1.4e+00 / 1.6e+00 / 7.8e-01 / 3.3e-01 | 2.9e-13 (1.7e-14) | 2.9e-13 (1.6e-14) |
| df polynomial | 127 | 7.4e-01 (1.4e-01) | 5.4e-01 / 7.4e-01 / 3.9e-01 / 1.7e-01 | 5.3e-13 (5.8e-15) | 5.3e-13 (5.3e-15) |
| textbook table, floor |  | 1.6e+01 (5.2e-01) | 1.6e+01 / 1.2e+00 / 3.2e-02 / 2.6e-02 | 2.2e-02 (1.6e-02) | 1.5e-03 (8.9e-04) |
| textbook table + inf, interp. in 1/df |  | 1.6e+00 (8.7e-02) | 1.6e+00 / 9.4e-02 / 1.7e-04 / 1.2e-04 | 8.3e-05 (3.3e-05) | 1.7e-06 (1.3e-06) |
| integer table, floor |  | 1.6e+01 (5.2e-01) | 1.6e+01 / 1.2e+00 / 3.2e-02 / 3.4e-03 | 4.0e-05 (5.8e-06) | 2.9e-06 (7.4e-07) |
| integer table, round |  | 6.0e+00 (1.4e-01) | 6.0e+00 / 6.3e-01 / 1.6e-02 / 1.7e-03 | 2.0e-05 (2.9e-06) | 1.4e-06 (3.7e-07) |
| integer table, interp. in 1/df |  | 1.6e+00 (8.7e-02) | 1.6e+00 / 9.4e-02 / 1.7e-04 / 2.1e-06 | 3.0e-10 (3.0e-11) | 1.5e-12 (3.3e-13) |

## alpha = 0.05

| method | degree | global MaxRE (MRE) | global MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 | Insurance MaxRE (MRE) | Adult MaxRE (MRE) |
|---|---|---|---|---|---|
| 1/df polynomial | 7 | 1.0e-03 (5.6e-04) | 2.7e-04 / 7.4e-04 / 9.4e-04 / 1.0e-03 | 2.7e-15 (7.7e-16) | 2.9e-15 (7.7e-16) |
| 1/df polynomial | 15 | 9.3e-09 (1.1e-09) | 9.8e-11 / 9.3e-09 / 3.6e-09 / 5.5e-09 | 1.9e-14 (5.5e-15) | 1.9e-14 (5.2e-15) |
| 1/df polynomial | 31 | 8.3e-09 (3.9e-10) | 1.9e-11 / 5.8e-09 / 2.4e-10 / 8.3e-09 | 5.9e-14 (6.9e-15) | 5.9e-14 (6.3e-15) |
| 1/df polynomial | 63 | 9.9e-09 (6.0e-10) | 2.7e-11 / 6.7e-09 / 1.8e-09 / 9.9e-09 | 2.9e-13 (1.7e-14) | 2.9e-13 (1.6e-14) |
| 1/df polynomial | 127 | 8.4e-09 (2.0e-10) | 1.4e-11 / 8.4e-09 / 6.5e-10 / 5.1e-09 | 5.3e-13 (5.9e-15) | 5.3e-13 (5.4e-15) |
| df polynomial | 7 | 9.7e-01 (1.6e-01) | 9.7e-01 / 7.6e-01 / 1.1e-01 / 4.1e-02 | 7.5e-06 (2.3e-06) | 3.6e-09 (1.6e-09) |
| df polynomial | 15 | 9.6e-01 (1.4e-01) | 9.6e-01 / 6.5e-01 / 1.7e-01 / 7.1e-02 | 3.1e-09 (9.9e-10) | 2.1e-14 (4.7e-15) |
| df polynomial | 31 | 8.9e-01 (1.5e-01) | 8.9e-01 / 6.0e-01 / 1.9e-01 / 8.9e-02 | 5.9e-14 (6.9e-15) | 5.9e-14 (6.2e-15) |
| df polynomial | 63 | 8.7e-01 (1.5e-01) | 8.3e-01 / 8.7e-01 / 3.3e-01 / 1.3e-01 | 2.9e-13 (1.7e-14) | 2.9e-13 (1.6e-14) |
| df polynomial | 127 | 2.8e-01 (5.1e-02) | 2.4e-01 / 2.8e-01 / 1.2e-01 / 4.9e-02 | 5.3e-13 (5.7e-15) | 5.3e-13 (5.3e-15) |
| textbook table, floor |  | 7.7e+00 (3.0e-01) | 7.7e+00 / 8.3e-01 / 2.5e-02 / 2.1e-02 | 1.8e-02 (1.3e-02) | 1.2e-03 (7.2e-04) |
| textbook table + inf, interp. in 1/df |  | 8.3e-01 (4.6e-02) | 8.3e-01 / 5.6e-02 / 1.1e-04 / 7.9e-05 | 5.5e-05 (2.1e-05) | 1.1e-06 (8.2e-07) |
| integer table, floor |  | 7.7e+00 (2.9e-01) | 7.7e+00 / 8.3e-01 / 2.5e-02 / 2.7e-03 | 3.2e-05 (4.7e-06) | 2.3e-06 (6.0e-07) |
| integer table, round |  | 3.5e+00 (9.3e-02) | 3.5e+00 / 4.5e-01 / 1.3e-02 / 1.3e-03 | 1.6e-05 (2.3e-06) | 1.1e-06 (3.0e-07) |
| integer table, interp. in 1/df |  | 8.3e-01 (4.6e-02) | 8.3e-01 / 5.6e-02 / 1.1e-04 / 1.3e-06 | 2.0e-10 (2.0e-11) | 9.9e-13 (2.1e-13) |

## alpha = 0.1

| method | degree | global MaxRE (MRE) | global MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 | Insurance MaxRE (MRE) | Adult MaxRE (MRE) |
|---|---|---|---|---|---|
| 1/df polynomial | 7 | 6.7e-05 (3.9e-05) | 2.6e-05 / 5.3e-05 / 6.3e-05 / 6.7e-05 | 3.3e-15 (8.1e-16) | 2.6e-15 (8.8e-16) |
| 1/df polynomial | 15 | 1.1e-08 (5.8e-10) | 1.7e-10 / 7.3e-09 / 1.1e-08 / 1.0e-08 | 1.9e-14 (5.4e-15) | 1.9e-14 (5.2e-15) |
| 1/df polynomial | 31 | 1.5e-08 (8.0e-10) | 6.8e-11 / 8.2e-09 / 1.5e-08 / 1.4e-08 | 5.9e-14 (6.8e-15) | 5.9e-14 (6.3e-15) |
| 1/df polynomial | 63 | 9.7e-09 (4.9e-10) | 9.6e-11 / 9.4e-09 / 9.7e-09 / 9.3e-09 | 2.9e-13 (1.7e-14) | 2.9e-13 (1.6e-14) |
| 1/df polynomial | 127 | 1.5e-08 (4.6e-10) | 3.0e-11 / 7.0e-09 / 6.2e-09 / 1.5e-08 | 5.3e-13 (5.9e-15) | 5.3e-13 (5.3e-15) |
| df polynomial | 7 | 9.2e-01 (1.4e-01) | 9.2e-01 / 6.4e-01 / 8.6e-02 / 3.0e-02 | 5.7e-06 (1.8e-06) | 2.7e-09 (1.3e-09) |
| df polynomial | 15 | 9.0e-01 (1.1e-01) | 9.0e-01 / 5.3e-01 / 1.2e-01 / 4.9e-02 | 2.4e-09 (7.5e-10) | 2.0e-14 (4.8e-15) |
| df polynomial | 31 | 7.8e-01 (9.8e-02) | 7.8e-01 / 3.6e-01 / 1.1e-01 / 5.0e-02 | 5.9e-14 (6.8e-15) | 5.9e-14 (6.2e-15) |
| df polynomial | 63 | 4.5e-01 (7.3e-02) | 4.5e-01 / 4.5e-01 / 1.3e-01 / 5.4e-02 | 2.9e-13 (1.7e-14) | 2.9e-13 (1.6e-14) |
| df polynomial | 127 | 1.1e-01 (1.8e-02) | 1.0e-01 / 1.1e-01 / 3.4e-02 / 1.4e-02 | 5.3e-13 (5.7e-15) | 5.3e-13 (5.3e-15) |
| textbook table, floor |  | 3.7e+00 (1.6e-01) | 3.7e+00 / 5.4e-01 / 1.9e-02 / 1.6e-02 | 1.4e-02 (9.6e-03) | 9.5e-04 (5.5e-04) |
| textbook table + inf, interp. in 1/df |  | 3.9e-01 (2.2e-02) | 3.9e-01 / 2.9e-02 / 6.1e-05 / 4.7e-05 | 3.2e-05 (1.3e-05) | 6.5e-07 (4.9e-07) |
| integer table, floor |  | 3.7e+00 (1.6e-01) | 3.7e+00 / 5.4e-01 / 1.9e-02 / 2.0e-03 | 2.5e-05 (3.6e-06) | 1.8e-06 (4.6e-07) |
| integer table, round |  | 1.9e+00 (5.8e-02) | 1.9e+00 / 3.0e-01 / 9.7e-03 / 1.0e-03 | 1.2e-05 (1.8e-06) | 8.8e-07 (2.3e-07) |
| integer table, interp. in 1/df |  | 3.9e-01 (2.2e-02) | 3.9e-01 / 2.9e-02 / 6.1e-05 / 7.8e-07 | 1.2e-10 (1.2e-11) | 5.9e-13 (1.3e-13) |
## Observations

- Global range [1, 2000]: the 1/df polynomial reaches MaxRE <= 1.1e-8 from degree 15 for alpha >= 0.01 (alpha = 0.001: 9.4e-4 at
  degree 15, 3.0e-8 at degree 31). The df polynomial stays at MaxRE >= 0.28 up to degree 127 (it does not converge on this range:
  c^2 changes by 4e5 / 11 (alpha = 0.001) or 161 / 3.8 (alpha = 0.05) between df = 1 and df = 2000).
- Case ranges: both inputs are accurate. At degree 7 the 1/df polynomial is at ~3e-15 and the df polynomial at 7.5e-6 (Insurance,
  alpha = 0.05) / 3.6e-9 (Adult); from degree 31 both are at the double-precision level (~6e-14).
- Lookup tables: on the global range the error is dominated by df < 2 for every mapping (interpolation in 1/df: 0.83 at
  alpha = 0.05); restricted to df >= 30, integer-table interpolation gives 1.3e-6 (alpha = 0.05). On the case ranges, integer-table
  interpolation is as accurate as the polynomials (2e-10 Insurance, 1e-12 Adult at alpha = 0.05), floor / round give 1e-6 to 3e-5,
  and the textbook grid 1e-3 to 2e-2 (floor) or 1e-6 to 6e-5 (interpolation with df = infinity).
- The 1/df polynomial at alpha = 0.001 is not monotone in the degree on the global range (degree 31: 3.0e-8, degree 63: 3.8e-7,
  degree 127: 1.0e-7).

## HE evaluation

`tcrit_experiment.py --part he`, raw results `results/he_raw.csv`: `HEHypothesisTesting._critical_value_from_inv_df` (stored 1/df
coefficients, u = 1/df in [1/2000, 1]) on encrypted u for df log-spaced in [1, 2000] (one ciphertext, 32768 slots), input
level 12, 10 repetitions. MaxRE against the exact c^2; time excludes encryption and decryption.

| degree | levels | time (s) | MaxRE alpha = 0.001 | 0.01 | 0.025 | 0.05 | 0.1 |
|---|---|---|---|---|---|---|---|
| 7 | 12 -> 8 | 0.226 +- 0.023 | 4.8e+01 | 1.5e-01 | 1.0e-02 | 1.0e-03 | 6.7e-05 |
| 15 | 12 -> 7 | 0.329 +- 0.023 | 9.5e-04 | 1.4e-06 | 1.1e-06 | 9.4e-07 | 6.0e-07 |
| 31 | 12 -> 6 | 0.598 +- 0.088 | 3.0e-05 | 1.3e-06 | 1.0e-06 | 8.4e-07 | 5.9e-07 |
| 63 | 12 -> 5 | 0.999 +- 0.239 | 2.6e-05 | 1.3e-06 | 1.1e-06 | 8.2e-07 | 5.5e-07 |
| 127 | 12 -> 4 | 1.589 +- 0.586 | 3.0e-05 | 1.4e-06 | 9.5e-07 | 7.0e-07 | 4.9e-07 |

From degree 15 the HE error is at the CKKS noise level (5e-7 to 1.5e-6) except alpha = 0.001, where it equals the
plaintext polynomial error (9.4e-4 at degree 15, ~3e-5 from degree 31).
