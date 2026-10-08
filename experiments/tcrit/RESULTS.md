# Critical-value approximation: method comparison

Target c^2 = t_{1-alpha/2, df}^2 (two-sided), reference SciPy `t.ppf`; relative error. Script:
`experiments/tcrit/tcrit_experiment.py --part plain`, raw results `results/plaintext.csv` (all degrees, alphas, ranges, df bands).

Conditions (every method is fitted/tabulated and evaluated on the same df range; evaluation grid = 20001 log-spaced df plus every
integer df in the range; the global range has 21,999 points):
- global: df in [1, 2000];
- case range: df in [min(n1, n2) - 1, n1 + n2 - 2] for each dataset comparison (public; it always contains the Welch df):
  Insurance [273, 1336], Adult [1023, 2046], Heart [138, 301], Diabetes [11356, 66219], Credit [6635, 29998], Bank [5288, 45209], Wine [1598, 6495] (Adult: both comparisons).

Methods: Chebyshev interpolation with input 1/df or df (degrees 7..127; the global 1/df fit is the stored table
`coeffs/t_critical_coeffs.py`, case-range fits are refitted on the case range); lookup tables with the textbook df grid
(1..30, 40, 50, 60, 80, 100, 120, 1000; interpolation adds df = infinity) or every integer df up to the largest range, mapped
by floor, round or linear interpolation of c^2 in 1/df.

## alpha = 0.001

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 |
|---|---|---|---|
| 1/df polynomial | 7 | 4.8e+01 (2.3e+01) | 8.5e-01 / 2.0e+01 / 3.9e+01 / 4.8e+01 |
| 1/df polynomial | 15 | 9.4e-04 (4.2e-04) | 1.4e-05 / 3.7e-04 / 7.6e-04 / 9.4e-04 |
| 1/df polynomial | 31 | 3.0e-08 (1.2e-08) | 9.6e-09 / 1.7e-08 / 2.6e-08 / 3.0e-08 |
| 1/df polynomial | 63 | 3.8e-07 (1.7e-07) | 8.6e-09 / 2.1e-07 / 3.1e-07 / 3.8e-07 |
| 1/df polynomial | 127 | 1.0e-07 (4.7e-08) | 9.2e-09 / 5.6e-08 / 8.7e-08 / 1.0e-07 |
| df polynomial | 7 | 1.0e+00 (2.5e-01) | 1.0e+00 / 9.8e-01 / 2.8e-01 / 1.1e-01 |
| df polynomial | 15 | 1.0e+00 (2.8e-01) | 1.0e+00 / 9.6e-01 / 6.5e-01 / 3.0e-01 |
| df polynomial | 31 | 6.2e+00 (1.2e+00) | 1.0e+00 / 6.2e+00 / 3.3e+00 / 1.7e+00 |
| df polynomial | 63 | 5.3e+01 (1.0e+01) | 1.1e+01 / 5.3e+01 / 4.5e+01 / 2.3e+01 |
| df polynomial | 127 | 8.7e+01 (1.7e+01) | 1.4e+01 / 8.7e+01 / 7.6e+01 / 3.9e+01 |
| textbook table, floor |  | 4.0e+02 (7.4e+00) | 4.0e+02 / 5.0e+00 / 6.9e-02 / 5.4e-02 |
| textbook table + inf, interp. in 1/df |  | 2.5e+01 (1.0e+00) | 2.5e+01 / 5.0e-01 / 7.0e-04 / 4.8e-04 |
| integer table, floor |  | 4.0e+02 (7.3e+00) | 4.0e+02 / 5.0e+00 / 6.9e-02 / 6.8e-03 |
| integer table, round |  | 5.8e+01 (7.9e-01) | 5.8e+01 / 2.0e+00 / 3.6e-02 / 3.4e-03 |
| integer table, interp. in 1/df |  | 2.5e+01 (1.0e+00) | 2.5e+01 / 5.0e-01 / 7.0e-04 / 8.0e-06 |

Case ranges (MaxRE):

| method | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank | Wine |
|---|---|---|---|---|---|---|---|---|
| 1/df polynomial | 7 | 3.1e-15 | 2.9e-15 | 2.8e-15 | 2.6e-15 | 3.0e-15 | 3.4e-15 | 3.1e-15 |
| 1/df polynomial | 15 | 1.9e-14 | 1.9e-14 | 1.9e-14 | 1.9e-14 | 1.9e-14 | 1.9e-14 | 1.9e-14 |
| 1/df polynomial | 31 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 |
| 1/df polynomial | 63 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 |
| 1/df polynomial | 127 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 |
| df polynomial | 7 | 1.9e-05 | 8.9e-09 | 1.9e-07 | 9.1e-07 | 5.1e-07 | 7.5e-06 | 1.2e-06 |
| df polynomial | 15 | 8.1e-09 | 2.6e-14 | 4.3e-13 | 7.9e-10 | 1.4e-10 | 2.5e-08 | 2.1e-10 |
| df polynomial | 31 | 5.8e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 2.5e-13 | 5.9e-14 |
| df polynomial | 63 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 |
| df polynomial | 127 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 |
| textbook table, floor |  | 4.5e-02 | 3.0e-03 | 3.1e-02 | 5.8e-03 | 5.7e-03 | 5.8e-03 | 5.0e-03 |
| textbook table + inf, interp. in 1/df |  | 3.1e-04 | 6.2e-06 | 3.4e-04 | 2.0e-06 | 3.2e-06 | 3.8e-06 | 6.2e-06 |
| integer table, floor |  | 8.0e-05 | 5.6e-06 | 3.1e-04 | 4.4e-08 | 1.3e-07 | 2.1e-07 | 2.3e-06 |
| integer table, round |  | 4.0e-05 | 2.8e-06 | 1.6e-04 | 2.3e-08 | 6.7e-08 | 1.0e-07 | 1.2e-06 |
| integer table, interp. in 1/df |  | 1.1e-09 | 5.6e-12 | 1.7e-08 | 8.2e-16 | 3.4e-15 | 8.2e-15 | 9.4e-13 |

## alpha = 0.01

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 |
|---|---|---|---|
| 1/df polynomial | 7 | 1.5e-01 (7.9e-02) | 1.5e-02 / 9.0e-02 / 1.3e-01 / 1.5e-01 |
| 1/df polynomial | 15 | 1.4e-07 (6.6e-08) | 1.3e-08 / 7.8e-08 / 1.2e-07 / 1.4e-07 |
| 1/df polynomial | 31 | 1.6e-08 (6.4e-09) | 9.7e-09 / 9.8e-09 / 1.3e-08 / 1.6e-08 |
| 1/df polynomial | 63 | 1.2e-08 (3.5e-09) | 5.2e-09 / 8.5e-09 / 6.4e-09 / 1.2e-08 |
| 1/df polynomial | 127 | 1.1e-08 (2.6e-09) | 5.9e-09 / 6.7e-09 / 5.3e-09 / 1.1e-08 |
| df polynomial | 7 | 1.0e+00 (2.1e-01) | 1.0e+00 / 9.1e-01 / 1.8e-01 / 6.7e-02 |
| df polynomial | 15 | 1.0e+00 (1.9e-01) | 1.0e+00 / 8.4e-01 / 3.2e-01 / 1.4e-01 |
| df polynomial | 31 | 1.6e+00 (3.5e-01) | 9.8e-01 / 1.6e+00 / 6.4e-01 / 3.1e-01 |
| df polynomial | 63 | 3.7e+00 (7.8e-01) | 2.8e+00 / 3.7e+00 / 2.5e+00 / 1.1e+00 |
| df polynomial | 127 | 2.8e+00 (5.3e-01) | 1.5e+00 / 2.8e+00 / 1.8e+00 / 8.0e-01 |
| textbook table, floor |  | 4.0e+01 (1.1e+00) | 4.0e+01 / 1.9e+00 / 4.1e-02 / 3.4e-02 |
| textbook table + inf, interp. in 1/df |  | 3.6e+00 (1.8e-01) | 3.6e+00 / 1.7e-01 / 2.7e-04 / 2.0e-04 |
| integer table, floor |  | 4.0e+01 (1.1e+00) | 4.0e+01 / 1.9e+00 / 4.1e-02 / 4.3e-03 |
| integer table, round |  | 1.2e+01 (2.3e-01) | 1.2e+01 / 9.2e-01 / 2.1e-02 / 2.2e-03 |
| integer table, interp. in 1/df |  | 3.6e+00 (1.8e-01) | 3.6e+00 / 1.7e-01 / 2.7e-04 / 3.3e-06 |

Case ranges (MaxRE):

| method | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank | Wine |
|---|---|---|---|---|---|---|---|---|
| 1/df polynomial | 7 | 3.4e-15 | 2.9e-15 | 2.8e-15 | 3.3e-15 | 2.7e-15 | 3.1e-15 | 2.8e-15 |
| 1/df polynomial | 15 | 1.8e-14 | 1.8e-14 | 1.9e-14 | 1.9e-14 | 1.8e-14 | 1.9e-14 | 1.8e-14 |
| 1/df polynomial | 31 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 |
| 1/df polynomial | 63 | 2.9e-13 | 2.9e-13 | 3.0e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 |
| 1/df polynomial | 127 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 |
| df polynomial | 7 | 1.2e-05 | 5.7e-09 | 1.2e-07 | 5.8e-07 | 3.3e-07 | 4.8e-06 | 8.0e-07 |
| df polynomial | 15 | 5.1e-09 | 2.3e-14 | 2.6e-13 | 5.1e-10 | 9.3e-11 | 1.6e-08 | 1.3e-10 |
| df polynomial | 31 | 6.0e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 1.6e-13 | 5.9e-14 |
| df polynomial | 63 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 |
| df polynomial | 127 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 |
| textbook table, floor |  | 2.9e-02 | 2.0e-03 | 2.0e-02 | 3.8e-03 | 3.7e-03 | 3.7e-03 | 3.2e-03 |
| textbook table + inf, interp. in 1/df |  | 1.3e-04 | 2.6e-06 | 1.4e-04 | 8.5e-07 | 1.3e-06 | 1.6e-06 | 2.6e-06 |
| integer table, floor |  | 5.1e-05 | 3.6e-06 | 2.0e-04 | 2.9e-08 | 8.4e-08 | 1.3e-07 | 1.5e-06 |
| integer table, round |  | 2.6e-05 | 1.8e-06 | 1.0e-04 | 1.5e-08 | 4.3e-08 | 6.8e-08 | 7.4e-07 |
| integer table, interp. in 1/df |  | 4.7e-10 | 2.4e-12 | 7.3e-09 | 8.0e-16 | 1.7e-15 | 3.7e-15 | 4.0e-13 |

## alpha = 0.025

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 |
|---|---|---|---|
| 1/df polynomial | 7 | 1.0e-02 (5.6e-03) | 1.8e-03 / 6.9e-03 / 9.3e-03 / 1.0e-02 |
| 1/df polynomial | 15 | 9.8e-09 (8.8e-10) | 9.8e-09 / 5.0e-09 / 1.2e-09 / 6.7e-09 |
| 1/df polynomial | 31 | 9.8e-09 (3.0e-09) | 8.5e-09 / 7.1e-09 / 6.1e-09 / 9.8e-09 |
| 1/df polynomial | 63 | 1.0e-08 (1.9e-09) | 9.8e-09 / 1.0e-08 / 3.8e-09 / 9.5e-09 |
| 1/df polynomial | 127 | 9.6e-09 (6.2e-10) | 7.3e-09 / 8.7e-09 / 1.3e-09 / 9.6e-09 |
| df polynomial | 7 | 9.9e-01 (1.8e-01) | 9.9e-01 / 8.4e-01 / 1.4e-01 / 5.2e-02 |
| df polynomial | 15 | 9.8e-01 (1.6e-01) | 9.8e-01 / 7.5e-01 / 2.2e-01 / 9.6e-02 |
| df polynomial | 31 | 9.4e-01 (2.1e-01) | 9.4e-01 / 9.4e-01 / 3.3e-01 / 1.5e-01 |
| df polynomial | 63 | 1.6e+00 (3.0e-01) | 1.4e+00 / 1.6e+00 / 7.8e-01 / 3.3e-01 |
| df polynomial | 127 | 7.4e-01 (1.4e-01) | 5.4e-01 / 7.4e-01 / 3.9e-01 / 1.7e-01 |
| textbook table, floor |  | 1.6e+01 (5.2e-01) | 1.6e+01 / 1.2e+00 / 3.2e-02 / 2.6e-02 |
| textbook table + inf, interp. in 1/df |  | 1.6e+00 (8.7e-02) | 1.6e+00 / 9.4e-02 / 1.7e-04 / 1.2e-04 |
| integer table, floor |  | 1.6e+01 (5.2e-01) | 1.6e+01 / 1.2e+00 / 3.2e-02 / 3.4e-03 |
| integer table, round |  | 6.0e+00 (1.4e-01) | 6.0e+00 / 6.3e-01 / 1.6e-02 / 1.7e-03 |
| integer table, interp. in 1/df |  | 1.6e+00 (8.7e-02) | 1.6e+00 / 9.4e-02 / 1.7e-04 / 2.1e-06 |

Case ranges (MaxRE):

| method | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank | Wine |
|---|---|---|---|---|---|---|---|---|
| 1/df polynomial | 7 | 3.4e-15 | 3.0e-15 | 3.2e-15 | 3.0e-15 | 3.0e-15 | 2.8e-15 | 3.2e-15 |
| 1/df polynomial | 15 | 1.9e-14 | 1.9e-14 | 1.8e-14 | 1.9e-14 | 1.9e-14 | 1.9e-14 | 1.9e-14 |
| 1/df polynomial | 31 | 6.0e-14 | 5.9e-14 | 6.0e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 |
| 1/df polynomial | 63 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 |
| 1/df polynomial | 127 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 |
| df polynomial | 7 | 9.3e-06 | 4.5e-09 | 8.9e-08 | 4.6e-07 | 2.6e-07 | 3.8e-06 | 6.3e-07 |
| df polynomial | 15 | 3.9e-09 | 2.3e-14 | 2.0e-13 | 4.0e-10 | 7.3e-11 | 1.3e-08 | 1.0e-10 |
| df polynomial | 31 | 5.9e-14 | 5.9e-14 | 6.0e-14 | 5.9e-14 | 5.9e-14 | 1.3e-13 | 5.9e-14 |
| df polynomial | 63 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 |
| df polynomial | 127 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 |
| textbook table, floor |  | 2.2e-02 | 1.5e-03 | 1.5e-02 | 3.0e-03 | 2.9e-03 | 3.0e-03 | 2.6e-03 |
| textbook table + inf, interp. in 1/df |  | 8.3e-05 | 1.7e-06 | 9.0e-05 | 5.3e-07 | 8.5e-07 | 1.0e-06 | 1.7e-06 |
| integer table, floor |  | 4.0e-05 | 2.9e-06 | 1.6e-04 | 2.3e-08 | 6.6e-08 | 1.0e-07 | 1.2e-06 |
| integer table, round |  | 2.0e-05 | 1.4e-06 | 7.9e-05 | 1.1e-08 | 3.4e-08 | 5.3e-08 | 5.9e-07 |
| integer table, interp. in 1/df |  | 3.0e-10 | 1.5e-12 | 4.6e-09 | 1.1e-15 | 1.4e-15 | 2.5e-15 | 2.5e-13 |

## alpha = 0.05

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 |
|---|---|---|---|
| 1/df polynomial | 7 | 1.0e-03 (5.6e-04) | 2.7e-04 / 7.4e-04 / 9.4e-04 / 1.0e-03 |
| 1/df polynomial | 15 | 9.3e-09 (1.1e-09) | 9.8e-11 / 9.3e-09 / 3.6e-09 / 5.5e-09 |
| 1/df polynomial | 31 | 8.3e-09 (3.9e-10) | 1.9e-11 / 5.8e-09 / 2.4e-10 / 8.3e-09 |
| 1/df polynomial | 63 | 9.9e-09 (6.0e-10) | 2.7e-11 / 6.7e-09 / 1.8e-09 / 9.9e-09 |
| 1/df polynomial | 127 | 8.4e-09 (2.0e-10) | 1.4e-11 / 8.4e-09 / 6.5e-10 / 5.1e-09 |
| df polynomial | 7 | 9.7e-01 (1.6e-01) | 9.7e-01 / 7.6e-01 / 1.1e-01 / 4.1e-02 |
| df polynomial | 15 | 9.6e-01 (1.4e-01) | 9.6e-01 / 6.5e-01 / 1.7e-01 / 7.1e-02 |
| df polynomial | 31 | 8.9e-01 (1.5e-01) | 8.9e-01 / 6.0e-01 / 1.9e-01 / 8.9e-02 |
| df polynomial | 63 | 8.7e-01 (1.5e-01) | 8.3e-01 / 8.7e-01 / 3.3e-01 / 1.3e-01 |
| df polynomial | 127 | 2.8e-01 (5.1e-02) | 2.4e-01 / 2.8e-01 / 1.2e-01 / 4.9e-02 |
| textbook table, floor |  | 7.7e+00 (3.0e-01) | 7.7e+00 / 8.3e-01 / 2.5e-02 / 2.1e-02 |
| textbook table + inf, interp. in 1/df |  | 8.3e-01 (4.6e-02) | 8.3e-01 / 5.6e-02 / 1.1e-04 / 7.9e-05 |
| integer table, floor |  | 7.7e+00 (2.9e-01) | 7.7e+00 / 8.3e-01 / 2.5e-02 / 2.7e-03 |
| integer table, round |  | 3.5e+00 (9.3e-02) | 3.5e+00 / 4.5e-01 / 1.3e-02 / 1.3e-03 |
| integer table, interp. in 1/df |  | 8.3e-01 (4.6e-02) | 8.3e-01 / 5.6e-02 / 1.1e-04 / 1.3e-06 |

Case ranges (MaxRE):

| method | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank | Wine |
|---|---|---|---|---|---|---|---|---|
| 1/df polynomial | 7 | 2.7e-15 | 2.9e-15 | 2.7e-15 | 2.9e-15 | 3.0e-15 | 3.1e-15 | 2.9e-15 |
| 1/df polynomial | 15 | 1.9e-14 | 1.9e-14 | 1.9e-14 | 1.9e-14 | 1.9e-14 | 1.9e-14 | 1.8e-14 |
| 1/df polynomial | 31 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 6.0e-14 |
| 1/df polynomial | 63 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 |
| 1/df polynomial | 127 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 |
| df polynomial | 7 | 7.5e-06 | 3.6e-09 | 7.1e-08 | 3.7e-07 | 2.1e-07 | 3.1e-06 | 5.0e-07 |
| df polynomial | 15 | 3.1e-09 | 2.1e-14 | 1.6e-13 | 3.2e-10 | 5.9e-11 | 1.0e-08 | 8.4e-11 |
| df polynomial | 31 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 1.0e-13 | 5.9e-14 |
| df polynomial | 63 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 |
| df polynomial | 127 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 |
| textbook table, floor |  | 1.8e-02 | 1.2e-03 | 1.2e-02 | 2.4e-03 | 2.3e-03 | 2.4e-03 | 2.1e-03 |
| textbook table + inf, interp. in 1/df |  | 5.5e-05 | 1.1e-06 | 5.9e-05 | 3.5e-07 | 5.6e-07 | 6.7e-07 | 1.1e-06 |
| integer table, floor |  | 3.2e-05 | 2.3e-06 | 1.3e-04 | 1.8e-08 | 5.3e-08 | 8.4e-08 | 9.3e-07 |
| integer table, round |  | 1.6e-05 | 1.1e-06 | 6.3e-05 | 9.2e-09 | 2.7e-08 | 4.3e-08 | 4.7e-07 |
| integer table, interp. in 1/df |  | 2.0e-10 | 9.9e-13 | 3.0e-09 | 5.8e-16 | 9.2e-16 | 1.6e-15 | 1.7e-13 |

## alpha = 0.1

Global range:

| method | degree | MaxRE (MRE) | MaxRE df 1-2 / 2-10 / 10-30 / 30-2000 |
|---|---|---|---|
| 1/df polynomial | 7 | 6.7e-05 (3.9e-05) | 2.6e-05 / 5.3e-05 / 6.3e-05 / 6.7e-05 |
| 1/df polynomial | 15 | 1.1e-08 (5.8e-10) | 1.7e-10 / 7.3e-09 / 1.1e-08 / 1.0e-08 |
| 1/df polynomial | 31 | 1.5e-08 (8.0e-10) | 6.8e-11 / 8.2e-09 / 1.5e-08 / 1.4e-08 |
| 1/df polynomial | 63 | 9.7e-09 (4.9e-10) | 9.6e-11 / 9.4e-09 / 9.7e-09 / 9.3e-09 |
| 1/df polynomial | 127 | 1.5e-08 (4.6e-10) | 3.0e-11 / 7.0e-09 / 6.2e-09 / 1.5e-08 |
| df polynomial | 7 | 9.2e-01 (1.4e-01) | 9.2e-01 / 6.4e-01 / 8.6e-02 / 3.0e-02 |
| df polynomial | 15 | 9.0e-01 (1.1e-01) | 9.0e-01 / 5.3e-01 / 1.2e-01 / 4.9e-02 |
| df polynomial | 31 | 7.8e-01 (9.8e-02) | 7.8e-01 / 3.6e-01 / 1.1e-01 / 5.0e-02 |
| df polynomial | 63 | 4.5e-01 (7.3e-02) | 4.5e-01 / 4.5e-01 / 1.3e-01 / 5.4e-02 |
| df polynomial | 127 | 1.1e-01 (1.8e-02) | 1.0e-01 / 1.1e-01 / 3.4e-02 / 1.4e-02 |
| textbook table, floor |  | 3.7e+00 (1.6e-01) | 3.7e+00 / 5.4e-01 / 1.9e-02 / 1.6e-02 |
| textbook table + inf, interp. in 1/df |  | 3.9e-01 (2.2e-02) | 3.9e-01 / 2.9e-02 / 6.1e-05 / 4.7e-05 |
| integer table, floor |  | 3.7e+00 (1.6e-01) | 3.7e+00 / 5.4e-01 / 1.9e-02 / 2.0e-03 |
| integer table, round |  | 1.9e+00 (5.8e-02) | 1.9e+00 / 3.0e-01 / 9.7e-03 / 1.0e-03 |
| integer table, interp. in 1/df |  | 3.9e-01 (2.2e-02) | 3.9e-01 / 2.9e-02 / 6.1e-05 / 7.8e-07 |

Case ranges (MaxRE):

| method | degree | Insurance | Adult | Heart | Diabetes | Credit | Bank | Wine |
|---|---|---|---|---|---|---|---|---|
| 1/df polynomial | 7 | 3.3e-15 | 2.6e-15 | 2.8e-15 | 3.0e-15 | 3.3e-15 | 3.6e-15 | 3.1e-15 |
| 1/df polynomial | 15 | 1.9e-14 | 1.9e-14 | 1.9e-14 | 1.9e-14 | 1.8e-14 | 1.9e-14 | 1.9e-14 |
| 1/df polynomial | 31 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 |
| 1/df polynomial | 63 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 |
| 1/df polynomial | 127 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 |
| df polynomial | 7 | 5.7e-06 | 2.7e-09 | 5.3e-08 | 2.8e-07 | 1.6e-07 | 2.3e-06 | 3.9e-07 |
| df polynomial | 15 | 2.4e-09 | 2.0e-14 | 1.3e-13 | 2.5e-10 | 4.5e-11 | 7.8e-09 | 6.4e-11 |
| df polynomial | 31 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 5.9e-14 | 8.0e-14 | 5.9e-14 |
| df polynomial | 63 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 | 2.9e-13 |
| df polynomial | 127 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 | 5.3e-13 |
| textbook table, floor |  | 1.4e-02 | 9.5e-04 | 9.4e-03 | 1.8e-03 | 1.8e-03 | 1.8e-03 | 1.6e-03 |
| textbook table + inf, interp. in 1/df |  | 3.2e-05 | 6.5e-07 | 3.5e-05 | 2.1e-07 | 3.3e-07 | 4.0e-07 | 6.5e-07 |
| integer table, floor |  | 2.5e-05 | 1.8e-06 | 9.7e-05 | 1.4e-08 | 4.1e-08 | 6.4e-08 | 7.1e-07 |
| integer table, round |  | 1.2e-05 | 8.8e-07 | 4.8e-05 | 7.1e-09 | 2.1e-08 | 3.3e-08 | 3.6e-07 |
| integer table, interp. in 1/df |  | 1.2e-10 | 5.9e-13 | 1.8e-09 | 6.6e-16 | 8.2e-16 | 1.3e-15 | 9.8e-14 |

## Observations

- Global range [1, 2000]: the 1/df polynomial reaches MaxRE <= 1.4e-7 from degree 15 for alpha >= 0.01 (<= 1.1e-8 for
  alpha >= 0.025; alpha = 0.001: 9.4e-4 at degree 15, 3.0e-8 at degree 31). The df polynomial does not converge: MaxRE >= 0.11
  for every alpha and degree up to 127 (alpha = 0.05: 0.28 at degree 127; alpha = 0.001: 87).
- Case ranges (all alphas): the 1/df polynomial is at 3e-15 from degree 7 for every range. The df polynomial gives 2.7e-9 to
  1.9e-5 at degree 7 (largest for Insurance) and <= 2.5e-8 at degree 15; from degree 31 both inputs are at double-precision
  level (<= 3e-13).
- Lookup tables: on the global range the error is dominated by df < 2 for every mapping (interpolation in 1/df: 0.83 at
  alpha = 0.05); restricted to df >= 30, integer-table interpolation gives 1.3e-6 (alpha = 0.05). On the case ranges,
  integer-table interpolation gives 6e-16 to 1.7e-8, round 7e-9 to 1.6e-4, and the textbook grid 9.5e-4 to 4.5e-2 (floor) or
  2e-7 to 3.4e-4 (interpolation with df = infinity). The integer table covers df up to 66,219 (largest case range, Diabetes).
- The 1/df polynomial at alpha = 0.001 is not monotone in the degree on the global range (degree 31: 3.0e-8, degree 63: 3.8e-7,
  degree 127: 1.0e-7).

## HE evaluation

`tcrit_experiment.py --part he` (about 5.5 h on the CPU container), raw results `results/he_raw.csv` (polynomials) and
`results/he_lookup_raw.csv` (lookup tables). Every method starts from the same encrypted u = 1/df at level 12 (fresh
encryption, absolute error 1.8e-7) and is evaluated under the same conditions as the plaintext part (global range and one
condition per distinct case range; the Welch df of the cases are added to the evaluated points), 10 repetitions. MaxRE is
taken against the exact c^2 over all repetitions; time excludes encryption and decryption.

- 1/df (stored): `HEHypothesisTesting._critical_value_from_inv_df`, the coefficients used by PP-TEST (u in [1/2000, 1]);
  one ciphertext, every slot one df (global: 32768 log-spaced df; case ranges: 32767 log-spaced df plus the Welch df).
- 1/df (range fit): 1/df polynomial refitted on the case range [1/df_hi, 1/df_lo].
- df (HE division): df = 1/u computed in HE (invSqrt of u * df_lo on [df_lo / df_hi, 1], squared; degree 63 with 4 Newton
  iterations on the global range, degree 31 with 1 iteration on the case ranges), then the df polynomial fitted on the
  condition range. The time includes the division stage.
- Lookup tables without division: df >= d_j is evaluated as u <= 1/d_j, i.e. one composite sign of (u - 1/d_j) / width
  for all entries at once (packed in the slots; 32768 entries per ciphertext). Floor and linear interpolation in 1/df are
  public weighted slot sums of the same step vector (floor) or of (u - 1/d_j) * step (interpolation). One df per
  ciphertext: 20 interior points of the global range, 3 interior points of each case range and the Welch df. The textbook
  table has 38 entries (with df = infinity), the integer table every integer df of the range.

### Polynomials, global range [1, 2000]

| input | degree | time (s) | levels | BTS | MaxRE alpha = 0.001 | 0.01 | 0.025 | 0.05 | 0.1 |
|---|---|---|---|---|---|---|---|---|---|
| 1/df (stored) | 7 | 0.19 +- 0.02 | 12 -> 8 | 0 | 4.8e+01 | 1.5e-01 | 1.0e-02 | 1.0e-03 | 6.7e-05 |
| 1/df (stored) | 15 | 0.33 +- 0.03 | 12 -> 7 | 0 | 9.5e-04 | 1.4e-06 | 9.6e-07 | 7.2e-07 | 5.0e-07 |
| 1/df (stored) | 31 | 0.59 +- 0.09 | 12 -> 6 | 0 | 2.8e-05 | 1.4e-06 | 1.2e-06 | 9.0e-07 | 6.3e-07 |
| 1/df (stored) | 63 | 0.92 +- 0.23 | 12 -> 5 | 0 | 2.8e-05 | 1.6e-06 | 1.0e-06 | 8.3e-07 | 5.1e-07 |
| 1/df (stored) | 127 | 1.60 +- 0.60 | 12 -> 4 | 0 | 3.1e-05 | 1.2e-06 | 1.1e-06 | 9.5e-07 | 6.0e-07 |
| df (HE division) | 7 | 8.18 +- 0.38 | 12 -> 4 | 2 | 1.0e+00 | 1.0e+00 | 9.9e-01 | 9.7e-01 | 9.2e-01 |
| df (HE division) | 15 | 8.31 +- 0.38 | 12 -> 3 | 2 | 1.0e+00 | 1.0e+00 | 9.8e-01 | 9.6e-01 | 9.0e-01 |
| df (HE division) | 31 | 10.61 +- 0.42 | 12 -> 7 | 3 | 6.2e+00 | 1.6e+00 | 9.4e-01 | 8.9e-01 | 7.8e-01 |
| df (HE division) | 63 | 11.12 +- 0.49 | 12 -> 6 | 3 | 5.3e+01 | 3.7e+00 | 1.6e+00 | 8.7e-01 | 4.5e-01 |
| df (HE division) | 127 | 12.09 +- 0.53 | 12 -> 5 | 3 | 8.7e+01 | 2.8e+00 | 7.4e-01 | 2.8e-01 | 1.1e-01 |

Division stage of the df input (df = 1/u, invSqrt of u * df_lo then squared): time, bootstraps, MaxRE of df.

| condition | time (s) | BTS | MaxRE of df |
|---|---|---|---|
| global | 8.02 +- 0.39 | 2 | 3.0e-04 |
| Insurance | 5.31 +- 0.35 | 1 | 2.1e-04 |
| Adult | 0.53 +- 0.02 | 0 | 3.5e-04 |
| Heart | 5.00 +- 0.36 | 1 | 5.9e-05 |
| Diabetes | 5.02 +- 0.49 | 1 | 1.3e-02 |
| Credit | 5.30 +- 0.32 | 1 | 4.7e-03 |
| Bank | 5.12 +- 0.42 | 1 | 7.2e-03 |
| Wine | 5.46 +- 0.26 | 1 | 1.3e-03 |

### Polynomials, case ranges, alpha = 0.05 (MaxRE over the range / at the Welch df)

| input | degree | time (s) | Insurance | Adult | Heart | Diabetes | Credit | Bank | Wine |
|---|---|---|---|---|---|---|---|---|---|
| 1/df (stored) | 7 | 0.20 | 9.8e-04 / 6.8e-04 | 1.0e-03 / 1.0e-03 | 6.6e-04 / 6.6e-04 | 1.1e-03 / 1.1e-03 | 1.1e-03 / 1.1e-03 | 1.1e-03 / 1.1e-03 | 1.1e-03 / 1.0e-03 |
| 1/df (stored) | 15 | 0.34 | 5.9e-07 / 1.1e-07 | 4.6e-07 / 1.7e-07 | 5.5e-07 / 1.6e-07 | 5.3e-07 / 2.0e-07 | 6.2e-07 / 1.7e-07 | 5.2e-07 / 2.5e-07 | 5.4e-07 / 1.9e-07 |
| 1/df (stored) | 31 | 0.57 | 5.6e-07 / 2.0e-07 | 5.4e-07 / 2.0e-07 | 5.3e-07 / 2.9e-07 | 6.0e-07 / 3.0e-07 | 4.6e-07 / 1.5e-07 | 6.3e-07 / 2.2e-07 | 6.0e-07 / 2.1e-07 |
| 1/df (range fit) | 7 | 0.18 | 7.2e-07 / 1.6e-07 | 5.1e-07 / 1.8e-07 | 5.6e-07 / 2.0e-07 | 5.3e-07 / 2.0e-07 | 5.5e-07 / 1.3e-07 | 4.6e-07 / 2.6e-07 | 4.7e-07 / 7.8e-08 |
| 1/df (range fit) | 15 | 0.25 | 5.2e-07 / 3.0e-07 | 5.7e-07 / 1.7e-07 | 5.6e-07 / 2.0e-07 | 6.0e-07 / 2.5e-07 | 5.1e-07 / 1.9e-07 | 5.2e-07 / 2.0e-07 | 5.3e-07 / 1.4e-07 |
| 1/df (range fit) | 31 | 0.38 | 5.5e-07 / 1.7e-07 | 5.4e-07 / 2.5e-07 | 5.6e-07 / 1.5e-07 | 5.1e-07 / 1.7e-07 | 5.4e-07 / 1.7e-07 | 5.9e-07 / 1.9e-07 | 5.4e-07 / 1.1e-07 |
| df (HE division) | 7 | 4.99 | 7.5e-06 / 6.6e-06 | 5.2e-07 / 2.0e-07 | 5.4e-07 / 1.8e-07 | 6.8e-07 / 2.4e-07 | 5.7e-07 / 3.1e-07 | 3.2e-06 / 1.9e-06 | 8.2e-07 / 3.6e-07 |
| df (HE division) | 15 | 5.13 | 4.8e-07 / 2.2e-07 | 5.2e-07 / 2.0e-07 | 5.1e-07 / 1.8e-07 | 5.6e-07 / 1.5e-07 | 4.7e-07 / 3.3e-07 | 4.9e-07 / 1.3e-07 | 5.0e-07 / 1.6e-07 |
| df (HE division) | 31 | 5.26 | 4.8e-07 / 2.2e-07 | 5.2e-07 / 2.0e-07 | 5.1e-07 / 1.8e-07 | 5.6e-07 / 1.5e-07 | 4.7e-07 / 3.3e-07 | 4.9e-07 / 1.3e-07 | 5.0e-07 / 1.6e-07 |

### Polynomials, case ranges, alpha = 0.001 (MaxRE over the range / at the Welch df)

| input | degree | time (s) | Insurance | Adult | Heart | Diabetes | Credit | Bank | Wine |
|---|---|---|---|---|---|---|---|---|---|
| 1/df (stored) | 7 | 0.20 | 4.6e+01 / 3.2e+01 | 4.8e+01 / 4.8e+01 | 3.1e+01 / 3.1e+01 | 5.1e+01 / 5.1e+01 | 5.1e+01 / 5.0e+01 | 5.1e+01 / 5.0e+01 | 5.0e+01 / 4.9e+01 |
| 1/df (stored) | 15 | 0.34 | 8.4e-04 / 9.0e-05 | 9.6e-04 / 9.5e-04 | 8.0e-04 / 1.3e-04 | 1.2e-03 / 1.2e-03 | 1.2e-03 / 1.1e-03 | 1.2e-03 / 1.1e-03 | 1.1e-03 / 1.0e-03 |
| 1/df (stored) | 31 | 0.57 | 3.7e-05 / 2.7e-06 | 3.2e-05 / 5.9e-06 | 3.9e-05 / 4.4e-06 | 3.2e-05 / 5.9e-06 | 3.3e-05 / 5.0e-06 | 3.4e-05 / 3.9e-06 | 4.2e-05 / 7.7e-06 |
| 1/df (range fit) | 7 | 0.18 | 1.2e-06 / 4.2e-07 | 1.2e-06 / 5.4e-07 | 1.6e-06 / 4.5e-07 | 1.3e-06 / 3.2e-07 | 1.6e-06 / 3.9e-07 | 1.4e-06 / 6.2e-07 | 1.3e-06 / 3.8e-07 |
| 1/df (range fit) | 15 | 0.25 | 1.3e-06 / 2.6e-07 | 1.5e-06 / 4.5e-07 | 1.2e-06 / 5.0e-07 | 1.2e-06 / 4.3e-07 | 1.2e-06 / 6.6e-07 | 1.4e-06 / 4.3e-07 | 1.3e-06 / 2.7e-07 |
| 1/df (range fit) | 31 | 0.38 | 1.1e-06 / 4.8e-07 | 1.6e-06 / 5.2e-07 | 1.2e-06 / 3.7e-07 | 1.2e-06 / 3.9e-07 | 1.4e-06 / 3.9e-07 | 1.4e-06 / 4.2e-07 | 1.2e-06 / 4.3e-07 |
| df (HE division) | 7 | 4.99 | 1.9e-05 / 1.7e-05 | 1.3e-06 / 4.8e-07 | 1.3e-06 / 4.6e-07 | 1.7e-06 / 5.9e-07 | 1.4e-06 / 7.7e-07 | 7.8e-06 / 4.6e-06 | 2.0e-06 / 8.9e-07 |
| df (HE division) | 15 | 5.13 | 1.2e-06 / 5.4e-07 | 1.3e-06 / 4.8e-07 | 1.3e-06 / 4.3e-07 | 1.4e-06 / 3.6e-07 | 1.2e-06 / 8.1e-07 | 1.2e-06 / 3.2e-07 | 1.2e-06 / 4.0e-07 |
| df (HE division) | 31 | 5.26 | 1.2e-06 / 5.4e-07 | 1.3e-06 / 4.8e-07 | 1.3e-06 / 4.3e-07 | 1.4e-06 / 3.6e-07 | 1.2e-06 / 8.1e-07 | 1.2e-06 / 3.2e-07 | 1.2e-06 / 4.0e-07 |

### Lookup tables in HE (one df per ciphertext; MaxRE over the points, alpha = 0.05 / 0.001)

| condition | points | table | entries | ciphertexts | compare (s) | BTS | floor MaxRE | interp MaxRE | max rel. diff HE vs plain |
|---|---|---|---|---|---|---|---|---|---|
| global | 20 | textbook | 38 | 1 | 10.26 +- 0.63 | 4 | 2.9e+00 / 4.0e+01 | 8.2e-01 / 1.5e+01 | 3.3e-05 |
| global | 20 | integer | 2000 | 1 | 10.20 +- 0.61 | 4 | 2.9e+00 / 4.0e+01 | 8.2e-01 / 1.5e+01 | 9.3e-05 |
| Insurance | 4 | textbook | 38 | 1 | 9.83 +- 0.77 | 4 | 1.8e-02 / 4.4e-02 | 5.0e-05 / 2.9e-04 | 3.1e-05 |
| Insurance | 4 | integer | 1064 | 1 | 10.16 +- 0.74 | 4 | 2.1e-05 / 5.2e-05 | 2.2e-07 / 5.5e-07 | 5.5e-07 |
| Adult | 5 | textbook | 38 | 1 | 10.13 +- 0.61 | 4 | 1.2e-03 / 3.0e-03 | 1.2e-06 / 3.8e-05 | 3.3e-05 |
| Adult | 5 | integer | 1024 | 1 | 10.44 +- 0.46 | 4 | 2.0e-06 / 4.9e-06 | 2.2e-07 / 4.7e-07 | 3.9e-06 |
| Heart | 4 | textbook | 38 | 1 | 10.25 +- 0.63 | 4 | 1.2e-02 / 3.0e-02 | 5.9e-05 / 3.5e-04 | 2.5e-05 |
| Heart | 4 | integer | 164 | 1 | 10.29 +- 0.70 | 4 | 6.2e-05 / 1.5e-04 | 1.6e-07 / 3.5e-07 | 3.4e-07 |
| Diabetes | 4 | textbook | 38 | 1 | 10.30 +- 0.78 | 4 | 2.4e-03 / 5.8e-03 | 4.3e-07 / 2.7e-05 | 2.8e-05 |
| Diabetes | 4 | integer | 54864 | 2 | 19.46 +- 1.16 | 4 | 8.8e-07 / 2.1e-06 | 2.0e-07 / 4.9e-07 | 2.1e-06 |
| Credit | 4 | textbook | 38 | 1 | 10.35 +- 0.65 | 4 | 2.3e-03 / 5.7e-03 | 5.0e-07 / 2.7e-05 | 2.5e-05 |
| Credit | 4 | integer | 23364 | 1 | 10.56 +- 0.79 | 4 | 9.0e-07 / 2.2e-06 | 1.5e-07 / 3.8e-07 | 2.2e-06 |
| Bank | 4 | textbook | 38 | 1 | 10.07 +- 1.39 | 4 | 2.3e-03 / 5.7e-03 | 6.7e-07 / 3.0e-05 | 3.0e-05 |
| Bank | 4 | integer | 39922 | 2 | 19.93 +- 1.08 | 4 | 1.1e-06 / 2.8e-06 | 2.1e-07 / 5.0e-07 | 2.8e-06 |
| Wine | 4 | textbook | 38 | 1 | 10.49 +- 0.65 | 4 | 1.9e-03 / 4.7e-03 | 1.3e-06 / 3.6e-05 | 3.1e-05 |
| Wine | 4 | integer | 4898 | 1 | 10.47 +- 0.57 | 4 | 9.8e-07 / 2.3e-06 | 2.3e-07 / 5.1e-07 | 3.4e-06 |

Weighted slot sum per output (one alpha, one mapping): 0.16 s on average.

### Observations

- The stored 1/df polynomial (degree 15) takes 0.33 s with no bootstrap and stays at the CKKS noise level (MaxRE <= 1.4e-6
  for alpha >= 0.01) on the global range and on every case range, including df > 2000 (Adult, Diabetes, Credit, Bank, Wine:
  <= 6.2e-7 at alpha = 0.05). At alpha = 0.001 its error is the polynomial error (9.5e-4 to 1.2e-3; degree 31: 2.8e-5 to
  4.2e-5).
- The df input needs a division in HE: 5 to 8 s and 1 to 2 bootstraps, about 25 times the cost of the 1/df polynomial. The
  computed df carries a relative error up to 1.3e-2 for large df (the absolute input noise 1.8e-7 is large relative to
  u = 1/df), which the flat c^2 at large df absorbs. On the global range the df polynomial fails as in plaintext (MaxRE
  0.11 to 87 up to degree 127). On the case ranges it reaches the noise level from degree 15 (<= 5.6e-7 at alpha = 0.05,
  <= 1.4e-6 at alpha = 0.001); at degree 7 it is less accurate (Insurance 7.5e-6, Bank 3.2e-6).
- On the same case range, the 1/df polynomial reaches the noise level from degree 7 (<= 7.2e-7 at alpha = 0.05, <= 1.6e-6
  at alpha = 0.001) in 0.18 s, i.e. the same accuracy as the df input without the division.
- Lookup tables cost one composite sign per table: 10 s and 4 bootstraps (about 20 s for the integer tables of Diabetes
  and Bank, which need 2 ciphertexts); the weighted sum adds 0.16 s per output. The HE values agree with the plaintext
  lookup (relative difference <= 9.3e-5), so the accuracy is that of the table: on the global range floor and
  interpolation fail at df < 2 (MaxRE 2.9 and 0.82 at alpha = 0.05; df >= 30: integer interpolation 6.6e-7); on the case
  ranges integer-table interpolation reaches the noise level (<= 2.3e-7 at alpha = 0.05), integer floor <= 6.2e-5,
  textbook interpolation <= 5.9e-5 and textbook floor 4.3e-4 to 1.8e-2.
