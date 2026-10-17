# Experiment 2: Critical-Boundary Approximation

Section 5.5 and Table 3 of the paper.

## Task

Compare Chebyshev approximations of the squared two-sided t critical value t_{1-alpha/2, df}^2 built over the degrees of freedom df (Direct-DF) and over the inverse degrees of freedom u = 1/df (Inverse-DF, ApproxTCrit).

## Setting

- Reference: SciPy `t.ppf`. Error: relative error, reported as MaxRE and MRE.
- 65,536 df values log-spaced over [3, 10^4], packed into two ciphertexts.
- Degrees 7, 15, 31, 63 and 127.
- Inverse-DF uses the stored coefficients `coeffs/t_critical_coeffs.py`, fitted on u in [0, 1]. Direct-DF fits a polynomial on df in [3, 10^4] and computes df = 1/u with the inverse square root.
- Encrypted evaluation: alpha = 0.05, ten repetitions. The plaintext part covers all four significance levels.

## Results (V100, alpha = 0.05, ten repetitions)

Mean (standard deviation). PP-TEST uses degree 15.

| Degree | Direct-DF MaxRE | Direct-DF MRE | Inverse-DF MaxRE | Inverse-DF MRE |
| --- | --- | --- | --- | --- |
| 7 | 6.09e-01 (1.45e-09) | 7.68e-02 (1.51e-09) | 1.00e-03 (3.35e-08) | 6.61e-04 (1.39e-10) |
| **15** | 5.76e-01 (2.25e-08) | 6.16e-02 (2.60e-09) | **5.33e-07** (1.18e-07) | **4.34e-08** (1.23e-10) |
| 31 | 4.67e-01 (3.05e-07) | 4.09e-02 (5.03e-09) | 5.30e-07 (7.39e-08) | 4.33e-08 (1.43e-10) |
| 63 | 2.37e-01 (2.79e-06) | 2.11e-02 (4.77e-09) | 5.42e-07 (6.27e-08) | 4.35e-08 (1.41e-10) |
| 127 | 4.58e-02 (1.21e-05) | 4.46e-03 (1.27e-08) | 5.50e-07 (5.87e-08) | 4.33e-08 (1.66e-10) |

The plaintext MaxRE for every significance level is in `results/SUMMARY.md`.

## Run

```
python3 experiments/experiment2/run_critical_boundary.py --part plain
python3 experiments/experiment2/run_critical_boundary.py --part he --reps 10 --alphas 0.05
```

## Output

| File | Contents |
| --- | --- |
| `plaintext.csv` | Plaintext MaxRE and MRE per method, degree and significance level, over [3, 10^4] and over the df range of each workload, with MaxRE per df band |
| `he_raw.csv` | Encrypted evaluation per condition, method, degree and repetition: time, bootstrapping count, levels, MaxRE and MRE |

The results of the paper are in `results/experiment2/`. Table 3 uses the rows with `condition = global`.
