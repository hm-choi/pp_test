# Boundary-sensitive decisions across degrees of freedom

Script `run_boundary_decision.py`; raw results `results/boundary_decision_raw.csv`.

Synthetic groups whose Welch statistic is (t_crit + margin) SE next to the two-sided critical value, margins
+-0.001, +-0.002, +-0.005, +-0.01, 10 repetitions; degree-15 critical-value polynomial and step function.
df >= 2: two groups of n = df/2 + 1 samples, 10 + 5 linspace(-1, 1, n) +- gap/2 (equal variances, df = 2n - 2);
df = 1: n1 = n2 = 2 with a constant second group (df = n1 - 1 = 1). [L, U] = minimum and maximum of the two groups.

Each cell: matching decisions / trials (8 margins x repetitions). * = inputs outside the invSqrt domain
(V / V_max < V_BAR_MIN = 0.01).

| alpha | df 1 | df 2 | df 10 | df 50 | df 100 | df 2000 | df 10000 |
|---|---|---|---|---|---|---|---|
| 0.001 | 42/80* | 40/80* | 80/80 | 80/80 | 70/80 | 70/80 | 70/80 |
| 0.01 | 40/80* | 78/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 |
| 0.025 | 40/80* | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 79/80 |
| 0.05 | 78/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 79/80 |
| 0.1 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 |

Mismatches inside the domain (df, alpha, margin: count): (1, 0.05, 0.001: 2), (2, 0.01, -0.001: 1), (2, 0.01, 0.001: 1), (100, 0.001, 0.001: 10), (2000, 0.001, -0.001: 10), (10000, 0.001, -0.001: 10), (10000, 0.025, 0.002: 1), (10000, 0.05, -0.001: 1).

Smallest plaintext |s| / B per df (over alphas and margins): df 1: 3.14e-06, df 2: 1.56e-07, df 10: 3.62e-06, df 50: 6.70e-07, df 100: 3.32e-07, df 2000: 1.65e-08, df 10000: 3.29e-09.
