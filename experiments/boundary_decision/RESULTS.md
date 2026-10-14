# Boundary-sensitive decisions

Script `run_boundary_decision.py`; raw results `results/boundary_decision_raw.csv`.

HEaaN FGb (GPU, NVIDIA Tesla V100), log_slots = 15. invSqrt domain V in [0.01, 1000] (degree-63 Chebyshev + 7 Newton iterations); ApproxTCrit degree 15. Synthetic groups of n = ceil(df / 2) + 1 samples whose Welch t statistic is t_crit + margin,
margins +-0.001, +-0.002, +-0.005, +-0.01, 10 repetitions (80 trials per cell). The variance ratio of the groups
sets the Welch-Satterthwaite df exactly; both groups are scaled so that V >= 0.02.

| alpha | df 3 | df 4 | df 5 | df 6 | df 7 | df 8 | df 9 | df 10 | df 50 | df 100 | df 2000 | df 10000 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.01 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 |
| 0.025 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 |
| 0.05 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 |
| 0.1 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 | 80/80 |

Total 3840/3840; in-domain inputs 3840/3840, errors 0.
