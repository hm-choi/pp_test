# Experiment 1 — Encrypted Welch's t-test

Configuration: HEaaN FGb on CPU, `log_slots=15`, two-sided \(\alpha=0.05\),
and degree-15 Chebyshev approximation for the encrypted t critical value.

## Real-data validation

| Dataset and comparison | Group sizes | \(|T-\hat T|\) | \(|\mathrm{InvDF}-\widehat{\mathrm{InvDF}}|\) | \(|s-\hat s|\) | Step match | Runtime (s) | Bootstraps |
|---|---:|---:|---:|---:|---|---:|---:|
| Insurance charges: smoker yes vs no | 274, 1064 | \(4.11\times10^{-8}\) | \(2.18\times10^{-6}\) | \(8.97\times10^{-3}\) | Yes | 87.87 | 11 |
| Adult educational-num: >50K vs <=50K | 1024, 1024 | \(2.44\times10^{-7}\) | \(3.54\times10^{-7}\) | \(3.51\times10^{-3}\) | Yes | 80.90 | 11 |
| Adult age: Female vs Male | 1024, 1024 | \(2.37\times10^{-8}\) | \(6.15\times10^{-7}\) | \(1.25\times10^{-4}\) | Yes | 76.69 | 11 |

Here \(s=T^2-c_\alpha(\mathrm{df})^2\). The encrypted step is evaluated on
a public normalization bound for \(s\).

## Boundary and degree stress test

Synthetic groups have 128 samples each and Welch df \(=254\). The stress
range uses public score bound \(0.05\), with target t margins from
\(-0.01\) to \(+0.01\) around the two-sided critical value.

| Critical-value degree | Score-sign match rate | Degree-15 integrated step match | Max critical-value error | Max score error | Mean critical evaluation (s) | Minimum score level |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 0.875 | N/A | \(1.191\times10^{-3}\) | \(4.713\times10^{-3}\) | 0.72 | 6 |
| 15 | 1.000 | 1.000 | \(2.0\times10^{-6}\) | \(1.1\times10^{-5}\) | 1.43 | 5 |
| 31 | 1.000 | N/A | \(2.0\times10^{-6}\) | \(1.1\times10^{-5}\) | 2.25 | 4 |
| 63 | 1.000 | N/A | \(2.0\times10^{-6}\) | \(1.1\times10^{-5}\) | 4.41 | 3 |
| 127 | 1.000 | N/A | \(2.0\times10^{-6}\) | \(1.1\times10^{-5}\) | 8.58 | 2 |

Degree 7 fails at the \(+0.001\) target margin. Degree 15 is used as the
default because higher degrees did not improve observed accuracy, while they
increased evaluation time and consumed more levels.
