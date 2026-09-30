# Experiment 2 — Encrypted Z-test and F-test

Configuration: HEaaN FGb on CPU, `log_slots=15`, two-sided \(\alpha=0.05\).
The Adult dataset uses 1,024 Female and 1,024 Male age records.

| Test | Statistic | Plaintext statistic | Encrypted statistic | Absolute error | Plaintext two-sided p-value | Decision match | Runtime (s) | Bootstraps |
|---|---|---:|---:|---:|---:|---|---:|---:|
| Z | Age mean, Female vs Male | -3.381726 | -3.381726 | \(8.30\times10^{-9}\) | 0.000720 | Yes (encrypted step) | 37.12 | 4 |
| F | Age variance, Female vs Male | 1.159098 | 1.159098 | \(1.54\times10^{-8}\) | 0.018301 | Yes (encrypted step) | 62.88 | 8 |

For the Z-test, the input variances \(207.398670\) and \(178.931092\) are
explicitly treated as public known variances. The F critical values are public
constants derived from the public group sizes and \(\alpha\); the final
two-sided decision is encrypted.

The Z-test additionally evaluates the encrypted score
\(Z^2-z_{1-\alpha/2}^2\). Its observed score error is
\(5.61\times10^{-8}\), and the encrypted step output is 1.0. The F-test
evaluates \((F-f_L)(F-f_U)\); its observed score error is
\(8.28\times10^{-7}\), and its encrypted step output is also 1.0.
