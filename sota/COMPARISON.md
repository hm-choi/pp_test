# SoTA comparison scope

This document separates an end-to-end protocol comparison from a same-engine
arithmetic ablation. They answer different questions and must not be merged
into a single runtime ranking.

## Protocol-level comparison

| Aspect | This work | Annamalai, Jin, and Aung (2022) | Zhou et al. (2026) |
|---|---|---|---|
| Cryptographic setting | Single-key CKKS prototype in this repository | MHE with an MHE-to-SMPC conversion protocol | Hybrid MHE + MPC platform |
| Welch input layout | Two encrypted groups | Encrypted normalized values and encrypted binary labels; vertical partitioning | Genomic values and phenotype/group information at separate sites |
| Public parameters | Group sizes, \(\alpha\), input/score bounds; F critical values where applicable | Group sizes, normalized \([0,1]\) domain, denominator upper bounds | Total sample/group counts; protocol-specific public configuration |
| Nonlinear arithmetic for Welch | HE inverse square root plus Newton refinement | MHE Goldschmidt division | MPC division after HE-to-MPC conversion |
| Welch output before inference | Encrypted \(T\) and encrypted \(1/\nu\) | Encrypted \(T^2\) and \(\nu\) | Secure t-statistic and df in the hybrid protocol |
| p-value mechanism | Not evaluated; avoids Student-t CDF by comparing against encrypted critical value | Private table lookup after conversion to secret sharing | Private lookup table for t-statistic/df |
| Final hypothesis decision | Encrypted score and encrypted step for Welch, Z, and F | Hybrid p-value protocol; not the same HE-only step interface | Hybrid MHE/MPC protocol; lookup-derived inference |
| Interaction during nonlinear inference | No MPC conversion in the current prototype | Collective MHE operations plus HE-to-SMPC conversion/lookup | MHE/MPC conversion and MPC nonlinear stages |
| Direct end-to-end runtime comparison with this repository | Not applicable | Not valid: different key model, interaction, batching, and MPC stages | Not valid: different platform, workloads, and MHE/MPC stages |

The proposed implementation keeps the decision ciphertext-encrypted, but its
single-key prototype is **not** a replacement for a multi-party key-management
or distributed-decryption protocol. Conversely, the cited MHE/MPC systems
address federation and communication, but their interactive stages must be
included in any end-to-end performance claim.

## Same-engine arithmetic ablation

`goldschmidt_welch_he_baseline.py` is retained only for this narrower test:
under the same HEaaN engine and ciphertext representation, compare the
inverse-square-root realization against a residual Goldschmidt quotient.

| Item | Proposed arithmetic | Goldschmidt baseline |
|---|---|---|
| Division realization | \(1/\sqrt{x}\) Chebyshev seed + Newton refinement | Residual Goldschmidt recurrence |
| Welch representation | Signed \(T\) and \(1/\nu\) | \(T^2\) and \(\nu\) |
| Bounds | Public upper range bound for normalization | Requires public positive lower and upper denominator bounds |
| Encrypted critical value / decision | Yes, in the proposed API | No |
| Permitted claim | End-to-end result for this repository | HE arithmetic diagnostic only |

The oracle-bound insurance script (`test_sota_comparison.py`) is deliberately
labelled a diagnostic. Its bounds are selected from the plaintext instance and
must never be cited as a privacy-valid comparison. The paper-layout script
(`test_annamalai2022_comparison.py`) uses the 2022 paper's encrypted
value/label layout and public upper bounds, but remains a single-key arithmetic
trace: it does not reproduce collective MHE bootstrapping, MHE-to-SMPC
conversion, private lookup, or packed batch division.

On the insurance diagnostic with 12 direct-HDiv rounds, the trace produced a
large numerical error (including an invalid negative df), whereas the proposed
native two-group circuit retained its expected accuracy under its public
`charges / 1000` scaling. This is **not** evidence of a speed or accuracy win
over the paper: it only establishes that the current single-key trace is not a
numerically faithful replacement for the paper's complete MHE realization.
Do not include the trace's runtime or error as a SoTA baseline result.

## Reporting rule

Report the following as separate tables or subsections:

1. The real-data and boundary accuracy/runtime results for this repository.
2. The protocol-level table above, including interaction and output scope.
3. The optional HE-only inverse-square-root versus Goldschmidt ablation, with
   its public-bound assumptions printed next to every runtime.

Do not claim an end-to-end speedup over either MHE/MPC paper from the
single-key arithmetic traces.

## References

1. M. S. M. S. Annamalai, C. Jin, and K. M. M. Aung, “Communication-Efficient
   Secure Federated Statistical Tests from Multiparty Homomorphic Encryption,”
   *Applied Sciences*, 12(22):11462, 2022.
   https://doi.org/10.3390/app122211462
2. W. Zhou et al., “Secure bioinformatics: privacy-preserving federated
   analytics using homomorphic encryption,” *Bioinformatics*, 42(5):btag081,
   2026. https://doi.org/10.1093/bioinformatics/btag081
