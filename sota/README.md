# Welch baseline naming and comparison scope

For the paper-ready protocol table and the rules distinguishing it from the
same-engine arithmetic ablation, see [COMPARISON.md](COMPARISON.md).

`goldschmidt_welch_he_baseline.py` is the single comparison baseline. It uses
the same encrypted sufficient statistics and reports

\[
T^2=\frac{(\bar x_1-\bar x_2)^2}{V},\qquad
\nu=\frac{V^2}{D},
\]

with normalized Goldschmidt division. It is appropriate for a same-engine HE
arithmetic comparison with the proposed inverse-square-root Welch method.

The baseline uses the residual recurrence
\(q_{i+1}=q_i(1+r_i),\ r_{i+1}=r_i^2\).  Its two ciphertext multiplications
are independent, so each iteration has multiplicative depth one.  The code
keeps only the residual state and does not bootstrap the redundant
``w_i = 1 + r_i`` ciphertext separately.

It is **not** a full MHE/MPC-paper implementation: it has no distributed
decryption, secret sharing, MHE-to-MPC conversion, private t-table lookup,
encrypted comparison, or final encrypted decision. These must be reported as
separate end-to-end protocol stages rather than represented as HE operations.

Relevant protocol-level papers are Annamalai, Jin, and Aung, *Communication-
Efficient Secure Federated Statistical Tests from Multiparty Homomorphic
Encryption*, Applied Sciences 12(22), 2022, doi:10.3390/app122211462; and Zhou
et al., *Secure bioinformatics: privacy-preserving federated analytics using
homomorphic encryption*, Bioinformatics, 2026,
doi:10.1093/bioinformatics/btag081.

`PublicWelchBounds` is mandatory. Its interval and output bounds are public
protocol parameters agreed before encryption. Inferring them from plaintext
reference statistics is permitted only for post-run evaluation, never as input
to the encrypted algorithm.

The current insurance ``test_sota_comparison.py`` uses oracle-selected bounds
only as a clearly labelled best-case diagnostic.  Its timing must not be used
as a privacy-valid or paper-faithful comparison.  The 2022 MHE paper bounds its
normalized Welch denominators from above, but Goldschmidt initialization also
requires a strictly positive lower bound; the paper text does not provide the
insurance-specific lower bounds used by this diagnostic.

`annamalai2022_welch_mhe_arithmetic.py` and
`experiments/experiment1/test_annamalai2022_comparison.py` instead follow the
2022 paper's published arithmetic setup: encrypted values and encrypted binary
labels, values normalized to ``[0, 1]``, public group sizes, direct HDiv with
the paper's upper bounds, and lower bound zero.  They are still an arithmetic
trace rather than a full replication because this repository does not provide
MHE collective operations, a packed HDiv batch, or the paper's MHE-to-SMPC
p-value protocol.
