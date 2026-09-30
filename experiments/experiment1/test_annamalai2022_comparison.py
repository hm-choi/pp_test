"""Compare the proposed arithmetic with the 2022-paper Welch setup.

This follows the published value/class input layout and denominator upper
bounds.  It is not a wall-clock reproduction of the paper's MHE benchmark:
this repository uses single-key HEaaN and has no collective bootstrap,
collective decryption, MHE-to-SMPC conversion, or packed batch-division API.
"""

from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT_ROOT))

from pp_test.engine.HEengine import HEEngine
from operators.hypothesis_testing import HEHypothesisTesting
from operators.operator import HEOperator
from sota.annamalai2022_welch_mhe_arithmetic import (
    Annamalai2022WelchArithmetic,
)


# A public domain cap, chosen independently of the observed maximum charge.
# Dividing by it maps each input value into the paper's [0, 1] domain.
PUBLIC_CHARGE_UPPER_DOLLARS = 65_000.0
# The proposed sufficient-statistics circuit needs a numerically well-scaled
# public representation before the Q - S^2/n cancellation.  This scale is
# public and leaves Welch T and df invariant.
PROPOSED_CHARGE_SCALE = 1_000.0
PROPOSED_RANGE_BOUND = 65.0
# The paper specifies Goldschmidt HDiv but does not state a single universal T.
# Twelve rounds are reported as an explicit experiment parameter, not attributed
# to the paper; this is needed for the paper's loose [0, b_u] D bound here.
PAPER_HDIV_ITERATIONS = 12


def plaintext_welch(values: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    one = values[labels == 1.0]
    zero = values[labels == 0.0]
    var_one = float(np.var(one, ddof=1))
    var_zero = float(np.var(zero, ddof=1))
    v = var_one / len(one) + var_zero / len(zero)
    d = ((var_one / len(one)) ** 2) / (len(one) - 1) + (
        (var_zero / len(zero)) ** 2
    ) / (len(zero) - 1)
    t = float((np.mean(one) - np.mean(zero)) / np.sqrt(v))
    return {"T": t, "T2": t * t, "DF": v * v / d, "V": v, "D": d}


def main() -> None:
    insurance = pd.read_csv(PROJECT_ROOT / "datasets" / "insurance.csv")
    raw_charges = insurance["charges"].to_numpy(dtype=np.float64)
    if float(np.max(raw_charges)) > PUBLIC_CHARGE_UPPER_DOLLARS:
        raise ValueError("The public charge cap does not bound the input data.")

    # This is the paper's vertically partitioned representation: one owner has
    # values, another has class labels.  Both vectors are encrypted below.
    values = raw_charges / PUBLIC_CHARGE_UPPER_DOLLARS
    labels = (insurance["smoker"] == "yes").to_numpy(dtype=np.float64)
    n_one = int(np.sum(labels))
    n_zero = int(len(labels) - n_one)
    reference = plaintext_welch(values, labels)

    # Warm-up is outside both measured regions, avoiding first-bootstrap bias.
    engine = HEEngine(warmup_bootstrap=True)
    operator = HEOperator(engine)
    proposed = HEHypothesisTesting(engine)
    paper_trace = Annamalai2022WelchArithmetic(engine, verbose=True)
    encrypted_values = operator.encrypt(values)
    encrypted_labels = operator.encrypt(labels)

    # The proposed API is natively a two-group interface.  It uses the same
    # raw observations but its public /1000 representation avoids catastrophic
    # cancellation in Q-S^2/n that occurs for this CKKS implementation when
    # all values are compressed into [0, 1].  Welch T and df are invariant to
    # this common positive rescaling.
    proposed_values = raw_charges / PROPOSED_CHARGE_SCALE
    encrypted_one = operator.encrypt(proposed_values[labels == 1.0])
    encrypted_zero = operator.encrypt(proposed_values[labels == 0.0])

    proposed.reset_bootstrap_count()
    started_at = perf_counter()
    proposed_result = proposed.HE_Welch_T_Test(
        encrypted_one,
        encrypted_zero,
        R=PROPOSED_RANGE_BOUND,
        alpha=0.05,
        critical_degree=15,
        score_bound=1_200.0,
        debug=True,
    )
    proposed_seconds = perf_counter() - started_at
    proposed_bootstraps = proposed.bootstrap_count()
    proposed_t = float(operator.decrypt(proposed_result["T"])[0])
    proposed_inv_df = float(operator.decrypt(proposed_result["INV_DF"])[0])

    paper_trace.reset_bootstrap_count()
    started_at = perf_counter()
    paper_result = paper_trace.welch_t2_df(
        encrypted_values,
        encrypted_labels,
        n_one,
        n_zero,
        PAPER_HDIV_ITERATIONS,
        return_debug=True,
    )
    paper_seconds = perf_counter() - started_at
    paper_bootstraps = paper_trace.bootstrap_count()
    paper_t2 = float(operator.decrypt(paper_result["T2"])[0])
    paper_df = float(operator.decrypt(paper_result["DF"])[0])

    print("========== Annamalai et al. (2022) Welch arithmetic setup ==========")
    print("input layout            : encrypted normalized values + encrypted binary labels")
    print("input domain            : values in [0, 1]")
    print(f"public charge cap       : {PUBLIC_CHARGE_UPPER_DOLLARS}")
    print(f"public group sizes      : label=1: {n_one}, label=0: {n_zero}")
    print(f"paper V upper           : {paper_result['V_upper']}")
    print(f"paper D upper           : {paper_result['D_upper']}")
    print(f"paper HDiv iterations   : {PAPER_HDIV_ITERATIONS} (experiment parameter)")
    print("paper lower bound       : 0, as stated; valid inputs require V,D > 0")
    print("\nplaintext T2, df       : " f"{reference['T2']}, {reference['DF']}")
    print("\n[Proposed InvSqrt arithmetic]")
    print("input layout            : native encrypted group vectors")
    print("underlying observations  : same raw charges as paper trace")
    print(f"public input scale       : charges / {PROPOSED_CHARGE_SCALE}")
    print(f"public range bound       : {PROPOSED_RANGE_BOUND}")
    print(f"encrypted T             : {proposed_t}")
    print(f"reported T2             : {proposed_t * proposed_t}")
    print(f"encrypted InvDF         : {proposed_inv_df}")
    print(f"reported df (=1/InvDF)  : {1.0 / proposed_inv_df}")
    print(f"T absolute error        : {abs(proposed_t - reference['T'])}")
    print(f"InvDF absolute error    : {abs(proposed_inv_df - 1.0 / reference['DF'])}")
    print(f"runtime (seconds)       : {proposed_seconds}")
    print(f"bootstraps              : {proposed_bootstraps}")
    print("\n[Annamalai et al. 2022 direct Goldschmidt arithmetic trace]")
    print(f"encrypted T2            : {paper_t2}")
    print(f"encrypted df            : {paper_df}")
    print(f"T2 absolute error       : {abs(paper_t2 - reference['T2'])}")
    print(f"df absolute error       : {abs(paper_df - reference['DF'])}")
    print(f"runtime (seconds)       : {paper_seconds}")
    print(f"bootstraps              : {paper_bootstraps}")
    print("\n========== Scope required for interpreting runtime ==========")
    print("This run uses single-key HEaaN, not the paper's MHE collective bootstrap.")
    print("The two HDiv calls are not packed into one batch ciphertext here.")
    print("The paper's MHE-to-SMPC p-value lookup is not part of this arithmetic trace.")
    print("A large direct-HDiv error indicates that this trace is not an accuracy")
    print("baseline outside the paper's full MHE range/bootstrapping realization.")


if __name__ == "__main__":
    main()
