"""Boundary-sensitive decisions of Welch's t-test across degrees of freedom.

Synthetic groups whose Welch t statistic lies next to the two-sided critical
value t_crit, (t_crit + margin) * SE, for every df in DFS, alpha in ALPHAS
and margin in MARGINS: two groups of n = ceil(df / 2) + 1 samples,
10 + k * 5 * linspace(-1, 1, n) and 10 + k * r * 5 * linspace(-1, 1, n),
shifted by +/- mean_gap / 2. The variance ratio r^2 sets the
Welch-Satterthwaite df exactly: with w_i = A_i / V,
1 / df = (w1^2 + w2^2) / (n - 1), so r = 1 (equal variances) for an even df
and r < 1 for an odd df. The spread factor k >= 1 raises
V = s1^2 / n1 + s2^2 / n2 to at least V_TARGET inside the invSqrt domain;
T and df do not depend on k.
The public interval [L, U] is the minimum and maximum of the two groups.
The encrypted decision uses the degree-15 critical-value polynomial and the
step function. Inputs whose V lies outside the invSqrt domain are still
evaluated and marked with in_domain = 0.

Run (inside the HEaaN container, cwd = project root):
    python3 experiments/experiment3/run_boundary_decision.py
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(PROJECT_ROOT))
sys.path.append(str(PROJECT_ROOT / "experiments"))

from common.he_utils import completed, first, make_engine, open_append, timed
from coeffs.invSqrt_coeffs import _INV_SQRT_MIN, _INV_SQRT_MAX


ALPHAS = [0.01, 0.025, 0.05, 0.1]
DFS = (3, 4, 5, 6, 7, 8, 9, 10, 50, 100, 2000, 10000)
MARGINS = (-0.01, -0.005, -0.002, -0.001, 0.001, 0.002, 0.005, 0.01)
CRITICAL_LOG_DEGREE = 4  # degree 15
REPS = 2

# Smallest V of the synthetic groups (twice the lower end of the invSqrt
# domain). The public score bound grows with R^2, so the spread is enlarged
# only as far as needed.
V_TARGET = 2.0 * _INV_SQRT_MIN

RESULT_DIR = Path(__file__).resolve().parent / "results"

RAW_FIELDS = [
    "df", "n1", "n2", "v", "in_domain", "error", "alpha", "margin", "rep",
    "stats_time", "decision_time", "stats_bootstraps", "decision_bootstraps",
    "plain_t2", "he_t2", "t2_abs_err", "plain_inv_df", "he_inv_df", "inv_df_abs_err",
    "plain_critical", "he_critical", "critical_abs_err",
    "plain_score", "he_score", "score_abs_err", "score_bound",
    "he_step", "plain_decision", "score_decision", "step_decision",
]


def groups(alpha, margin, df):

    critical_target = float(stats.t.ppf(1.0 - alpha / 2.0, df))

    n = math.ceil(df / 2) + 1
    noise = 5.0 * np.linspace(-1.0, 1.0, n)

    # w1^2 + w2^2 = (n - 1) / df with w1 + w2 = 1 and w2 / w1 = r^2.
    w1 = (1.0 + math.sqrt(2.0 * (n - 1) / df - 1.0)) / 2.0
    r = math.sqrt((1.0 - w1) / w1)

    v = float(np.var(noise, ddof=1)) / n * (1.0 + r ** 2)
    noise = max(1.0, math.sqrt(V_TARGET / v)) * noise

    standard_error = math.sqrt(float(np.var(noise, ddof=1)) / n * (1.0 + r ** 2))
    mean_gap = (critical_target + margin) * standard_error

    return 10.0 + noise + mean_gap / 2.0, 10.0 + r * noise - mean_gap / 2.0


def plain_values(x1, x2, alpha):

    n1, n2 = len(x1), len(x2)
    a1, a2 = np.var(x1, ddof=1) / n1, np.var(x2, ddof=1) / n2
    v = a1 + a2
    t2 = (np.mean(x1) - np.mean(x2)) ** 2 / v
    df = v ** 2 / (a1 ** 2 / (n1 - 1) + a2 ** 2 / (n2 - 1))
    crit = stats.t.ppf(1.0 - alpha / 2.0, df)

    return t2, 1.0 / df, crit, t2 - crit ** 2


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=REPS)
    parser.add_argument("--dfs", default="all")
    parser.add_argument("--alphas", default="all")
    parser.add_argument("--margins", default="all")
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    args = parser.parse_args()

    from engine.HEdata import Message
    from operators.HEhypo_test import HEHypothesisTesting

    dfs = DFS if args.dfs == "all" else tuple(int(d) for d in args.dfs.split(","))
    alphas = ALPHAS if args.alphas == "all" else [float(a) for a in args.alphas.split(",")]
    margins = MARGINS if args.margins == "all" else tuple(float(m) for m in args.margins.split(","))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.out_dir / "boundary_decision_raw.csv"
    done = completed(raw_path, ("df", "alpha", "margin", "rep"), 1)

    engine = make_engine()
    ht = HEHypothesisTesting(engine)
    file, writer = open_append(raw_path, RAW_FIELDS)

    for df, alpha, margin in [(d, a, m) for d in dfs for a in alphas for m in margins]:
        x1, x2 = groups(alpha, margin, df)
        n1, n2 = len(x1), len(x2)
        R = float(max(np.max(x1), np.max(x2)) - min(np.min(x1), np.min(x2)))
        v = np.var(x1, ddof=1) / n1 + np.var(x2, ddof=1) / n2
        in_domain = int(_INV_SQRT_MIN <= v <= _INV_SQRT_MAX)

        p_t2, p_inv_df, p_crit, p_score = plain_values(x1, x2, alpha)
        base = {"df": df, "n1": n1, "n2": n2, "v": v, "in_domain": in_domain,
                "alpha": alpha, "margin": margin, "plain_t2": p_t2, "plain_inv_df": p_inv_df,
                "plain_critical": p_crit, "plain_score": p_score, "plain_decision": int(p_score > 0.0)}

        for rep in range(args.reps):
            if (str(df), str(alpha), str(margin), str(rep)) in done:
                continue

            c1 = engine.enc(Message(x1, engine.log_slots))
            c2 = engine.enc(Message(x2, engine.log_slots))
            tr, btr = {}, {}

            try:
                st_, st, sb = timed(ht, lambda: ht.HE_Welch_statistics(c1, c2, n1, n2, R, tr))
                step, dt, db = timed(ht, lambda: ht.HE_Welch_decision(st_, alpha, CRITICAL_LOG_DEGREE, btr))
            except RuntimeError as error:
                # Out-of-domain inputs can exceed the bootstrap range.
                writer.writerow({**base, "rep": rep, "error": str(error)[:120]})
                file.flush()
                print(f"df={df} alpha={alpha} margin={margin:+.3f} rep={rep} error: {error}", flush=True)
                continue

            # Reporting only (not timed).
            B = btr["score_bound"]
            he_t2 = first(engine, btr["t_sq"])
            he_inv_df = first(engine, tr["inv_df"])
            he_crit = math.sqrt(max(first(engine, btr["critical_sq"]), 0.0))
            he_score = first(engine, btr["score_normalized"]) * B
            he_step = first(engine, step)

            writer.writerow({
                **base, "rep": rep,
                "stats_time": st, "decision_time": dt, "stats_bootstraps": sb, "decision_bootstraps": db,
                "he_t2": he_t2, "t2_abs_err": abs(he_t2 - p_t2),
                "he_inv_df": he_inv_df, "inv_df_abs_err": abs(he_inv_df - p_inv_df),
                "he_critical": he_crit, "critical_abs_err": abs(he_crit - p_crit),
                "he_score": he_score, "score_abs_err": abs(he_score - p_score), "score_bound": B,
                "he_step": he_step, "score_decision": int(he_score > 0.0), "step_decision": int(he_step > 0.5),
            })
            file.flush()
            print(f"df={df} v={v:.3g} alpha={alpha} margin={margin:+.3f} rep={rep} "
                  f"step={he_step:.4f} plain={base['plain_decision']}", flush=True)

    file.close()
    print(f"saved {raw_path}")


if __name__ == "__main__":
    main()
