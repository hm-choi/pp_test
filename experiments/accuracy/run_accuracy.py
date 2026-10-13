"""Accuracy and decision agreement of the encrypted hypothesis tests.

For every dataset comparison (experiments/common/cases.py), test and
repetition, the encrypted statistics are computed once and every significance
level is evaluated as a separate decision. Copies of the intermediate
ciphertexts are decrypted afterwards (for the measurement only) and compared
with the plaintext (SciPy) values.

Tests:
- welch : s = T^2 - c^2 with the public bound max(R^2 / V_min, c_max^2),
          V_min the lower end of the invSqrt domain; reports T^2, 1/df, critical value, score, sign, step and
          the p-value from the decrypted T^2 and df;
- f     : s = (s1^2 - F_L s2^2)(s1^2 - F_U s2^2) with the public bound
          m_L m_U; reports score, sign and step;
- z     : s = Z^2 - z^2 with the public bound max(R^2 / V_public, z^2)
          (sample variances used as the public variances).

invSqrt (welch): V in [_INV_SQRT_MIN, _INV_SQRT_MAX] (coeffs/invSqrt_coeffs.py)
with the default degree and Newton iterations of HEHypothesisTesting.

Run (inside the HEaaN container, cwd = project root):
    python3 experiments/accuracy/run_accuracy.py
"""

import argparse
import csv
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(PROJECT_ROOT))
sys.path.append(str(PROJECT_ROOT / "experiments"))

from common.cases import load_cases
from common.he_utils import completed, first, make_engine, open_append
from common.plain import f_plain, welch_plain, z_plain


ALPHAS = [0.01, 0.025, 0.05, 0.1]
CRITICAL_LOG_DEGREE = 4  # degree 15
REPS = 2
TESTS = ["welch", "f", "z"]

RESULT_DIR = Path(__file__).resolve().parent / "results"

RAW_FIELDS = [
    "case", "test", "rep", "alpha",
    "plain_stat", "he_stat", "stat_abs_err", "stat_rel_err",
    "plain_inv_df", "he_inv_df", "inv_df_abs_err",
    "plain_critical", "he_critical", "critical_abs_err",
    "plain_score", "he_score", "score_abs_err", "score_bound",
    "he_sign", "he_step", "plain_decision", "step_decision", "score_decision",
    "plain_p", "he_p", "p_abs_err", "p_decision",
]


def run_welch(engine, ht, case, c1, c2):

    tr = {}
    st_ = ht.HE_Welch_statistics(c1, c2, case["n1"], case["n2"], case["R"], tr)
    branches = []

    for alpha in ALPHAS:
        btr = {}
        step = ht.HE_Welch_decision(st_, alpha, CRITICAL_LOG_DEGREE, btr)
        branches.append((alpha, step, btr))

    he_t2 = first(engine, branches[0][2]["t_sq"])
    he_inv_df = first(engine, tr["inv_df"])
    rows = []

    for alpha, step, btr in branches:
        B = btr["score_bound"]
        rows.append(dict(alpha=alpha, he_stat=he_t2, he_inv_df=he_inv_df,
                         he_critical=math.sqrt(max(first(engine, btr["critical_sq"]), 0.0)),
                         he_score=first(engine, btr["score_normalized"]) * B, score_bound=B,
                         he_sign=first(engine, btr["sign"]), he_step=first(engine, step),
                         he_p=2.0 * stats.t.sf(math.sqrt(max(he_t2, 0.0)), 1.0 / he_inv_df)))

    return rows


def run_f(engine, ht, case, c1, c2, plain):

    st_ = ht.HE_F_statistics(c1, c2, case["n1"], case["n2"], case["R"])
    rows = []

    for alpha in ALPHAS:
        btr = {}
        step = ht.HE_F_decision(st_, plain[alpha]["fl"], plain[alpha]["fu"], btr)
        B = btr["score_bound"]
        rows.append(dict(alpha=alpha, he_score=first(engine, btr["score_normalized"]) * B, score_bound=B,
                         he_sign=first(engine, btr["sign"]), he_step=first(engine, step)))

    return rows


def run_z(engine, ht, case, c1, c2):

    s1, s2 = np.var(case["x1"], ddof=1), np.var(case["x2"], ddof=1)
    tr = {}
    st_ = ht.HE_Z_statistic(c1, c2, case["n1"], case["n2"], s1, s2, case["R"], 0.0, tr)
    he_z = first(engine, tr["z"])
    rows = []

    for alpha in ALPHAS:
        btr = {}
        step = ht.HE_Z_decision(st_, alpha, btr)
        rows.append(dict(alpha=alpha, he_stat=he_z, he_critical=stats.norm.ppf(1.0 - alpha / 2.0),
                         he_score=first(engine, btr["score"]), score_bound=btr["score_bound"],
                         he_sign=first(engine, btr["sign"]), he_step=first(engine, step),
                         he_p=2.0 * stats.norm.sf(abs(he_z))))

    return rows


def complete(row, p, alpha):

    row["plain_score"] = p["score"]
    row["score_abs_err"] = abs(row["he_score"] - p["score"])

    if "stat" in p:
        row["plain_stat"] = p["stat"]
        row["plain_p"] = p["p"]
        row["stat_abs_err"] = abs(row["he_stat"] - p["stat"])
        row["stat_rel_err"] = row["stat_abs_err"] / abs(p["stat"])
        row["p_abs_err"] = abs(row["he_p"] - p["p"])
        row["p_decision"] = int(row["he_p"] < alpha)

    if "inv_df" in p:
        row["plain_inv_df"] = p["inv_df"]
        row["inv_df_abs_err"] = abs(row["he_inv_df"] - p["inv_df"])

    if "critical" in p:
        row["plain_critical"] = p["critical"]
        row["critical_abs_err"] = abs(row["he_critical"] - p["critical"])

    row["plain_decision"] = int(p["score"] > 0.0)
    row["step_decision"] = int(row["he_step"] > 0.5)
    row["score_decision"] = int(row["he_score"] > 0.0)

    return row


def summarize(raw_path, summary_path):

    groups = {}

    with raw_path.open() as file:
        for row in csv.DictReader(file):
            groups.setdefault((row["case"], row["test"], row["alpha"]), []).append(row)

    fields = ["case", "test", "alpha", "reps", "plain_stat", "stat_rel_err_max", "inv_df_abs_err_max",
              "critical_abs_err_max", "plain_score", "score_abs_err_max", "score_bound", "plain_p", "p_abs_err_max",
              "plain_decision", "step_match_rate", "score_decision_match_rate", "p_decision_match_rate",
              "step_min", "step_max"]

    def col(rs, name):
        vals = [float(r[name]) for r in rs if r[name] not in ("", None)]
        return np.array(vals) if vals else np.array([np.nan])

    rows = []

    for key, rs in groups.items():
        dec = int(rs[0]["plain_decision"])
        p_rows = [r for r in rs if r["p_decision"] != ""]
        rows.append({
            **dict(zip(fields[:3], key)), "reps": len(rs),
            "plain_stat": rs[0]["plain_stat"],
            "stat_rel_err_max": col(rs, "stat_rel_err").max(),
            "inv_df_abs_err_max": col(rs, "inv_df_abs_err").max(),
            "critical_abs_err_max": col(rs, "critical_abs_err").max(),
            "plain_score": rs[0]["plain_score"],
            "score_abs_err_max": col(rs, "score_abs_err").max(),
            "score_bound": rs[0]["score_bound"],
            "plain_p": rs[0]["plain_p"],
            "p_abs_err_max": col(rs, "p_abs_err").max(),
            "plain_decision": dec,
            "step_match_rate": np.mean([int(r["step_decision"]) == dec for r in rs]),
            "score_decision_match_rate": np.mean([int(r["score_decision"]) == dec for r in rs]),
            "p_decision_match_rate": np.mean([int(r["p_decision"]) == dec for r in p_rows]) if p_rows else "",
            "step_min": col(rs, "he_step").min(),
            "step_max": col(rs, "he_step").max(),
        })

    with summary_path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    return rows


def main():

    global ALPHAS

    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=REPS)
    parser.add_argument("--cases", default="all")
    parser.add_argument("--tests", default=",".join(TESTS))
    parser.add_argument("--alphas", default="all")
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()

    if args.alphas != "all":
        ALPHAS = [float(a) for a in args.alphas.split(",")]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.out_dir / "accuracy_raw.csv"

    if not args.summary_only:
        from engine.HEdata import Message
        from operators.HEhypo_test import HEHypothesisTesting

        engine = make_engine()
        cases = load_cases()

        if args.cases != "all":
            cases = {k: v for k, v in cases.items() if k in args.cases.split(",")}

        done = completed(raw_path, ("case", "test", "rep"), len(ALPHAS))
        file, writer = open_append(raw_path, RAW_FIELDS)

        for case_name, case in cases.items():
            x1, x2 = case["x1"], case["x2"]
            s1, s2 = np.var(x1, ddof=1), np.var(x2, ddof=1)

            plain = {
                "welch": {a: welch_plain(x1, x2, a) for a in ALPHAS},
                "f": {a: f_plain(x1, x2, a) for a in ALPHAS},
                "z": {a: z_plain(x1, x2, a, s1, s2) for a in ALPHAS},
            }

            for test in args.tests.split(","):
                ht = HEHypothesisTesting(engine)

                for rep in range(args.reps):
                    if (case_name, test, str(rep)) in done:
                        continue

                    c1 = engine.enc(Message(x1, engine.log_slots))
                    c2 = engine.enc(Message(x2, engine.log_slots))

                    if test == "welch":
                        rows = run_welch(engine, ht, case, c1, c2)
                    elif test == "f":
                        rows = run_f(engine, ht, case, c1, c2, plain[test])
                    else:
                        rows = run_z(engine, ht, case, c1, c2)

                    for row in rows:
                        row = complete(row, plain[test][row["alpha"]], row["alpha"])
                        writer.writerow({"case": case_name, "test": test, "rep": rep,
                                         **{k: row.get(k, "") for k in RAW_FIELDS[3:]}})

                    file.flush()
                    print(f"{case_name} {test} rep={rep} "
                          f"score_abs_err={max(r['score_abs_err'] for r in rows):.2e} "
                          f"match={sum(r['step_decision'] == r['plain_decision'] for r in rows)}/{len(rows)}",
                          flush=True)

        file.close()

    rows = summarize(raw_path, args.out_dir / "accuracy_summary.csv")
    print(f"{len(rows)} summary rows")


if __name__ == "__main__":
    main()
