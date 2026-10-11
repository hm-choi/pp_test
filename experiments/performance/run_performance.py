"""Execution time and bootstrapping counts of the encrypted hypothesis tests.

For every dataset comparison (experiments/common/cases.py), test and
repetition: the statistics part and the decision part are timed separately
(encryption and decryption excluded), together with their bootstrap counts.
For Welch's t-test, ApproxTCrit alone is timed once more on the same
encrypted 1/df (it is part of the decision). The decision is decrypted after
timing only to check it against the plaintext decision.

Run (inside the HEaaN container, cwd = project root):
    python3 experiments/performance/run_performance.py
"""

import argparse
import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(PROJECT_ROOT))
sys.path.append(str(PROJECT_ROOT / "experiments"))

from common.cases import load_cases
from common.he_utils import completed, first, invsqrt_config, make_engine, open_append, timed
from common.plain import f_plain, welch_plain, z_plain


ALPHAS = [0.05]
CRITICAL_LOG_DEGREE = 4  # degree 15
REPS = 2
TESTS = ["welch", "f", "z"]

RESULT_DIR = Path(__file__).resolve().parent / "results"

RAW_FIELDS = [
    "case", "test", "rep", "alpha",
    "stats_time", "decision_time", "total_time", "stats_bootstraps", "decision_bootstraps",
    "critical_time", "plain_decision", "step_decision",
]


def run(ht, test, case, c1, c2, plain):

    n1, n2, R = case["n1"], case["n2"], case["R"]

    if test == "welch":
        st_, st, sb = timed(ht, lambda: ht.HE_Welch_statistics(c1, c2, n1, n2, R))
    elif test == "f":
        st_, st, sb = timed(ht, lambda: ht.HE_F_statistics(c1, c2, n1, n2, R))
    else:
        s1, s2 = np.var(case["x1"], ddof=1), np.var(case["x2"], ddof=1)
        st_, st, sb = timed(ht, lambda: ht.HE_Z_statistic(c1, c2, n1, n2, s1, s2, R))

    rows = []

    for alpha in ALPHAS:
        p = plain[alpha]

        if test == "welch":
            step, dt, db = timed(ht, lambda: ht.HE_Welch_decision(st_, alpha, CRITICAL_LOG_DEGREE))
            _, ct, _ = timed(ht, lambda: ht._critical_value_from_inv_df(
                st_["inv_df"], alpha, CRITICAL_LOG_DEGREE, min_output_level=5))
        elif test == "f":
            step, dt, db = timed(ht, lambda: ht.HE_F_decision(st_, p["fl"], p["fu"]))
            ct = ""
        else:
            step, dt, db = timed(ht, lambda: ht.HE_Z_decision(st_, alpha))
            ct = ""

        rows.append({"alpha": alpha, "stats_time": st, "decision_time": dt, "total_time": st + dt,
                     "stats_bootstraps": sb, "decision_bootstraps": db, "critical_time": ct,
                     "plain_decision": int(p["score"] > 0.0), "step": step})

    return rows


def main():

    global ALPHAS

    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=REPS)
    parser.add_argument("--cases", default="all")
    parser.add_argument("--tests", default=",".join(TESTS))
    parser.add_argument("--alphas", default=",".join(str(a) for a in ALPHAS))
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    args = parser.parse_args()

    ALPHAS = [float(a) for a in args.alphas.split(",")]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.out_dir / "performance_raw.csv"

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
        HEHypothesisTesting.check_inv_sqrt_domain(x1, x2, case["R"])

        plain = {
            "welch": {a: welch_plain(x1, x2, a) for a in ALPHAS},
            "f": {a: f_plain(x1, x2, a) for a in ALPHAS},
            "z": {a: z_plain(x1, x2, a, s1, s2) for a in ALPHAS},
        }

        for test in args.tests.split(","):
            ht = HEHypothesisTesting(engine, invsqrt_config() if test == "welch" else None)

            for rep in range(args.reps):
                if (case_name, test, str(rep)) in done:
                    continue

                c1 = engine.enc(Message(x1, engine.log_slots))
                c2 = engine.enc(Message(x2, engine.log_slots))
                rows = run(ht, test, case, c1, c2, plain[test])

                for row in rows:
                    row["step_decision"] = int(first(engine, row.pop("step")) > 0.5)
                    writer.writerow({"case": case_name, "test": test, "rep": rep,
                                     **{k: row.get(k, "") for k in RAW_FIELDS[3:]}})

                file.flush()
                r0 = rows[0]
                print(f"{case_name} {test} rep={rep} total={r0['total_time']:.1f}s "
                      f"bts={r0['stats_bootstraps'] + r0['decision_bootstraps']} "
                      f"match={int(r0['step_decision'] == r0['plain_decision'])}", flush=True)

    file.close()
    print(f"saved {raw_path}")


if __name__ == "__main__":
    main()
