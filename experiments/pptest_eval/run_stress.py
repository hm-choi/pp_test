"""Boundary-sensitive decision (stress) test for Welch's t-test.

Follows experiments/experiment1/test3.py: two synthetic groups of 128 samples
(10 + 5 * linspace(-1, 1, 128), shifted by +/- mean_gap / 2) whose Welch
t statistic is (t_crit + margin) * SE for margins in MARGINS around the
two-sided critical value t_crit at each alpha in ALPHAS.

For every alpha, margin and repetition, the encrypted statistics are computed once;
then for every critical-value degree 2^3-1 .. 2^7-1 the encrypted score
s = (mean1 - mean2)^2 - c^2 V (public bound max(R^2, c_max^2 V_max)) is
computed, and for degree 15 the full decision (sign / step) is evaluated.

invSqrt configuration: degree 63, 7 iterations (the HE-DAP tables are built
for the dataset domains only).

Run (inside the HEaaN container, cwd = project root):
    python3 experiments/pptest_eval/run_stress.py
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
sys.path.append(str(PROJECT_ROOT / "experiments" / "invsqrt"))
sys.path.append(str(Path(__file__).resolve().parent))

from run_eval import invsqrt_config, first, timed


ALPHAS = [0.001, 0.01, 0.025, 0.05, 0.1]
MARGINS = (-0.01, -0.005, -0.002, -0.001, 0.001, 0.002, 0.005, 0.01)
GROUP_SIZE = 128
LOG_DEGREES = range(3, 8)
STEP_LOG_DEGREE = 4
REPS = 10
CONFIGS = ["default"]

RESULT_DIR = Path(__file__).resolve().parent / "results"

RAW_FIELDS = [
    "alpha", "margin", "config", "rep", "degree", "full_decision",
    "stats_time", "branch_time", "stats_bootstraps", "branch_bootstraps",
    "plain_t2", "he_t2", "t2_abs_err", "plain_inv_df", "he_inv_df", "inv_df_abs_err",
    "plain_critical", "he_critical", "critical_abs_err",
    "plain_score", "he_score", "score_abs_err", "score_level",
    "he_step", "plain_decision", "score_decision", "step_decision",
    "score_bound",
]


def groups(alpha, margin):

    noise = 5.0 * np.linspace(-1.0, 1.0, GROUP_SIZE)
    variance = float(np.var(noise, ddof=1))
    df_target = float(2 * GROUP_SIZE - 2)
    critical_target = float(stats.t.ppf(1.0 - alpha / 2.0, df_target))
    standard_error = np.sqrt(2.0 * variance / GROUP_SIZE)
    mean_gap = (critical_target + margin) * standard_error

    return 10.0 + noise + mean_gap / 2.0, 10.0 + noise - mean_gap / 2.0


def plain_values(x1, x2, alpha):

    n1, n2 = len(x1), len(x2)
    a1, a2 = np.var(x1, ddof=1) / n1, np.var(x2, ddof=1) / n2
    v = a1 + a2
    t2 = (np.mean(x1) - np.mean(x2)) ** 2 / v
    df = v ** 2 / (a1 ** 2 / (n1 - 1) + a2 ** 2 / (n2 - 1))
    crit = stats.t.ppf(1.0 - alpha / 2.0, df)

    # s = (mean1 - mean2)^2 - c^2 V = V (T^2 - c^2).
    return t2, 1.0 / df, crit, (t2 - crit ** 2) * v


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--reps", type=int, default=REPS)
    parser.add_argument("--alphas", default="all")
    parser.add_argument("--margins", default="all")
    parser.add_argument("--out-dir", type=Path, default=RESULT_DIR)
    args = parser.parse_args()

    from engine.HEengine import HEengine
    from engine.HEdata import Message
    from operators.HEhypo_test import HEHypothesisTesting
    from hedap_optimizer import open_append

    args.out_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.out_dir / "stress_raw.csv"

    done = set()

    if raw_path.is_file():
        counts = {}
        with raw_path.open() as file:
            for row in csv.DictReader(file):
                k = (row["alpha"], row["margin"], row["config"], row["rep"])
                counts[k] = counts.get(k, 0) + 1
        per_run = len(LOG_DEGREES) + 1
        done = {k for k, n in counts.items() if n == per_run}

    engine = HEengine(device_type="cpu", setting_root="/heaan_setting/")
    cfgs = {name: invsqrt_config(name, None) for name in CONFIGS}

    alphas = ALPHAS if args.alphas == "all" else [float(a) for a in args.alphas.split(",")]
    margins = MARGINS if args.margins == "all" else tuple(float(m) for m in args.margins.split(","))

    file, writer = open_append(raw_path, RAW_FIELDS)

    for alpha, margin in [(a, m) for a in alphas for m in margins]:
        x1, x2 = groups(alpha, margin)
        n1, n2 = len(x1), len(x2)
        R = float(max(np.max(x1), np.max(x2)))
        HEHypothesisTesting.check_inv_sqrt_domain(x1, x2, R)
        p_t2, p_inv_df, p_crit, p_score = plain_values(x1, x2, alpha)

        for config_name, cfg in cfgs.items():
            ht = HEHypothesisTesting(engine, cfg)

            for rep in range(args.reps):
                if (str(alpha), str(margin), config_name, str(rep)) in done:
                    continue

                c1 = engine.enc(Message(x1, engine.log_slots))
                c2 = engine.enc(Message(x2, engine.log_slots))
                tr = {}

                st_, st, sb = timed(ht, lambda: ht.HE_Welch_statistics(c1, c2, n1, n2, R, tr))

                branches = []

                for ld in LOG_DEGREES:
                    # Critical value + score only (no sign).
                    def branch():
                        btr = {}
                        score = ht.welch_score(st_, alpha, ld, btr)
                        return btr["critical_scaled"], score, btr["score_bound"]

                    (crit, score, B), bt, bb = timed(ht, branch)
                    branches.append((2 ** ld - 1, 0, crit, score, None, bt, bb, B))

                btr = {}
                step, bt, bb = timed(ht, lambda: ht.HE_Welch_decision(st_, alpha, STEP_LOG_DEGREE, btr))
                branches.append((2 ** STEP_LOG_DEGREE - 1, 1,
                                 btr["critical_scaled"], btr["score_normalized"], step, bt, bb, btr["score_bound"]))

                # Reporting only (not timed).
                he_t2 = first(engine, ht.HE_Welch_t_squared(st_))
                he_inv_df = first(engine, tr["inv_df"])
                plain_decision = int(p_score > 0.0)

                for degree, full, crit, score, step, bt, bb, B in branches:
                    he_crit = math.sqrt(max(first(engine, crit) * B, 0.0))
                    he_score = first(engine, score) * B
                    he_step = first(engine, step) if step is not None else ""

                    writer.writerow({
                        "alpha": alpha, "margin": margin, "config": config_name, "rep": rep, "degree": degree,
                        "full_decision": full,
                        "stats_time": st, "branch_time": bt, "stats_bootstraps": sb, "branch_bootstraps": bb,
                        "plain_t2": p_t2, "he_t2": he_t2, "t2_abs_err": abs(he_t2 - p_t2),
                        "plain_inv_df": p_inv_df, "he_inv_df": he_inv_df, "inv_df_abs_err": abs(he_inv_df - p_inv_df),
                        "plain_critical": p_crit, "he_critical": he_crit, "critical_abs_err": abs(he_crit - p_crit),
                        "plain_score": p_score, "he_score": he_score, "score_abs_err": abs(he_score - p_score),
                        "score_bound": B,
                        "score_level": score.level(),
                        "he_step": he_step, "plain_decision": plain_decision,
                        "score_decision": int(he_score > 0.0),
                        "step_decision": int(he_step > 0.5) if step is not None else "",
                    })

                file.flush()

                print(f"alpha={alpha} margin={margin:+.3f} {config_name} rep={rep} stats={st:.1f}s "
                      f"step={first(engine, branches[-1][4]):.4f} plain={plain_decision}", flush=True)

    file.close()
    print(f"saved {raw_path}")


if __name__ == "__main__":
    main()
