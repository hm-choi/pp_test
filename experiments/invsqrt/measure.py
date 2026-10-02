"""Repeated measurement of HEApprox.invSqrt with HE-DAP selected parameters.

Configurations (per domain), taken from hedap_optimizer.py outputs:
- degree_sweep   : input level 12, no Pre-BTS, every degree 2^4-1 .. 2^9-1
                   with the HE-DAP (GetOptIter) iteration count;
- hedap_accuracy : HE-DAP u1 per input level (fastest under the MRE bound);
- hedap_speed    : HE-DAP u2 per input level (fastest regardless of MRE);
- hedap_accuracy_maxre_alg : supplementary u1 using the algorithmic MaxRE
                   (against 1/sqrt(dec(x))) instead of MRE;
- pptest_default : HEHypothesisTesting default (degree 63, 7 iterations).

Each configuration runs HEApprox.invSqrt (including the input mapping)
`--reps` times on the HE-DAP input set (linspace over the domain);
encryption / decryption are not timed. Errors are reported against the true
1/sqrt(x) and against 1/sqrt(dec(x)) (algorithmic error, excluding the input
encryption noise). One extra run per configuration evaluates a log-spaced
input set over the domain ("geom"), and a separate ciphertext evaluates
out-of-domain inputs in [lo / 10, lo) ("ood"). They are kept apart because
diverging Newton iterations on out-of-domain slots can overflow and corrupt
every slot.

Run (inside the HEaaN container, cwd = project root):
    python3 experiments/invsqrt/measure.py --domains raw
"""

import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[2]
sys.path.append(str(PROJECT_ROOT))
sys.path.append(str(Path(__file__).resolve().parent))

from engine.HEengine import HEengine
from engine.HEdata import Message
from operators.HEapprox import HEApprox

from domains import get_domains
from hedap_optimizer import (
    D_MIN, D_MAX, L_MAX, RESULT_DIR,
    hedap_data, decrypt, load_raw, select_iteration, open_append, map_input,
)


REPS = 10

RAW_FIELDS = [
    "domain", "method", "lo", "hi", "level", "degree", "pre_bts", "iteration", "tags",
    "input", "rep", "time", "bootstrap_count",
    "mre", "max_re", "mre_alg", "max_re_alg",
    "max_re_out", "neg_out",
]

KEY_FIELDS = ("domain", "level", "degree", "pre_bts", "iteration")


def config_key(c):

    return tuple(c[k] for k in KEY_FIELDS)


def build_configs(optimal, raw, domains):

    configs = {}

    def add(domain, level, degree, pre_bts, iteration, tag):
        c = {
            "domain": domain["name"], "method": domain["method"],
            "lo": domain["lo"], "hi": domain["hi"], "level": level,
            "degree": degree, "pre_bts": pre_bts, "iteration": iteration,
        }
        configs.setdefault(config_key(c), {**c, "tags": []})["tags"].append(tag)

    for domain in domains:
        name = domain["name"]

        for d in range(D_MIN, D_MAX + 1):
            rows = raw.get((name, L_MAX, 2 ** d - 1, 0))

            if rows is not None:
                i = select_iteration([r["mre"] for r in rows])

                if i is not None:
                    add(domain, L_MAX, 2 ** d - 1, 0, i, "degree_sweep")

        for level, entry in optimal.get(name, {}).get("levels", {}).items():
            for kind in ("accuracy", "speed", "accuracy_maxre_alg"):
                u = entry[kind]
                add(domain, int(level), u["degree"], u["pre_bts"], u["iteration"], f"hedap_{kind}")

        add(domain, L_MAX, 63, 0, 7, "pptest_default")

    return list(configs.values())


def run_once(engine, approx, config, x):

    dom = (config["lo"], config["hi"])
    log_degree = config["degree"].bit_length()

    ctxt = engine.enc(Message(x, engine.log_slots), level=config["level"])
    x_dec = decrypt(engine, ctxt, len(x))

    # One-off plaintext work (coefficients) outside the timer.
    approx.inv_sqrt_coeffs(log_degree, dom, config["method"])
    approx.reset_bootstrap_count()

    start = time.perf_counter()

    x_half, x_cheb = map_input(engine, ctxt, config)
    y = approx.invSqrt(
        x_half,
        x_cheb,
        log_degree=log_degree,
        iteration=config["iteration"],
        domain=dom,
        method=config["method"],
        pre_bts=bool(config["pre_bts"]),
    )

    elapsed = time.perf_counter() - start

    return elapsed, approx.bootstrap_count(), decrypt(engine, y, len(x)), x_dec


def errors(x, x_dec, y):

    re = np.abs(1.0 - y * np.sqrt(x))

    # Algorithmic error: reference computed from the decrypted (noisy) input.
    x_ref = np.maximum(x_dec, np.finfo(float).tiny)
    re_alg = np.abs(1.0 - y * np.sqrt(x_ref))

    return re, re_alg


def load_done(path):

    done = {}

    if not path.is_file():
        return done

    with path.open() as file:
        for row in csv.DictReader(file):
            key = (row["domain"], int(row["level"]), int(row["degree"]),
                   int(row["pre_bts"]), int(row["iteration"]))
            done.setdefault(key, set()).add((row["input"], int(row["rep"])))

    return done


def summarize(raw_path, summary_path):

    groups = {}

    with raw_path.open() as file:
        for row in csv.DictReader(file):
            key = tuple(row[k] for k in ("domain", "method", "lo", "hi") + KEY_FIELDS[1:])
            groups.setdefault(key, {"tags": row["tags"], "lin": [], "geom": [], "ood": []})[row["input"]].append(row)

    fields = [
        "domain", "method", "lo", "hi", "level", "degree", "pre_bts", "iteration", "tags", "reps",
        "time_mean", "time_std", "bootstrap_count",
        "mre", "max_re", "mre_alg", "max_re_alg",
        "geom_mre", "geom_max_re", "geom_max_re_alg", "out_max_re", "out_neg",
    ]

    rows = []

    for key, g in groups.items():
        lin = g["lin"]

        if not lin:
            continue

        f = lambda name, rs=lin: np.array([float(r[name]) for r in rs])
        geom = g["geom"][0] if g["geom"] else None
        ood = g["ood"][0] if g["ood"] else None

        rows.append({
            **dict(zip(fields[:8], key)),
            "tags": g["tags"],
            "reps": len(lin),
            "time_mean": f("time").mean(),
            "time_std": f("time").std(ddof=1) if len(lin) > 1 else 0.0,
            "bootstrap_count": int(np.median(f("bootstrap_count"))),
            "mre": f("mre").mean(),
            "max_re": f("max_re").max(),
            "mre_alg": f("mre_alg").mean(),
            "max_re_alg": f("max_re_alg").max(),
            "geom_mre": float(geom["mre"]) if geom else "",
            "geom_max_re": float(geom["max_re"]) if geom else "",
            "geom_max_re_alg": float(geom["max_re_alg"]) if geom else "",
            "out_max_re": float(ood["max_re_out"]) if ood else "",
            "out_neg": int(ood["neg_out"]) if ood else "",
        })

    rows.sort(key=lambda r: (r["domain"], -int(r["level"]), int(r["degree"]), int(r["pre_bts"])))

    with summary_path.open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    return rows


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--domains", default="raw", help="normalized, raw, all or comma-separated names")
    parser.add_argument("--reps", type=int, default=REPS)
    parser.add_argument("--dir", type=Path, default=None, help="directory of hedap_optimizer.py outputs")
    parser.add_argument("--summary-only", action="store_true")
    args = parser.parse_args()

    out_dir = args.dir or RESULT_DIR / args.domains.replace(",", "_")
    raw_path = out_dir / "measure_raw.csv"
    summary_path = out_dir / "measure_summary.csv"

    if not args.summary_only:
        domains = get_domains(args.domains)
        optimal = json.loads((out_dir / "hedap_optimal.json").read_text())
        configs = build_configs(optimal, load_raw(out_dir / "hedap_raw.csv"), domains)

        (out_dir / "measure_configs.json").write_text(json.dumps(configs, indent=2))
        print(f"{len(configs)} configurations x {args.reps} reps", flush=True)

        engine = HEengine(device_type="cpu", setting_root="/heaan_setting/")
        approx = HEApprox(engine)
        done = load_done(raw_path)

        file, writer = open_append(raw_path, RAW_FIELDS)

        for config in configs:
            lo, hi = config["lo"], config["hi"]
            key_done = done.get(config_key(config), set())

            inputs = {
                "lin": hedap_data(config, engine.num_slots),
                "geom": np.geomspace(lo, hi, engine.num_slots),
                "ood": np.geomspace(lo / 10, lo, engine.num_slots, endpoint=False),
            }

            plan = [("lin", r) for r in range(args.reps)] + [("geom", 0), ("ood", 0)]

            for input_name, rep in plan:
                if (input_name, rep) in key_done:
                    continue

                x = inputs[input_name]
                elapsed, bts, y, x_dec = run_once(engine, approx, config, x)
                re, re_alg = errors(x, x_dec, y)

                row = {
                    **{k: config[k] for k in ("domain", "method", "lo", "hi") + KEY_FIELDS[1:]},
                    "tags": "|".join(config["tags"]),
                    "input": input_name, "rep": rep,
                    "time": elapsed, "bootstrap_count": bts,
                    "mre": re.mean(), "max_re": re.max(),
                    "mre_alg": re_alg.mean(), "max_re_alg": re_alg.max(),
                    "max_re_out": "", "neg_out": "",
                }

                if input_name == "ood":
                    row.update(max_re_out=re.max(), neg_out=int(np.sum(y <= 0)))

                writer.writerow(row)
                file.flush()

            print(
                f"{config['domain']} l={config['level']:2d} deg={config['degree']:3d} "
                f"c={config['pre_bts']} iter={config['iteration']:2d} "
                f"[{'|'.join(config['tags'])}] done",
                flush=True,
            )

        file.close()

    rows = summarize(raw_path, summary_path)

    for r in rows:
        print(
            f"{r['domain']} l={r['level']:>2} deg={r['degree']:>3} c={r['pre_bts']} "
            f"iter={r['iteration']:>2} time={r['time_mean']:.2f}+-{r['time_std']:.2f}s "
            f"BTS={r['bootstrap_count']} MRE={r['mre']:.2e} MaxRE={r['max_re']:.2e} "
            f"MaxRE_alg={r['max_re_alg']:.2e} [{r['tags']}]"
        )

    print(f"saved {summary_path}")


if __name__ == "__main__":
    main()
