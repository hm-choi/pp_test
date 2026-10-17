"""Helpers shared by the HE experiments."""

import csv
import os
import time

import numpy as np

# HEaaN device ("cpu" or "gpu") and key directory, overridable per server.
HEAAN_DEVICE = os.environ.get("PPTEST_DEVICE", "cpu")
HEAAN_SETTING = os.environ.get("PPTEST_HEAAN_SETTING", "/heaan_setting/")


def make_engine():

    from engine.HEengine import HEengine

    return HEengine(device_type=HEAAN_DEVICE, setting_root=HEAAN_SETTING)


def first(engine, ctxt):

    # First slot of a decrypted ciphertext.
    return float(np.array(engine.dec(ctxt)[0], dtype=np.complex128).real[0])


def timed(ht, func):

    # Returns (result, elapsed seconds, bootstraps).
    ht.reset_bootstrap_count()
    start = time.perf_counter()
    out = func()
    return out, time.perf_counter() - start, ht.bootstrap_count()


def open_append(path, fields):

    # Append-mode CSV writer; drops a partial last line left by a hard kill.
    if path.is_file():
        data = path.read_bytes()

        if data and not data.endswith(b"\n"):
            path.write_bytes(data[:data.rfind(b"\n") + 1])

    new_file = not path.is_file() or path.stat().st_size == 0
    file = path.open("a", newline="")
    writer = csv.DictWriter(file, fieldnames=fields)

    if new_file:
        writer.writeheader()

    return file, writer


def completed(path, key_fields, per_key):

    # Keys (tuples of strings) that already have `per_key` rows in the CSV.
    if not path.is_file():
        return set()

    counts = {}

    with path.open() as file:
        for row in csv.DictReader(file):
            key = tuple(row[k] for k in key_fields)
            counts[key] = counts.get(key, 0) + 1

    return {k for k, n in counts.items() if n == per_key}
