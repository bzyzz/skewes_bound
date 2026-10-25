#!/usr/bin/env python3
"""
Extract the 2,325,572 certified zeta zeros (gamma_1 .. gamma_2325572, up to T0 = 1.3e6)
from $MERTENS_ZERO_CACHE into a fast float64 numpy array
for exploratory scans, while verifying contiguity and strict monotonicity.
"""
import os
import glob
from multiprocessing import Pool
import numpy as np
ROOT = os.path.dirname(os.path.abspath(__file__))

CACHE_DIR = os.environ.get("MERTENS_ZERO_CACHE", os.path.join(ROOT, "..", "mertens_bound", "certification", "zero_cache"))
OUT_PATH = os.path.join(ROOT, "gammas_2325572.npy")

def read_batch(args):
    n0, num, path = args
    gammas = []
    with open(path) as f:
        for line in f:
            gs = line.split("\t")[0]
            if gs.startswith("["):
                mid = gs[1:gs.index(" +/-")]
            else:
                mid = gs
            gammas.append(float(mid))
    assert len(gammas) == num, f"{path}: got {len(gammas)}, expected {num}"
    return (n0, np.array(gammas, dtype=np.float64))

def main():
    # Prefer 45-digit files when available, otherwise 20-digit
    by_start = {}
    for path in glob.glob(os.path.join(CACHE_DIR, "*.txt")):
        base = os.path.basename(path)[:-4]
        parts = base.split("_")
        if len(parts) != 3:
            continue
        n0, num, dps = int(parts[0]), int(parts[1]), int(parts[2])
        if (n0 - 1) % 250 != 0:
            continue
        if n0 < 2325501 and num != 250:
            continue
        if n0 not in by_start or (num, dps) > (by_start[n0][0], by_start[n0][1]):
            by_start[n0] = (num, dps, path)

    starts = sorted(by_start.keys())
    tasks = []
    expected_n0 = 1
    for n0 in starts:
        num, dps, path = by_start[n0]
        assert n0 == expected_n0, f"Gap at {expected_n0} -> {n0}"
        tasks.append((n0, num, path))
        expected_n0 += num

    print(f"Found {len(tasks)} contiguous batches covering N = {expected_n0 - 1} zeros.")
    with Pool(24) as p:
        res = p.map(read_batch, tasks)
    res.sort(key=lambda x: x[0])
    all_gammas = np.concatenate([r[1] for r in res])
    assert np.all(np.diff(all_gammas) > 0), "Not strictly increasing!"
    np.save(OUT_PATH, all_gammas)
    print(f"Saved {len(all_gammas)} zeros to {OUT_PATH}")
    print(f"gamma_1 = {all_gammas[0]:.12f}, gamma_{len(all_gammas)} = {all_gammas[-1]:.12f}")

if __name__ == "__main__":
    main()
