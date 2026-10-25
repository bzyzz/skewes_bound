#!/usr/bin/env python3
"""
Compute and cache certified Arb enclosures of zeta zero ordinates gamma_n
for n = 2,325,573 .. N_MAX (extending the existing 2,325,572 certified zeros in mertens_bound).
Each batch of 2000 zeros is verified to lie on Re(rho) = 1/2 and be strictly increasing,
and saved with rigorous outward-rounded Arb ball strings (radius=True).
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
import sys
import time
from multiprocessing import Pool
import flint
from flint import acb, arb
ROOT = os.path.dirname(os.path.abspath(__file__))

CACHE_DIR = os.path.join(ROOT, "zero_cache_ext")
DPS = 16
BATCH = 2000
N_START = 2325573
N_END = int(sys.argv[1]) if len(sys.argv) > 1 else 13955181
WORKERS = int(sys.argv[2]) if len(sys.argv) > 2 else 56

def compute_batch(args):
    n0, num = args
    path = os.path.join(CACHE_DIR, f"{n0:09d}_{num}_{DPS}.txt")
    if os.path.exists(path):
        return (n0, num, True)
    flint.ctx.dps = DPS
    zs = acb.zeta_zeros(n0, num)
    lines = []
    prev_g = None
    for rho in zs:
        if not (rho.real == arb("0.5")):
            raise RuntimeError(f"Zero near {n0} not certified on critical line!")
        g = rho.imag
        if prev_g is not None and not (prev_g < g):
            raise RuntimeError(f"Zeros near {n0} not strictly increasing!")
        prev_g = g
        lines.append(g.str(DPS + 4, radius=True) + "\n")
    tmp = f"{path}.tmp.{os.getpid()}"
    with open(tmp, "w") as f:
        f.writelines(lines)
    os.replace(tmp, path)
    return (n0, num, False)

def main():
    os.makedirs(CACHE_DIR, exist_ok=True)
    tasks = []
    for n0 in range(N_START, N_END + 1, BATCH):
        num = min(BATCH, N_END + 1 - n0)
        tasks.append((n0, num))
    print(f"Building certified zero cache for n = {N_START} .. {N_END} ({len(tasks)} batches of {BATCH}) with {WORKERS} workers at dps={DPS}...", flush=True)
    t0 = time.time()
    with Pool(WORKERS) as p:
        for i, (n0, num, cached) in enumerate(p.imap_unordered(compute_batch, tasks), 1):
            if i % 100 == 0 or i == len(tasks):
                dt = time.time() - t0
                rate = (i * BATCH) / max(dt, 1e-3)
                print(f"  [{i:4d}/{len(tasks)}] n0={n0} | {dt:6.1f}s | {rate:,.0f} zeros/s", flush=True)
    print(f"DONE in {time.time()-t0:.1f}s!", flush=True)

if __name__ == "__main__":
    main()
