#!/usr/bin/env python3
"""Certify zeros of L(1/2+it, chi) in Arb, seeded by float_zeros.py.

For a located zero g (a float64, hence an exact binary number) we evaluate, in ball arithmetic,
Z(t) = Re(exp(i theta(t)) L(1/2+it, chi)) at t = g - eps and t = g + eps. If the two balls have
certified opposite signs, Z has a zero of odd order in the open interval (g - eps, g + eps).
Z is real for real primitive chi (W_chi = 1), and we check that the imaginary part of the
ball contains 0. The intervals of consecutive zeros are checked to be disjoint, so every
certified interval contains a distinct zero.

On failure we retry with eps = 1e-6, then scan an Arb grid on [g - 1e-4, g + 1e-4] for a sign
change. Anything still failing is listed in the output and stops the Turing step.

usage: arb_certify.py q n float_zeros.npy out_dir nproc [eps=1e-8] [dps=20]
Output: out_dir/chunk_XXXXX.npz with arrays lo, hi (float64 interval endpoints, exact) and
fail (indices into the input array). Rerunning skips chunks that already exist.
"""
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

CHUNK = 2000


def make_Z(q, n, dps):
    import flint
    from flint import acb, arb, dirichlet_char
    flint.ctx.dps = dps
    chi = dirichlet_char(q, n)
    delta = int(chi.parity())
    lq = (arb(q)/arb.pi()).log()
    shift = (arb(1)/2 + delta)/2

    def Z(t):
        """Ball for Z(t) at the exact float64 t; returns (Re ball, Im ball)."""
        ta = arb(float(t))
        th = ta/2*lq + acb(shift, ta/2).lgamma().imag
        v = acb(0, th).exp()*acb.dirichlet_l(acb(arb(1)/2, ta), chi)
        return v.real, v.imag

    return Z


def sign(ball):
    if ball > 0:
        return 1
    if ball < 0:
        return -1
    return 0


def certify_one(Z, g, eps):
    zl, il = Z(g - eps)
    zr, ir = Z(g + eps)
    if not (il.contains(0) and ir.contains(0)):
        raise RuntimeError(f"Im Z ball excludes 0 near t={g}")
    sl, sr = sign(zl), sign(zr)
    return sl != 0 and sr != 0 and sl != sr


def work(args):
    q, n, gs, i0, out_dir, eps, dps = args
    path = os.path.join(out_dir, f"chunk_{i0:08d}.npz")
    if os.path.exists(path):
        return path, 0
    Z = make_Z(q, n, dps)
    lo = np.empty(len(gs))
    hi = np.empty(len(gs))
    fail = []
    for k, g in enumerate(gs):
        g = float(g)
        ok = False
        for e in (eps, 1e-6):
            if certify_one(Z, g, e):
                lo[k], hi[k], ok = g - e, g + e, True
                break
        if not ok:
            grid = g + np.linspace(-1e-4, 1e-4, 201)
            sg = [sign(Z(t)[0]) for t in grid]
            ch = [j for j in range(200) if sg[j] != 0 and sg[j + 1] != 0 and sg[j] != sg[j + 1]]
            if len(ch) == 1:
                lo[k], hi[k], ok = grid[ch[0]], grid[ch[0] + 1], True
        if not ok:
            lo[k], hi[k] = np.nan, np.nan
            fail.append(i0 + k)
    tmp = path + ".tmp.npz"
    np.savez(tmp, lo=lo, hi=hi, fail=np.array(fail, dtype=np.int64), i0=i0)
    os.replace(tmp, path)
    return path, len(fail)


def main():
    q, n, src, out_dir, nproc = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], int(sys.argv[5])
    eps = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-8
    dps = int(sys.argv[7]) if len(sys.argv) > 7 else 20
    os.makedirs(out_dir, exist_ok=True)
    g = np.load(src)
    assert np.all(np.diff(g) > 2e-4), "located zeros closer than 2e-4: certify that pair by hand"
    jobs = [(q, n, g[i:i + CHUNK], i, out_dir, eps, dps) for i in range(0, len(g), CHUNK)]
    jobs.sort(key=lambda j: -j[2][-1])        # tallest (slowest) chunks first
    t0 = time.time()
    nfail = 0
    with Pool(nproc) as p:
        for k, (path, nf) in enumerate(p.imap_unordered(work, jobs)):
            nfail += nf
            if k % 20 == 0:
                print(f"  {k + 1}/{len(jobs)} chunks, {nfail} failures, {time.time() - t0:.0f}s", flush=True)
    # assemble and check disjointness
    lo, hi, fail = [], [], []
    for i in range(0, len(g), CHUNK):
        d = np.load(os.path.join(out_dir, f"chunk_{i:08d}.npz"))
        lo.append(d["lo"]); hi.append(d["hi"]); fail.extend(d["fail"].tolist())
    lo, hi = np.concatenate(lo), np.concatenate(hi)
    good = ~np.isnan(lo)
    assert np.all(lo[good] < hi[good])
    assert np.all(hi[good][:-1] < lo[good][1:]), "certified intervals overlap"
    np.savez(os.path.join(out_dir, "certified.npz"), lo=lo, hi=hi, fail=np.array(fail, dtype=np.int64),
             q=q, n=n, eps=eps, dps=dps)
    print(f"q={q} n={n}: {good.sum()} certified zeros, {len(fail)} failures, max width {np.nanmax(hi - lo):.1e}, "
          f"{time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
