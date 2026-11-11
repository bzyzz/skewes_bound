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
fail (indices into the input array). Rerunning reuses only chunks whose complete identity matches the request.
"""
import hashlib
import json
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


def chunk_identity(q, n, gs, i0, eps, dps):
    """Bind a chunk to its character, exact seeds, precision and implementation."""
    seeds = np.ascontiguousarray(gs, dtype="<f8")
    with open(__file__, "rb") as handle:
        code_hash = hashlib.sha256(handle.read()).hexdigest()
    return dict(schema=1, q=int(q), n=int(n), i0=int(i0), count=len(seeds),
                seeds_sha256=hashlib.sha256(seeds.tobytes()).hexdigest(),
                eps_hex=float(eps).hex(), dps=int(dps), code_sha256=code_hash)


def read_chunk(path, expected):
    """Reject stale, partial, foreign-character or unbound cached chunks."""
    with np.load(path, allow_pickle=False) as data:
        if "identity" not in data or json.loads(str(data["identity"])) != expected:
            raise ValueError(f"Cache identity mismatch: {path}; use a fresh output directory")
        lo, hi, fail = data["lo"].copy(), data["hi"].copy(), data["fail"].copy()
        if lo.shape != (expected["count"],) or hi.shape != lo.shape:
            raise ValueError(f"Cache length mismatch: {path}")
        good = np.isfinite(lo) & np.isfinite(hi)
        failed = np.nonzero(~good)[0] + expected["i0"]
        if not np.array_equal(failed, fail) or np.any(lo[good] >= hi[good]):
            raise ValueError(f"Invalid cached intervals/failure list: {path}")
        return lo, hi, fail


def work(args):
    q, n, gs, i0, out_dir, eps, dps = args
    path = os.path.join(out_dir, f"chunk_{i0:08d}.npz")
    identity = chunk_identity(q, n, gs, i0, eps, dps)
    if os.path.exists(path):
        _, _, failures = read_chunk(path, identity)
        return path, len(failures)
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
            values = [Z(t) for t in grid]
            if not all(imag.contains(0) for _, imag in values):
                raise ArithmeticError(f"Im Z excludes zero in fallback grid near {g}")
            sg = [sign(real) for real, _ in values]
            ch = [j for j in range(200) if sg[j] != 0 and sg[j + 1] != 0 and sg[j] != sg[j + 1]]
            if len(ch) == 1:
                lo[k], hi[k], ok = grid[ch[0]], grid[ch[0] + 1], True
        if not ok:
            lo[k], hi[k] = np.nan, np.nan
            fail.append(i0 + k)
    tmp = path + ".tmp.npz"
    np.savez(tmp, lo=lo, hi=hi, fail=np.array(fail, dtype=np.int64), i0=i0,
             identity=json.dumps(identity, sort_keys=True))
    os.replace(tmp, path)
    return path, len(fail)


def main():
    q, n, src, out_dir, nproc = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], int(sys.argv[5])
    eps = float(sys.argv[6]) if len(sys.argv) > 6 else 1e-8
    dps = int(sys.argv[7]) if len(sys.argv) > 7 else 20
    os.makedirs(out_dir, exist_ok=True)
    if nproc < 1 or dps < 15 or not np.isfinite(eps) or eps <= 0:
        raise ValueError("Require nproc >= 1, dps >= 15 and positive finite eps")
    g = np.load(src, allow_pickle=False)
    if g.ndim != 1 or len(g) == 0 or not np.all(np.isfinite(g)) or np.any(g <= 0):
        raise ValueError("Seeds must be a nonempty, positive, finite one-dimensional array")
    if not np.all(np.diff(g) > 2e-4):
        raise ValueError("Located zeros closer than 2e-4: certify that pair separately")
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
        identity = chunk_identity(q, n, g[i:i + CHUNK], i, eps, dps)
        a, b, failures = read_chunk(os.path.join(out_dir, f"chunk_{i:08d}.npz"), identity)
        lo.append(a); hi.append(b); fail.extend(failures.tolist())
    lo, hi = np.concatenate(lo), np.concatenate(hi)
    good = ~np.isnan(lo)
    if len(lo) != len(g) or not np.all(lo[good] < hi[good]):
        raise ValueError("Assembled interval count or width is invalid")
    if not np.all(hi[good][:-1] < lo[good][1:]):
        raise ValueError("Certified intervals overlap")
    np.savez(os.path.join(out_dir, "certified.npz"), lo=lo, hi=hi, fail=np.array(fail, dtype=np.int64),
             q=q, n=n, eps=eps, dps=dps,
             identity=json.dumps(chunk_identity(q, n, g, 0, eps, dps), sort_keys=True))
    print(f"q={q} n={n}: {good.sum()} certified zeros, {len(fail)} failures, max width {np.nanmax(hi - lo):.1e}, "
          f"{time.time() - t0:.0f}s")
    if fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
