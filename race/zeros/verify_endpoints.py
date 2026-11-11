#!/usr/bin/env python3
"""Check every stored zero interval of L(1/2+it, chi) again, from the stored file alone.

arb_certify.py writes certified.npz but keeps no transcript of its sign evaluations. This script is the
consumer-side check: it reads only the stored float64 endpoints lo[i] < hi[i] and evaluates, in Arb ball
arithmetic, Z(t) = Re(e^{i theta(t)} L(1/2+it, chi)) at t = lo[i] and t = hi[i]. It requires, for every i,
  - the imaginary part of each ball contains 0 (Z is real for real primitive chi),
  - the two real balls have certified, opposite signs,
and, for the whole file, that the endpoints are finite, positive, and strictly ordered
(lo[0] < hi[0] < lo[1] < hi[1] < ...). Each interval then contains a zero of odd order, and the zeros are
distinct. Nothing from the float seeds or from the producer's run is trusted.

The working precision defaults to 30 digits; the producer used 20. The result is written next to the input
as endpoints_<q>_<n>.json, with the SHA-256 of the input file and of this script. As a consistency
diagnostic it also records whether sign Z(hi[i]) = sign Z(lo[i+1]) for all i, which must hold when the list
is complete and all zeros are simple (that completeness is what turing.py proves).

usage: verify_endpoints.py certified.npz nproc [dps=30]
"""
import hashlib
import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

CHUNK = 2000


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def make_Z(q, n, dps):
    import flint
    from flint import acb, arb, dirichlet_char
    flint.ctx.dps = dps
    chi = dirichlet_char(q, n)
    delta = int(chi.parity())
    lq = (arb(q)/arb.pi()).log()
    shift = (arb(1)/2 + delta)/2

    def Z(t):
        ta = arb(float(t))                    # a float64 is an exact binary number
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


def work(args):
    q, n, lo, hi, i0, dps, cache = args
    path = os.path.join(cache, f"chunk_{i0:08d}.npz")
    if os.path.exists(path):
        with np.load(path, allow_pickle=False) as d:
            if int(d["count"]) == len(lo):
                return i0, d["sl"].copy(), d["sr"].copy(), d["im_ok"].copy(), float(d["min_abs"])
    Z = make_Z(q, n, dps)
    sl = np.zeros(len(lo), dtype=np.int8)
    sr = np.zeros(len(lo), dtype=np.int8)
    im_ok = np.zeros(len(lo), dtype=bool)
    min_abs = float("inf")
    for k in range(len(lo)):
        zl, il = Z(lo[k])
        zr, ir = Z(hi[k])
        im_ok[k] = bool(il.contains(0) and ir.contains(0))
        sl[k], sr[k] = sign(zl), sign(zr)
        if sl[k] != 0 and sr[k] != 0:
            min_abs = min(min_abs, float(abs(zl).lower()), float(abs(zr).lower()))
    tmp = path + ".tmp.npz"
    np.savez(tmp, sl=sl, sr=sr, im_ok=im_ok, min_abs=min_abs, count=len(lo))
    os.replace(tmp, path)
    return i0, sl, sr, im_ok, min_abs


def main():
    path, nproc = sys.argv[1], int(sys.argv[2])
    dps = int(sys.argv[3]) if len(sys.argv) > 3 else 30
    if nproc < 1 or dps < 15:
        raise ValueError("require nproc >= 1 and dps >= 15")
    with np.load(path, allow_pickle=False) as d:
        q, n = int(d["q"]), int(d["n"])
        lo, hi = d["lo"].copy(), d["hi"].copy()
        nfail = len(d["fail"])
    digest = sha256_file(path)
    code = sha256_file(os.path.abspath(__file__))
    shape_ok = bool(lo.ndim == 1 and hi.shape == lo.shape and len(lo) > 0)
    order_ok = bool(shape_ok and nfail == 0 and np.all(np.isfinite(lo)) and np.all(np.isfinite(hi))
                    and np.all(lo > 0) and np.all(lo < hi) and np.all(hi[:-1] < lo[1:]))
    if not order_ok:
        raise SystemExit(f"{path}: endpoints are not finite, positive and strictly ordered, or failures are listed")
    # the cache is bound to the input file, this script and the precision
    cache = os.path.join(os.path.dirname(os.path.abspath(path)), f"endpoint_chunks_{digest[:12]}_{code[:12]}_{dps}")
    os.makedirs(cache, exist_ok=True)
    jobs = [(q, n, lo[i:i + CHUNK], hi[i:i + CHUNK], i, dps, cache) for i in range(0, len(lo), CHUNK)]
    jobs.sort(key=lambda j: -j[3][-1])        # tallest (slowest) chunks first
    SL = np.zeros(len(lo), dtype=np.int8)
    SR = np.zeros(len(lo), dtype=np.int8)
    IM = np.zeros(len(lo), dtype=bool)
    min_abs = float("inf")
    t0 = time.time()
    with Pool(nproc) as p:
        for k, (i0, sl, sr, im_ok, ma) in enumerate(p.imap_unordered(work, jobs)):
            SL[i0:i0 + len(sl)], SR[i0:i0 + len(sr)], IM[i0:i0 + len(sl)] = sl, sr, im_ok
            min_abs = min(min_abs, ma)
            if k % 20 == 0:
                print(f"  {q}.{n}: {k + 1}/{len(jobs)} chunks, {time.time() - t0:.0f}s", flush=True)
    good = IM & (SL != 0) & (SR != 0) & (SL != SR)
    bad = np.nonzero(~good)[0]
    import flint
    res = dict(schema=1, q=q, n=n, count=int(len(lo)), dps=dps, input_sha256=digest, code_sha256=code,
               python_flint=flint.__version__, t_min=float(lo[0]), t_max=float(hi[-1]),
               max_width=float(np.max(hi - lo)), order_ok=order_ok, sign_changes=int(good.sum()),
               failures=[int(i) for i in bad[:1000]], n_failures=int(len(bad)),
               min_abs_Z_at_endpoints=min_abs,
               signs_sha256=hashlib.sha256(SL.tobytes() + SR.tobytes()).hexdigest(),
               consecutive_signs_agree=bool(np.all(SR[:-1] == SL[1:])),
               all_ok=bool(len(bad) == 0), seconds=round(time.time() - t0))
    out = os.path.join(os.path.dirname(os.path.abspath(path)), f"endpoints_{q}_{n}.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    print(f"{q}.{n}: {res['sign_changes']}/{res['count']} sign changes certified, {len(bad)} failures, "
          f"min |Z| at endpoints {min_abs:.2e}, consecutive signs agree: {res['consecutive_signs_agree']}, "
          f"{res['seconds']}s -> {out}", flush=True)
    if not res["all_ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
