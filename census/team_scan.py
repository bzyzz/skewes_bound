#!/usr/bin/env python3
"""Heuristic scan for the first sign change of the residue vs non-residue ("team") race modulo q.

Residues R = squares of units mod q. avg_R pi(x;q,.) - avg_NR pi(x;q,.) has the sign of sum_{chi real, chi != chi0} pi(x, chi),
and E(u) := u e^{-u/2} sum_chi pi(e^u, chi) ~ -(r - 1) Phi(u) - sum_chi 2 Re sum_{gamma > 0} e^{i gamma u} (1/rho + 1/(rho^2 u)),
where r = number of real characters mod q (the index of the squares). Residues lead iff E > 0.
Stage 1: Jurkat-Peyerimhoff weights at height T1 on a grid of step 0.25/T1. Stage 2: blocks whose JP maximum exceeds
-MARGIN are rescanned with full weights and all zeros <= T2 (step 0.02/T2). Floating point; proves nothing.
usage: team_scan.py q UMAX [T1 T2 MARGIN]
"""
import ctypes, glob, json, math, os, sys
import numpy as np
from flint import dirichlet_char

HERE = os.path.dirname(os.path.abspath(__file__))
_lib = ctypes.CDLL(os.path.join(HERE, "libscan.so"))
_dp = ctypes.POINTER(ctypes.c_double)
_lib.scan.argtypes = [_dp]*5 + [ctypes.c_long, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_long, ctypes.c_long, _dp, _dp, _dp]


def real_prims(q):
    """Labels (d, n) of the real primitive characters with conductor d > 1 dividing q (one per real char mod q)."""
    out = []
    for d in range(3, q + 1):
        if q % d:
            continue
        for n in range(1, d + 1):
            if math.gcd(n, d) == 1:
                c = dirichlet_char(d, n)
                if c.is_primitive() and c.is_real():
                    out.append((d, n))
    return out


def load_zeros(d, n, T):
    cands = sorted(glob.glob(os.path.join(HERE, "zeros", f"z_{d}_{n}_T*.npy")), key=lambda p: -float(p.split("_T")[1][:-4]))
    for p in cands:
        z = np.load(p)
        if z[-1] >= 0.98*T:
            return z[z <= T]
    raise FileNotFoundError(f"no zeros for {d}.{n} up to {T}")


def weights(g, kind, T):
    x = g/T
    if kind == "full":
        return (x <= 1).astype(float)
    x = np.minimum(x, 1.0)
    return (1 - x)*np.cos(np.pi*x) + np.sin(np.pi*x)/np.pi


def coeffs(q, kind, T):
    gs, A, B = [], [], []
    for (d, n) in real_prims(q):
        g = load_zeros(d, n, T)
        w = weights(g, kind, T)
        rho = 0.5 + 1j*g
        gs.append(g); A.append(2*w/rho); B.append(2*w/rho**2)   # coefficients of -E
    return np.concatenate(gs), np.concatenate(A), np.concatenate(B)


def scan_negE(q, kind, T, u0, u1, h, blk=2048, full=False):
    """min over blocks of -E(u) (so E > 0 <=> -E < 0)."""
    g, A, B = coeffs(q, kind, T)
    c0 = float(len(real_prims(q)))   # (r - 1): r - 1 nonprincipal real characters
    npts = int(math.ceil((u1 - u0)/h)) + 1
    nblk = (npts + blk - 1)//blk
    arrs = [np.ascontiguousarray(x, dtype=np.float64) for x in (g, A.real, A.imag, B.real, B.imag)]
    bmin = np.empty(nblk); bu = np.empty(nblk); out = np.empty(npts) if full else None
    p = lambda x: x.ctypes.data_as(_dp)
    _lib.scan(*[p(x) for x in arrs], len(g), c0, float(u0), float(h), npts, blk, p(bmin), p(bu), p(out) if full else None)
    return bmin, bu, out


def main():
    q, UMAX = int(sys.argv[1]), float(sys.argv[2])
    T1 = float(sys.argv[3]) if len(sys.argv) > 3 else 2000
    T2 = float(sys.argv[4]) if len(sys.argv) > 4 else 1e4
    MARGIN = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3
    U0 = 20.0
    h1 = 0.25/T1
    bmin, bu, _ = scan_negE(q, "jp", T1, U0, UMAX, h1)
    width = 2048*h1
    cand = np.nonzero(bmin < MARGIN)[0]          # increasing u
    rows = []; extra = None
    for b in cand:
        if extra is not None:
            extra -= 1
            if extra < 0:
                break
        lo = U0 + b*width - width; hi = lo + 3*width
        h2 = 0.02/T2
        m2, u2, full = scan_negE(q, "full", T2, lo, hi, h2, full=True)
        k = int(np.argmin(full)); neg = full < 0
        rows.append(dict(u_jp=float(bu[b]), jp_max=float(-bmin[b]), u=float(lo + k*h2), full_max=float(-full[k]),
                         width_pos=float(neg.sum()*h2)))
        if extra is None and rows[-1]["full_max"] > 0:
            extra = int(os.environ.get("AFTER_HIT", "5"))   # refine a few more candidates after the first hit
    rows.sort(key=lambda r: r["u"])
    hits = [r for r in rows if r["full_max"] > 0]
    res = dict(stopped_after_hit=extra is not None, n_jp_candidates=int(len(cand)), q=q, chars=[f"{d}.{n}" for d, n in real_prims(q)], r=len(real_prims(q)) + 1, UMAX=UMAX, T1=T1, T2=T2,
               n_candidates=len(rows), first=hits[0] if hits else None, candidates=rows[:40],
               best_jp=dict(u=float(bu[np.argmin(bmin)]), jp_max=float(-bmin.min())))
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(HERE, "scans", f"team_q{q}_U{int(UMAX)}.json"), "w"), indent=1)


if __name__ == "__main__":
    os.makedirs(os.path.join(HERE, "scans"), exist_ok=True)
    main()
