#!/usr/bin/env python3
"""Heuristic scan for the first sign change in Mertens' third theorem over K = Q(sqrt D) (and over Q for comparison).

sign Delta_K(x) = sign f_K(u) + O(u e^{-u/2}),  f_K(u) = 1 - 2/u - 2 Re sum_{gamma>0} x^{i gamma}/(rho - 1),
the sum running over the zeros of zeta_K = zeta * L(s, chi_D) (for K = Q: zeta only). Negative f => Delta_K < 0.
usage: field.py q n   (q.n = Conrey label of chi_D; 'Q' for the rationals)
"""
import math, os, sys
import numpy as np
from wide import lib, dp, G, jp

def zeros(lab):
    zz = [G]
    if lab != "Q":
        d = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "zeros", f"arb_{lab}", "certified.npz"))
        zz.append(0.5*(d["lo"] + d["hi"]))
    return np.sort(np.concatenate(zz))

def f_scan(g, T, kind, u0, u1, h, blk, full=False):
    g = g[g <= T]; rho = 0.5 + 1j*g
    w = jp(g/T) if kind == "jp" else np.ones_like(g)
    A = -2*w/(rho - 1)
    npts = int(math.ceil((u1 - u0)/h)) + 1; nb = (npts + blk - 1)//blk
    arrs = [np.ascontiguousarray(v, dtype=np.float64) for v in (g, A.real, A.imag, np.zeros_like(g), np.zeros_like(g))]
    bmin, bu = np.empty(nb), np.empty(nb); out = np.empty(npts) if full else None
    p = lambda v: v.ctypes.data_as(dp)
    lib.scan(*[p(v) for v in arrs], len(g), 0.0, float(u0), float(h), npts, blk, p(bmin), p(bu), p(out) if full else None)
    if full:
        u = u0 + h*np.arange(npts); return u, 1 - 2/u + out
    return bu, 1 - 2/bu + bmin   # block minima (2/u varies negligibly within a block)

lab = sys.argv[1]; U1 = float(sys.argv[2]) if len(sys.argv) > 2 else 1500
g = zeros(lab); Tmax = 2e4
bu, bm = f_scan(g, 2000, "jp", 10, U1, 5e-5, 2048)
order = np.argsort(bm); cands = []
for k in order:
    if all(abs(bu[k] - c) > 0.5 for c in cands): cands.append(bu[k])
    if len(cands) == 40: break
rows = []
for uc in cands:
    u, f = f_scan(g, Tmax, "jp", uc - 0.01, uc + 0.01, 0.1/Tmax, 256, full=True)
    k = f.argmin(); u1 = u[k]
    u, f = f_scan(g, Tmax, "full", u1 - 5e-4, u1 + 5e-4, 0.1/Tmax, 128, full=True)
    k = f.argmin(); neg = f < 0
    uj, fj = f_scan(g, Tmax, "jp", u1 - 5e-4, u1 + 5e-4, 0.1/Tmax, 128, full=True)
    rows.append((u[k], f[k], neg.sum()*0.1/Tmax, fj.min(), (fj < 0).sum()*0.1/Tmax))
rows.sort()
print(f"K = {'Q' if lab == 'Q' else 'Q(sqrt ' + lab.split('_')[0] + ')'}: {len(g[g <= Tmax])} zeros up to {Tmax:g}; candidate regions in increasing u")
for u, fv, wd, fj, wj in rows:
    flag = "  <-- negative" if fv < 0 else ""
    print(f"  u={u:10.5f} (10^{u/math.log(10):7.2f})  full f={fv:+.4f} (neg width {wd:.1e})  JP f={fj:+.4f} (neg width {wj:.1e}){flag}")
