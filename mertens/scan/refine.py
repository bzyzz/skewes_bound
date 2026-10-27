#!/usr/bin/env python3
"""Refine wide-scan candidates: per candidate, minima of f1 with JP weights at T = 1e4, 1e5 and full weight at 1.3e6,
and the width of the region where the T = 1e5 JP value is negative (a proxy for certifiability with T0 ~ 1e5-1e6)."""
import math, sys
import numpy as np
from wide import scan, G, jp, lib, dp

d = np.load(sys.argv[1]); ncand = int(sys.argv[2])
order = np.argsort(d["bmin"]); cands = []
for k in order:
    if all(abs(d["bu"][k] - c) > 0.5 for c in cands):
        cands.append(float(d["bu"][k]))
    if len(cands) == ncand:
        break
rows = []
for uc in cands:
    b, u, _ = scan(1e4, uc - 0.02, uc + 0.02, 1e-5); k = b.argmin(); u1 = u[k]; v1 = b[k]
    b, u, f = scan(1e5, u1 - 1e-3, u1 + 1e-3, 1e-6, blk=64, full=True)
    uu = u1 - 1e-3 + 1e-6*np.arange(len(f)); k = f.argmin(); u2 = uu[k]; v2 = f[k]; w2 = (f < 0).sum()*1e-6
    # full weight at 1.3e6 (JP with T beyond the data ~ full weight up to 1.3e6)
    g = G; rho = 0.5 + 1j*g; A = -2/(rho - 1)
    import ctypes
    h = 0.2/1.3e6; half = 2e-4; npts = int(2*half/h) + 1; blk = 128; nb = (npts + blk - 1)//blk
    arrs = [np.ascontiguousarray(v, dtype=np.float64) for v in (g, A.real, A.imag, np.zeros_like(g), np.zeros_like(g))]
    bmin, bu, full = np.empty(nb), np.empty(nb), np.empty(npts); p = lambda v: v.ctypes.data_as(dp)
    lib.scan(*[p(v) for v in arrs], len(g), 0.0, float(u2 - half), float(h), npts, blk, p(bmin), p(bu), p(full))
    full = 1 + full; k = full.argmin(); u3 = u2 - half + h*k; v3 = full[k]; w3 = (full < 0).sum()*h
    rows.append((v3, uc, u3, v1, v2, w2, w3))
    print(f"u={u3:.7f} (10^{u3/math.log(10):.2f}): JP1e4 {v1:+.4f} | JP1e5 {v2:+.4f} (neg width {w2:.1e}) | full1.3e6 {v3:+.5f} (neg width {w3:.1e}) | f2=f1-2/u: {v3-2/u3:+.5f}", flush=True)
rows.sort()
print("\nRanked by full-weight min at 1.3e6:")
for r in rows[:12]:
    print(f"  u={r[2]:.7f} (10^{r[2]/math.log(10):.2f}) full={r[0]:+.5f} width={r[6]:.1e}  JP1e5={r[4]:+.4f} width={r[5]:.1e}")
