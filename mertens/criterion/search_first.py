#!/usr/bin/env python3
"""Float parameter search for the first-theorem criterion near u = 1859.1291847 (the final U is recomputed in Arb).

U(omega) = 1 + 3 eps + S_M(omega) + R4 + R5 + R_small,  S_M = -2 sum_{0<gamma<=T0} hatK(gamma) Re(e^{i gamma omega}/(rho - 1)).
usage: search_first.py T0 [margin]
"""
import math, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scan"))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "certification"))
from wide import G, lib, dp
from flint import arb
import error_bounds as eb

def khat(x, m):
    x = np.asarray(x, float)
    for k in range(1, m + 1): x = np.where(np.abs(np.abs(x) - k*np.pi) < 1e-7, x + 2e-7, x)
    r = np.sinc(x/np.pi)
    for k in range(1, m + 1): r = r*(k*np.pi)**2/((k*np.pi)**2 - x*x)
    return r

T0 = float(sys.argv[1]); margin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
uc, delta, M = 1859.1291847, 2e-9, 494
g = G[G <= T0]; rho = 0.5 + 1j*g
best = None
for c in [12, 14, 16, 20]:
    for lam in [(4, 52, 44), (0, 0, 100), (14, 72, 14)]:
        eta1 = c/T0
        w = sum(l/100*khat(eta1*g, m) for l, m in zip(lam, (1, 2, 3)) if l)*khat(delta*g, M)
        A = -2*w/(rho - 1)
        half = 5e-4; h = eta1/40; npts = int(2*half/h) + 1; blk = 256; nb = (npts + blk - 1)//blk
        arrs = [np.ascontiguousarray(v, dtype=np.float64) for v in (g, A.real, A.imag, np.zeros_like(g), np.zeros_like(g))]
        bmin, bu, full = np.empty(nb), np.empty(nb), np.empty(npts); p = lambda v: v.ctypes.data_as(dp)
        lib.scan(*[p(v) for v in arrs], len(g), 0.0, float(uc - half), float(h), npts, blk, p(bmin), p(bu), p(full))
        la = [arb(l)/100 for l in lam]
        r4 = min(float(eb.bound_R4_envelope(arb(uc), arb(eta1), arb(T0), *la).mid()), float(eb.bound_R4_stieltjes(arb(uc), arb(eta1), arb(T0), *la).mid()))
        U = 1 + 3*1.5423e-9 + full + r4
        u = uc - half + h*np.arange(npts); ok = np.nonzero(U < -margin)[0]
        if not len(ok): print(f"c={c} lam={lam}: min U={U.min():+.4f}"); continue
        om = u[ok[0]]; endp = om + eta1 + delta
        print(f"c={c} lam={lam}: min U={U.min():+.4f}; first omega with U<-{margin}: {om:.9f} -> endpoint {endp:.9f} (10^{endp/math.log(10):.6f}); R4={r4:.2e}", flush=True)
        if best is None or endp < best[0]: best = (endp, om, c, lam, eta1)
print("BEST:", best)
