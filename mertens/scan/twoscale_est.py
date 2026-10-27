# Float estimate of a two-scale certificate for Mertens' first theorem near u = 1859.129 (stand-in error bounds from Skewes).
import os, sys, math, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "certification"))
from flint import arb
import error_bounds as eb
from wide import G, lib, dp
def khat(x, m):
    x = np.asarray(x, float)
    for k in range(1, m + 1): x = np.where(np.abs(np.abs(x) - k*np.pi) < 1e-7, x + 2e-7, x)
    r = np.sinc(x/np.pi)
    for k in range(1, m + 1): r = r*(k*np.pi)**2/((k*np.pi)**2 - x*x)
    return r
uc = 1859.1291847
for T0 in [2e4, 5e4, 1e5, 3e5, 1.3e6]:
    best = None
    for c in [12, 16, 24]:
        for lam in [(0.04, 0.52, 0.44), (0, 0, 1)]:
            eta1 = c/T0
            # mollifier: smallest delta with R5 < 1e-8 (M chosen near the optimum)
            M = int(uc/4) + 30
            delta = None
            for d in np.geomspace(5e-10, 1e-7, 60):
                try:
                    if eb.bound_R5(arb(uc), arb(eta1 + d), arb(d), M_moll=M) < arb("1e-8"): delta = d; break
                except Exception: pass
            g = G[G <= T0]; rho = 0.5 + 1j*g
            w = sum(l*khat(eta1*g, m) for l, m in zip(lam, (1, 2, 3)) if l)*khat(delta*g, M)
            A = -2*w/(rho - 1)
            half = 4e-4; h = min(eta1, 1e-5)/20; npts = int(2*half/h) + 1; blk = 256; nb = (npts + blk - 1)//blk
            arrs = [np.ascontiguousarray(v, dtype=np.float64) for v in (g, A.real, A.imag, np.zeros_like(g), np.zeros_like(g))]
            bmin, bu, full = np.empty(nb), np.empty(nb), np.empty(npts); p = lambda v: v.ctypes.data_as(dp)
            lib.scan(*[p(v) for v in arrs], len(g), 0.0, float(uc - half), float(h), npts, blk, p(bmin), p(bu), p(full))
            S = 1 + full
            lamA = [arb(str(l)) for l in lam]
            R4 = float(eb.bound_R4(arb(uc), arb(eta1), arb(T0), *lamA).mid())
            crit = S + R4 + 1e-8
            k = crit.argmin(); neg = np.nonzero(crit < 0)[0]
            u = uc - half + h*np.arange(npts); eta0 = eta1 + delta
            endp = u[neg[0]] + eta0 if len(neg) else math.inf
            rec = (endp, crit[k], c, lam, delta, M, R4)
            if best is None or rec < best: best = rec
    e = best[0]
    print(f"T0={T0:.1e} (N~{(G<=T0).sum() if T0<=1.3e6 else -1}): best min U={best[1]:+.4f} (c={best[2]}, lam={best[3]}, delta={best[4]:.1e}, M={best[5]}, R4={best[6]:.1e})"
          + (f" -> x < e^{e:.7f} = 10^{e/math.log(10):.4f}" if math.isfinite(e) else " -> not certified"), flush=True)
