#!/usr/bin/env python3
"""Float parameter search for the team race mod q near u = uc (final L recomputed in Arb by team_certify.py).
usage: team_search.py q T0 uc"""
import ctypes, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "criterion")); sys.path.insert(0, HERE)
from flint import arb
import race_criterion as rc, team_criterion as tc
lib = ctypes.CDLL(os.path.join(ROOT, "criterion", "libscan.so")); dp = ctypes.POINTER(ctypes.c_double)
lib.scan.argtypes = [dp]*5 + [ctypes.c_long, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_long, ctypes.c_long, dp, dp, dp]
q, T0, uc = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
cfg = tc.MODULI[q]; chars = cfg["chars"]
zs = []
for (cq, cn, d, A, _) in chars:
    c = np.load(os.path.join(ROOT, "data", f"arb_{cq}_{cn}", "certified.npz")); g = 0.5*(c["lo"] + c["hi"]); zs.append(g[g <= T0])
s2, s3 = [], []
for (cq, cn, d, A, _), g in zip(chars, zs):
    s2.append(arb(2*np.sum(g**-2.0)) + rc.sum_inv_pow_tail(arb(T0), cq, 2, 2*len(g)))
    s3.append(arb(2*np.sum(g**-3.0)) + rc.sum_inv_pow_tail(arb(T0), cq, 3, 2*len(g)))
g_all = np.concatenate(zs); rho = 0.5 + 1j*g_all
best = None
rc_phi = lambda v: 1 + 2/v + 8/v**2 + 58.56/v**3      # Phi_D, decreasing
for M in (40, 60):
    for c in (12, 14, 16, 20, 24):
        eta1 = c/T0
        delta = tc.min_delta(chars, uc - 0.01, uc + 0.01, M, arb("1e-6"))
        for lam in [(4, 52, 44), (0, 0, 100), (0, 100, 0), (14, 72, 14)]:
            w = rc.Khat_f(g_all, eta1, delta, M, tuple(l/100 for l in lam))
            A_ = 2*w/rho; B_ = 2*w/rho**2
            half = 4e-3; h = min(eta1, delta)/20; npts = int(2*half/h) + 1; blk = 256; nb = (npts + blk - 1)//blk
            arrs = [np.ascontiguousarray(v, dtype=np.float64) for v in (g_all, A_.real, A_.imag, B_.real, B_.imag)]
            bmin, bu, full = np.empty(nb), np.empty(nb), np.empty(npts); p = lambda v: v.ctypes.data_as(dp)
            lib.scan(*[p(v) for v in arrs], len(g_all), 0.0, float(uc - half), float(h), npts, blk, p(bmin), p(bu), p(full))
            err = tc.errors(q, uc, eta1, delta, M, T0, [arb(l)/100 for l in lam], s2, s3)
            u = uc - half + h*np.arange(npts); eta0 = eta1 + delta
            other = float(err["total"].mid()) - float(err["bias"].mid())
            crit = full + other + (cfg["phi"] - 1)*rc_phi(u - eta0)
            k = crit.argmin(); ok = np.nonzero(crit < -0.01)[0]
            endp = u[ok[0]] + eta0 if len(ok) else math.inf
            line = f"M={M} c={c:2d} lam={lam} delta={delta:.2e}: min crit={crit[k]:+.4f} at {u[k]:.6f}; R4={float(err['R4'].mid()):.1e} R5={float(err['R5'].mid()):.1e}"
            if len(ok): line += f"; first omega with crit<-0.01 -> endpoint {endp:.7f} (10^{endp/math.log(10):.5f})"
            print(line, flush=True)
            rec = (endp, u[ok[0]] if len(ok) else None, M, c, lam, delta, crit[k])
            if best is None or rec[0] < best[0]: best = rec
print("BEST:", best)
