#!/usr/bin/env python3
"""Choose (omega, eta1, lam) for the certificate near u = 495.7028331 from the binned zero sums:
our zeros <= 8e7 (../sim/bins_495.npz) plus LMFDB zeros in (8e7, T0] (bins/*.npz). Floating point.
usage: search_215.py T0 margin"""
import glob, math, os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "certification"))
from flint import arb
import error_bounds as eb
HERE = os.path.dirname(os.path.abspath(__file__))
U_REF = 495.7028331345
T0 = float(sys.argv[1]); margin = float(sys.argv[2])
d = np.load(os.path.join(HERE, "..", "sim", "bins_495.npz"))
parts = [(d["gb"], d["M0"], d["M1"], d["M2"])]
for p in sorted(glob.glob(os.path.join(HERE, "bins", "*.npz"))):
    e = np.load(p)
    if len(e["gb"]): parts.append((e["gb"], e["M0"], e["M1"], e["M2"]))
gb = np.concatenate([q[0] for q in parts]); M0 = np.concatenate([q[1] for q in parts])
M1 = np.concatenate([q[2] for q in parts]); M2 = np.concatenate([q[3] for q in parts])
# merge duplicate bins (a 2000-wide bin may be split between our data and a file, or between two files)
ub, inv = np.unique(gb, return_inverse=True)
A0 = np.zeros(len(ub), complex); A1 = np.zeros(len(ub), complex); A2 = np.zeros(len(ub), complex)
np.add.at(A0, inv, M0); np.add.at(A1, inv, M1); np.add.at(A2, inv, M2)
gb, M0, M1, M2 = ub, A0, A1, A2
sel = gb <= T0 + 1000; gb, M0, M1, M2 = gb[sel], M0[sel], M1[sel], M2[sel]
print(f"{len(gb)} bins up to {gb.max():.4g}")
def kap(xi, m):
    xi = np.asarray(xi, float)
    for k in range(1, m + 1): xi = np.where(np.abs(np.abs(xi) - k*np.pi) < 1e-7, xi + 2e-7, xi)
    r = np.sinc(xi/np.pi)
    for k in range(1, m + 1): r = r*(k*np.pi)**2/((k*np.pi)**2 - xi*xi)
    return r
DELTA, MM = 5e-10, 135
dy = np.arange(-6e-8, 6e-8, 2e-11); u = U_REF + dy
# full-weight function (all zeros <= T0)
def evaluate(w0, w1, w2):
    out = np.zeros(len(dy))
    for b0 in range(0, len(gb), 20000):
        sl = slice(b0, b0 + 20000); E = np.exp(1j*np.outer(dy, gb[sl]))
        out += np.real(E @ (w0[sl]*M0[sl] + w1[sl]*M1[sl] + 0.5*w2[sl]*M2[sl]) + 1j*dy*(E @ (w0[sl]*M1[sl] + w1[sl]*M2[sl])) - 0.5*dy**2*(E @ (w0[sl]*M2[sl])))
    return out
one = np.ones(len(gb)); zero = np.zeros(len(gb))
F = 1 + evaluate(one, zero, zero)
neg = F < 0
print(f"full weight to {T0:g}: min F = {F.min():+.6f} at u = {u[F.argmin()]:.10f}; negative on {neg.sum()*2e-11:.2e} total"
      + (f", from {u[neg][0]:.10f} to {u[neg][-1]:.10f}" if neg.any() else ""))
best = None
for c in [12, 14, 16, 20, 24]:
    for lam in [(4, 52, 44), (0, 0, 100), (14, 72, 14)]:
        eta1 = c/T0
        f = lambda g: sum(l/100*kap(eta1*g, m) for l, m in zip(lam, (1, 2, 3)) if l)*kap(DELTA*g, MM)*(g <= T0)
        w0 = f(gb); w1 = (f(gb + 1) - f(gb - 1))/2; w2 = f(gb + 1) - 2*w0 + f(gb - 1)
        S = evaluate(w0, w1, w2)
        la = [arb(l)/100 for l in lam]
        r4 = min(float(eb.bound_R4_envelope(arb(U_REF), arb(eta1), arb(T0), *la).mid()), float(eb.bound_R4_stieltjes(arb(U_REF), arb(eta1), arb(T0), *la).mid()))
        U = 1 + S + r4
        ok = np.nonzero(U < -margin)[0]
        line = f"c={c} lam={lam}: min U={U.min():+.6f} at u={u[U.argmin()]:.10f}; R4={r4:.2e}"
        if len(ok):
            om = u[ok[0]]; endp = om + eta1 + DELTA
            line += f"; first omega with U<-{margin}: {om:.11f} -> endpoint {endp:.11f} (10^{endp/math.log(10):.7f})"
            if best is None or endp < best[0]: best = (endp, om, c, lam, eta1, U.min())
        print(line, flush=True)
print("BEST:", best)
