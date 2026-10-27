#!/usr/bin/env python3
"""What would it take to certify the region u ~ 495.7028331 (x ~ 10^215.28) for Mertens' first theorem?

F(u) = 1 - 2 Re sum_{gamma>0} x^{i gamma}/(rho - 1) (bias 1; lower-order terms < 1e-40 here).
Zeros <= 8e7: all 195,565,267 Arb zeros, binned (2000-wide bins, 2nd-order Taylor moments around U_REF).
Zeros in (8e7, 1.2e11]: random-phase Gaussian process with variance s^2 * 2/gamma^2 per zero (density log(g/2pi)/2pi).
For T0 in a list, criterion U(omega) = 1 + S_M(omega, T0) + R4(c, T0) (+ negligible R5): certifiable if min U < 0.
usage: sim215.py s n_real seed
"""
import glob, json, math, os, sys, time
os.environ.setdefault("OPENBLAS_NUM_THREADS", "16")
import numpy as np
REPO = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
sys.path.insert(0, os.path.join(REPO, "certification"))
from flint import arb
import error_bounds as eb

HERE = os.path.dirname(os.path.abspath(__file__))
U_REF = 495.7028331345
s_cal = float(sys.argv[1]); n_real = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
BIN = os.path.join(HERE, "bins_495.npz")

if not os.path.exists(BIN):
    t0 = time.time()
    low = np.load(os.path.join(REPO, "gammas_2325572.npy"))
    ext = []
    for p in sorted(glob.glob(os.path.join(REPO, "zero_cache_ext", "*.txt"))):
        with open(p) as f:
            ext.append(np.array([float(l[1:l.index(" ")]) for l in f]))
    ext = np.concatenate(ext); ext = ext[ext > low[-1]]
    hi_lo = [np.load(p) for p in sorted(glob.glob(os.path.join(REPO, "zero_cache_arb", "z_*.npz")))]
    hi = np.concatenate([d["hi"] for d in hi_lo]); lo = np.concatenate([d["lo"] for d in hi_lo])
    g_hi = np.concatenate([low, ext, hi]); g_lo = np.concatenate([np.zeros(len(low) + len(ext)), lo])
    assert np.all(np.diff(g_hi) >= 0)
    print(f"loaded {len(g_hi)} zeros in {time.time()-t0:.0f}s", flush=True)
    edges = np.arange(0, 8e7 + 2000, 2000.0)
    idx = np.searchsorted(g_hi, edges)
    gb = 0.5*(edges[:-1] + edges[1:])
    Z0 = np.zeros(len(gb), complex); Z1 = np.zeros(len(gb), complex); Z2 = np.zeros(len(gb), complex)
    Uh = np.longdouble(U_REF)
    for b in range(len(gb)):
        i0, i1 = idx[b], idx[b + 1]
        if i1 <= i0: continue
        gh = g_hi[i0:i1].astype(np.longdouble) + g_lo[i0:i1].astype(np.longdouble)
        ph = np.fmod(gh*Uh, np.longdouble(2*np.pi)).astype(np.float64)
        g = gh.astype(np.float64); dg = g - gb[b]
        z = -2*np.exp(1j*ph)/(-0.5 + 1j*g)
        Z0[b] = z.sum(); Z1[b] = (z*dg).sum(); Z2[b] = (z*dg*dg).sum()
    np.savez(BIN, gb=gb, M0=Z0, M1=Z1, M2=Z2)
    print(f"bins built in {time.time()-t0:.0f}s", flush=True)
d = np.load(BIN); gb, M0, M1, M2 = d["gb"], d["M0"], d["M1"], d["M2"]

# grid
P = 4e-6; N = 131072; h = P/N; dGam = 2*np.pi/P
dy0 = -1e-6; nj = int(2e-6/h); dy = dy0 + h*np.arange(nj); u = U_REF + dy

def kap(xi, m):
    xi = np.asarray(xi, float)
    for k in range(1, m + 1): xi = np.where(np.abs(np.abs(xi) - k*np.pi) < 1e-7, xi + 2e-7, xi)
    r = np.sinc(xi/np.pi)
    for k in range(1, m + 1): r = r*(k*np.pi)**2/((k*np.pi)**2 - xi*xi)
    return r
LAM = (0.04, 0.52, 0.44); DELTA, MM = 6e-10, 175
def Khat(g, eta1): return sum(l*kap(eta1*g, m) for l, m in zip(LAM, (1, 2, 3)))*kap(DELTA*g, MM)
T0S = [8e7, 2e8, 5e8, 1e9, 3e9, 1e10, 3e10]; CS = [12, 16, 24]
configs = [(T0, c) for T0 in T0S for c in CS]
R4 = {}
for T0, c in configs:
    R4[(T0, c)] = float(eb.bound_R4(arb(U_REF), arb(c)/arb(T0), arb(T0), *[arb(x) for x in ("0.04", "0.52", "0.44")]).mid())

cache = os.path.join(HERE, "meas_495.npz")
if os.path.exists(cache):
    meas = np.load(cache)["meas"]
else:
    cols = []
    for T0, c in configs:
        eta1 = c/T0; f = lambda g: Khat(g, eta1)*(g <= T0)
        K0 = f(gb); K1 = (f(gb + 1) - f(gb - 1))/2; K2 = f(gb + 1) - 2*K0 + f(gb - 1)
        cols.append((K0*M0 + K1*M1 + 0.5*K2*M2, K0*M1 + K1*M2, K0*M2))
    cols.append((M0, M1, M2))
    Aw = np.stack([c_[0] for c_ in cols], 1); Bw = np.stack([c_[1] for c_ in cols], 1); Cw = np.stack([c_[2] for c_ in cols], 1)
    meas = np.zeros((nj, len(cols)))
    for b0 in range(0, len(gb), 1000):
        sl = slice(b0, b0 + 1000); E = np.exp(1j*np.outer(dy, gb[sl]))
        meas += np.real(E @ Aw[sl] + 1j*dy[:, None]*(E @ Bw[sl]) - 0.5*dy[:, None]**2*(E @ Cw[sl]))
    np.savez(cache, meas=meas)
    print("measured part done", flush=True)

TMAX = 1.2e11; k = np.arange(int((TMAX - 8e7)/dGam)); gk = 8e7 + (k + 0.5)*dGam
var_k = 4*dGam*(np.log(gk/(2*np.pi))/(2*np.pi))/gk**2*s_cal**2
phase0 = np.exp(1j*k*dGam*dy0)
def synth(coef):
    full = np.zeros(N, complex); full[:len(coef)] = coef*phase0
    return np.real(np.exp(1j*(8e7 + 0.5*dGam)*dy)*(np.fft.ifft(full)[:nj]*N))
Kt = [Khat(gk, c/T0)*(gk <= T0) for T0, c in configs]
res = {f"{T0:g}": {"cert": 0, "endp": []} for T0 in T0S}; true_neg = 0; true_min = []
for r in range(n_real):
    cc = (rng.standard_normal(len(k)) + 1j*rng.standard_normal(len(k)))*np.sqrt(var_k/2)
    Ft = 1 + meas[:, -1] + synth(cc); true_min.append(Ft.min()); true_neg += Ft.min() < 0
    best = {f"{T0:g}": math.inf for T0 in T0S}
    for i, (T0, c) in enumerate(configs):
        U = 1 + meas[:, i] + synth(cc*Kt[i]) + R4[(T0, c)]
        ok = np.nonzero(U < 0)[0]
        if len(ok): best[f"{T0:g}"] = min(best[f"{T0:g}"], u[ok[0]] + c/T0 + DELTA)
    for kk, v in best.items():
        if math.isfinite(v): res[kk]["cert"] += 1; res[kk]["endp"].append(v)
print(f"s={s_cal}: P(true F < 0 somewhere in window) = {true_neg/n_real:.3f}; median true min {np.median(true_min):+.5f}")
for T0 in T0S:
    rr = res[f"{T0:g}"]; Nz = T0/(2*np.pi)*math.log(T0/(2*np.pi*math.e))
    med = np.median(rr["endp"]) if rr["endp"] else float("nan")
    print(f"  T0={T0:.0e} (N~{Nz:.2e} zeros): P(certifiable) = {rr['cert']/n_real:.3f}" + (f", median endpoint u={med:.9f} (10^{med/math.log(10):.6f})" if rr["endp"] else "")
          + f"; R4(c=12/16/24) = {R4[(T0,12)]:.1e}/{R4[(T0,16)]:.1e}/{R4[(T0,24)]:.1e}")
json.dump({"s": s_cal, "res": {k_: {"cert": v["cert"], "endp": v["endp"]} for k_, v in res.items()}, "true_min": true_min},
          open(os.path.join(HERE, f"sim215_s{s_cal}.json"), "w"))
