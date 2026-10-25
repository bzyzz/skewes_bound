#!/usr/bin/env python3
"""Float parameter search for the race criterion; the final U(omega) is recomputed in Arb by certify.py.

usage: search.py T0 zero_source [u0 u1]
  zero_source: 'float5e4' (heuristics/zeros/*_5e4.npy, from heuristics/zeros_L.py) or 'cert' (data/arb_*/certified.npz midpoints)
"""
import ctypes, os, sys, math
import numpy as np
from flint import arb
import race_criterion as rc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
lib = ctypes.CDLL(os.path.join(HERE, "libscan.so"))
dp = ctypes.POINTER(ctypes.c_double)
lib.scan.argtypes = [dp]*5 + [ctypes.c_long, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                              ctypes.c_long, ctypes.c_long, dp, dp, dp]


def load_zeros(src):
    out = []
    for (q, n, d, A, _) in rc.CHARS:
        if src == "float5e4":
            g = np.load(os.path.join(ROOT, "heuristics/zeros", f"z_{q}_{n}_5e4.npy"))
        else:
            c = np.load(os.path.join(ROOT, "data", f"arb_{q}_{n}", "certified.npz"))
            g = 0.5*(c["lo"] + c["hi"])
        out.append(g)
    return out


def scan(g, A, B, u0, u1, h):
    npts = int(math.ceil((u1 - u0)/h)) + 1
    blk = 64
    nblk = (npts + blk - 1)//blk
    arrs = [np.ascontiguousarray(x, dtype=np.float64) for x in (g, A.real, A.imag, B.real, B.imag)]
    bmin, bu, full = np.empty(nblk), np.empty(nblk), np.empty(npts)
    p = lambda x: x.ctypes.data_as(dp)
    lib.scan(*[p(x) for x in arrs], len(g), 0.0, float(u0), float(h), npts, blk, p(bmin), p(bu), p(full))
    return u0 + h*np.arange(npts), full


def zero_sums(zs, T0):
    s2, s3 = [], []
    for (q, n, d, A, _), g in zip(rc.CHARS, zs):
        g = g[g <= T0]
        NT0 = 2*len(g)
        T = arb(T0)
        s2.append(arb(2*np.sum(g**-2.0)) + rc.sum_inv_pow_tail(T, q, 2, NT0))
        s3.append(arb(2*np.sum(g**-3.0)) + rc.sum_inv_pow_tail(T, q, 3, NT0))
    return s2, s3


def main():
    T0 = float(sys.argv[1]); src = sys.argv[2]
    u0, u1 = (float(sys.argv[3]), float(sys.argv[4])) if len(sys.argv) > 4 else (111.538, 111.554)
    zs = load_zeros(src)
    s2, s3 = zero_sums(zs, T0)
    g_all = np.concatenate([g[g <= T0] for g in zs])
    rho = 0.5 + 1j*g_all
    best = None
    LAMS = [(0, 0, 1), (0, 1, 0), (0.04, 0.52, 0.44), (0.2, 0.4, 0.4)]
    for M in (40,):
        for c in (10, 11, 12, 13, 14, 16, 20):
            eta1 = c/T0
            u_lo_c, u_hi_c = u0 - 3*eta1, u1 + 3*eta1
            delta = rc.min_delta(u_lo_c, u_hi_c, M, arb("1e-5"))
            for lam in LAMS:
                w = rc.Khat_f(g_all, eta1, delta, M, lam)
                A_ = 4*w/rho; B_ = 4*w/rho**2
                h = min(eta1, delta)/30
                u, S = scan(g_all, A_, B_, u0, u1, h)     # S = 2 sum_chi S_chi
                err = rc.errors(0.5*(u0 + u1), eta1, delta, M, T0, lam, s2, s3)
                eta0 = eta1 + delta
                # main term uses u_lo = omega - eta0; Phi_D is decreasing, so evaluating at u0 - eta0 is safe
                fixed = float(err["total"].mid()) - float(err["main"].mid()) + 8*float(rc.Phi_D(arb(u0 - eta0)).mid())
                crit = S + fixed
                k = crit.argmin()
                neg = np.nonzero(crit < 0)[0]
                endp = u[neg[0]] + eta0 if len(neg) else np.inf
                rec = (endp, crit[k], u[k], c, lam, delta, M, float(err["R4"].mid()), float(err["R3"].mid()), float(err["R5"].mid()))
                print(f"c={c:2d} lam={lam} delta={delta:.2e}: min U={crit[k]:+.5f} at {u[k]:.6f}; S_min={S[k]:+.4f}; "
                      f"R4={rec[7]:.2e} R3={rec[8]:.1e} R5={rec[9]:.1e}; endpoint {'%.6f' % endp if len(neg) else '-'}", flush=True)
                if best is None or rec < best:
                    best = rec
    print("BEST:", best, "-> 10^%.4f" % (best[0]/math.log(10)) if np.isfinite(best[0]) else "")


if __name__ == "__main__":
    main()
