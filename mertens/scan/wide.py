#!/usr/bin/env python3
"""Heuristic wide scan for sign changes in Mertens' first and second theorems over Q.

sqrt(x) * Delta_1(x) ~ 1 - 2 Re sum_{gamma>0} x^{i gamma}/(rho - 1)                (first theorem, exact bias 1)
sqrt(x) log x * Delta_2(x) ~ 1 - 2/u - 2 Re sum_{gamma>0} x^{i gamma}/(rho - 1)     (second theorem, to O(1/u^2))
Negative values indicate sign changes. Zeros are weighted by the Jurkat-Peyerimhoff taper at height T (smoothing).

usage: wide.py T u0 u1 out.npz
"""
import ctypes, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
lib = ctypes.CDLL(os.path.join(HERE, "libscan.so"))
dp = ctypes.POINTER(ctypes.c_double)
lib.scan.argtypes = [dp]*5 + [ctypes.c_long, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_long, ctypes.c_long, dp, dp, dp]
G = np.load(os.path.join(HERE, "..", "..", "gammas_2325572.npy"))   # first 2,325,572 zeros (see ../../README.md, Data)


def jp(x):
    x = np.minimum(x, 1.0)
    return (1 - x)*np.cos(np.pi*x) + np.sin(np.pi*x)/np.pi


def scan(T, u0, u1, h, blk=2048, full=False):
    g = G[G <= T]; rho = 0.5 + 1j*g
    A = -2*jp(g/T)/(rho - 1)
    npts = int(math.ceil((u1 - u0)/h)) + 1; nb = (npts + blk - 1)//blk
    arrs = [np.ascontiguousarray(v, dtype=np.float64) for v in (g, A.real, A.imag, np.zeros_like(g), np.zeros_like(g))]
    bmin, bu = np.empty(nb), np.empty(nb); out = np.empty(npts) if full else None
    p = lambda v: v.ctypes.data_as(dp)
    lib.scan(*[p(v) for v in arrs], len(g), 0.0, float(u0), float(h), npts, blk, p(bmin), p(bu), p(out) if full else None)
    return 1.0 + bmin, bu, (1.0 + out if full else None)


if __name__ == "__main__":
    T, u0, u1, out = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    bmin, bu, _ = scan(T, u0, u1, 0.1/T)
    np.savez(out, bmin=bmin, bu=bu, T=T)
    order = np.argsort(bmin); seen = []
    print(f"T={T:g}, u in [{u0},{u1}]: deepest separated minima of f1 (first theorem); f2 = f1 - 2/u")
    for k in order:
        if all(abs(bu[k] - s) > 0.5 for s in seen):
            seen.append(bu[k])
            print(f"  u={bu[k]:11.5f} (10^{bu[k]/math.log(10):8.2f})  f1={bmin[k]:+.4f}  f2={bmin[k]-2/bu[k]:+.4f}")
        if len(seen) == 25:
            break
