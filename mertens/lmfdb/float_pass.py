#!/usr/bin/env python3
"""Floating-point pass over LMFDB files: 2000-wide bin moments of z = -2 e^{i gamma U}/(rho - 1) at U_REF for zeros
with index > N_ARB (= N(8e7), our Arb zeros) and gamma <= TMAX. Also checks the overlap with our Arb zero batches.
usage: float_pass.py filelist TMAX nproc"""
import os, sys, glob, json
from multiprocessing import Pool
import numpy as np
import platt
HERE = os.path.dirname(os.path.abspath(__file__))
U_REF = np.longdouble("495.7028331345")
N_ARB = 195565267
ARB_DIR = os.environ.get("ARB_ZERO_DIR", os.path.join(HERE, "..", "..", "zero_cache_arb"))
DATA = os.environ.get("LMFDB_ZETA_DIR", os.path.join(HERE, "data"))
TWO_PI = np.longdouble(2)*np.arccos(np.longdouble(-1))

def overlap_check(d):
    """Compare LMFDB zeros with index <= N_ARB against our Arb balls (hi + lo +- rad)."""
    idx = d["idx"].astype(np.int64); sel = np.nonzero(idx <= N_ARB)[0]
    if not len(sel): return None
    g = platt.gamma_ld(d)[sel]; worst = 0.0
    for p0 in sorted(set(((idx[sel] - 13955182)//2000)*2000 + 13955182)):
        if p0 < 13955182: continue
        z = np.load(os.path.join(ARB_DIR, f"z_{p0:09d}.npz")); m = len(z["hi"])
        ii = np.nonzero((idx[sel] >= p0) & (idx[sel] < p0 + m))[0]
        ours = z["hi"].astype(np.longdouble)[idx[sel][ii] - p0] + z["lo"].astype(np.longdouble)[idx[sel][ii] - p0]
        dev = np.abs(g[ii] - ours) - z["rad"][idx[sel][ii] - p0]
        worst = max(worst, float(dev.max()))
    return worst

def work(fname):
    out = os.path.join(HERE, "bins", fname.replace(".dat", ".npz"))
    if os.path.exists(out): return out, None
    d = platt.decode(os.path.join(DATA, fname))
    ov = overlap_check(d)
    g = platt.gamma_ld(d); idx = d["idx"].astype(np.int64)
    keep = (idx > N_ARB) & (g <= TMAX)
    g = g[keep]
    edges = np.arange(0, 3.1e10, 2000.0)
    b = np.floor(g.astype(np.float64)/2000.0).astype(np.int64)
    ub = np.unique(b); gb = (ub + 0.5)*2000.0
    ph = np.fmod(g*U_REF, TWO_PI).astype(np.float64)
    gf = g.astype(np.float64); z = -2*np.exp(1j*ph)/(-0.5 + 1j*gf)
    pos = np.searchsorted(ub, b); dg = gf - gb[pos]
    M0 = np.zeros(len(ub), complex); M1 = np.zeros(len(ub), complex); M2 = np.zeros(len(ub), complex)
    np.add.at(M0, pos, z); np.add.at(M1, pos, z*dg); np.add.at(M2, pos, z*dg*dg)
    np.savez(out, gb=gb, M0=M0, M1=M1, M2=M2, n=len(g), idx_min=int(idx[keep].min()) if keep.any() else -1,
             idx_max=int(idx[keep].max()) if keep.any() else -1, gmax=float(g.max()) if len(g) else 0.0)
    return out, ov

if __name__ == "__main__":
    files = [l.strip() for l in open(sys.argv[1]) if l.strip()]
    TMAX = np.longdouble(sys.argv[2]); nproc = int(sys.argv[3])
    os.makedirs(os.path.join(HERE, "bins"), exist_ok=True)
    with Pool(nproc) as p:
        for out, ov in p.imap_unordered(work, files):
            if ov is not None:
                print(f"overlap check {os.path.basename(out)}: max(|LMFDB - Arb| - rad_Arb) = {ov:.3e} (must be <= 2^-102 = {2.0**-102:.1e})", flush=True)
    print("FLOAT PASS DONE", flush=True)
