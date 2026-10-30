#!/usr/bin/env python3
"""Floating-point zeros of L(1/2+it, chi) for any primitive chi (real or complex), 0 < t <= T (heuristic).

Generalises ../race/zeros/float_zeros.py: Z(t) = eps^{-1/2} exp(i theta(t)) L(1/2+it, chi) is real, where
eps = tau(chi)/(i^delta sqrt q) is the root number (Lambda(s, chi) = eps Lambda(1-s, conj chi)). For complex chi
the zeros with t < 0 are the conjugates of the zeros of L(s, conj chi) with t > 0.

The output only seeds arb_certify.py; nothing here is rigorous. Completeness is decided
later by the Turing check in turing.py.

L is
evaluated by Euler-Maclaurin on each Hurwitz zeta(s, a/q). The grid scan (step ZSTEP,
default 0.04) uses incremental phases. Sign changes are refined by Illinois iteration.
Every sampled local minimum of |Z| that is not next to a sign change is rescanned on a
grid 50 times finer, to catch close pairs that fall inside one step.

usage: float_zeros_any.py q n T nproc out.npy
"""
import math
import os
import sys
from multiprocessing import Pool

import numpy as np

# Bernoulli numbers B_2k / (2k)!
B2K = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330,
       854513/138, -236364091/2730, 8553103/6, -23749461029/870]
COEF = [B2K[k]/math.factorial(2*k + 2) for k in range(len(B2K))]
FINE = 50


_CV = {}


def char_vals(q, n):
    """(values chi(a) as complex array scaled by eps^{-1/2}, parity). The root-number phase is folded into
    the values so that the summation code below is unchanged."""
    if (q, n) in _CV:
        return _CV[(q, n)]
    from flint import dirichlet_char
    c = dirichlet_char(q, n)
    assert c.is_primitive(), "character must be primitive"
    v = np.zeros(q, dtype=complex)
    for a in range(q):
        z = c(a)
        v[a] = 0.0 if z == 0 else complex(float(z.real), float(z.imag))
    par = c.parity()
    tau = sum(v[a]*np.exp(2j*np.pi*a/q) for a in range(q))
    eps = tau/((1j)**par*math.sqrt(q))
    assert abs(abs(eps) - 1) < 1e-9, eps
    _CV[(q, n)] = (v/np.sqrt(eps), par)
    return _CV[(q, n)]


def lngamma_imag(z):
    """Im log Gamma(z) for complex array z with Re z > 0 (Stirling after shifting by 10)."""
    z = np.asarray(z, dtype=complex)
    shift = np.zeros_like(z)
    w = z.copy()
    for _ in range(10):
        shift += np.log(w)
        w = w + 1
    s = (w - 0.5)*np.log(w) - w + 0.5*np.log(2*np.pi)
    w2 = w*w
    s += 1/(12*w) - 1/(360*w*w2) + 1/(1260*w*w2*w2) - 1/(1680*w*w2**3)
    return (s - shift).imag


def theta(t, q, par):
    return 0.5*np.asarray(t)*np.log(q/np.pi) + lngamma_imag((0.5 + par)/2 + 0.5j*np.asarray(t))


def em_tail(s, a, q, M):
    """Euler-Maclaurin tail of sum_{m >= M} (m + a/q)^{-s}."""
    X = M + a/q
    XmS = np.exp(-s*math.log(X))
    part = XmS*X/(s - 1) + 0.5*XmS
    poch = s.copy()
    Xp = XmS/X
    for k in range(12):
        part += COEF[k]*poch*Xp
        poch = poch*(s + 2*k + 1)*(s + 2*k + 2)
        Xp = Xp/(X*X)
    return part


def Z_direct(t, q, chi, par, blk=512):
    """Z(t) for an array t (points of similar height), direct summation."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.empty(len(t))
    M = int(np.max(t)/np.pi) + 20
    m = np.arange(M, dtype=float)
    for i0 in range(0, len(t), blk):
        tt = t[i0:i0 + blk]
        s = 0.5 + 1j*tt
        tot = np.zeros(len(tt), dtype=complex)
        for a in range(1, q):
            if chi[a] == 0:
                continue
            lx = np.log(m + a/q)
            part = np.exp(-np.outer(s, lx)).sum(axis=1) + em_tail(s, a, q, M)
            tot += chi[a]*part
        out[i0:i0 + blk] = (np.exp(1j*theta(tt, q, par))*tot*np.exp(-s*math.log(q))).real
    return out


def illinois(lo, hi, flo, fhi, q, chi, par):
    for _ in range(60):
        if len(lo) == 0 or np.max(hi - lo) < 1e-12:
            break
        mid = hi - fhi*(hi - lo)/(fhi - flo)
        mid = np.where((mid <= lo) | (mid >= hi), 0.5*(lo + hi), mid)
        fm = Z_direct(mid, q, chi, par)
        left = np.sign(fm) == np.sign(flo)
        lo, hi, flo, fhi = (np.where(left, mid, lo), np.where(left, hi, mid),
                            np.where(left, fm, flo*0.5), np.where(left, fhi*0.5, fm))
    return 0.5*(lo + hi)


def scan_block(args):
    q, n, t0, t1, h = args
    chi, par = char_vals(q, n)
    nsteps = int(math.ceil((t1 - t0)/h))
    tg = t0 + h*np.arange(nsteps + 1)
    M = int(t1/np.pi) + 20
    m = np.arange(M, dtype=float)
    res = [(a, m + a/q) for a in range(1, q) if chi[a] != 0]
    logs = [np.log(x) for (_, x) in res]
    amp = [np.exp(-0.5*lx) for lx in logs]
    step = [np.exp(-1j*h*lx) for lx in logs]
    Zg = np.empty(nsteps + 1)
    ANCH = 400
    i = 0
    while i <= nsteps:
        j = min(nsteps, i + ANCH - 1)
        tt = tg[i:j + 1]
        s = 0.5 + 1j*tt
        tot = np.zeros(len(tt), dtype=complex)
        for k, (a, x) in enumerate(res):
            v = amp[k]*np.exp(-1j*tt[0]*logs[k])
            acc = np.empty(len(tt), dtype=complex)
            for r in range(len(tt)):
                acc[r] = v.sum()
                v *= step[k]
            tot += chi[a]*(acc + em_tail(s, a, q, M))
        Zg[i:j + 1] = (np.exp(1j*theta(tt, q, par))*tot*np.exp(-s*math.log(q))).real
        i = j + 1
    sc = np.sign(Zg[:-1])*np.sign(Zg[1:]) < 0
    idx = np.nonzero(sc)[0]
    zs = [illinois(tg[idx].copy(), tg[idx + 1].copy(), Zg[idx].copy(), Zg[idx + 1].copy(), q, chi, par)]
    # near misses: local minima of |Z| at samples with no sign change on either side
    a_ = np.abs(Zg)
    lm = np.nonzero((a_[1:-1] < a_[:-2]) & (a_[1:-1] < a_[2:]) & ~sc[:-1] & ~sc[1:])[0] + 1
    nfound = 0
    for c0 in range(0, len(lm), 40):
        cl = lm[c0:c0 + 40]
        fine = (tg[cl - 1][:, None] + (2*h/FINE)*np.arange(FINE + 1)[None, :]).ravel()
        zf = Z_direct(fine, q, chi, par).reshape(len(cl), FINE + 1)
        tf = fine.reshape(len(cl), FINE + 1)
        r, c = np.nonzero(np.sign(zf[:, :-1])*np.sign(zf[:, 1:]) < 0)
        if len(r):
            nfound += len(r)
            zs.append(illinois(tf[r, c].copy(), tf[r, c + 1].copy(), zf[r, c].copy(), zf[r, c + 1].copy(), q, chi, par))
    z = np.sort(np.concatenate(zs))
    z = z[(z > t0) & (z <= t1)]
    return z, len(lm), nfound


def main():
    q, n, T, nproc, out = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    h = float(os.environ.get("ZSTEP", "0.04"))
    B = 8*nproc
    edges = [0.5] + [T*math.sqrt(k/B) for k in range(1, B + 1)]   # equal work (cost ~ t)
    jobs = [(q, n, edges[k], edges[k + 1], h) for k in range(B)]
    with Pool(nproc) as p:
        res = p.map(scan_block, jobs)
    z = np.unique(np.concatenate([r[0] for r in res]))
    chi, par = char_vals(q, n)
    np.save(out, z)
    print(f"q={q} n={n} par={par} zeros={len(z)} min gap={np.diff(z).min():.4g} "
          f"near-miss windows={sum(r[1] for r in res)} extra zeros from them={sum(r[2] for r in res)}")
    for TT in np.linspace(T/10, T, 10):
        nf = int(np.sum(z <= TT))
        th = float(theta(np.array([TT]), q, par)[0]/np.pi)
        print(f"  T={TT:.0f} found={nf} theta/pi={th:.2f} diff={nf - th:+.2f}")


if __name__ == "__main__":
    main()
