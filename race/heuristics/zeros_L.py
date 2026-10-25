#!/usr/bin/env python3
"""Floating-point zeros of L(s, chi) for real primitive Dirichlet characters (heuristic only).

Z(t) = exp(i theta(t)) L(1/2+it, chi) is real (root number 1 for real primitive chi).
L is evaluated by Euler-Maclaurin applied to each Hurwitz zeta(s, a/q).
Grid scan with incremental phase updates, sign changes refined by Illinois iteration.

usage: zeros_L.py q n T nproc out.npy
"""
import sys, math, os
import numpy as np
from multiprocessing import Pool

# Bernoulli numbers B_2k / (2k)!
B2K = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510, 43867/798, -174611/330,
       854513/138, -236364091/2730, 8553103/6, -23749461029/870]
FACT = [math.factorial(2*k) for k in range(1, len(B2K)+1)]
COEF = [B2K[k]/FACT[k] for k in range(len(B2K))]


def char_vals(q, n):
    from flint import dirichlet_char
    c = dirichlet_char(q, n)
    v = np.zeros(q)
    for a in range(q):
        z = c(a)
        v[a] = 0.0 if z == 0 else round(float(z.real))
    return v, c.parity()


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
    return 0.5*t*np.log(q/np.pi) + lngamma_imag((0.5+par)/2 + 0.5j*np.asarray(t))


def L_direct(t, q, chi, M=None, K=12):
    """L(1/2+it) for array t (all in a narrow range), Euler-Maclaurin per residue."""
    t = np.atleast_1d(np.asarray(t, dtype=float))
    s = 0.5 + 1j*t
    tm = float(np.max(np.abs(t)))
    if M is None:
        M = int(tm/np.pi) + 20
    out = np.zeros(t.shape, dtype=complex)
    m = np.arange(M, dtype=float)
    for a in range(1, q):
        if chi[a] == 0:
            continue
        al = a/q
        x = m + al                        # (m+alpha), m < M
        lx = np.log(x)
        part = np.exp(-np.outer(s, lx)).sum(axis=1)
        X = M + al
        lX = math.log(X)
        XmS = np.exp(-s*lX)
        part += XmS*X/(s-1) + 0.5*XmS
        # EM corrections: sum_k B2k/(2k)! * (s)_{2k-1} X^{-s-2k+1}
        poch = s.copy()                   # (s)_1
        Xp = XmS/X                        # X^{-s-1}
        for k in range(K):
            part += COEF[k]*poch*Xp
            poch = poch*(s+2*k+1)*(s+2*k+2)
            Xp = Xp/(X*X)
        out += chi[a]*part
    return out*np.exp(-s*math.log(q))


def Z_direct(t, q, chi, par):
    t = np.atleast_1d(np.asarray(t, dtype=float))
    return (np.exp(1j*theta(t, q, par))*L_direct(t, q, chi)).real


def scan_block(args):
    q, n, t0, t1, h = args
    chi, par = char_vals(q, n)
    # grid t0..t1 inclusive, incremental phases, re-anchored every 400 steps
    nsteps = int(math.ceil((t1-t0)/h))
    tg = t0 + h*np.arange(nsteps+1)
    M = int(t1/np.pi) + 20
    m = np.arange(M, dtype=float)
    res = []
    for a in range(1, q):
        if chi[a] == 0:
            continue
        res.append((a, m + a/q))
    Zg = np.empty(nsteps+1)
    # precompute static pieces
    logs = [np.log(x) for (_, x) in res]
    amp = [np.exp(-0.5*lx) for lx in logs]
    step = [np.exp(-1j*h*lx) for lx in logs]
    ANCH = 400
    i = 0
    while i <= nsteps:
        j = min(nsteps, i+ANCH-1)
        tt = tg[i:j+1]
        s = 0.5 + 1j*tt
        tot = np.zeros(len(tt), dtype=complex)
        for k, (a, x) in enumerate(res):
            v = amp[k]*np.exp(-1j*tt[0]*logs[k])
            acc = np.empty(len(tt), dtype=complex)
            for r in range(len(tt)):
                acc[r] = v.sum()
                v *= step[k]
            al = a/q
            X = M + al
            lX = math.log(X)
            XmS = np.exp(-s*lX)
            part = acc + XmS*X/(s-1) + 0.5*XmS
            poch = s.copy()
            Xp = XmS/X
            for kk in range(12):
                part += COEF[kk]*poch*Xp
                poch = poch*(s+2*kk+1)*(s+2*kk+2)
                Xp = Xp/(X*X)
            tot += chi[a]*part
        Lv = tot*np.exp(-s*math.log(q))
        Zc = np.exp(1j*theta(tt, q, par))*Lv
        Zg[i:j+1] = Zc.real
        i = j+1
    # sign changes
    idx = np.nonzero(np.sign(Zg[:-1])*np.sign(Zg[1:]) < 0)[0]
    lo = tg[idx].copy(); hi = tg[idx+1].copy()
    flo = Zg[idx].copy(); fhi = Zg[idx+1].copy()
    # Illinois refinement, vectorised
    for it in range(40):
        if len(lo) == 0:
            break
        mid = hi - fhi*(hi-lo)/(fhi-flo)
        fm = Z_direct(mid, q, chi, par)
        left = np.sign(fm) == np.sign(flo)
        # replace lo where same sign as lo, else hi
        lo_new = np.where(left, mid, lo); flo_new = np.where(left, fm, flo*0.5)
        hi_new = np.where(left, hi, mid); fhi_new = np.where(left, fhi*0.5, fm)
        # Illinois: halve the retained endpoint's f
        flo = np.where(left, fm, flo*0.5); fhi = np.where(left, fhi*0.5, fm)
        lo, hi = lo_new, hi_new
        if np.max(hi-lo) < 1e-11:
            break
    zeros = 0.5*(lo+hi)
    zeros = zeros[(zeros > t0) & (zeros <= t1)]
    return zeros, float(np.max(np.abs(Zg)))


def main():
    q, n, T, nproc, out = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    h = float(os.environ.get("ZSTEP", "0.04"))
    # blocks with roughly equal work (cost ~ t * len)
    edges = [0.5]
    B = 4*nproc
    for k in range(1, B+1):
        edges.append(T*math.sqrt(k/B))
    jobs = [(q, n, edges[k], edges[k+1], h) for k in range(B)]
    with Pool(nproc) as p:
        res = p.map(scan_block, jobs)
    z = np.sort(np.concatenate([r[0] for r in res]))
    chi, par = char_vals(q, n)
    # count check against theta/pi (S(T) small on average)
    chk = []
    for TT in np.linspace(T/10, T, 10):
        nf = int(np.sum(z <= TT))
        chk.append((TT, nf, float(theta(np.array([TT]), q, par)[0]/np.pi)))
    np.save(out, z)
    gaps = np.diff(z)
    print(f"q={q} n={n} par={par} zeros={len(z)} first={z[:4]} min gap={gaps.min():.4g}")
    for TT, nf, th in chk:
        print(f"  T={TT:.0f} found={nf} theta/pi={th:.2f} diff={nf-th:+.2f}")


if __name__ == "__main__":
    main()
