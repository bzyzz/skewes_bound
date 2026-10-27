#!/usr/bin/env python3
"""Direct (floating-point) evaluation of the sign of Delta_K(x) = prod_{N p <= x} (1 - 1/N p)^{-1} - e^gamma kappa_K log x
for K = Q, Q(sqrt 5), Q(sqrt 13), x up to X. Reports the first x with Delta_K(x) < 0 and the minimum of the normalized error."""
import math, sys
import numpy as np
from flint import acb, dirichlet_char
X = int(float(sys.argv[1]))
sieve = np.ones(X + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(X**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
P = np.nonzero(sieve)[0].astype(np.int64)
EG = 0.5772156649015329
for name, q, n in [("Q", 1, 0), ("Q(sqrt 5)", 5, 4), ("Q(sqrt 13)", 13, 12)]:
    if q == 1:
        norms, mult = P.astype(np.float64), np.ones(len(P)); kappa = 1.0
    else:
        c = dirichlet_char(q, n); kappa = float(acb.dirichlet_l(acb(1), c).real.mid())
        chi = np.array([int(round(float(c(a).real))) for a in range(q)])[P % q]
        split, inert, ram = P[chi == 1], P[chi == -1], P[chi == 0]
        inert = inert[inert.astype(np.float64)**2 <= X]
        norms = np.concatenate([split, inert.astype(np.float64)**2, ram]).astype(np.float64)
        mult = np.concatenate([2*np.ones(len(split)), np.ones(len(inert)), np.ones(len(ram))])
    o = np.argsort(norms, kind="stable"); norms, mult = norms[o], mult[o]
    terms = -mult*np.log1p(-1/norms)
    S = np.cumsum(terms)                     # log prod over N p <= norms[k]
    nxt = np.append(norms[1:], X)            # E3 is decreasing until the next norm: test just below it
    E3 = S - EG - math.log(kappa) - np.log(np.log(nxt))
    neg = np.nonzero(E3 < 0)[0]
    norm_err = E3*np.sqrt(nxt)*np.log(nxt)
    k = norm_err[norms > 1e3].argmin() + np.searchsorted(norms, 1e3, side="right")
    print(f"{name:11s} kappa={kappa:.6f}: " + (f"FIRST NEGATIVE just below x = {nxt[neg[0]]:.0f} (E3 = {E3[neg[0]]:.3e}); {len(neg)} negative intervals up to {X:.0e}"
          if len(neg) else f"no negative value up to {X:.0e}") + f"; min sqrt(x)log(x)*E3 = {norm_err[k]:+.4f} at x~{nxt[k]:.4g}", flush=True)
