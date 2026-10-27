#!/usr/bin/env python3
"""Survey: first integer x <= X with Delta_K(x) < 0 for quadratic fields K = Q(sqrt D), D fundamental, |D| <= Dmax.
chi_D(p) from the periodic Kronecker symbol; kappa_K = L(1, chi_D)*(residue of zeta) = L(1, chi_D)."""
import math, sys
import numpy as np
from flint import arb
X = int(float(sys.argv[1])); Dmax = int(sys.argv[2]); EG = float(arb.const_euler().mid())
sieve = np.ones(X + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(X**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
P = np.nonzero(sieve)[0].astype(np.int64); Pf = P.astype(np.float64)

def fundamental(D):
    if D in (0, 1): return False
    if D % 4 == 1:
        m = abs(D)
    elif D % 4 == 0:
        m = abs(D)//4
        if (D//4) % 4 not in (2, 3): return False
    else:
        return False
    return all(m % (p*p) for p in range(2, int(math.isqrt(m)) + 1))

def kron(D, n):
    """Kronecker symbol (D/n) for n >= 1."""
    r = 1
    while n % 2 == 0:
        n //= 2
        if D % 2 == 0: return 0
        r *= 1 if D % 8 in (1, 7) else -1
    a, b = D % n if n > 1 else 0, n
    if n == 1: return r
    # Jacobi (a/b)
    j = 1
    while a:
        while a % 2 == 0:
            a //= 2
            if b % 8 in (3, 5): j = -j
        a, b = b, a
        if a % 4 == 3 and b % 4 == 3: j = -j
        a %= b
    return r*j if b == 1 else 0

out = []
for D in [d for d in range(-Dmax, Dmax + 1) if fundamental(d)]:
    m = abs(D)
    tab = np.array([kron(D, a) if a else 0 for a in range(m)])
    L1 = float((-sum((int(tab[a])*(arb(a)/m).digamma() for a in range(1, m)), arb(0))/m).mid())
    chi = tab[P % m]
    split, inert, ram = Pf[chi == 1], Pf[chi == -1], Pf[chi == 0]
    inert = inert[inert**2 <= X]
    norms = np.concatenate([split, inert**2, ram]); mult = np.concatenate([2*np.ones(len(split)), np.ones(len(inert)), np.ones(len(ram))])
    o = np.argsort(norms, kind="stable"); norms, mult = norms[o], mult[o]
    S = np.cumsum(-mult*np.log1p(-1/norms))
    thr = np.exp(np.exp(S - EG - math.log(L1)))           # within [norm_k, norm_{k+1}), Delta < 0 iff x > thr
    nxt = np.append(norms[1:], X + 1.0)
    negk = np.nonzero(thr < nxt - 1e-9)[0]
    # first integer x >= max(norm_k, ceil(thr)) and < next norm
    first = None
    for k in negk:
        xi = max(math.ceil(thr[k]), int(norms[k]))
        if xi < nxt[k]: first = xi; break
    lo = np.maximum(thr, norms); dens = np.sum(np.where(nxt > lo, np.log(np.minimum(nxt, X)) - np.log(lo), 0))/math.log(X/2)
    out.append((D, first, dens, L1))
    print(f"D={D:5d} kappa={L1:.4f} first negative integer x: {first if first else '> %.0e' % X:>10}  negativity log-density {dens:.4f}", flush=True)
none = [r for r in out if r[1] is None]
print(f"\n{len(out)} fields; no negative value up to {X:.0e}: {[r[0] for r in none]}")
big = sorted([r for r in out if r[1]], key=lambda r: -r[1])[:10]
print("latest first sign changes:", [(r[0], r[1]) for r in big])
