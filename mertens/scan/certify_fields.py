#!/usr/bin/env python3
"""Rigorous (Arb) first integer x with Delta_K(x) < 0 for K = Q(sqrt D), D fundamental, |D| <= Dmax.
Delta_K(n) = prod_{N p <= n} (1 - 1/N p)^{-1} - e^gamma kappa_K log n, kappa_K = L(1, chi_D) (Arb, via digamma).
Checks Delta_K(n) > 0 for 2 <= n < x0 and Delta_K(x0) < 0 with certified signs; writes fields_first_negative.json."""
import json, math, sys
import flint
from flint import arb
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

flint.ctx.dps = 40
Dmax = int(sys.argv[1]); NMAX = 400
primes = [p for p in range(2, NMAX + 1) if all(p % q for q in range(2, int(p**0.5) + 1))]
eg = arb.const_euler().exp()
out = {}
for D in [d for d in range(-Dmax, Dmax + 1) if fundamental(d)]:
    m = abs(D)
    kappa = -sum((kron(D, a)*(arb(a)/m).digamma() for a in range(1, m) if kron(D, a)), arb(0))/m
    norms = []
    for p in primes:
        c = kron(D, p)
        if c == 1: norms += [p, p]
        elif c == 0: norms.append(p)
        elif p*p <= NMAX: norms.append(p*p)
    first = None
    for n in range(2, NMAX + 1):
        prod = arb(1)
        for N in norms:
            if N <= n: prod *= 1/(1 - arb(1)/N)
        val = prod - eg*kappa*arb(n).log()
        if val < 0: first = n; break
        assert val > 0, f"sign of Delta_K({n}) not certified for D={D}"
    assert first is not None, f"no negative value up to {NMAX} for D={D}"
    out[D] = first
json.dump(out, open("fields_first_negative.json", "w"), indent=0)
vals = sorted(out.items(), key=lambda kv: -kv[1])
print(f"{len(out)} fields, all certified in Arb; max first negative x = {vals[0][1]} (D = {[d for d, v in vals if v == vals[0][1]]})")
print("D=5:", out[5], " D=13:", out[13], " D=-4:", out[-4], " D=-3:", out[-3], " D=8:", out[8])
