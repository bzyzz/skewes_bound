#!/usr/bin/env python3
"""Certify (Arb, adaptive bisection) that Axler's bound implies the Dusart-form bound used in Phi(u).

Axler, Integers 18 (2018) #A52, Theorem 2 (open access): for every x >= 49,
  pi(x) < x / (log x - 1 - 1/log x - 3.15/log^2 x - 12.85/log^3 x - 71.3/log^4 x - 463.2275/log^5 x - 4585/log^6 x).
With L = log y, this is pi(y) < (y/L) A(L), A(L) = 1/(1 - 1/L - 1/L^2 - 3.15/L^3 - ... - 4585/L^7).
Dusart (2018, Thm 5.1) gives pi(y) < (y/L) D(L), D(L) = 1 + 1/L + 2/L^2 + 7.32/L^3 for y >= 4e9.
We certify A(L) < D(L) for 30 <= L <= 400, i.e. for u = 2L in [60, 800]: on that range every Phi(u) upper bound
used in this repository (Skewes u ~ 728, race 24;11 u ~ 111.5, team races u ~ 64.8, 158.6, 339.6) also follows
from Axler's theorem, independently of Dusart's constant 7.32.
"""
import sys
from flint import arb
sys.setrecursionlimit(10000)

def A(L):
    return 1/(1 - 1/L - 1/L**2 - arb("3.15")/L**3 - arb("12.85")/L**4 - arb("71.3")/L**5
              - arb("463.2275")/L**6 - arb(4585)/L**7)

def D(L):
    return 1 + 1/L + 2/L**2 + arb("7.32")/L**3

def cert(lo, hi, depth=0):
    Lb = arb(((lo + hi)/2).mid(), ((hi - lo)/2).mid()*1.000001)
    if (D(Lb) - A(Lb)) > 0:
        return True
    if depth > 40:
        return False
    mid = (lo + hi)/2
    return cert(lo, mid, depth + 1) and cert(mid, hi, depth + 1)

if __name__ == "__main__":
    ok = all(cert(arb(30) + k, arb(31) + k) for k in range(370))
    print("A(L) < D(L) for all L in [30, 400]:", ok)
    sys.exit(0 if ok else 1)
