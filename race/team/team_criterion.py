#!/usr/bin/env python3
"""Two-scale criterion for the residue-vs-nonresidue ("team") race modulo q in {8, 12}.

For q = 8 and q = 12 the squares in (Z/qZ)^x are {1}. Put G(x) := phi(q) pi(x;q,1) - sum_{(b,q)=1} pi(x;q,b),
so pi(x;q,1) exceeds the average of the other classes iff G(x) > 0. By orthogonality, the corrected explicit
formula (../THEOREM.md section 2) and p^2 = 1 (mod q) for p not dividing q,

  u e^{-u/2} G(e^u) >= -(phi-1) Phi_D(u) - sum_{chi in X_q} sum_rho u e^{-u/2} li(e^{rho u}) - eps0(u),

where X_q are the three real nonprincipal characters (coefficient 1 each); the term +(phi-1) omega(q)/2 is dropped
(lower bound). Integrating against K gives

  I(omega) >= L := -(phi-1) Phi_D(u_lo) - sum_chi S_chi(omega, T0) - R3 - R4 - R5 - eps0,

and L > 0 implies pi(x;q,1) > average of the other classes for some log x in [omega - eta0, omega + eta0].
R3, R4, R5 are the race bounds (race_criterion.py) with coefficient 1 instead of 2.
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "criterion"))
from flint import arb
import race_criterion as rc
from race_criterion import Phi_D, Phi2m, C_M, bmor_E, bmor_dE, C_const, require

# (q, n, parity delta, A_chi (Platt 2016, rounded down), n_imprimitive = #{p | q : p does not divide cond})
MODULI = {
    8:  dict(phi=4, chars=[(4, 3, 1, arb("2.5e7"), 0), (8, 3, 1, arb("1.25e7"), 0), (8, 5, 0, arb("1.25e7"), 0)]),
    12: dict(phi=4, chars=[(3, 2, 1, arb("3.3e7"), 1), (4, 3, 1, arb("2.5e7"), 1), (12, 11, 0, arb("8.3e6"), 0)]),
}
COEF = 1     # each character enters with coefficient 1


def R3_bound(u_lo, eta0, sums2, sums3):
    tot = arb(0)
    for s2, s3 in zip(sums2, sums3):
        tot += 2*s3 + eta0/2*s2
    return COEF*tot/(u_lo*u_lo)


def R4_bound(chars, omega, eta1, T0, lam, n_grid=6000):
    pi = arb.pi()
    require(eta1*T0 > 3*pi, "eta1 T0 must exceed 3 pi")
    H = lambda xi: sum(l*Phi2m(xi, m) for l, m in zip(lam, (1, 2, 3)) if l != 0)
    phi = lambda t: H(eta1*t)*(1/t + 1/(t*t*omega))
    tot = arb(0)
    for (q, n, d, A, _) in chars:
        require(T0 < A, f"T0 must be below the GRH height A_chi for {q}.{n}")
        r = (A/T0)**(arb(1)/n_grid); t = T0; integ = arb(0)
        for k in range(n_grid):
            t2 = t*r
            integ += phi(t)*((q*t2/(2*pi)).log()/(2*pi) + bmor_dE(t, q)/2)*(t2 - t)
            t = t2
        tot += 2*COEF*(integ + bmor_E(T0, q)*phi(T0))
    return tot.upper()


def R5_bound(chars, u_lo, u_hi, delta, M):
    pi = arb.pi(); N = 2*M; R = u_lo/2; tot = arb(0)
    for (q, n, d, A, _) in chars:
        rate = (M*pi/delta + N/R)/A
        require(rate < 1, "R5 rate must be < 1")
        tail = (N + 1)/pi*((q*A/(2*pi)).log()/N + 1/arb(N)**2)
        tot += COEF*(u_hi/2).exp()*4*C_M(M)*rate**N*tail
    return tot.upper()


def eps0_bound(chars, phi_q, u_lo, u_hi):
    tot = arb(0)
    for (q, n, d, A, nim) in chars:
        tot += COEF*((-u_lo).exp() + (1 - d)*u_hi.log() + C_const(q, n)) + COEF*nim*(u_hi.log() + arb("1.37"))
    # prime powers m >= 3: |-phi A1 + A_all| <= phi * 1.3 x^{1/3}/log x  (FH (2.12)), times u e^{-u/2}
    return (u_lo*(-u_lo/2).exp()*tot + phi_q*arb("1.3")*(-u_lo/6).exp()).upper()


def errors(q, omega, eta1, delta, M, T0, lam, sums2, sums3):
    cfg = MODULI[q]; chars = cfg["chars"]
    omega, eta1, delta, T0 = arb(omega), arb(eta1), arb(delta), arb(T0)
    eta0 = eta1 + delta; u_lo, u_hi = omega - eta0, omega + eta0
    # u_lo >= 45 suffices: Dusart (Thm 5.1) needs e^{u/2} >= 4e9 (u >= 44.22); Ford-Hudson (2.12) needs x >= e^30;
    # (2.10) needs x > 10; the R5 disk argument (radius u_lo/2, |arg z| <= pi/6) and the li expansions hold for any u_lo > 0.
    require(u_lo >= 45, "need u_lo >= 45")
    out = dict(bias=((cfg["phi"] - 1)*Phi_D(u_lo)).upper(), R3=R3_bound(u_lo, eta0, sums2, sums3).upper(),
               R4=R4_bound(chars, omega, eta1, T0, lam), R5=R5_bound(chars, u_lo, u_hi, delta, M),
               eps0=eps0_bound(chars, cfg["phi"], u_lo, u_hi))
    out["total"] = (out["bias"] + out["R3"] + out["R4"] + out["R5"] + out["eps0"]).upper()
    return out


def min_delta(chars, u_lo, u_hi, M, target):
    import math
    lo, hi = 1e-7, 1e-1
    for _ in range(60):
        mid = math.sqrt(lo*hi)
        try:
            ok = R5_bound(chars, arb(u_lo), arb(u_hi), arb(mid), M) <= target
        except ValueError:
            ok = False
        lo, hi = (lo, mid) if ok else (mid, hi)
    return hi
