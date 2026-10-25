#!/usr/bin/env python3
"""Two-scale criterion for a negative value of Delta(x) = pi(x;24,11) - pi(x;24,1).

Statement (proof in THEOREM.md). Let X = {3.2, 4.3, 8.5, 24.11} (primitive real characters with chi(11) = -1).
Let K = K_core * psi_delta (the kernel of ../THEOREM.md, section 1), with supp K = [-eta0, eta0] and eta0 = eta1 + delta.
Put u_lo = omega - eta0 and u_hi = omega + eta0. Then

  I(omega) := int K(u - omega) 8 u e^{-u/2} Delta(e^u) du  <=  U(omega),
  U = 8 Phi_D(u_lo) + 2 sum_{chi in X} S_chi(omega, T0) + R3 + R4 + R5 + eps0,
  S_chi = sum_{0<|gamma|<=T0} hatK(gamma) e^{i gamma omega} (1/rho + 1/(rho^2 omega))   (real).

If U < 0, then Delta(e^u) < 0 for some u in [u_lo, u_hi].

Inputs:
- Ford-Hudson (2001) (2.10)-(2.13);
- Dusart (2018) Thm 5.1;
- BMOR (2021) Thm 1.1 for N(T, chi);
- Platt (2016) GRH heights A_chi;
- our Turing-certified zeros up to T0.
"""
import math

import numpy as np
import flint
from flint import acb, arb, dirichlet_char

flint.ctx.dps = 40

# (q, n, parity delta, A_chi, n_imprimitive): A_chi is rounded down from Platt (2016), max(1e8/q, ...);
# n_imprimitive = #{p | 24 : p does not divide q}.
CHARS = [(3, 2, 1, arb("3.3e7"), 1), (4, 3, 1, arb("2.5e7"), 1), (8, 5, 0, arb("1.25e7"), 1),
         (24, 11, 0, arb("4.1e6"), 0)]


# ---------------- kernel ----------------
def kap_hat(xi, m):
    """hat kappa_{2m}(xi) in Arb, xi a real ball (removable singularities handled by Arb sinc)."""
    pi = arb.pi()
    r = xi.sinc()
    for k in range(1, m + 1):
        kp2 = (k*pi)**2
        r = r*kp2/(kp2 - xi*xi)
    return r


def Phi2m(xi, m):
    """Nonincreasing majorant of |hat kappa_{2m}(xi)| for xi > m pi."""
    pi = arb.pi()
    r = 1/xi
    for k in range(1, m + 1):
        kp2 = (k*pi)**2
        r = r*kp2/(xi*xi - kp2)
    return r


def C_M(M):
    return arb(4**M*math.factorial(M)**2)/arb(2*math.factorial(2*M))


# ---------------- error terms ----------------
def Phi_D(u):
    return 1 + 2/u + 8/(u*u) + arb("58.56")/(u*u*u)


def bmor_E(t, q):
    """BMOR Thm 1.1: |N(t,chi) - (t/pi) log(qt/(2 pi e)) + chi(-1)/4| <= E(t), for |gamma| <= t counts."""
    ell = (q*(t + 2)/(2*arb.pi())).log()
    assert ell > arb("1.567")
    return arb("0.22737")*ell + 2*(1 + ell).log() - arb("0.5")


def bmor_dE(t, q):
    ell = (q*(t + 2)/(2*arb.pi())).log()
    return (arb("0.22737") + 2/(1 + ell))/(t + 2)


def sum_inv_pow_tail(T0, q, n, NT0):
    """Upper bound for sum over zeros with |gamma| > T0 of |gamma|^{-n} (both signs), given N(T0) = NT0 exactly.

    Stieltjes: sum = int_{T0}^inf t^{-n} dN <= T0^{-n} (Mn(T0) + E(T0) - NT0) + int_{T0}^inf t^{-n} (Mn' + E') dt,
    with Mn(t) = (t/pi) log(qt/(2 pi e)) - chi(-1)/4 (we bound chi(-1)/4 by 1/4), Mn' = (1/pi) log(qt/(2 pi)),
    and E' <= 1/t (for ell > 1.567).
    """
    pi = arb.pi()
    L = (q*T0/(2*pi)).log()
    Mn = T0/pi*(L - 1) + arb(1)/4
    first = T0**(-n)*(Mn + bmor_E(T0, q) - NT0)
    # int t^{-n} (1/pi) log(qt/2pi) dt = T0^{1-n} (L/(n-1) + 1/(n-1)^2)/pi ; int t^{-n-1} dt = T0^{-n}/n
    integ = T0**(1 - n)*(L/(n - 1) + 1/arb(n - 1)**2)/pi + T0**(-n)/n
    return first + integ


def R3_bound(u_lo, eta0, sums2, sums3):
    """(2/u_lo^2) sum_chi [2 sum |gamma|^{-3} + (eta0/2) sum gamma^{-2}], sums over all zeros of each L."""
    tot = arb(0)
    for s2, s3 in zip(sums2, sums3):
        tot += 2*s3 + eta0/2*s2
    return 2*tot/(u_lo*u_lo)


def R4_bound(omega, eta1, T0, lam, n_grid=6000):
    """4 sum_chi [int_{T0}^{A} phi (1/2pi) log(qt/2pi) dt + 2 e(T0) phi(T0) + int_{T0}^{A} phi e' dt],
    e = E/2, phi(t) = H(eta1 t)(1/t + 1/(t^2 omega)), H = sum lam_m Phi_{2m}: nonincreasing for eta1 t > 3 pi.
    Upper Riemann sums with left endpoints on a geometric grid (both integrands are nonincreasing)."""
    pi = arb.pi()
    require(eta1*T0 > 3*pi, "eta1 T0 must exceed 3 pi")
    H = lambda xi: sum(l*Phi2m(xi, m) for l, m in zip(lam, (1, 2, 3)) if l != 0)
    phi = lambda t: H(eta1*t)*(1/t + 1/(t*t*omega))
    tot = arb(0)
    for (q, n, d, A, _) in CHARS:
        require(T0 < A, f"T0 must be below the GRH height A_chi for {q}.{n}")
        r = (A/T0)**(arb(1)/n_grid)
        t = T0
        integ = arb(0)
        for k in range(n_grid):
            t2 = t*r
            f = phi(t)*((q*t2/(2*pi)).log()/(2*pi) + bmor_dE(t, q)/2)   # log factor at the right end: upper bound
            integ += f*(t2 - t)
            t = t2
        tot += 4*(integ + bmor_E(T0, q)*phi(T0))
    return tot.upper()


def R5_bound(u_lo, u_hi, delta, M):
    """2 sum_chi e^{u_hi/2} 4 C_M (M pi/delta + 2M/R)^{2M} sum_{|gamma|>A_chi} |gamma|^{-2M-1}, with R = u_lo/2,
    sum_{|gamma|>A} |gamma|^{-N-1} <= ((N+1)/pi) A^{-N} (log(qA/2pi)/N + 1/N^2)."""
    pi = arb.pi()
    N = 2*M
    R = u_lo/2
    tot = arb(0)
    for (q, n, d, A, _) in CHARS:
        rate = (M*pi/delta + N/R)/A
        require(rate < 1, "R5 rate must be < 1")
        tail = (N + 1)/pi*((q*A/(2*pi)).log()/N + 1/arb(N)**2)
        tot += 2*(u_hi/2).exp()*4*C_M(M)*rate**N*tail
    return tot.upper()


def C_const(q, n):
    """|C_chi| in the explicit formula (THEOREM.md §2, corrected FH (2.10)) for real primitive chi:
    C_chi = log L(0, chi) (odd), log L'(0, chi) - gamma_E (even). Arb Taylor series of L at s = 0."""
    from flint import acb_series
    chi = dirichlet_char(q, n)
    c = acb_series.dirichlet_l(acb_series([0, 1], prec=2), chi).coeffs()
    if int(chi.parity()) == 1:
        require(c[0].real > 0, "L(0, chi) > 0")
        return abs(c[0].real.log())
    require(c[1].real > 0, "L'(0, chi) > 0")
    return abs(c[1].real.log() - arb.const_euler())


def eps0_bound(u_lo, u_hi):
    tot = arb(0)
    for (q, n, d, A, nim) in CHARS:
        # imprimitive correction (corrected FH (2.13)): sum_{m <= log x/log p} 1/m <= 1 + loglog x - loglog p
        # <= loglog x + 1.37 for p = 2, 3
        tot += 2*((-u_lo).exp() + (1 - d)*u_hi.log() + C_const(q, n)) + 2*nim*(u_hi.log() + arb("1.37"))
    return (u_lo*(-u_lo/2).exp()*tot + arb("10.4")*(-u_lo/6).exp()).upper()


def require(cond, msg):
    if not cond:
        raise ValueError(msg)


def errors(omega, eta1, delta, M, T0, lam, sums2, sums3):
    omega, eta1, delta, T0 = arb(omega), arb(eta1), arb(delta), arb(T0)
    eta0 = eta1 + delta
    u_lo, u_hi = omega - eta0, omega + eta0
    require(u_lo > 100, "need u_lo > 100")
    out = dict(main=(8*Phi_D(u_lo)).upper(), R3=R3_bound(u_lo, eta0, sums2, sums3).upper(),
               R4=R4_bound(omega, eta1, T0, lam), R5=R5_bound(u_lo, u_hi, delta, M), eps0=eps0_bound(u_lo, u_hi))
    out["total"] = (out["main"] + out["R3"] + out["R4"] + out["R5"] + out["eps0"]).upper()
    return out


# ---------------- float helpers for the parameter search ----------------
def khat_f(x, m):
    x = np.asarray(x, dtype=float)
    for k in range(1, m + 1):     # step off the removable singularities at k pi
        x = np.where(np.abs(np.abs(x) - k*np.pi) < 1e-7, x + 2e-7, x)
    out = np.sinc(x/np.pi)
    for k in range(1, m + 1):
        out = out*(k*np.pi)**2/((k*np.pi)**2 - x*x)
    return out


def Khat_f(g, eta1, delta, M, lam):
    core = sum(l*khat_f(eta1*g, m) for l, m in zip(lam, (1, 2, 3)) if l != 0)
    return core*khat_f(delta*g, M)


def min_delta(u_lo, u_hi, M, target):
    """Smallest delta (float search, then checked in Arb) with R5 <= target for mollifier order M."""
    lo, hi = 1e-7, 1e-2
    for _ in range(60):
        mid = math.sqrt(lo*hi)
        try:
            ok = R5_bound(arb(u_lo), arb(u_hi), arb(mid), M) <= target
        except ValueError:
            ok = False
        lo, hi = (lo, mid) if ok else (mid, hi)
    return hi
