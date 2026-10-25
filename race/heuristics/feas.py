"""Heuristic certifiability of a race region with the two-scale kernel (and the Gaussian kernel for comparison).

Two-scale kernel K = K_core * psi_delta, K_core = kappa_2m(./eta1)/eta1, psi = kappa_2M(./delta)/delta,
hat K(g) = khat_m(eta1 g) khat_M(delta g),  khat_m(x) = sin x / x prod_{k<=m} (k pi)^2/((k pi)^2 - x^2).

Zeros with |gamma| > A (possibly off the line) are bounded by 2M-fold integration by parts; we require
  2M log(A delta/(M pi)) >= omega/2 + L,   L = 30 (a crude stand-in for log(1/eps) + log factors).
Optimum: A delta/(M pi) = e, 2M = omega/2 + L  =>  delta_min = e pi (omega/2 + L)/(2A).
Gaussian (Lehman / Ford-Hudson): A^2/(2 alpha) >= omega/2 + L and eta >= 2A/alpha  =>  eta >= 4(omega/2+L)/A,
and the computed zeros must reach T0 ~ sqrt(2 alpha * 20).
Platt (2016) GRH heights: t0(q) = max(1e8/q, 7.5e7/q + 200) (even q), max(1e8/q, 3.75e7/q + 200) (odd q).
"""
import math
import numpy as np
import race

L_SAFETY = 30.0
COND = {'-3': 3, '-4': 4, '-8': 8, '8': 8, '12': 12, '-24': 24, '24': 24}


def platt_height(q):
    return max(1e8/q, (7.5e7 if q % 2 == 0 else 3.75e7)/q + 200)


def race_A(q, a, b):
    labs = [lab for lab in race.MODCHARS[q] if race.chival(lab, a) != race.chival(lab, b)]
    return min(platt_height(COND[l]) for l in labs), labs


def khat(x, m):
    x = np.asarray(x, dtype=float)
    out = np.where(np.abs(x) < 1e-12, 1.0, np.sin(x)/np.where(x == 0, 1, x))
    for k in range(1, m+1):
        d = (k*np.pi)**2 - x**2
        d = np.where(np.abs(d) < 1e-9, 1e-9, d)
        out = out*(k*np.pi)**2/d
    return out


def two_scale_params(omega, A):
    M = max(1, int(math.ceil((omega/2 + L_SAFETY)/2)))
    delta = math.e*math.pi*M/A
    return M, delta


def R4_bound(labs_c, eta1, m, T0, A, M, delta):
    """sum_chi |c_chi| * 2 * int_{T0}^{A} |hatK(t)| (1/t) dN_chi(t), dN = log(q t/2pi)/(2 pi) dt (smooth density)."""
    t = np.geomspace(T0, A, 200000)
    w = np.abs(khat(eta1*t, m))*np.abs(khat(delta*t, M))
    tot = 0.0
    for lab, c in labs_c:
        dens = np.log(COND[lab]*t/(2*np.pi))/(2*np.pi)
        f = w*2*abs(c)/t*dens
        tot += np.trapezoid(f, t) if hasattr(np, 'trapezoid') else np.trapz(f, t)
    return tot


def gaussian_T0(omega, eta):
    """Smallest zero height for a Gaussian kernel with window eta at omega (A chosen minimal)."""
    A = 4*(omega/2 + L_SAFETY)/eta
    alpha = 2*A/eta
    return math.sqrt(40*alpha), A
