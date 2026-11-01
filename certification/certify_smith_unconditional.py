#!/usr/bin/env python3
"""
Rigorous Arb ball-arithmetic verification of the unconditional refinement of C. Smith (2016):

Using:
  1) The asymmetric smooth-mollified Gaussian kernel
       K(y) := G_alpha(y) * w_{eta_L, eta_R, delta}(y),
       w_{eta_L, eta_R, delta}(y) := (1_{[-eta_L, eta_R]} * psi_delta)(y),
     supported on [-eta_{L,0}, eta_{R,0}] = [-(eta_L + delta), eta_R + delta],
     which equals G_alpha(y) identically on [-(eta_L - delta), eta_R - delta] and vanishes
     to order 2M = 380 at the endpoints -(eta_L + delta) and +(eta_R + delta);
  2) The unconditional explicit-formula lower bound for -(1/2)*pi(x^{1/2}) via Pi(x^{1/2}),
     together with Platt-Trudgian (2021) (RH verified up to A = 3*10^12);
  3) Our one-sided complex integration-by-parts bound on the Gaussian spatial tail:
       |int_{eta}^\infty G_alpha(x) e^{i gamma x} dx| <= (1 + 1/(alpha*eta^2)) * G_alpha(eta) / sqrt((alpha*eta)^2 + gamma^2),
     combined with Fubini's identity 1 - w(y) = int_{-delta}^{delta} psi_delta(t) 1_{[eta_R + t, infty)}(y) dt;
  4) Exact binomial Fourier moments mu_{M, j-1} for ||psi_delta^{(j-1)}||_{L^infty} in the 380-fold
     integration-by-parts bound on off-line zeros |gamma| > A = 3*10^12;
  5) David Platt's Arb-computed value of the Gaussian zero sum at omega = 727.951335401,
     alpha = 4.5*10^19, T = 30,610,046,000 (N(T) = 103,800,788,359 zeros) from Table 4 (p. 54)
     and Section 8 (p. 46, |delta S| < 6*10^{-26}) of C. Smith (2016);

Repairs after review (2026-10-30; derivations in the paper, Section sec:smith, and paper/proof_audit.md):
  - C_0 is enclosed by a direct Arb evaluation of w e^{-w} Ei(w) over the whole window (the asymptotic remainder
    bound previously used is false at these parameters);
  - the half-frequency term li(e^{rho u/2}) costs a factor 2^{2M} in the integration by parts for |gamma| > A;
  - the two transition strips are bounded separately in the off-line bound;
  - the prime-power, secondary-term, spectral-tail and u-variation errors are derived directly
    (Backlund bound, absolute first moment) instead of transferring Smith's Gaussian error terms;
  - in the RH variant the spatial tails and the secondary term also include the zeros above A.

this script certifies UNCONDITIONALLY in Arb ball arithmetic that pi(x) > li(x) in
  [exp(727.95133539950), exp(727.95133540203)],
removing the Riemann Hypothesis assumption of Smith (2016) and sharpening Smith's Theorem 9.1
upper bound exp(727.95133540230) to exp(727.95133540203) < 1.39716616297 * 10^316.
"""
import json
import math
import os
import sys

import flint
from flint import arb

from certio import CertificationError, ball, fail, provenance, require

HERE = os.path.dirname(os.path.abspath(__file__))
CODE = [os.path.abspath(__file__), os.path.join(HERE, "certio.py")]
flint.ctx.dps = 60

A_PT = arb("3000000000000")          # Platt-Trudgian (2021): A = 3 * 10^12
SUM_INV_GAMMA2 = arb("0.0231050")    # STD (2015) Lemma 2.9
SUM_INV_GAMMA3 = arb("0.000729549")  # STD (2015) Lemma 2.9


def C0_window_arb(omega, half_width):
    """
    Arb enclosure of C_0(u) = u li(e^{u/2}) / (2 e^{u/2}) = w e^{-w} Ei(w) (w = u/2) for ALL u in
    [omega - half_width, omega + half_width]: the ball is evaluated directly (no monotonicity or asymptotic expansion).
    The bias term of the criterion is int K(y) C_0(omega + y) dy <= sup_window C_0, since K >= 0 and int K <= 1.
    """
    u = arb(omega, half_width.upper())
    w = u / 2
    require(w > 300, "w = u/2 must exceed 300")
    return w * w.ei() * (-w).exp()


def G_alpha_arb(eta, alpha):
    """G_alpha(eta) = sqrt(alpha / (2*pi)) * exp(-alpha * eta^2 / 2)."""
    pi = arb.pi()
    return (alpha / (2 * pi)).sqrt() * (-alpha * eta * eta / 2).exp()


def one_sided_gaussian_tail_bound(eta, alpha, omega, A=A_PT, n_sub=2048):
    """
    Rigorous upper bound on the one-sided hard-cutoff Gaussian tail at cut height eta > 0:
      T_hard(eta) := sum_{0 < |gamma| <= A} |int_{eta}^\infty G_alpha(x) e^{i*gamma*x} dx| * (1/|rho| + 1/(|rho|^2*(omega-eta))).
    By complex integration by parts on [eta, infty):
      |int_{eta}^\infty G_alpha(x) e^{i*gamma*x} dx| <= (1 + 1/(alpha*eta^2)) * G_alpha(eta) / sqrt((alpha*eta)^2 + gamma^2).
    """
    pi = arb.pi()
    require(alpha * eta * eta > 4, "alpha * eta^2 must exceed 4")
    G_eta = G_alpha_arb(eta, alpha)
    fac = (1 + 1 / (alpha * eta * eta)) * G_eta
    a_eta = alpha * eta

    gamma1 = arb("14.134725141734693790457251983562470270784257115699")
    t1 = 2 * pi * arb(1).exp()
    f_g1 = 1 / (gamma1 * (a_eta * a_eta + gamma1 * gamma1).sqrt())
    f_t1 = 1 / (t1 * (a_eta * a_eta + t1 * t1).sqrt())
    rem_backlund = 4 * f_t1 * t1.log() + 2 / (a_eta * t1)

    u_max = (A / t1).log()
    du = u_max / n_sub
    int_main = arb(0)
    for j in range(n_sub):
        ua = j * du
        ub = (j + 1) * du
        iv_u = arb((ua + ub) / 2, (du / 2).upper())
        t_iv = t1 * iv_u.exp()
        ig = (t_iv / (2 * pi)).log() / (a_eta * a_eta + t_iv * t_iv).sqrt()
        int_main += (ig * du).upper()
    int_main = int_main / (2 * pi)

    sum_1_rho = 2 * fac * (f_g1 + int_main + rem_backlund)
    sum_1_rho2 = 2 * (G_eta / (a_eta * (omega - eta))) * SUM_INV_GAMMA2
    return (sum_1_rho + sum_1_rho2).upper()


def mollified_spatial_tail_bound(eta, delta, alpha, omega, M_moll=190, K_sub=64):
    """
    Rigorous upper bound on the mollified spatial tail on [eta - delta, infty):
      By Fubini's theorem, 1 - w(y) = int_{-delta}^{delta} psi_delta(t) 1_{[eta + t, infty)}(y) dt,
      so the mollified tail is bounded by int_{-delta}^{delta} psi_delta(t) T_hard(eta + t) dt.
      Since T_hard(y) is decreasing in y, on each subinterval [ta, tb] subset [-delta, delta],
      T_hard(eta + t) <= T_hard(eta + ta).
    """
    pi = arb.pi()
    N = 2 * M_moll
    C_M = arb(4**M_moll * (math.factorial(M_moll) ** 2)) / arb(2 * math.factorial(N))
    h = (2 * delta) / K_sub
    total = arb(0)
    t_hard_left_edge = None
    for k in range(K_sub):
        ta = -delta + k * h
        tb = -delta + (k + 1) * h
        iv_t = arb((ta + tb) / 2, (h / 2).upper())
        psi_iv = (C_M / delta) * ((pi * iv_t / (2 * delta)).cos() ** N)
        weight = (psi_iv * h).upper()
        if weight > arb("1e-11"):
            total += weight * one_sided_gaussian_tail_bound(eta + ta, alpha, omega, n_sub=1024)
        else:
            if t_hard_left_edge is None:
                t_hard_left_edge = one_sided_gaussian_tail_bound(eta - delta, alpha, omega, n_sub=256)
            total += weight * t_hard_left_edge
    return total.upper()


def bound_unconditional_off_line(omega, eta_L, eta_R, delta, alpha, M_moll=190, A=A_PT):
    """
    Rigorous upper bound on the off-line zero sum |gamma| > A = 3*10^12 for BOTH li(e^{rho*u}) and
    (1/2)*li(e^{rho*u/2}) integrated against K(y) = G_alpha(y) * w_{eta_L, eta_R, delta}(y), via N = 2*M_moll
    integrations by parts on [omega - eta_L - delta, omega + eta_R + delta].
    - j = 0 (no derivative on w): the interior, where |w| <= 1.
    - j >= 1: the transition strips [eta - delta, eta + delta] for eta = eta_L and eta = eta_R, bounded SEPARATELY:
      ||psi_delta^{(j-1)}||_inf <= (C_M/delta) (pi/delta)^{j-1} mu[j-1] (exact binomial moments), and the m-th
      derivative of the Gaussian factor by Cauchy's estimate on a circle of radius r_m. The Cauchy exponent
      E(y) = -alpha y^2/2 + alpha y r_m + alpha r_m^2/2 is concave with maximum at y = r_m; its integral over the strip
      is <= e^{E(a)}/(alpha (a - r_m)) if r_m < a (a = eta - delta), else <= 2 delta e^{E(min(r_m, eta + delta))}.
    - The secondary term li(e^{rho u/2}) has frequency gamma/2, so N integrations by parts gain only (gamma/2)^{-N}:
      the factor (e^{u_hi/2} + 1) of the previous version becomes (e^{u_hi/2} + 2^N).
    """
    N = 2 * M_moll
    pi = arb.pi()
    u_hi = omega + eta_R + delta
    y_min = min(eta_L, eta_R) - delta
    require(y_min > 0, "y_min = min(eta_L, eta_R) - delta must be positive")
    C_M = arb(4**M_moll * (math.factorial(M_moll) ** 2)) / arb(2 * math.factorial(N))
    # Exact rational binomial moments mu[r] = 2^{-2M} * sum_{k=-M}^M comb(2M, M+k) * |k|^r
    binom_w = [arb(math.comb(N, M_moll + k)) / (arb(2) ** N) for k in range(0, M_moll + 1)]
    mu = [binom_w[0] * (0**r) + 2 * sum(binom_w[k] * (arb(k)**r) for k in range(1, M_moll + 1))
          for r in range(N + 1)]
    g0 = 2 * (alpha / (2 * pi)).sqrt()        # 2 * sup G_alpha, with |h_rho| <= 2/|gamma|

    total = (eta_L + eta_R + 2 * delta) * g0 * arb(math.factorial(N)) * ((alpha * arb(1).exp() / N) ** M_moll) / (A ** N)
    for j in range(1, N + 1):
        m = N - j
        psi_sup = (C_M / delta) * ((pi / delta) ** (j - 1)) * mu[j - 1]
        if m == 0:
            r_m, fac_m = arb(0), arb(1)
        else:
            r_m = (-alpha * y_min + (alpha * alpha * y_min * y_min + 4 * alpha * m).sqrt()) / (2 * alpha)
            fac_m = arb(math.factorial(m)) / (r_m ** m)
        strips = arb(0)
        for eta in (eta_L, eta_R):
            a, b = eta - delta, eta + delta
            slope = alpha * (a - r_m)
            length = 2 * delta
            if slope > 0:
                y_star = a
                if 1 / slope < length:
                    length = 1 / slope
            else:
                y_star = r_m if r_m < b else b
            strips += length * (-alpha * y_star * y_star / 2 + alpha * y_star * r_m + alpha * r_m * r_m / 2).exp()
        total += arb(math.comb(N, j)) * psi_sup * g0 * fac_m * strips / (A ** N)
    # sum_{|gamma| > A} |gamma|^{-(N+1)} < 2 A^{-N} log A; the two explicit-formula terms give e^{u_hi/2} + 2^N.
    return (2 * A.log() * ((u_hi / 2).exp() + arb(2) ** N) * total).upper()


def gaussian_spectral_tail(T, alpha, omega):
    """
    Zeros with |gamma| > T in the Gaussian sum: sum_{|gamma| > T} e^{-gamma^2/(2 alpha)} (1/|rho| + 1/(|rho|^2 omega)).
    With e = exp(-T^2/(2 alpha)) and the Backlund bound for N(t) (Lemma lem:backlund):
      <= e { (alpha/(pi T^2)) log(T/2pi) (1 + 1/(omega T)) + 8 (1/T + 1/(omega T^2)) log T + (4 alpha/T^3)(1 + 1/(omega T)) }.
    """
    pi = arb.pi()
    e = (-T * T / (2 * alpha)).exp()
    a = 1 / T + 1 / (omega * T * T)
    return (e * (alpha / pi * (T / (2 * pi)).log() / T**2 * (1 + 1 / (omega * T)) + 8 * a * T.log()
                 + 4 * alpha / T**3 * (1 + 1 / (omega * T)))).upper()


def variation_bound(alpha, u_lo):
    """u-variation of 1/rho + 1/(rho^2 u) across the kernel and the li remainder, with int |y| G_alpha(y) dy =
    sqrt(2/pi)/sqrt(alpha): 2 (sum gamma^-2) sqrt(2/pi)/(sqrt(alpha) u_lo^2) + 4 (sum |gamma|^-3)/u_lo^2."""
    pi = arb.pi()
    return (2 * SUM_INV_GAMMA2 * (2 / pi).sqrt() / alpha.sqrt() / u_lo**2 + 4 * SUM_INV_GAMMA3 / u_lo**2).upper()


def secondary_bound(u_lo, A=A_PT):
    """The term (1/2) li(e^{rho u/2}) for 0 < |gamma| <= A (on the critical line): |li(e^{rho u/2})| e^{-u/4} u <= ...;
    bounded directly by 2 e^{-u_lo/4} (sum_{0<gamma<=A} 1/gamma + 2 sum gamma^-2 / u_lo), with
    sum_{0<gamma<=A} 1/gamma <= log^2(A/2pi)/(4 pi) + 0.9321."""
    pi = arb.pi()
    sum1 = (A / (2 * pi)).log() ** 2 / (4 * pi) + arb("0.9321")
    return (2 * (-u_lo / 4).exp() * (sum1 + 2 * SUM_INV_GAMMA2 / u_lo)).upper()


def small_terms(u_lo, u_hi):
    """Prime powers p^m, m >= 3, the constant and the trivial-zero integral: u_hi ((log 2 + 1) e^{-u_lo/2} + (2/log 2) e^{-u_lo/6})."""
    l2 = arb(2).log()
    return (u_hi * ((l2 + 1) * (-u_lo / 2).exp() + 2 / l2 * (-u_lo / 6).exp())).upper()


def rh_spatial_above_A(c_list, alpha, omega, A=A_PT):
    """RH variant (hard cutoff): the one-sided spatial tails also get contributions from the zeros with |gamma| > A
    (on the line under RH): 2 sum_c (1 + 1/(alpha c^2)) G_alpha(c) log(A)/A (1 + 1/(omega A))."""
    fac = sum(((1 + 1 / (alpha * c * c)) * G_alpha_arb(c, alpha) for c in c_list), arb(0))
    return (2 * fac * A.log() / A * (1 + 1 / (omega * A))).upper()


def rh_secondary_above_A(u_lo, alpha, A=A_PT):
    """RH variant: the secondary term li(e^{rho u/2}) for |gamma| > A: 8 e^{-u_lo/4} (4 sqrt(alpha/2pi) + 1/4 + 1/300) log(A)/A."""
    pi = arb.pi()
    g0 = (alpha / (2 * pi)).sqrt()
    return (8 * (-u_lo / 4).exp() * (4 * g0 + arb(1) / 4 + arb(1) / 300) * A.log() / A).upper()


def main():
    omega = arb("727.951335401")
    alpha = arb("4.5e19")
    T = arb("30610046000")
    eta_L = arb("1.00e-9")
    eta_R = arb("0.725e-9")
    delta = arb("0.305e-9")
    M_moll = 190

    eta_L0 = eta_L + delta
    eta_R0 = eta_R + delta
    u_lo = omega - eta_L0
    u_hi = omega + eta_R0

    # Assumption (the theorem's hypothesis): Platt's value of the Gaussian zero sum, Smith (2016) Table 4.
    S12_platt = arb("1.00288049783641096204", arb("1e-20"))

    # Unconditional criterion:  I >= S12 - C - (spatial + spectral + variation + secondary + small + off-line).
    C = C0_window_arb(omega, eta_L0 if eta_L0 > eta_R0 else eta_R0)
    spatial_L = mollified_spatial_tail_bound(eta_L, delta, alpha, omega, M_moll=M_moll)
    spatial_R = mollified_spatial_tail_bound(eta_R, delta, alpha, omega, M_moll=M_moll)
    spectral = gaussian_spectral_tail(T, alpha, omega)
    variation = variation_bound(alpha, u_lo)
    secondary = secondary_bound(u_lo)
    small = small_terms(u_lo, u_hi)
    R_off = bound_unconditional_off_line(omega, eta_L, eta_R, delta, alpha, M_moll=M_moll)
    R_total = (spatial_L + spatial_R + spectral + variation + secondary + small + R_off).upper()
    I = S12_platt - C - R_total
    I_lower = I.lower()
    require(I_lower > 0, f"Unconditional margin I_lower = {I.str(10)} is NOT positive!")

    # Some x in [e^u_lo, e^u_hi] has pi(x) - li(x) >= Q := I_lower e^{u_lo/2}/u_hi; the integer rounding of
    # check_theorem.py gives more than (Q - 1) u_lo - 1 consecutive integers.
    peak_diff_lower = (I_lower * (u_lo / 2).exp() / u_hi).lower()
    consec_ints_lower = ((peak_diff_lower - 1) * u_lo - 1).lower()

    # RH-conditional variant with hard asymmetric cutoff [-eta_L_rh, +eta_R_rh] (delta = 0; all zeros on the line):
    eta_L_rh = arb("1.30e-9")
    eta_R_rh = arb("0.708e-9")
    u_lo_rh = omega - eta_L_rh
    u_hi_rh = omega + eta_R_rh
    C_rh = C0_window_arb(omega, eta_L_rh)
    spatial_rh = (one_sided_gaussian_tail_bound(eta_L_rh, alpha, omega, n_sub=2048)
                  + one_sided_gaussian_tail_bound(eta_R_rh, alpha, omega, n_sub=2048)).upper()
    spatial_above_A = rh_spatial_above_A([eta_L_rh, eta_R_rh], alpha, omega)
    secondary_above_A = rh_secondary_above_A(u_lo_rh, alpha)
    R_total_rh = (spatial_rh + spatial_above_A + spectral + variation + secondary + secondary_above_A + small).upper()
    I_rh = S12_platt - C_rh - R_total_rh
    I_lower_rh = I_rh.lower()
    require(I_lower_rh > 0, f"RH-conditional margin I_lower_rh = {I_rh.str(10)} is NOT positive!")

    params = {
        "omega": "727.951335401", "alpha": "4.5e19", "T": "30610046000",
        "eta_L": "1.00e-9", "eta_R": "0.725e-9", "delta": "0.305e-9", "M_moll": M_moll,
    }
    cert = provenance(params, CODE)
    cert.update({
        "assumption": "S1 + S2 = 1.00288049783641096204 +- 1e-20 (Platt, reported in Smith 2016, Table 4)",
        "omega": "727.951335401",
        "eta_L0": ball(eta_L0), "eta_R0": ball(eta_R0), "u_lo": ball(u_lo), "u_hi": ball(u_hi),
        "S12_platt": ball(S12_platt), "C_window": ball(C),
        "spatial_left": ball(spatial_L), "spatial_right": ball(spatial_R), "spectral_tail_T": ball(spectral),
        "variation": ball(variation), "secondary": ball(secondary), "small": ball(small),
        "off_line_above_A": ball(R_off), "R_total": ball(R_total),
        "I": ball(I), "I_lower": ball(I_lower),
        "Xi_upper_bound": ball(u_hi.exp().upper(), 20),
        "consecutive_integers_lower": ball(consec_ints_lower, 15),
        "rh_conditional": {
            "eta_L_rh": ball(eta_L_rh), "eta_R_rh": ball(eta_R_rh), "u_lo_rh": ball(u_lo_rh), "u_hi_rh": ball(u_hi_rh),
            "C_window": ball(C_rh), "spatial": ball(spatial_rh), "spatial_above_A": ball(spatial_above_A),
            "secondary_above_A": ball(secondary_above_A), "R_total_rh": ball(R_total_rh),
            "I_rh": ball(I_rh), "I_lower_rh": ball(I_lower_rh),
            "Xi_upper_bound_rh": ball(u_hi_rh.exp().upper(), 20),
        },
    })
    out_path = os.path.join(HERE, "smith_unconditional_certificate.json")
    with open(out_path, "w") as f:
        json.dump(cert, f, indent=2)

    print("======================================================================")
    print("REFINEMENT OF SMITH (2016), CONDITIONAL ON PLATT'S GAUSSIAN SUM")
    print("======================================================================")
    print(f"S1 + S2 (assumed)   : {S12_platt.str(22)}")
    print(f"C on window         : {C.str(15)}")
    print(f"  spatial (L, R)    <= {spatial_L.str(6, radius=False)}, {spatial_R.str(10, radius=False)}")
    print(f"  spectral |g| > T  <= {spectral.str(10, radius=False)}")
    print(f"  variation         <= {variation.str(6, radius=False)}")
    print(f"  secondary, small  <= {secondary.str(4, radius=False)}, {small.str(4, radius=False)}")
    print(f"  zeros |g| > A     <= {R_off.str(6, radius=False)}")
    print(f"UNCONDITIONAL MARGIN: I = {I.str(12)} (lower {I_lower.str(10, radius=False)}) > 0")
    print(f"Window              : u in [{u_lo.lower().str(16, radius=False)}, {u_hi.upper().str(16, radius=False)}]")
    print(f"Consecutive integers: > {consec_ints_lower.str(8, radius=False)}")
    print(f"RH: spatial {spatial_rh.str(8, radius=False)}, above A {spatial_above_A.str(4, radius=False)} + {secondary_above_A.str(4, radius=False)}")
    print(f"RH MARGIN           : I = {I_rh.str(12)} (lower {I_lower_rh.str(10, radius=False)}) > 0")
    print(f"RH window           : u in [{u_lo_rh.lower().str(16, radius=False)}, {u_hi_rh.upper().str(16, radius=False)}]")
    print("======================================================================")


if __name__ == "__main__":
    try:
        main()
    except CertificationError as e:
        fail(str(e))
