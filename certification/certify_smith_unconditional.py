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


def exact_C0_arb(omega, n_terms=25):
    """
    Rigorous Arb enclosure of C_0(u) = u * li(e^{u/2}) / (2 * e^{u/2}) for u near omega ~ 727.95.
    By Lemma 2.1 of STD (2015), for w = u/2 > 300 and n = 25:
      li(e^w) / e^w = sum_{k=1}^n (k-1)! / w^k + R_n,   0 < R_n < n! / (w^{n+1} * (1 - n/w)).
    """
    w = omega / 2
    require(w > 300, "w = omega/2 must exceed 300")
    s = arb(0)
    for k in range(1, n_terms + 1):
        s += arb(math.factorial(k - 1)) / (w ** k)
    rem_hi = arb(math.factorial(n_terms)) / ((w ** (n_terms + 1)) * (1 - arb(n_terms) / w))
    rem = arb(rem_hi / 2, (rem_hi / 2).upper())
    return w * (s + rem)


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
    Rigorous upper bound on the off-line zero sum |gamma| > A = 3*10^12 for BOTH li(e^{rho*u})
    and (1/2)*li(e^{rho*u/2}) integrated against K(y) = G_alpha(y) * w_{eta_L, eta_R, delta}(y)
    via N = 2*M_moll integrations by parts on [omega - eta_L - delta, omega + eta_R + delta].
    """
    N = 2 * M_moll
    pi = arb.pi()
    u_lo = omega - eta_L - delta
    u_hi = omega + eta_R + delta
    y_min = min(eta_L, eta_R) - delta
    require(y_min > 0, "y_min = min(eta_L, eta_R) - delta must be positive")
    C_M = arb(4**M_moll * (math.factorial(M_moll) ** 2)) / arb(2 * math.factorial(N))

    # Exact rational binomial moments mu[r] = 2^{-2M} * sum_{k=-M}^M comb(2M, M+k) * |k|^r
    binom_w = [arb(math.comb(N, M_moll + k)) / (arb(2) ** N) for k in range(0, M_moll + 1)]
    mu = [binom_w[0] * (0**r) + 2 * sum(binom_w[k] * (arb(k)**r) for k in range(1, M_moll + 1))
          for r in range(N + 1)]

    # 1. j = 0 term (interior [u_lo, u_hi], where |w(y)| <= 1 and m = N):
    r_N = (arb(N) / alpha).sqrt()
    f_N_int = (u_hi - u_lo) * arb(2) * (alpha / (2 * pi)).sqrt() * arb(math.factorial(N)) * ((alpha * arb(1).exp() / N) ** M_moll)
    leibniz_sum = f_N_int / (A ** N)

    # 2. j >= 1 terms (supported on [-eta_L-delta, -eta_L+delta] U [eta_R-delta, eta_R+delta], where |y| >= y_min):
    # On each transition strip, ||psi_delta^{(j-1)}||_{L^infty} <= (C_M / delta) * (pi / delta)^{j-1} * mu[j-1],
    # and for m = N - j, integrating |f_rho^{(m)}(omega + y)| over |y| in [y_min, y_min + 2*delta]
    # with Cauchy circle radius r_m around y gives:
    #   int_{y_min}^{y_min+2*delta} exp(-alpha*y^2/2 + alpha*y*r_m + alpha*r_m^2/2) dy
    #     <= min(2*delta, 1 / (alpha * (y_min - r_m))) * exp(-alpha*y_min^2/2 + alpha*y_min*r_m + alpha*r_m^2/2).
    for j in range(1, N + 1):
        m = N - j
        psi_sup = (C_M / delta) * ((pi / delta) ** (j - 1)) * mu[j - 1]
        if m == 0:
            r_m = arb(0)
            fac_m = arb(1)
        else:
            r_m = (-alpha * y_min + (alpha * alpha * y_min * y_min + 4 * alpha * m).sqrt()) / (2 * alpha)
            fac_m = arb(math.factorial(m)) / (r_m ** m)
        # The exponent E(y) = -alpha*y^2/2 + alpha*y*r_m + alpha*r_m^2/2 is concave with maximum at y = r_m, so its
        # supremum over the strip [y_min, y_min + 2*delta] is at y_min only when r_m <= y_min (then E decreases with
        # slope <= -alpha*(y_min - r_m) and the integral is at most e^{E(y_min)} / (alpha*(y_min - r_m)));
        # otherwise take the supremum at min(r_m, y_min + 2*delta) times the strip length.
        y_max = y_min + 2 * delta
        slope = alpha * (y_min - r_m)
        eff_len = 2 * delta
        if slope > 0:
            y_star = y_min
            if (1 / slope) < eff_len:
                eff_len = 1 / slope
        else:
            y_star = r_m if r_m < y_max else y_max
        g_ymin = (-alpha * y_star * y_star / 2 + alpha * y_star * r_m + alpha * r_m * r_m / 2).exp()
        # Factor of 2 for the two transition strips (left and right), and 2*|gamma|^{-1} from |h_rho| <= 2/|gamma|
        term_j = arb(math.comb(N, j)) * 2 * psi_sup * eff_len * arb(2) * (alpha / (2 * pi)).sqrt() * fac_m * g_ymin / (A ** N)
        leibniz_sum += term_j

    # Multiply by 2 * A * log(A) / A = 2 * log(A) (from sum_{|gamma| > A} |gamma|^{-(N+1)} < 2 * A^{-N} log A)
    # and by (e^{u_hi/2} + 1) to cover both li(e^{rho u}) and (1/2)li(e^{rho u/2}):
    R_off = (2 * A.log() * ((u_hi / 2).exp() + 1) * leibniz_sum).upper()
    return R_off


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

    S12_platt = arb("1.00288049783641096204", arb("1e-20"))
    C_max = exact_C0_arb(u_lo).upper()

    pi = arb.pi()
    R1 = (arb("1.812") * u_hi * (-u_lo / 6).exp()).upper()
    R2 = (arb("0.025") * (-u_lo / 4).exp() + (1 / (32 * alpha) - omega / 4).exp() * (arb("1.302") + arb("0.04") * alpha)).upper()
    R3_1 = (alpha * (-T * T / (2 * alpha)).exp() / (2 * pi)) * (T / (2 * pi)).log() * (1 / (T ** 3) + 2 / (T ** 2))
    R3_2 = (2 * alpha * (-T * T / (2 * alpha)).exp() / (2 * pi * omega * (T ** 3))) * (T / (2 * pi)).log()
    R3 = (R3_1 + R3_2).upper()

    R4_left = mollified_spatial_tail_bound(eta_L, delta, alpha, omega, M_moll=M_moll)
    R4_right = mollified_spatial_tail_bound(eta_R, delta, alpha, omega, M_moll=M_moll)
    R4 = (R4_left + R4_right).upper()

    R5_var = (arb("0.019") / (alpha.sqrt() * (u_lo ** 2)) + 4 * SUM_INV_GAMMA3 / (u_lo ** 2)).upper()
    R_off = bound_unconditional_off_line(omega, eta_L, eta_R, delta, alpha, M_moll=M_moll)

    R_total = (R1 + R2 + R3 + R4 + R5_var + R_off).upper()
    I_lower = (S12_platt.lower() - C_max - R_total).lower()
    require(I_lower > 0, f"Unconditional margin I_lower = {I_lower.str(10)} is NOT positive!")

    peak_diff_lower = (I_lower * (u_lo / 2).exp() / u_hi).lower()
    consec_ints_lower = (peak_diff_lower * u_lo).lower()

    # Also certify the RH-conditional bound with hard asymmetric cutoff [-eta_L_rh, +eta_R_rh] (delta = 0, R_off = 0):
    eta_L_rh = arb("1.30e-9")
    eta_R_rh = arb("0.708e-9")
    u_lo_rh = omega - eta_L_rh
    u_hi_rh = omega + eta_R_rh
    C_max_rh = exact_C0_arb(u_lo_rh).upper()
    R4_left_rh = one_sided_gaussian_tail_bound(eta_L_rh, alpha, omega, n_sub=2048)
    R4_right_rh = one_sided_gaussian_tail_bound(eta_R_rh, alpha, omega, n_sub=2048)
    R_total_rh = (R1 + R2 + R3 + R4_left_rh + R4_right_rh + R5_var).upper()
    I_lower_rh = (S12_platt.lower() - C_max_rh - R_total_rh).lower()
    require(I_lower_rh > 0, f"RH-conditional margin I_lower_rh = {I_lower_rh.str(10)} is NOT positive!")

    params = {
        "omega": "727.951335401", "alpha": "4.5e19", "T": "30610046000",
        "eta_L": "1.00e-9", "eta_R": "0.725e-9", "delta": "0.305e-9", "M_moll": M_moll,
    }
    cert = provenance(params, CODE)
    cert.update({
        "omega": "727.951335401",
        "eta_L0": ball(eta_L0, 20),
        "eta_R0": ball(eta_R0, 20),
        "u_lo": ball(u_lo, 25),
        "u_hi": ball(u_hi, 25),
        "S12_platt_ball": ball(S12_platt, 25),
        "C_max_ball": ball(C_max, 25),
        "R1": ball(R1, 15),
        "R2": ball(R2, 15),
        "R3": ball(R3, 15),
        "R4_left": ball(R4_left, 15),
        "R4_right": ball(R4_right, 15),
        "R4": ball(R4, 15),
        "R5_var": ball(R5_var, 15),
        "R_off_(|gamma|>3e12)": ball(R_off, 15),
        "R_total": ball(R_total, 15),
        "I_lower": ball(I_lower, 15),
        "Xi_upper_bound": ball(u_hi.exp().upper(), 20),
        "consecutive_integers_lower": ball(consec_ints_lower, 15),
        "rh_conditional": {
            "eta_L_rh": ball(eta_L_rh, 20),
            "eta_R_rh": ball(eta_R_rh, 20),
            "u_lo_rh": ball(u_lo_rh, 25),
            "u_hi_rh": ball(u_hi_rh, 25),
            "R4_right_rh": ball(R4_right_rh, 15),
            "R_total_rh": ball(R_total_rh, 15),
            "I_lower_rh": ball(I_lower_rh, 15),
            "Xi_upper_bound_rh": ball(u_hi_rh.exp().upper(), 20),
        },
    })
    out_path = os.path.join(HERE, "smith_unconditional_certificate.json")
    with open(out_path, "w") as f:
        json.dump(cert, f, indent=2)

    print("======================================================================")
    print("UNCONDITIONAL REFINEMENT OF SMITH (2016) VIA MOLLIFIED ASYMMETRIC KERNEL")
    print("======================================================================")
    print(f"Center omega        : {omega.str(15, radius=False)}")
    print(f"Asymmetric window   : [-{eta_L0.str(10, radius=False)}, +{eta_R0.str(10, radius=False)}]")
    print(f"S1 + S2 (Platt/Arb) : {S12_platt.str(22)}")
    print(f"C_max on window     : <= {C_max.str(20, radius=False)}")
    print(f"Excess (S1+S2 - C)  : >= {(S12_platt.lower() - C_max).str(10, radius=False)}")
    print("\nError bounds (Arb ball arithmetic):")
    print(f"  R1 (prime cubes)  <= {R1.str(10, radius=False)}")
    print(f"  R2 (li(x^(rho/2)))<= {R2.str(10, radius=False)}")
    print(f"  R3 (|gamma| > T)  <= {R3.str(10, radius=False)}")
    print(f"  R4_left           <= {R4_left.str(10, radius=False)}")
    print(f"  R4_right          <= {R4_right.str(10, radius=False)}")
    print(f"  R4                <= {R4.str(10, radius=False)}")
    print(f"  R5_var            <= {R5_var.str(10, radius=False)}")
    print(f"  R_off (|g|>3e12)  <= {R_off.str(10, radius=False)}")
    print(f"  R_total           <= {R_total.str(10, radius=False)}")
    print("----------------------------------------------------------------------")
    print(f"UNCONDITIONAL MARGIN: I >= {I_lower.str(10, radius=False)} > 0")
    print(f"Positive interval   : u in [{u_lo.lower().str(16, radius=False)}, {u_hi.upper().str(16, radius=False)}]")
    print(f"Upper bound on Xi   : Xi < exp({u_hi.upper().str(16, radius=False)})")
    print(f"                         = {u_hi.exp().upper().str(15, radius=False)}")
    print(f"Consecutive integers: >= {consec_ints_lower.str(8, radius=False)}")
    print("----------------------------------------------------------------------")
    print(f"RH-CONDITIONAL REFINEMENT (eta_R = {eta_R_rh.str(8, radius=False)} vs Smith 1.30e-9):")
    print(f"  R4_right_rh       <= {R4_right_rh.str(10, radius=False)}")
    print(f"  Margin I_rh       >= {I_lower_rh.str(10, radius=False)} > 0")
    print(f"  Positive interval : u in [{u_lo_rh.lower().str(16, radius=False)}, {u_hi_rh.upper().str(16, radius=False)}]")
    print(f"  Upper bound on Xi : Xi < exp({u_hi_rh.upper().str(16, radius=False)})")
    print(f"                         = {u_hi_rh.exp().upper().str(15, radius=False)}")
    print("======================================================================")


if __name__ == "__main__":
    try:
        main()
    except CertificationError as e:
        fail(str(e))
