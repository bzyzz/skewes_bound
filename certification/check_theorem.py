#!/usr/bin/env python3
"""
Combine the certified zero-sum enclosure in s_certificate.json with the rigorous Arb
error bounds in error_bounds.py to verify:

    I(omega, eta0) >= J_K(omega, T0) - (R_1 + R_2 + R_3 + R_4 + R_5) > 0,

and deduce the unconditional and RH-conditional upper bounds on Skewes's number Xi,
the positive-region interval [exp(omega - eta0), exp(omega + eta0)], and the certified
number of consecutive integers x with pi(x) > li(x).
"""
import json
import os
import sys

import flint
from flint import arb

from certio import CertificationError, ball, fail, require
from error_bounds import phi_weight, total_error

HERE = os.path.dirname(os.path.abspath(__file__))
flint.ctx.dps = 60


def kappa_core_tail_integral_arb(x0, lam1, lam2, lam3):
    """
    Exact closed-form Arb evaluation of \int_{x0}^1 kappa_core(x) dx for x0 in [0, 1]:
      \int_{x0}^1 kappa_2(x) dx = (1 - x0)/2 - sin(pi*x0)/(2*pi)
      \int_{x0}^1 kappa_4(x) dx = (1 - x0)/2 - 2*sin(pi*x0)/(3*pi) - sin(2*pi*x0)/(12*pi)
      \int_{x0}^1 kappa_6(x) dx = (1 - x0)/2 - 3*sin(pi*x0)/(4*pi) - 3*sin(2*pi*x0)/(16*pi) - sin(3*pi*x0)/(48*pi)
    """
    require(x0 >= 0 and x0 <= 1, "x0 must lie in [0, 1]")
    pi = arb.pi()
    s1 = (pi * x0).sin() / pi
    s2 = (2 * pi * x0).sin() / pi
    s3 = (3 * pi * x0).sin() / pi
    half_rem = (1 - x0) / 2
    i2 = half_rem - s1 / 2
    i4 = half_rem - (2 * s1) / 3 - s2 / 12
    i6 = half_rem - (3 * s1) / 4 - (3 * s2) / 16 - s3 / 48
    return lam1 * i2 + lam2 * i4 + lam3 * i6


def rh_sharpened_eta(omega, eta1, delta, lam1, lam2, lam3, I_lower):
    """
    Under RH, Schoenfeld (1976) [STD 2015 Theorem 9.2] gives |pi(x) - li(x)| < sqrt(x) log(x) / (8*pi)
    for x >= 2657. Thus for 0 < eta_rh < eta0 = eta1 + delta:
      |I(omega, eta0) - I(omega, eta_rh)|
        <= (( (omega - eta_rh)^2 + (omega + eta0)^2 ) / (8 * pi * Phi(omega + eta0))) * \int_{(eta_rh - delta)/eta1}^1 kappa_core(x) dx.
    Also for one-sided right-only trimming on [omega - eta0, omega + eta_plus]:
      |I(omega, eta0) - I_one_sided(omega, -eta0, eta_plus)|
        <= ((omega + eta0)^2 / (8 * pi * Phi(omega + eta0))) * \int_{(eta_plus - delta)/eta1}^1 kappa_core(x) dx.
    """
    pi = arb.pi()
    eta0 = eta1 + delta
    u_hi = omega + eta0
    phi_hi = phi_weight(u_hi)
    # 1. Symmetric RH sharpening: binary search for smallest x0 in [0, 1] with tail_err < I_lower * 0.98
    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        x0 = arb(mid)
        eta_cand = x0 * eta1 + delta
        fac = ((omega - eta_cand) ** 2 + u_hi ** 2) / (8 * pi * phi_hi)
        err = (fac * kappa_core_tail_integral_arb(x0, lam1, lam2, lam3)).upper()
        if err < I_lower * arb("0.98"):
            hi = mid
        else:
            lo = mid
    x0_sym = arb(hi)
    eta_rh_sym = x0_sym * eta1 + delta
    fac_sym = ((omega - eta_rh_sym) ** 2 + u_hi ** 2) / (8 * pi * phi_hi)
    err_sym = (fac_sym * kappa_core_tail_integral_arb(x0_sym, lam1, lam2, lam3)).upper()
    I_rh_sym = (I_lower - err_sym).lower()
    require(I_rh_sym > 0, "RH symmetric sharpened integral not certified positive")

    # 2. One-sided right-only RH sharpening:
    lo, hi = 0.0, 1.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        x0 = arb(mid)
        fac = (u_hi ** 2) / (8 * pi * phi_hi)
        err = (fac * kappa_core_tail_integral_arb(x0, lam1, lam2, lam3)).upper()
        if err < I_lower * arb("0.98"):
            hi = mid
        else:
            lo = mid
    x0_one = arb(hi)
    eta_rh_one = x0_one * eta1 + delta
    fac_one = (u_hi ** 2) / (8 * pi * phi_hi)
    err_one = (fac_one * kappa_core_tail_integral_arb(x0_one, lam1, lam2, lam3)).upper()
    I_rh_one = (I_lower - err_one).lower()
    require(I_rh_one > 0, "RH one-sided sharpened integral not certified positive")

    return (eta_rh_sym, I_rh_sym, err_sym), (eta_rh_one, I_rh_one, err_one)


def main():
    cert_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "s_certificate.json")
    with open(cert_path) as f:
        cert = json.load(f)

    require(cert.get("schema") == 2, "expected schema-2 certificate")
    omega = arb(cert["omega"])
    eta1 = arb(cert["eta1"])
    delta = arb(cert["delta"])
    eta0 = eta1 + delta
    T0 = arb(cert["T0"])
    M_moll = int(cert["M_moll"])
    parts_lam = [int(x) for x in cert["lam"].split(",")]
    q_lam = sum(parts_lam)
    lam1, lam2, lam3 = arb(parts_lam[0]) / q_lam, arb(parts_lam[1]) / q_lam, arb(parts_lam[2]) / q_lam

    S_K = arb(cert["S_K"])
    Phi_om = phi_weight(omega)
    J_K = -1 + S_K / Phi_om
    require(J_K.overlaps(arb(cert["J_K"])), "recomputed J_K does not overlap certificate J_K")
    J_K_lo = arb(cert["J_K_lower_ball"]).lower()
    require(J_K.overlaps(arb(cert["J_K_lower_ball"])), "J_K_lower_ball does not overlap J_K")

    print("======================================================================")
    print("RIGOROUS VERIFICATION OF SKEWES NUMBER UPPER BOUND")
    print("======================================================================")
    print(f"Certificate : {cert_path}")
    print(f"Parameters  : omega = {cert['omega']}, eta1 = {cert['eta1']}, delta = {cert['delta']}")
    print(f"              T0 = {cert['T0']} (N(T0) = {cert['N_T0']:,}), lam = {cert['lam']}, M_moll = {M_moll}")
    print(f"Phi(omega)  : {Phi_om.str(20)}")
    print(f"S_K(omega)  : {S_K.str(20)}")
    print(f"J_K(omega)  : {arb(cert['J_K']).str(15)}  (lower bound >= {J_K_lo.str(12, radius=False)})")
    print("\nError bounds (Arb ball arithmetic):")
    R_tot, r_parts = total_error(omega, eta1, delta, T0, lam1, lam2, lam3, M_moll=M_moll, verbose=True)

    I_lower = (J_K_lo - R_tot).lower()
    R_tot_env = (r_parts["R2"] + r_parts["R3"] + r_parts["R4_envelope"] + r_parts["R5"]).upper()
    I_lower_env = (J_K_lo - R_tot_env).lower()
    require(I_lower > 0, f"Margin J_K - R_tot = {I_lower.str(10)} is NOT positive!")

    u_lo = (omega - eta0).lower()
    u_hi = (omega + eta0).upper()
    # By Lemma 9.4 & 9.5 of STD (2015), max_{u in [u_lo, u_hi]} (pi(e^u) - li(e^u)) >= I_lower * e^{u_lo/2} * Phi(u_hi) / u_hi
    use_I = I_lower_env if I_lower_env > 0 else I_lower
    peak_diff_lower = (use_I * (u_lo / 2).exp() * phi_weight(u_hi) / u_hi).lower()
    consec_ints_lower = (peak_diff_lower * u_lo).lower()

    (eta_rh_sym, I_rh_sym, err_sym), (eta_rh_one, I_rh_one, err_one) = rh_sharpened_eta(
        omega, eta1, delta, lam1, lam2, lam3, use_I
    )
    u_rh_sym_lo = (omega - eta_rh_sym).lower()
    u_rh_sym_hi = (omega + eta_rh_sym).upper()
    u_rh_one_hi = (omega + eta_rh_one).upper()

    print("\n----------------------------------------------------------------------")
    print("UNCONDITIONAL CONCLUSION (using Platt-Trudgian A = 3*10^12):")
    print(f"  Margin (envelope)   : I(omega, eta0) >= {I_lower_env.str(10, radius=False)}")
    print(f"  Margin (Stieltjes)  : I(omega, eta0) >= {I_lower.str(10, radius=False)}")
    print(f"  eta0 = eta1 + delta  = {eta0.str(12, radius=False)}")
    print(f"  Positive interval    : u in [{u_lo.str(15, radius=False)}, {u_hi.str(15, radius=False)}]")
    print(f"  Upper bound on Xi    : Xi < exp({u_hi.str(15, radius=False)})")
    print(f"                              = {u_hi.exp().upper().str(12, radius=False)}")
    print(f"  Max pi(x) - li(x)   >= {peak_diff_lower.str(8, radius=False)}")
    print(f"  Consecutive integers>= {consec_ints_lower.str(8, radius=False)}")
    print("----------------------------------------------------------------------")
    print("CONDITIONAL CONCLUSION (assuming RH, via Schoenfeld Theorem 9.2):")
    print(f"  Symmetric interval   : u in [{u_rh_sym_lo.str(15, radius=False)}, {u_rh_sym_hi.str(15, radius=False)}]")
    print(f"  One-sided upper bnd  : Xi < exp({u_rh_one_hi.str(15, radius=False)})")
    print("======================================================================")
    print("VERIFICATION PASSED.")


if __name__ == "__main__":
    try:
        main()
    except CertificationError as e:
        fail(str(e))
