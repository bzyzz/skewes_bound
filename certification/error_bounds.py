#!/usr/bin/env python3
"""
Rigorous Arb ball-arithmetic evaluation of the kernel Fourier transform hat{K}(gamma)
and all error terms R_1, R_2, R_3, R_4, R_5 for the Mollified Cosine-Power criterion
on pi(x) - li(x).

Kernel definition:
  K(y) = (K_core * psi_delta)(y),   supp(K) = [-eta0, eta0],   eta0 = eta1 + delta,
where:
  K_core(y) = (1 / eta1) * sum_{m=1}^3 lam_m * kappa_{2m}(y / eta1),
  psi_delta(y) = (1 / delta) * kappa_{2M}(y / delta),
  kappa_{2m}(x) = (4^m (m!)^2 / (2 (2m)!)) * cos^{2m}(pi * x / 2) * 1_{|x| <= 1}.

Fourier transform (entire function of gamma):
  hat{K}(gamma) = hat{kappa}_core(eta1 * gamma) * hat{kappa}_{2M}(delta * gamma),
  hat{kappa}_{2m}(xi) = (sin(xi) / xi) * prod_{k=1}^m ( (k*pi)^2 / ((k*pi)^2 - xi^2) ).
"""
import math
import flint
from flint import arb

from certio import require

flint.ctx.dps = 50

A_PT = arb("3000000000000")          # Platt-Trudgian (2021) RH verification height: 3 * 10^12
SUM_INV_GAMMA2 = arb("0.0231050")    # STD (2015) Lemma 2.9: sum_{gamma > 0} 1/gamma^2 < 2.31050e-2
SUM_INV_GAMMA3 = arb("0.000729549")  # STD (2015) Lemma 2.9: sum_{gamma > 0} 1/gamma^3 < 7.29549e-4


def phi_weight(u):
    """Phi(u) = 1 + 2/u + 8/u^2 + 58.56/u^3 from Dusart (2018) Theorem 5.1 (valid for u >= 44.22)."""
    return 1 + 2 / u + 8 / (u * u) + arb("58.56") / (u * u * u)


def kappa_2m_hat_arb(xi, m):
    """
    Rigorous Arb evaluation of hat{kappa}_{2m}(xi) = sinc(xi/pi) * prod_{k=1}^m ((k*pi)^2 / ((k*pi)^2 - xi^2)).
    Handles the removable singularities at xi = k*pi (1 <= k <= m) safely.
    """
    pi = arb.pi()
    xi2 = xi * xi
    for k in range(1, m + 1):
        kpi = k * pi
        if (xi - kpi).contains(0) or (xi + kpi).contains(0):
            u = (xi - kpi) / pi if (xi - kpi).contains(0) else (xi + kpi) / pi
            sinc_part = ((-1) ** (k - 1)) * u.sinc_pi() / (xi * (kpi + abs(xi)))
            res = sinc_part * (kpi * kpi)
            for j in range(1, m + 1):
                if j != k:
                    jp2 = (j * pi) * (j * pi)
                    res = res * (jp2 / (jp2 - xi2))
            return res
    res = (xi / pi).sinc_pi()
    for k in range(1, m + 1):
        kp2 = (k * pi) * (k * pi)
        res = res * (kp2 / (kp2 - xi2))
    return res


def kappa_core_hat_arb(xi, lam1, lam2, lam3):
    """hat{kappa}_core(xi) = lam1 * hat{kappa}_2(xi) + lam2 * hat{kappa}_4(xi) + lam3 * hat{kappa}_6(xi)."""
    return (lam1 * kappa_2m_hat_arb(xi, 1)
            + lam2 * kappa_2m_hat_arb(xi, 2)
            + lam3 * kappa_2m_hat_arb(xi, 3))


def kappa_core_rational_factor_arb(xi, lam1, lam2, lam3):
    """
    For xi > 3*pi, hat{kappa}_core(xi) = sin(xi) * R_core(xi), where
      R_core(xi) = -pi^2 * (lam1*xi^4 - pi^2*(13*lam1 + 4*lam2)*xi^2 + 36*pi^4)
                   / (xi * (xi^2 - pi^2) * (xi^2 - 4*pi^2) * (xi^2 - 9*pi^2)).
    When lam1 > 0, we complete the square in xi^2 in the numerator so xi^2 appears only once,
    avoiding interval-arithmetic dependency blowup on wide subintervals.
    """
    pi = arb.pi()
    pi2 = pi * pi
    pi4 = pi2 * pi2
    xi2 = xi * xi
    require(xi > 3 * pi, "kappa_core_rational_factor_arb requires xi > 3*pi")
    if lam1 > 0:
        c_half = pi2 * (13 * lam1 + 4 * lam2) / (2 * lam1)
        quad = lam1 * (xi2 - c_half) ** 2 + (36 * pi4 - lam1 * c_half * c_half)
        num = pi2 * abs(quad)
    else:
        num = pi2 * abs(-pi2 * (4 * lam2) * xi2 + 36 * pi4)
    den = xi * (xi2 - pi2) * (xi2 - 4 * pi2) * (xi2 - 9 * pi2)
    return num / den


def kernel_hat_arb(gamma, eta1, delta, lam1, lam2, lam3, M_moll=60):
    """Full Fourier transform hat{K}(gamma) = hat{kappa}_core(eta1 * gamma) * hat{kappa}_{2M}(delta * gamma)."""
    return kappa_core_hat_arb(eta1 * gamma, lam1, lam2, lam3) * kappa_2m_hat_arb(delta * gamma, M_moll)


def phi_2m_majorant_arb(xi, m):
    """Strictly decreasing majorant (1/xi) * prod_{k=1}^m ((k*pi)^2 / (xi^2 - (k*pi)^2)) for xi > m*pi."""
    pi = arb.pi()
    require(xi > m * pi, f"phi_2m_majorant_arb requires xi > {m}*pi")
    xi2 = xi * xi
    res = 1 / xi
    for k in range(1, m + 1):
        kp2 = (k * pi) * (k * pi)
        res = res * (kp2 / (xi2 - kp2))
    return res


def phi_core_tail_majorant_arb(xi, lam1, lam2, lam3):
    """Strictly decreasing majorant of |R_core(xi)| (and |hat{kappa}_core(xi)|) for xi > 3*pi."""
    return (lam1 * phi_2m_majorant_arb(xi, 1)
            + lam2 * phi_2m_majorant_arb(xi, 2)
            + lam3 * phi_2m_majorant_arb(xi, 3))


def bound_R2(omega, eta0):
    """Prime-cube and log(2) error term R2."""
    u_lo = omega - eta0
    u_hi = omega + eta0
    require(u_lo >= arb("44.22"), "hypothesis omega - eta0 >= 44.22 fails")
    log2 = arb(2).log()
    return (u_hi * (log2 * (-u_lo / 2).exp() + (2 / log2) * (-u_lo / 6).exp())).upper()


def bound_R3(omega, eta0, A=A_PT):
    """
    On-line zero variation and li(e^{rho u}) asymptotic truncation error R3 for 0 < |gamma| <= A.
    Because supp(K) = [-eta0, eta0], int_{-eta0}^{eta0} K(y) e^{i*gamma*y} dy = hat{K}(gamma) is exact,
    so the Gaussian truncation error 0.19*K(eta) of STD (2015) is identically zero.
    """
    u_lo = omega - eta0
    require(u_lo >= arb("44.22"), "omega - eta0 >= 44.22 fails")
    require(eta0 <= omega / 100, "eta0 <= omega / 100 fails")
    pi = arb.pi()
    c_f = 2 + 16 / u_lo + arb("175.68") / (u_lo * u_lo)
    c_g = arb(1)
    mom1 = eta0 / 2
    sum_inv_gamma = (A / (2 * pi)).log() ** 2 / (4 * pi) + arb("0.9321")
    term_y = 2 * mom1 * (c_f * sum_inv_gamma + c_g * SUM_INV_GAMMA2) / (u_lo * u_lo)
    term_li = 4 * SUM_INV_GAMMA3 / (u_lo * u_lo)
    return (term_y + term_li).upper()


def bound_R4_envelope(omega, eta1, T0, lam1, lam2, lam3, n_sub=8192, xi_tail_val=160):
    """
    Conservative R4 bound: applies Lemma 2.7 of STD (2015) directly to the nonincreasing
    step-envelope majorant Psi(eta1 * t) of |hat{kappa}_core(eta1 * t)|.
    """
    pi = arb.pi()
    xi_T = eta1 * T0
    require(xi_T > 3 * pi, "xi_T = eta1 * T0 must exceed 3*pi")
    xi_tail = arb(xi_tail_val)
    require(xi_tail > xi_T, "xi_tail must exceed xi_T")

    h = (xi_tail - xi_T) / n_sub
    a_list = [xi_T + j * h for j in range(n_sub)]
    b_list = [xi_T + (j + 1) * h for j in range(n_sub)]
    u_bounds = []
    for a, b in zip(a_list, b_list):
        mid = (a + b) / 2
        rad = (b - a) / 2
        iv = arb(mid, rad.upper())
        r_val = kappa_core_rational_factor_arb(iv, lam1, lam2, lam3)
        val = (abs(iv.sin()) * r_val).upper()
        require(val.is_finite(), "non-finite kappa_core_hat enclosure on subinterval")
        u_bounds.append(val)

    tail_at_end = phi_core_tail_majorant_arb(xi_tail, lam1, lam2, lam3).upper()
    M_env = [None] * n_sub
    cur = tail_at_end
    for j in range(n_sub - 1, -1, -1):
        cur = cur if cur > u_bounds[j] else u_bounds[j]
        M_env[j] = cur

    int_main = arb(0)
    int_rem = arb(0)
    for j in range(n_sub):
        t_a = a_list[j] / eta1
        t_b = b_list[j] / eta1
        fac = M_env[j] * (1 + 1 / (t_a * omega))
        la = (t_a / (2 * pi)).log()
        lb = (t_b / (2 * pi)).log()
        int_main += fac * (lb * lb - la * la) / (2 * pi)
        int_rem += 4 * fac * (1 / t_a - 1 / t_b)

    T_tail = xi_tail / eta1
    fac_tail = tail_at_end * (1 + 1 / (T_tail * omega))
    lt = (T_tail / (2 * pi)).log()
    int_main += fac_tail * (3 * lt + 1) / (9 * pi)
    int_rem += 4 * fac_tail / (4 * T_tail)

    phi_T0 = M_env[0] * (1 / T0 + 1 / (T0 * T0 * omega))
    term_b = 8 * phi_T0 * T0.log()
    return (int_main + int_rem + term_b).upper()


def bound_R4_stieltjes(omega, eta1, T0, lam1, lam2, lam3, n_sub=8192, xi_tail_val=160):
    """
    Sharp R4 bound via Stieltjes integration by parts on phi(t) = |sin(eta1*t)| * g(t),
    where g(t) = H_env(eta1*t) * (1/t + 1/(t^2*omega)) is positive and nonincreasing, and
    H_env(xi) is the nonincreasing step-envelope majorant of |R_core(xi)|.
    """
    pi = arb.pi()
    xi_T = eta1 * T0
    require(xi_T > 3 * pi, "xi_T = eta1 * T0 must exceed 3*pi")
    xi_tail = arb(xi_tail_val)
    require(xi_tail > xi_T, "xi_tail must exceed xi_T")

    h = (xi_tail - xi_T) / n_sub
    a_list = [xi_T + j * h for j in range(n_sub)]
    b_list = [xi_T + (j + 1) * h for j in range(n_sub)]
    r_bounds = []
    kh_bounds = []
    for a, b in zip(a_list, b_list):
        mid = (a + b) / 2
        rad = (b - a) / 2
        iv = arb(mid, rad.upper())
        r_arb = kappa_core_rational_factor_arb(iv, lam1, lam2, lam3)
        r_val = r_arb.upper()
        kh_val = (abs(iv.sin()) * r_arb).upper()
        require(r_val.is_finite() and kh_val.is_finite(), "non-finite enclosure on subinterval")
        r_bounds.append(r_val)
        kh_bounds.append(kh_val)

    tail_at_end = phi_core_tail_majorant_arb(xi_tail, lam1, lam2, lam3).upper()
    H_env = [None] * n_sub
    cur = tail_at_end
    for j in range(n_sub - 1, -1, -1):
        cur = cur if cur > r_bounds[j] else r_bounds[j]
        H_env[j] = cur

    int_main = arb(0)
    int_rem = arb(0)
    int_tv_sin = arb(0)
    for j in range(n_sub):
        t_a = a_list[j] / eta1
        t_b = b_list[j] / eta1
        fac_kh = kh_bounds[j] * (1 + 1 / (t_a * omega))
        fac_g = H_env[j] * (1 + 1 / (t_a * omega))
        la = (t_a / (2 * pi)).log()
        lb = (t_b / (2 * pi)).log()
        int_main += fac_kh * (lb * lb - la * la) / (2 * pi)
        int_rem += 4 * fac_g * (1 / t_a - 1 / t_b)
        int_tv_sin += 4 * eta1 * fac_g * (t_b - t_a) * (t_b.log() / t_a)

    T_tail = xi_tail / eta1
    fac_tail = tail_at_end * (1 + 1 / (T_tail * omega))
    lt = (T_tail / (2 * pi)).log()
    int_main += fac_tail * (3 * lt + 1) / (9 * pi)
    int_rem += 4 * fac_tail / (4 * T_tail)
    # Tail of 4 * eta1 * \int_{T_tail}^\infty g(t) log(t) dt:
    # Since g(t) <= fac_tail * (T_tail^3 / t^4), \int_{T_tail}^\infty t^{-4} log(t) dt = T_tail^{-3} (3 log T_tail + 1) / 9
    int_tv_sin += 4 * eta1 * fac_tail * (3 * T_tail.log() + 1) / 9

    g_T0 = H_env[0] * (1 / T0 + 1 / (T0 * T0 * omega))
    term_b = 8 * g_T0 * T0.log()
    return (int_main + int_rem + int_tv_sin + term_b).upper()


def bound_R4(omega, eta1, T0, lam1, lam2, lam3, n_sub=8192, xi_tail_val=160):
    """Return min(bound_R4_envelope, bound_R4_stieltjes) -- both are rigorous upper bounds."""
    b_env = bound_R4_envelope(omega, eta1, T0, lam1, lam2, lam3, n_sub=n_sub, xi_tail_val=xi_tail_val)
    b_st = bound_R4_stieltjes(omega, eta1, T0, lam1, lam2, lam3, n_sub=n_sub, xi_tail_val=xi_tail_val)
    return b_st if b_st < b_env else b_env


def bound_R5(omega, eta0, delta, M_moll=60, A=A_PT):
    """
    Rigorous upper bound for R5 (off-line zeros with |gamma| > A = 3*10^12) via N = 2*M_moll
    integrations by parts on [omega - eta0, omega + eta0].
    """
    N = 2 * M_moll
    pi = arb.pi()
    u_lo = omega - eta0
    u_hi = omega + eta0
    R_cauchy = arb(300)
    require(u_lo - R_cauchy > 100, "u_lo - R_cauchy > 100 fails")
    C_M = arb(4**M_moll * (math.factorial(M_moll) ** 2)) / arb(2 * math.factorial(N))
    rate = (M_moll * pi) / (A * delta) + arb(N) / (A * R_cauchy)
    require(rate < 1, "R5 rate must be < 1")
    log_R5 = (arb(8) * C_M * A.log()).log() + u_hi / 2 + N * rate.log()
    return log_R5.exp().upper()


def total_error(omega, eta1, delta, T0, lam1, lam2, lam3, M_moll=60, verbose=True):
    """Compute all error bounds R2, R3, R4, R5 (with R1 = 0) and return (R_total, parts_dict)."""
    require(lam1 >= 0 and lam2 >= 0 and lam3 >= 0, "convex weights lam_i must be >= 0")
    require((lam1 + lam2 + lam3).contains(1), "lam1 + lam2 + lam3 must enclose 1")
    eta0 = eta1 + delta
    r2 = bound_R2(omega, eta0)
    r3 = bound_R3(omega, eta0)
    r4_env = bound_R4_envelope(omega, eta1, T0, lam1, lam2, lam3)
    r4_st = bound_R4_stieltjes(omega, eta1, T0, lam1, lam2, lam3)
    r4 = r4_st if r4_st < r4_env else r4_env
    r5 = bound_R5(omega, eta0, delta, M_moll=M_moll)
    r_tot = (r2 + r3 + r4 + r5).upper()
    parts = {
        "R1": arb(0), "R2": r2, "R3": r3,
        "R4_envelope": r4_env, "R4_stieltjes": r4_st, "R4": r4,
        "R5": r5, "R_total": r_tot,
    }
    if verbose:
        for k in ["R1", "R2", "R3", "R4_envelope", "R4_stieltjes", "R4", "R5", "R_total"]:
            print(f"  {k:14s} <= {parts[k].str(10, radius=False)}")
    return r_tot, parts


if __name__ == "__main__":
    omega = arb("727.95133580")
    eta1 = arb("1.55e-6")
    delta = arb("1e-8")
    T0 = arb("6800000")
    lam1, lam2, lam3 = arb(12) / 100, arb(70) / 100, arb(18) / 100
    print(f"Testing error bounds at omega={omega.str(12, radius=False)}, eta1={eta1.str(6, radius=False)}, T0={T0.str(10, radius=False)}:")
    total_error(omega, eta1, delta, T0, lam1, lam2, lam3)
