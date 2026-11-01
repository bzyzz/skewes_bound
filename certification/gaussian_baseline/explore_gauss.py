# Baseline (non-rigorous exploration): Platt-Trudgian 2016 Theorem 3 + Sec. 3.1.6 trimming, Gaussian kernel,
# on the same zeros (|gamma| <= T0). Reports the smallest u_hi = omega + eta0 with a positive trimmed margin.
import math, sys
import numpy as np
exec(open(__file__.replace("explore_gauss.py", "explore.py")).read().split("if __name__")[0])

def bins_gauss(T0, alpha):
    g = G[G <= T0]; n = len(g) // BIN * BIN
    gf = g.astype(float)
    w = np.exp(-gf * gf / (2 * alpha))
    ph = np.fmod(g * OM0, TWO_PI).astype(float)
    a = w * np.exp(1j * ph) / (0.5 + 1j * gf)
    tail = (g[n:], a[n:])
    gb = g[:n].reshape(-1, BIN); ab = a[:n].reshape(-1, BIN)
    c = gb.mean(axis=1); dl = (gb - c[:, None]).astype(float)
    B = np.stack([(ab * (1j * dl) ** k).sum(axis=1) / math.factorial(k) for k in range(KT)])
    return c.astype(float), B, tail

def errs(T, alpha, omega):
    # choose the smallest A (<= 3e12) with R4 <= 1e-10; then eta = 2A/alpha (PT16 Thm 3 conditions)
    A = math.sqrt(2 * alpha * (omega / 2 + 60))
    for _ in range(50):
        eta = 2 * A / alpha
        lr4 = math.log(A * math.log(A) * (4 / math.sqrt(alpha) + 15 * eta)) - A * A / (2 * alpha) + (omega + eta) / 2
        if lr4 < math.log(1e-10): break
        A *= 1.01
    assert A <= 3e12 and 4 * A / omega <= alpha <= A * A
    eta = 2 * A / alpha
    R1 = 1.5423e-9
    R2 = 0.08 * math.sqrt(alpha) * math.exp(-alpha * eta**2 / 2) + math.exp(-T * T / (2 * alpha)) * (
        alpha / (math.pi * T * T) * math.log(T / (2 * math.pi)) + 8 * math.log(T) / T + 4 * alpha / T**3)
    R3 = math.exp(-(omega - eta) / 2) * math.log(2 * math.pi) + 3 * math.exp(-(omega - eta) / 6)
    return eta, R1 + R2 + R3 + math.exp(lr4), A

def trim(J, omega, eta, alpha):
    # smallest eta0 in (0, eta] with 1.3082e-9 (eta - eta0) K(eta0) (e^{(omega+eta)/2} + e^{(omega-eta0)/2}) < J
    def lhs(e0):
        return math.log(1.3082e-9 * (eta - e0) * math.sqrt(alpha / (2 * math.pi))) - alpha * e0 * e0 / 2 + \
               math.log(math.exp((omega + eta) / 2 - omega / 2) + math.exp((omega - e0) / 2 - omega / 2)) + omega / 2
    lo, hi = 0.0, eta
    if lhs(1e-30) < math.log(J): return 0.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if lhs(mid) < math.log(J): hi = mid
        else: lo = mid
    return hi

if __name__ == "__main__":
    T0 = float(sys.argv[1])
    d = np.arange(-2.2e-6, 2.0e-6, 2e-10)
    best = None
    for c in (3.0, 3.25, 3.5, 3.75, 4.0, 4.25, 4.5):
        alpha = (T0 / c) ** 2
        eta, R, A = errs(T0, alpha, 727.95)
        bn = bins_gauss(T0, alpha)
        J = -1 + S_theta(bn, d) - R
        for i in np.where(J > 0)[0]:
            om = float(OM0) + d[i]
            e0 = trim(J[i], om, eta, alpha)
            u = om + e0
            if best is None or u < best[0]:
                best = (u, c, alpha, eta, e0, om, J[i], A)
        print(f"T0={T0:.1e} T/sqrt(alpha)={c}: eta={eta:.3e} (A={A:.3e}), R={R:.2e}, best so far {best}", flush=True)
    u, c, alpha, eta, e0, om, j, A = best
    print(f"RESULT T0={T0:.1e}: Gaussian/PT16 method u_hi={u:.10f} (omega={om:.10f}, alpha={alpha:.3e}, eta={eta:.3e}, "
          f"trimmed eta0={e0:.3e}, eta0*T0={e0*T0:.1f}, J={j:.2e})")
