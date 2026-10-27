# Validate u e^{-u/2} G(e^u) ~ -(phi-1) Phi(u) - sum_chi 2 Re sum_{gamma>0} K(gamma) e^{i gamma u}(1/rho + 1/(rho^2 u)) by smoothing.
import math, os, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CH = {8: ["4_3", "8_3", "8_5"], 12: ["3_2", "4_3", "12_11"]}
for q in (8, 12):
    for omega in (16.5, 17.3, 18.1):
        eta = 0.02; lo, hi = omega - eta, omega + eta; X = int(math.exp(hi)) + 2
        sv = np.ones(X + 1, bool); sv[:2] = False
        for i in range(2, int(X**0.5) + 1):
            if sv[i]: sv[i*i::i] = False
        P = np.nonzero(sv)[0]; P = P[np.gcd(P, q) == 1]
        w = np.where(P % q == 1, 3.0, -1.0)            # phi*1[p = 1] - 1  (phi = 4)
        lp = np.log(P.astype(float))
        G0 = w[lp <= lo].sum(); jl = lp[lp > lo]; jc = w[lp > lo]
        C3 = 4**3*math.factorial(3)**2/(2*math.factorial(6))
        grid = np.linspace(lo, hi, 400001); K = C3*np.cos(np.pi*(grid - omega)/(2*eta))**6/eta
        cs = np.concatenate([[0], np.cumsum(jc)]); G = G0 + cs[np.searchsorted(jl, grid, side="right")]
        lhs = np.trapezoid(K*grid*np.exp(-grid/2)*G, grid)
        rhs = -3*np.trapezoid(K*(1 + 2/grid + 8/grid**2 + 48/grid**3), grid)
        for c in CH[q]:
            d = np.load(os.path.join(ROOT, "data", f"arb_{c}", "certified.npz")); g = 0.5*(d["lo"] + d["hi"]); g = g[g <= 2e4]
            rho = 0.5 + 1j*g; xi = eta*g; kh = np.sinc(xi/np.pi)
            for k in (1, 2, 3): kh = kh*(k*np.pi)**2/((k*np.pi)**2 - xi*xi)
            rhs -= 2*np.sum((kh*np.exp(1j*g*omega)*(1/rho + 1/(rho**2*omega))).real)
        print(f"q={q} omega={omega}: smoothed u e^(-u/2) G from primes = {lhs:+.4f}; bias + zero formula = {rhs:+.4f}; diff {lhs - rhs:+.4f}", flush=True)
