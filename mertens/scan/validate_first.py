# Validate sqrt(x) Delta_1(x) ~ 1 - 2 Re sum_{gamma>0} x^{i gamma}/(rho - 1) by smoothed comparison with primes.
import math, os, sys
import numpy as np
E = -1.3325822757332208817658287760710277
G = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "gammas_2325572.npy"))
for omega in (16.5, 17.3, 18.1):
    eta = 0.02; lo, hi = omega - eta, omega + eta
    X = int(math.exp(hi)) + 2
    sv = np.ones(X + 1, bool); sv[:2] = False
    for i in range(2, int(X**0.5) + 1):
        if sv[i]: sv[i*i::i] = False
    P = np.nonzero(sv)[0].astype(np.float64)
    C3 = 4**3*math.factorial(3)**2/(2*math.factorial(6))
    grid = np.linspace(lo, hi, 400001)
    K = C3*np.cos(np.pi*(grid - omega)/(2*eta))**6/eta
    lp = np.log(P); terms = lp/P
    S0 = terms[lp <= lo].sum(); jl = lp[lp > lo]; jc = terms[lp > lo]
    # F(u) = sqrt(x)(S(x) - u - E), S(x) = S0 + sum of jumps below u
    cum = np.searchsorted(jl, grid, side="right"); cs = np.concatenate([[0], np.cumsum(jc)])
    F = np.exp(grid/2)*(S0 + cs[cum] - grid - E)
    lhs = np.trapezoid(K*F, grid)
    g = G[G <= 1.3e6]; rho = 0.5 + 1j*g; xi = eta*g
    kh = np.sinc(xi/np.pi)
    for k in (1, 2, 3): kh = kh*(k*np.pi)**2/((k*np.pi)**2 - xi*xi)
    rhs = 1 - 2*np.sum((kh*np.exp(1j*g*omega)/(rho - 1)).real)
    print(f"omega={omega}: smoothed sqrt(x)Delta_1 from primes = {lhs:+.5f}; explicit formula 1 - 2Re sum = {rhs:+.5f}; diff {lhs - rhs:+.5f}")
