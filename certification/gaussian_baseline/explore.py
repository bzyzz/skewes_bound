# Exploration (non-rigorous, long double): smallest certifiable endpoint for theta(x) > x with the two-scale kernel.
# J_theta(omega) = -1 - 1.5423e-9 + S_theta(omega) - R4 - R5,  S_theta = -sum_{|g|<=T0} Khat(g) e^{i g omega}/rho.
import math, numpy as np, sys, time
L = np.longdouble
PI = L("3.14159265358979323846264338327950288"); TWO_PI = 2 * PI
S = "/tmp/claude-28261/-data-wang-zhenyub/6e51b187-df41-49b4-8d3b-37362fb380c5/scratchpad"
z = np.load(f"{S}/indep/zeros.npz"); G = np.concatenate([z["gi"].astype(L) + z["gf"].astype(L), np.load(f"{S}/indep2/g_hi_ld.npy")])
OM0 = L("727.951334"); KT = 5; BIN = 4096
def kh(xi, m):
    r = np.sin(xi) / xi
    for k in range(1, m + 1): r = r * (k*np.pi)**2 / ((k*np.pi)**2 - xi*xi)
    return r
def bins(T0, eta1, lam, delta=2e-9, M=60):
    g = G[G <= T0]; n = len(g) // BIN * BIN
    w = sum(l/100 * kh(float(eta1) * g.astype(float), m) for l, m in zip(lam, (1, 2, 3))) * kh(delta * g.astype(float), M)
    ph = np.fmod(g * OM0, TWO_PI).astype(float)
    a = w * np.exp(1j * ph) / (0.5 + 1j * g.astype(float))          # Khat e^{i g om0} / rho
    tail = (g[n:], a[n:])
    gb = g[:n].reshape(-1, BIN); ab = a[:n].reshape(-1, BIN)
    c = gb.mean(axis=1); dl = (gb - c[:, None]).astype(float)
    B = np.stack([(ab * (1j * dl) ** k).sum(axis=1) / math.factorial(k) for k in range(KT)])
    return c.astype(float), B, tail
def S_theta(bn, d):
    c, B, (gt, at) = bn
    poly = sum(B[k][:, None] * d[None, :] ** k for k in range(KT))
    s = (np.exp(1j * np.outer(c, d)) * poly).sum(axis=0) + (at[:, None] * np.exp(1j * np.outer(gt.astype(float), d))).sum(axis=0)
    return -2 * s.real
def R4_st(eta1, T0, lam, h=5e-5, xmax=600.0):
    l = [x/100 for x in lam]; xi = np.arange(eta1*T0, xmax, h); t = xi/eta1
    kc = l[0]*kh(xi,1) + l[1]*kh(xi,2) + l[2]*kh(xi,3)
    l1, l2 = l[0], l[1]
    rc = np.pi**2*np.abs(l1*xi**4 - np.pi**2*(13*l1+4*l2)*xi**2 + 36*np.pi**4)/(xi*(xi**2-np.pi**2)*(xi**2-4*np.pi**2)*(xi**2-9*np.pi**2))
    H = np.maximum.accumulate(rc[::-1])[::-1]; g = H/t; tr = np.trapezoid
    return tr(np.abs(kc)/t*np.log(t/(2*np.pi)), t)/np.pi + 8*g[0]*np.log(T0) + 4*tr(g/t, t) + 4*eta1*tr(g*np.log(t), t)
if __name__ == "__main__":
    d = np.arange(-2.2e-6, 2.0e-6, 2e-10)
    best = {}
    for T0 in [float(sys.argv[1])]:
        for lam in ((4, 52, 44), (14, 72, 14), (0, 50, 50), (0, 20, 80), (0, 0, 100)):
            for c_eta in (10.0, 11.0, 12.0, 13.5, 15.0, 17.0):
                eta1 = c_eta / T0
                if eta1 * T0 <= 3*np.pi: continue
                bn = bins(T0, eta1, lam); r4 = R4_st(eta1, T0, lam) * 1.04   # 4% cushion: repo bound ~3% above this estimate
                J = -1 - 1.5423e-9 + S_theta(bn, d) - r4
                pos = np.where(J > 0)[0]
                if len(pos) == 0: continue
                i = pos[0]; uhi = float(OM0) + d[i] + eta1 + 2e-9
                key = T0
                if key not in best or uhi < best[key][0]:
                    best[key] = (uhi, lam, eta1, float(OM0) + d[i], J[i], r4)
            print(f"T0={T0:.1e} lam={lam} done; best so far: {best.get(T0)}", flush=True)
    print("\nRESULT (smallest certifiable u_hi = omega + eta0 near 727.95133):")
    for T0, (u, lam, e1, om, j, r4) in sorted(best.items()):
        print(f"  T0={T0:.1e}: u_hi={u:.9f}  (omega={om:.10f}, eta1={e1:.3e}, lam={lam}, R4~{r4:.2e})   vs PT16 727.951332668")
