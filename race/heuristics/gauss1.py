"""Gaussian (Lehman/Ford-Hudson) kernel at one region with zeros to T0: E_G min + truncation tail + window floor."""
import sys, math, numpy as np, race, feas
q, a, b = map(int, sys.argv[1:4]); suffix = sys.argv[4] if sys.argv[4] != '-' else ''; uc = float(sys.argv[5]); T0 = float(sys.argv[6])
Aver = float(sys.argv[7])
_, labs = feas.race_A(q, a, b)
labs_c = [(l, race.chival(l, a) - race.chival(l, b)) for l in labs]
c0 = -(race.nsqrt(q, a) - race.nsqrt(q, b))
for alpha in [1e7, 3e7, 1e8, 2e8, 4e8, 1e9]:
    # truncation tail of on-line zeros above T0 (FH R3-like): sum_{T0<g} 2|c| e^{-g^2/2a}/g dN
    t = np.linspace(T0, T0 + 12*math.sqrt(alpha), 200000)
    tail = sum(np.trapezoid(2*abs(c)*np.exp(-t**2/(2*alpha))/t*np.log(feas.COND[l]*t/(2*np.pi))/(2*np.pi), t) for l, c in labs_c)
    # window: need A_eff with A_eff^2/(2 alpha) >= uc/2 + 30 and eta >= 2 A_eff/alpha, A_eff <= Aver
    A_eff = math.sqrt(2*alpha*(uc/2 + 30))
    if A_eff > Aver:
        print(f"alpha={alpha:.1e}: needs A_eff={A_eff:.3g} > verified {Aver:.3g}"); continue
    eta = 2*A_eff/alpha
    gs, As_, Bs = [], [], []
    for lab, c in labs_c:
        g = race.zeros(lab, suffix); g = g[g <= T0]
        w = np.exp(-g**2/(2*alpha)); rho = 0.5 + 1j*g
        gs.append(g); As_.append(-c*2*w/rho); Bs.append(-c*2*w/rho**2)
    g = np.concatenate(gs); Ac = np.concatenate(As_); Bc = np.concatenate(Bs)
    bmin, bu = race.scan(c0, g, Ac, Bc, uc-0.001, uc+0.001, 1/(40*math.sqrt(alpha)), blk=64)
    k = bmin.argmin()
    print(f"alpha={alpha:.1e}: eta={eta:.2e} A_eff={A_eff:.3g}: E_G={bmin[k]:+.4f} tail={tail:.2e} crit={bmin[k]+tail:+.4f} at u={bu[k]:.6f}", flush=True)
