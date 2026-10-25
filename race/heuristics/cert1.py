"""Two-scale criterion at a single region for several hypothetical verified heights A (heuristic)."""
import sys, math, numpy as np, race, feas
q, a, b = map(int, sys.argv[1:4]); suffix = sys.argv[4] if sys.argv[4] != '-' else ''; uc = float(sys.argv[5]); T0 = float(sys.argv[6])
As = [float(x) for x in sys.argv[7].split(',')]
_, labs = feas.race_A(q, a, b)
labs_c = [(l, race.chival(l, a) - race.chival(l, b)) for l in labs]
c0 = -(race.nsqrt(q, a) - race.nsqrt(q, b))
for A in As:
    M, delta = feas.two_scale_params(uc, A)
    best = None
    for m in (1, 2, 3):
        for f in (6, 8, 10, 12, 14, 20):
            eta1 = f/T0
            R4 = feas.R4_bound(labs_c, eta1, m, T0, max(A, T0*1.01), M, delta)
            Rm = 1.3*math.exp(-uc/6) + abs(c0)*10.6/uc**3
            gs, As_, Bs = [], [], []
            for lab, c in labs_c:
                g = race.zeros(lab, suffix); g = g[g <= T0]
                w = feas.khat(eta1*g, m)*feas.khat(delta*g, M)
                rho = 0.5 + 1j*g
                gs.append(g); As_.append(-c*2*w/rho); Bs.append(-c*2*w/rho**2)
            g = np.concatenate(gs); Ac = np.concatenate(As_); Bc = np.concatenate(Bs)
            h = min(eta1, delta)/20
            bmin, bu = race.scan(c0, g, Ac, Bc, uc-0.001, uc+0.001, h, blk=64)
            k = bmin.argmin(); crit = bmin[k] + R4 + Rm
            if best is None or crit < best[0]:
                best = (crit, bmin[k], R4, Rm, bu[k], m, eta1)
    crit, EK, R4, Rm, om, m, eta1 = best
    print(f"A={A:.3g}: M={M} delta={delta:.2e} | best m={m} eta1={eta1:.2e} (eta1*T0={eta1*T0:.0f}): E_K={EK:+.4f} R4={R4:.2e} Rm={Rm:.1e} "
          f"crit={crit:+.4f} {'CERTIFIABLE x<exp(%.6f)=10^%.3f' % (om+eta1+delta, (om+eta1+delta)/math.log(10)) if crit < 0 else 'no'}", flush=True)
