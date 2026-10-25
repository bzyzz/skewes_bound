"""For candidate regions of a race, test (heuristically, floating point) whether the two-scale criterion
  E_K(omega) + R4 + R_misc < 0
holds with the zeros we have (T0 = 2e4) and the Platt verified height A, and report the bound exp(omega + eta0).
E_K uses weights hatK(gamma) = khat_m(eta1 gamma) khat_M(delta gamma), delta at the IBP floor.
R_misc = 1.3 e^{-omega/6} + N(1)*10.6/omega^3 (prime powers k >= 3, Dusart 7.32 vs asymptotic 6 in li(sqrt x)).
usage: certscan.py q a b json_from_targets [T0]
"""
import sys, json, math
import numpy as np
import race, feas

q, a, b = map(int, sys.argv[1:4])
recs = json.load(open(sys.argv[4]))[f"{q};{a},{b}"]["cands"]
T0 = float(sys.argv[5]) if len(sys.argv) > 5 else 2e4
A, labs = feas.race_A(q, a, b)
labs_c = [(l, race.chival(l, a) - race.chival(l, b)) for l in labs]
c0 = -(race.nsqrt(q, a) - race.nsqrt(q, b))
regions = sorted([r for r in recs if r['Efull'] < 0.5], key=lambda r: r['u_full'])
print(f"race {q};{a},{b}: A={A:.3g} chars={labs}; {len(regions)} regions with full-weight min < 0.5")
best = None
for r in regions:
    uc = r['u_full']
    M, delta = feas.two_scale_params(uc, A)
    res = []
    for m in (2, 3):
        for f in (10, 14, 20, 30):
            eta1 = f/T0
            R4 = feas.R4_bound(labs_c, eta1, m, T0, A, M, delta)
            Rm = 1.3*math.exp(-uc/6) + abs(c0)*10.6/uc**3
            gs, As, Bs = [], [], []
            for lab, c in labs_c:
                g = race.zeros(lab); g = g[g <= T0]
                w = feas.khat(eta1*g, m)*feas.khat(delta*g, M)
                rho = 0.5 + 1j*g
                gs.append(g); As.append(-c*2*w/rho); Bs.append(-c*2*w/rho**2)
            g = np.concatenate(gs); Ac = np.concatenate(As); Bc = np.concatenate(Bs)
            h = min(eta1, delta)/20
            bmin, bu = race.scan(c0, g, Ac, Bc, uc-0.02, uc+0.02, h, blk=64)
            k = bmin.argmin()
            crit = bmin[k] + R4 + Rm
            res.append((crit, bmin[k], R4, Rm, bu[k], m, eta1))
    res.sort()
    crit, EK, R4, Rm, om, m, eta1 = res[0]
    eta0 = eta1 + delta
    ok = crit < 0
    print(f"  u_full={uc:.5f} Efull={r['Efull']:+.3f} Ejp={r['Ejp']:+.3f} | best m={m} eta1={eta1:.2e} delta={delta:.2e} M={M}: "
          f"E_K={EK:+.4f} R4={R4:.2e} Rm={Rm:.1e} crit={crit:+.4f} {'CERTIFIABLE' if ok else 'no'}"
          + (f" -> x < exp({om+eta0:.6f}) = 10^{(om+eta0)/math.log(10):.4f}" if ok else ""), flush=True)
    if ok and best is None:
        best = (om+eta0, r)
if best:
    print(f"EARLIEST CERTIFIABLE (T0={T0:g}): x < exp({best[0]:.6f}) = 10^{best[0]/math.log(10):.4f}")
else:
    print("no certifiable region among candidates")
