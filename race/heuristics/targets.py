"""Coarse JP(T=Tc) scan of E(u) for u in [u0, umax], then full-weight T=2e4 refinement of candidate blocks."""
import sys, time, math, json
import numpy as np
import race

RACES = [tuple(map(int, r.split(','))) for r in sys.argv[1].split(';')]
umax = float(sys.argv[2]); Tc = float(sys.argv[3]); thr = float(sys.argv[4])
Tf = 2e4
NNEG = int(sys.argv[6]) if len(sys.argv) > 6 else 8
out = {}
for (q, a, b) in RACES:
    t0 = time.time()
    c0, g, A, B = race.coeffs(q, a, b, 'jp', Tc)
    h = 1/(8*Tc)
    bmin, bu = race.scan(c0, g, A, B, 5.0, umax, h, blk=2048)
    blen = 2048*h
    cand = np.nonzero(bmin < thr)[0]
    print(f"race {q};{a},{b} c0={c0}: coarse JP T={Tc:g}: global min {bmin.min():.3f} at u={bu[bmin.argmin()]:.3f}; "
          f"{len(cand)} candidate blocks below {thr} [{time.time()-t0:.0f}s]", flush=True)
    cf0, gf, Af, Bf = race.coeffs(q, a, b, 'full', Tf)
    cj0, gj, Aj, Bj = race.coeffs(q, a, b, 'jp', Tf)
    rec = []
    nneg = 0
    for k in cand:
        ua = 5.0 + k*blen - 0.002; ub = ua + blen + 0.004
        fm, fu = race.scan(cf0, gf, Af, Bf, ua, ub, 0.1/Tf, blk=256)
        jm, ju = race.scan(cj0, gj, Aj, Bj, ua, ub, 0.1/Tf, blk=256)
        rec.append(dict(u_coarse=float(bu[k]), Ecoarse=float(bmin[k]), Efull=float(fm.min()), u_full=float(fu[fm.argmin()]),
                        Ejp=float(jm.min()), u_jp=float(ju[jm.argmin()])))
        if rec[-1]['Efull'] < 0:
            nneg += 1
            print('   neg region', {kk: round(v, 5) for kk, v in rec[-1].items()}, flush=True)
        if nneg >= NNEG:
            break
    rec.sort(key=lambda r: r['u_coarse'])
    firstneg = next((r for r in rec if r['Efull'] < 0), None)
    firstjp = next((r for r in rec if r['Ejp'] < 0), None)
    print(f"   first full-weight T=2e4 negative: {firstneg}", flush=True)
    print(f"   first JP T=2e4 negative:          {firstjp}", flush=True)
    for r in rec[:12]:
        print("    ", {k: round(v, 5) for k, v in r.items()}, flush=True)
    out[f"{q};{a},{b}"] = dict(c0=c0, coarse_min=float(bmin.min()), coarse_argmin=float(bu[bmin.argmin()]), cands=rec)
    json.dump(out, open(f"targets_{sys.argv[5]}.json", "w"), indent=1)
