import sys, time, math
import numpy as np
import race
KNOWN = [((4, 3, 1), 26861), ((8, 5, 1), 588067889), ((3, 2, 1), 608981813029), ((24, 13, 1), 9.78e11)]
T = float(sys.argv[1]) if len(sys.argv) > 1 else 2e4
for (q, a, b), x in KNOWN:
    for kind in ['full', 'jp']:
        t0 = time.time()
        c0, g, A, B = race.coeffs(q, a, b, kind, T)
        h = 0.1/T
        bmin, bu = race.scan(c0, g, A, B, 5.0, 40.0, h, blk=1024)
        neg = np.nonzero(bmin < 0)[0]
        first = bu[neg[0]] if len(neg) else None
        print(f"race {q};{a},{b} kind={kind} T={T:g} nzeros={len(g)} c0={c0}: first E<0 at u={first} "
              f"(x={math.exp(first) if first else 0:.4g}); known u={math.log(x):.4f} (x={x:g}); "
              f"min E on [5,40] = {bmin.min():.3f} at u={bu[bmin.argmin()]:.4f}  [{time.time()-t0:.0f}s]", flush=True)
