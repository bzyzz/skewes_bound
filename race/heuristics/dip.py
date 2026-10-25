"""Depth/width of a candidate region vs zero cutoff T (full weight and JP)."""
import sys, numpy as np, race
q, a, b = map(int, sys.argv[1:4]); uc = float(sys.argv[5]); half = float(sys.argv[6]) if len(sys.argv) > 6 else 0.003
suffix = sys.argv[4] if sys.argv[4] != '-' else ''
Ts = [float(t) for t in sys.argv[7].split(',')] if len(sys.argv) > 7 else [2.5e3, 5e3, 1e4, 2e4]
for kind in ['full', 'jp']:
    for T in Ts:
        c0, g, A, B = race.coeffs(q, a, b, kind, T, suffix)
        h = 0.05/T
        _, _, f = race.scan(c0, g, A, B, uc-half, uc+half, h, blk=64, full=True)
        u = uc-half + h*np.arange(len(f))
        k = f.argmin(); neg = f < 0
        w = neg.sum()*h
        # contiguous width around argmin
        l = k
        while l > 0 and f[l-1] < 0: l -= 1
        r = k
        while r < len(f)-1 and f[r+1] < 0: r += 1
        cw = (r-l+1)*h if f[k] < 0 else 0
        print(f"{kind:4s} T={T:8.0f}: min E={f[k]:+.4f} at u={u[k]:.6f}; total neg measure={w:.2e}; contiguous neg width at min={cw:.2e}", flush=True)
