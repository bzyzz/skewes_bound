#!/usr/bin/env python3
"""Arb isolation of all zeros of zeta(s) up to T0 = 8e7 and per-batch Arb partial sums of S_K.

Low zeros (n <= 13,955,181) are read from the Arb zero caches (see README, Data). High zeros
(13,955,181 < n <= 195,565,267) are isolated here with acb.zeta_zeros, which certifies that each lies on the
critical line, and are saved as exact balls in zero_cache_arb/z_<n0>.npz (also used by theta_certify.py).
Each batch writes arb_batches/<l|h>_<n0>.json with its first/last zero and the Arb partial sums of S_K for
every certificate in SETS. The run is resumable. Combine with arb_combine.py.

usage: arb_recert.py [workers]     (about 15 hours with 64 workers; 4.2 GB of zero files)
Environment: MERTENS_ZERO_CACHE (Arb zeros n <= 2,325,572; default ../../mertens_bound/certification/zero_cache).
"""
import glob, json, os, sys, time
os.environ["OMP_NUM_THREADS"] = "1"
from multiprocessing import Pool
import flint
from flint import acb, arb
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(HERE, "arb_batches")
ZD = os.environ.get("ARB_ZERO_DIR", os.path.join(ROOT, "zero_cache_arb"))
CERT = HERE
MER = os.environ.get("MERTENS_ZERO_CACHE", os.path.join(ROOT, "..", "mertens_bound", "certification", "zero_cache"))
EXT = os.path.join(ROOT, "zero_cache_ext")
N_LOW, N_5E7, N_8E7 = 13955181, 118488123, 195565267
DPS = 20
SETS = {}
for c in ["s_cert_T5e7", "s_certificate", "s_cert_T8e7_st"]:
    d = json.load(open(f"{CERT}/{c}.json"))
    SETS[c] = dict(omega=d["omega"], eta1=d["eta1"], delta=d["delta"], lam=[int(x) for x in d["lam"].split(",")],
                   M=int(d["M_moll"]), nmax=d["N_T0"])
def khat(xi, m):
    # sin(xi)/xi * prod_k (k pi)^2/((k pi)^2 - xi^2), removable singularity at xi = k pi handled via sinc
    pi = arb.pi(); res = None
    for k in range(1, m + 1):
        if (xi - k * pi).contains(0):
            u = (xi - k * pi) / pi                      # sin(xi)/(k pi - xi) = (-1)^(k+1) sinc(pi u)
            res = (-1) ** (k + 1) * u.sinc_pi() * (k * pi) ** 2 / (xi * (k * pi + xi))
            for j in range(1, m + 1):
                if j != k: res *= (j * pi) ** 2 / ((j * pi) ** 2 - xi * xi)
            return res
    res = xi.sin() / xi
    for k in range(1, m + 1): res *= (k * pi) ** 2 / ((k * pi) ** 2 - xi * xi)
    return res
def batch_sum(gs, n0):
    flint.ctx.dps = DPS + 10
    out = {}
    for name, s in SETS.items():
        om, e1, de = arb(s["omega"]), arb(s["eta1"]), arb(s["delta"])
        l = [arb(x) / 100 for x in s["lam"]]
        tot = arb(0)
        for i, g in enumerate(gs):
            if n0 + i > s["nmax"]: break
            w = (l[0] * khat(e1 * g, 1) + l[1] * khat(e1 * g, 2) + l[2] * khat(e1 * g, 3)) * khat(de * g, s["M"])
            rho = acb(arb("0.5"), g)
            tot += -2 * w * ((1 / rho + 1 / (rho * rho * om)) * acb(0, g * om).exp()).real
        out[name] = tot.str(25, radius=True)
    return out
def save_zeros(n0, gs):
    # exact double-double midpoint + outward radius: ball(hi+lo, rad) contains the Arb ball
    hi = np.empty(len(gs)); lo = np.empty(len(gs)); rd = np.empty(len(gs))
    for i, g in enumerate(gs):
        m = g.mid(); h = float(m); l = float(m - arb(h))
        hi[i], lo[i] = h, l
        rd[i] = float((g.rad() + abs(m - arb(h) - arb(l))).upper()) * (1 + 1e-15) + 1e-30
    np.savez(f"{ZD}/z_{n0:09d}.npz", hi=hi, lo=lo, rad=rd)
def high(args):
    n0, k = args
    p = f"{OUT}/h_{n0:09d}.json"
    zp = f"{ZD}/z_{n0:09d}.npz"
    if os.path.exists(p) and os.path.exists(zp): return p
    flint.ctx.dps = DPS
    zs = acb.zeta_zeros(n0, k)
    gs = []
    for z in zs:
        assert z.real == arb("0.5"), n0
        if gs: assert gs[-1] < z.imag, n0
        gs.append(z.imag)
    save_zeros(n0, gs)
    if os.path.exists(p):
        old = json.load(open(p)); assert arb(old["first"]).overlaps(gs[0]) and arb(old["last"]).overlaps(gs[-1]); return p
    res = {"n0": n0, "k": k, "first": gs[0].str(30, radius=True), "last": gs[-1].str(30, radius=True),
           "maxrad": max(float(g.rad()) for g in gs), "S": batch_sum(gs, n0)}
    json.dump(res, open(p + ".tmp", "w")); os.replace(p + ".tmp", p); return p
def low(args):
    path, n0, k, fmt = args
    p = f"{OUT}/l_{n0:09d}.json"
    if os.path.exists(p): return p
    flint.ctx.dps = DPS + 10
    gs = []
    for i, line in enumerate(open(path)):
        if i >= k: break
        gs.append(arb(line.split("\t")[0].strip()))
    assert len(gs) == k and all(a < b for a, b in zip(gs, gs[1:]))
    res = {"n0": n0, "k": k, "first": gs[0].str(30, radius=True), "last": gs[-1].str(30, radius=True),
           "maxrad": max(float(g.rad()) for g in gs), "S": batch_sum(gs, n0)}
    json.dump(res, open(p + ".tmp", "w")); os.replace(p + ".tmp", p); return p
def low_tasks():
    best = {}
    for d, kind in ((MER, "m"), (EXT, "e")):
        for p in glob.glob(d + "/*.txt"):
            n0, k, dps = map(int, os.path.basename(p)[:-4].split("_"))
            if n0 not in best or (k, dps) > best[n0][:2]: best[n0] = (k, dps, p, kind)
    tasks, cur = [], 1
    while cur <= N_LOW:
        k, dps, p, kind = best[cur]; k = min(k, N_LOW + 1 - cur)
        tasks.append((p, cur, k, kind)); cur += k
    return tasks
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True); os.makedirs(ZD, exist_ok=True)
    W = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    lt = low_tasks()
    ht = [(n0, min(2000, N_8E7 + 1 - n0)) for n0 in range(N_LOW + 1, N_8E7 + 1, 2000)]
    print(f"{len(lt)} low batches, {len(ht)} high batches", flush=True)
    t = time.time()
    with Pool(W) as pool:
        for i, _ in enumerate(pool.imap_unordered(low, lt, chunksize=8), 1):
            if i % 5000 == 0: print(f"low {i}/{len(lt)} {time.time()-t:.0f}s", flush=True)
        todo = [t for t in ht if not (os.path.exists(f"{OUT}/h_{t[0]:09d}.json") and os.path.exists(f"{ZD}/z_{t[0]:09d}.npz"))]
        todo.sort()   # increasing n, so zero files become usable for growing cutoffs T0
        for i, _ in enumerate(pool.imap(high, todo), 1):
            if i % 500 == 0: print(f"high {i}/{len(ht)} {time.time()-t:.0f}s", flush=True)
    print("DONE", time.time() - t, flush=True)
