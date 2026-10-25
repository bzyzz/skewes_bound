#!/usr/bin/env python3
"""Certify U(omega) < 0 for the race criterion in Arb, using the Turing-certified zeros.

usage: certify.py T0 omega eta1 delta M lam1,lam2,lam3 nproc
Reads data/arb_<q>_<n>/certified.npz and data/arb_<q>_<n>/turing_<q>_<n>.json. Writes criterion/cert_<tag>.json.
"""
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
import flint
from flint import arb, acb

import race_criterion as rc

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CHUNK = 20000


def zero_balls(q, n, T0):
    d = np.load(os.path.join(ROOT, "data", f"arb_{q}_{n}", "certified.npz"))
    tj = json.load(open(os.path.join(ROOT, "data", f"arb_{q}_{n}", f"turing_{q}_{n}.json")))
    rc.require(len(d["fail"]) == 0, f"{q}.{n}: uncertified zeros present")
    rc.require(tj["L_positive_half_to_two"], f"{q}.{n}: L > 0 on [1/2,2] not certified")
    cp = [c for c in tj["checkpoints"] if c.get("ok") and abs(c["T"] - T0) < 1e-9]
    rc.require(len(cp) == 1, f"{q}.{n}: no Turing certificate at T0 = {T0}")
    lo, hi = d["lo"], d["hi"]
    sel = hi <= T0
    rc.require(int(sel.sum()) == cp[0]["F_T"], f"{q}.{n}: zero count mismatch")
    rc.require(not np.any((lo < T0) & (hi > T0)), "a zero interval straddles T0")
    return lo[sel], hi[sel], cp[0]["F_T"]


def work(args):
    lo, hi, omega_s, eta1_s, delta_s, M, lam_s = args
    flint.ctx.dps = 30
    omega, eta1, delta = arb(omega_s), arb(eta1_s), arb(delta_s)
    lam = [arb(x) for x in lam_s]
    S = arb(0)
    s2 = arb(0)
    s3 = arb(0)
    for a, b in zip(lo, hi):
        g = arb((a + b)/2, (b - a)/2)       # contains [a, b]: (a+b)/2 and (b-a)/2 are exact or rounded outward by arb
        g = g.union(arb(a)).union(arb(b))
        core = sum(l*rc.kap_hat(eta1*g, m) for l, m in zip(lam, (1, 2, 3)) if l != 0)
        Kh = core*rc.kap_hat(delta*g, M)
        rho = acb(arb(1)/2, g)
        z = Kh*acb(0, g*omega).exp()*(1/rho + 1/(rho*rho*omega))
        S += 2*z.real                        # zeros come in pairs gamma, -gamma
        s2 += 2/(g*g)
        s3 += 2/(g*g*g)
    # arb is not picklable; Arb's decimal ball string encloses the ball, and parsing it back stays enclosing
    return S.str(40), s2.str(40), s3.str(40)


def main():
    T0 = float(sys.argv[1]); omega_s, eta1_s, delta_s = sys.argv[2], sys.argv[3], sys.argv[4]
    M = int(sys.argv[5]); lam_s = sys.argv[6].split(","); nproc = int(sys.argv[7])
    lam = [arb(x) for x in lam_s]
    rc.require(all(l >= 0 for l in lam) and (sum(lam)).contains(1), "convex weights")
    t0 = time.time()
    S_tot = arb(0)
    s2l, s3l, per = [], [], {}
    for (q, n, d, A, _) in rc.CHARS:
        lo, hi, NT0 = zero_balls(q, n, T0)
        jobs = [(lo[i:i + CHUNK], hi[i:i + CHUNK], omega_s, eta1_s, delta_s, M, lam_s) for i in range(0, len(lo), CHUNK)]
        with Pool(nproc) as p:
            res = [tuple(arb(x) for x in r) for r in p.map(work, jobs)]
        S = sum((r[0] for r in res), arb(0))
        s2 = sum((r[1] for r in res), arb(0)) + rc.sum_inv_pow_tail(arb(T0), q, 2, arb(2*NT0))
        s3 = sum((r[2] for r in res), arb(0)) + rc.sum_inv_pow_tail(arb(T0), q, 3, arb(2*NT0))
        S_tot += 2*S
        s2l.append(s2); s3l.append(s3)
        per[f"{q}.{n}"] = dict(N_pos=NT0, S_chi=S.str(15), sum_g2=s2.str(10), sum_g3=s3.str(10))
        print(f"{q}.{n}: {NT0} zeros, S_chi = {S.str(12)} ({time.time() - t0:.0f}s)", flush=True)
    err = rc.errors(arb(omega_s), arb(eta1_s), arb(delta_s), M, arb(T0), lam, s2l, s3l)
    U = S_tot + err["total"]
    eta0 = arb(eta1_s) + arb(delta_s)
    u_hi = arb(omega_s) + eta0
    ok = bool(U < 0)
    out = dict(T0=T0, omega=omega_s, eta1=eta1_s, delta=delta_s, M=M, lam=lam_s, chars=per,
               two_sum_S=S_tot.str(15), errors={k: v.str(10) for k, v in err.items()}, U=U.str(15), certified=ok,
               u_lo=(arb(omega_s) - eta0).str(15), u_hi=u_hi.str(15),
               bound_log10=float((u_hi/arb(10).log()).upper()) if ok else None)
    tag = f"T{int(T0)}_w{omega_s}"
    json.dump(out, open(os.path.join(HERE, f"cert_{tag}.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
    if ok:
        print(f"CERTIFIED: pi(x;24,11) < pi(x;24,1) for some x with log x in [{out['u_lo']}, {out['u_hi']}]; "
              f"first sign change < 10^{out['bound_log10']:.4f}")


if __name__ == "__main__":
    main()
