#!/usr/bin/env python3
"""Arb certificate for the team race mod q: L = -(phi-1) Phi_D(u_lo) - sum_chi S_chi - R3 - R4 - R5 - eps0 > 0
implies pi(x;q,1) > (sum of pi(x;q,b) over the other classes)/(phi-1) for some log x in [omega - eta0, omega + eta0].
usage: team_certify.py q T0 omega eta1 delta M lam1,lam2,lam3 nproc"""
import json, os, sys, time
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "criterion")); sys.path.insert(0, HERE)
from flint import arb
import certify as rcert, race_criterion as rc, team_criterion as tc
from certify import CHUNK

def main():
    q, T0 = int(sys.argv[1]), float(sys.argv[2]); omega_s, eta1_s, delta_s = sys.argv[3:6]
    M = int(sys.argv[6]); lam_s = sys.argv[7].split(","); nproc = int(sys.argv[8])
    lam = [arb(x) for x in lam_s]
    rc.require(all(l >= 0 for l in lam) and sum(lam).contains(1), "convex weights")
    cfg = tc.MODULI[q]; t0 = time.time()
    S_tot = arb(0); s2l, s3l, per = [], [], {}
    for (cq, cn, d, A, _) in cfg["chars"]:
        lo, hi, NT0 = rcert.zero_balls(cq, cn, T0)
        jobs = [(lo[i:i + CHUNK], hi[i:i + CHUNK], omega_s, eta1_s, delta_s, M, lam_s) for i in range(0, len(lo), CHUNK)]
        with Pool(nproc) as p:
            res = [tuple(arb(x) for x in r) for r in p.map(rcert.work, jobs)]
        S = sum((r[0] for r in res), arb(0))
        s2 = sum((r[1] for r in res), arb(0)) + rc.sum_inv_pow_tail(arb(T0), cq, 2, arb(2*NT0))
        s3 = sum((r[2] for r in res), arb(0)) + rc.sum_inv_pow_tail(arb(T0), cq, 3, arb(2*NT0))
        S_tot += S; s2l.append(s2); s3l.append(s3)
        per[f"{cq}.{cn}"] = dict(N_pos=NT0, S_chi=S.str(15), sum_g2=s2.str(10), sum_g3=s3.str(10))
        print(f"{cq}.{cn}: {NT0} zeros, S_chi = {S.str(12)} ({time.time() - t0:.0f}s)", flush=True)
    err = tc.errors(q, arb(omega_s), arb(eta1_s), arb(delta_s), M, arb(T0), lam, s2l, s3l)
    L = -S_tot - err["total"]
    eta0 = arb(eta1_s) + arb(delta_s); u_lo, u_hi = arb(omega_s) - eta0, arb(omega_s) + eta0
    ok = bool(L > 0)
    out = dict(q=q, T0=T0, omega=omega_s, eta1=eta1_s, delta=delta_s, M=M, lam=lam_s, chars=per,
               sum_S=S_tot.str(15), errors={k: v.str(10) for k, v in err.items()}, L=L.str(15), certified=ok,
               u_lo=u_lo.str(15), u_hi=u_hi.str(15), bound_log10=float((u_hi/arb(10).log()).upper()) if ok else None)
    json.dump(out, open(os.path.join(HERE, f"team_cert_q{q}_T{int(T0)}_w{omega_s}.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))
    if ok:
        print(f"CERTIFIED: pi(x;{q},1) exceeds the average of the other classes for some x with log x in [{out['u_lo']}, {out['u_hi']}]; "
              f"first such x < 10^{out['bound_log10']:.4f}")

if __name__ == "__main__":
    main()
